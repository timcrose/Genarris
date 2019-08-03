# -*- coding: utf-8 -*-

import networkx as nx
import numpy as np 
import pandas as pd

from ibslib.io import read,write
from ibslib.molecules import MoleculeBonding


class BondNeighborhood():
    """
    Returns the bonding neighborhood of each atom for a structure. User is 
    allowed to define a radius that the algorithm traverses to build 
    the neighborhood for each atom. If the radius is 0, this would 
    correspond to just return the atoms im the system.
    
    Arguments
    ---------
    radius: int
        Radius is the number of edges the model is allowed to traverse on the 
        graph representation of the molecule. 
    mb_kwargs: dict
        Keyword arguments for the molecule bonding module. The default values
        are the recommended settings. A user may potentially want to decrease
        the natural_cutoff_mult. This value is multiplied by covalent bond 
        radii in the MoleculeBonding class. It's highly recommended that the
        skin value is kept at 0. 
        For more details see ibslib.molecules.MoleculeBonding
        
    """
    def __init__(self, radius=1, 
                       mb_kwargs={"natural_cutoff_mult": 1.2,
                                  "skin": 0}
                ):
        self.radius = radius
        self.mb_kwargs = mb_kwargs
        if radius != 1:
            raise Exception("Radius greater than 1 not implemented")
        
    
    def calc(self, struct):
        self.struct = struct
        self.ele = struct.geometry["element"]
        g = self._build_graph(struct)
        n = self._calc_neighbors(g)
        n = self._sort(g,n)
        fragments = self._construct_fragment_list(n)
        fragments,count = np.unique(fragments, return_counts=True)
        
        return fragments.tolist(),count
        
    
    def _build_graph(self, struct):
        """
        Builds networkx graph of a structures bonding.
        """
        
        mb = MoleculeBonding(struct, **self.mb_kwargs)
        g = nx.Graph()
        
        # Add node to graph for each atom in struct
        self.ele = struct.geometry["element"]
        g.add_nodes_from(range(len(self.ele)))
        
        # Add edges
        for i,bond_list in enumerate(mb.bonding_list):
            [g.add_edge(i,x) for x in bond_list]
        
        return g
    
    
    def _calc_neighbors(self, g):
        """
        Calculates neighbors for each node in the graph. Uses the radius 
        which was declared when BondNeighborhood was initialized. Ordering
        of the neighbors follows:
            1. Terminal atom alphabetical 
            2. Self
            3. Bonded atom alphabetical
            4. If continuing from bonded group, print terminal groups of the 
               bonded group. 
            5. If there's ambiguity about which should come next, place in 
               alphabetical order. When the radius is small, there's less
               ambiguity. If the radius becomes large there will be more. 
               Although, only a radius of 1 is currently implemented.
            
        """
        neighbors = [[[]] for x in g.nodes]
        
        for i,idx_list in enumerate(neighbors):
            neighbor_list = [x for x in g.adj[i]]
            neighbor_list.append(i)
            idx_list[0] += neighbor_list
        
        return neighbors
    
    
    def _sort(self, g, neighbor_list):
        """
        Sorts neighborlist according to definition in _calc_neighbors. Only 
        works for a radius of 1.
        
        Arguments
        ---------
        g: nx.Graph
        neighbor_list: list of int
            List of adjacent nodes plus the node itself as i
        
        """
        sorted_list_final = [[[]] for x in g.nodes]
        for i,temp in enumerate(neighbor_list):
            # Preparing things which aren't writting well 
            idx_list = temp[0]
            current_node = i
            
            terminal_groups = []
            bonded_groups = []
            for idx in idx_list:
                if g.degree(idx) == 1:
                    terminal_groups.append(idx)
                else:
                    bonded_groups.append(idx)
            
            terminal_ele = self.ele[terminal_groups]
            alphabet_idx = np.argsort(terminal_ele)
            terminal_groups = [terminal_groups[x] for x in alphabet_idx]
            
            sorted_list = terminal_groups
            if current_node not in terminal_groups:
                sorted_list.append(current_node)
                remove_idx = bonded_groups.index(current_node)
                del(bonded_groups[remove_idx])
            
            bonded_ele = self.ele[bonded_groups]
            alphabet_idx = np.argsort(bonded_ele)
            bonded_groups = [bonded_groups[x] for x in alphabet_idx]
            
            sorted_list += bonded_groups
            
            sorted_list_final[i][0] = sorted_list
        
        return sorted_list_final
    
    
    def _construct_fragment_list(self,n):
        fragment_list = [self.struct.geometry["element"][tuple(x)] for x in n]
        fragment_list = [''.join(x) for x in fragment_list]
        return fragment_list


class construct_bond_neighborhood_model():
    """
    Manages the construction of a regularized bond neighborhood model. Can 
    be used for any user defined target value. Current use case is for the 
    prediction of solid state volumes of molecules.
    
    Usage
    -----
    0. Initialize a BondNeighborhood class. Initialized the constructor
       class with the BondNeighborhood and the list of extra properties 
       the user would like to include in the model as add_prop and the
       target property as target_prop.
    1. Call set_train with a Structure dictionary to construct the 
       neighborhood_df and target_df
    2. Call regularize to reduce the number of neighborhoods in the model.
    3-1. Can call get_feature_vector(struct) to get the feature vector for the 
         input structure which is correctly ordered with respect to the 
         regularized neighborhood dataframe. 
    3-2. Alternatively, can call set_test to build a neighborhood_df for a 
         test dataset and the target_df_t for the test dataset. 
         This neighborhood_df will have the same columns as the neighborhood_df 
         constructed by set_train.
    4. Call fit_model with a regressor object. This will return a dataframe 
       for the training set and testing set which contains the target values 
       and the predicted values. 
    
    Arguments
    ---------
    bn: BondNeighborhood object
        Initialized BondNeighborhood object.
    add_prop: list of str
        List of properties stored in the structures which the user wants to 
        add to the feature set.
    target_prop: list of str
        List of properties stored in the structures which the user wants to 
        use as the target values.
        
        
    """
    def __init__(self, bn, add_prop=["MC_volume"], target_prop=["molecule_volume"]):
        self.bn = bn
        self.add_prop = add_prop
        self.target_prop = target_prop
        
        ##### Initialize internals
        self.neighborhood_df = pd.DataFrame()
        self.neighborhood_df_r = pd.DataFrame() # Regularized neighborhood df
        self.total_neighborhoods_observed = 0
        # neighborhood df built for testing dataset populated by set_test
        self.neighborhood_df_t = pd.DataFrame() 
        
        self.target_df = pd.DataFrame() # Target df for training set
        self.target_df_t = pd.DataFrame() # Target df to testing set
        
        
    def fit_model(self, regressor, r=True, test=True):
        """
        Fits and tests the model. Regressor object must have a regressor.fit 
        method. 
        
        Arguments
        ---------
        regressor: object
            Regressor object which has a regressor.fit method. This can be
            any regressor from sklearn.
        r: bool
            True: The regularized self.neighborhood_df_r will be used
            False: The nonregularized self.neighborhood_df will be used
        test: bool
            True: Use the testing dataset
            False: Do not use the testing dataset
            
        """
        
        self._check_r_setting(r,"fit_model")
        
        # Choose training df based on r argument
        if r == True:
            temp_df = self.neighborhood_df_r
        else:
            temp_df = self.neighborhood_df
        
        if test and temp_df.shape[1] != self.neighborhood_df_t.shape[1]:
            raise Exception("The number of features of the training dataset"+
                    "is {} ".format(temp_df.shape[1]) + 
                    "which does not match the number of features "+
                    "in the testing dataset, which is {}."
                    .format(self.neighborhood_df_t.shape[1]))
            
        regressor.fit(temp_df.values, self.target_df.values)
        train_pred = regressor.predict(temp_df.values)
        
        #### Construct dataframe to return for the training dataframe
        pred_columns = ["predicted_"+ x for x in 
                        constructor.target_df.columns.tolist()]
        df_columns = self.target_df.columns.tolist() + pred_columns
        train_df = pd.DataFrame(data=np.zeros((self.target_df.shape[0],
                                               self.target_df.shape[1]*2)),
                                index=self.target_df.index,
                                columns=df_columns)
        
        self.train_df = train_df
        
        train_df.iloc[:,0:self.target_df.shape[1]] = self.target_df.values
        train_df.iloc[:,self.target_df.shape[1]:] = train_pred
        
        #### If using a testing dataset, calculate and construct dataframe
        if test:
            test_pred = regressor.predict(self.neighborhood_df_t.values)
            
            # Constructing test dataframe
            test_df = pd.DataFrame(data=np.zeros((self.target_df_t.shape[0],
                                               self.target_df_t.shape[1]*2)),
                                index=self.target_df_t.index,
                                columns=df_columns)
            test_df.iloc[:,0:self.target_df_t.shape[1]] = self.target_df_t.values
            test_df.iloc[:,self.target_df_t.shape[1]:] = test_pred
            
            return train_df,test_df
        else:
            return train_df
         
    
    def set_train(self, struct_dict):
        """
        Set training struct dict and builds an internal neighborhood dataframe.
        First collects all fragments from the dataset. Then builds a dataframe
        using only the unique fragments observed from the dataset. 
        
        Arugments
        ---------
        struct_dict: StructDict
            Structure dictionary for which to train the model on.
        add_prop: list of str
            List of any additional properties to add to the features set.
            
        """
        # Reinitialized internals when set_train is called
        self.neighborhood_df = pd.DataFrame()
        self.total_neighborhoods_observed = 0
        self.neighborhood_df_r = pd.DataFrame()
        
        # First get fragment and counts for each structure
        temp_columns = ["fragments", "counts"]
        temp_columns += self.add_prop
        struct_keys = [x for x in struct_dict.keys()]
        temp_df = pd.DataFrame(data=np.zeros((len(struct_dict),2+len(self.add_prop)),
                                             dtype=object),
                               index=struct_keys,
                               columns=temp_columns)
        # Store properties added for future reference neighborhood_df construction
        temp_add_prop_dict = {}
        
        for struct_id,struct in struct_dict.items():
            f,c = self.bn.calc(struct)
            struct.set_property("bond_neighborhood_fragments", f)
            struct.set_property("bond_neighborhood_counts", c)
            
            values = [np.array(f),c]
            if len(self.add_prop) > 0:
                prop_values = self._get_add_prop(struct)
                temp_add_prop_dict[struct_id] = prop_values
                values += prop_values
            
            temp_df.loc[struct_id] = np.array(values, dtype=object)
        
        # Building final neighborhood dataframe from unique bond neighborhoods
        all_fragments = np.hstack(temp_df["fragments"].values)
        unique_fragments = np.sort(np.unique(all_fragments))
        temp_columns = np.array([x for x in self.add_prop], dtype=object)
        temp_columns = np.hstack((temp_columns, unique_fragments))
        
        self.neighborhood_df = pd.DataFrame(
                                data=np.zeros((len(struct_dict),
                                               len(temp_columns))),
                                index = struct_keys,
                                columns = temp_columns)
        
        # Populate neighborhood_df                       
        for struct_id,struct in struct_dict.items():
            f = struct.properties["bond_neighborhood_fragments"]
            c = struct.properties["bond_neighborhood_counts"]
            self.neighborhood_df.loc[struct_id, f] = c
            self.neighborhood_df.loc[struct_id,self.add_prop] = temp_add_prop_dict[struct_id]
            
            # Updating the neighborhoods observed 
            self.total_neighborhoods_observed += np.sum(c)
        
        # Get target df for training
        self.target_df = self._get_target_df(struct_dict)
    
    
    def set_test(self, struct_dict, r=True):
        """
        Creates a neighborhood dataframe for the struct_dict which has the same
        columns as the internal neighborhood_df or neighborhood_df_r (the
        regularized dataframe). 
        
        Arugments
        ---------
        struct_dict: StructDict
            Dictionary indexed by structure ids with Structure obejct as the 
            values. 
        r: bool
            True: self.neighborhood_df_r will be used
            False: self.neighborhood_df will be used
            
        """
        # Check if self.neighborhood_df_r has been constructed
        self._check_r_setting(r,"set_test")
                
        if r == False:
            temp_df = self.neighborhood_df
        else:
            temp_df = self.neighborhood_df_r
            
        
        self.neighborhood_df_t = pd.DataFrame(
                 data = np.zeros((len(struct_dict),
                                  temp_df.columns.shape[0])),
                 index = [x for x in struct_dict.keys()],
                 columns = temp_df.columns.values)
        
        # Populate testing dataframe
        for struct_id,struct in struct_dict.items():
            nv = self.get_neighborhood_vector(struct,r=r)
            self.neighborhood_df_t.loc[struct_id] = nv
        
        # Populate target dataframe for test set
        self.target_df_t = self._get_target_df(struct_dict)
        
            
    def regularize(self, tol=10):
        """
        Reduces the neighborhood_df to only those neighborhoods which have 
        been observed >= tol number of times in the training dataset.
        
        Arguments
        ---------
        tol: int or float
            If int, a neighborhood must be observed more than tol time in order 
            for it to be kept in the neighborhood_df_r.
            If float, a neighborhood must make up a this fraction of all 
            observed neighborhoods. This tolerance value should be extremely
            small such as 0.001 or 0.0001.
        
        """
        if len(self.neighborhood_df) == 0:
            raise Exception("Please call "+
                    "construct_bond_neighborhood_model.set_train before you "+
                    "call regularize.")
        drop_list = []
        for name in self.neighborhood_df.columns:
            # Find where columns are nonzero
            non_zero = np.sum(self.neighborhood_df[name] > 0) 
            
            if type(tol) == int:
                if non_zero < tol:
                    drop_list.append(name)
            elif type(tol) == float:
                if non_zero / self.total_neighborhoods_observed < tol:
                    drop_list.append(name)
        
        self.neighborhood_df_r = self.neighborhood_df.drop(drop_list,axis=1)
    
    
    def get_neighborhood_vector(self, struct, r=True):
        """
        Returns the neighborhood vector commensurate for either the regularized 
        neighborhood dataframe or the neighborhood dataframe from set_train.
        
        """
        # Check if self.neighborhood_df_r has been constructed
        self._check_r_setting(r,"get_neighborhood_vector")
                
        if r == False:
            temp_df = self.neighborhood_df
        else:
            temp_df = self.neighborhood_df_r
        
        fragment_list = temp_df.columns.values
        
        # Initialize feature vector for structure
        neighborhood_vector = np.zeros(fragment_list.shape)
        
        # Calc structure fragments and counts
        f,c = self.bn.calc(struct)
        # Check which fragments in struct are in fragment list and return
        # the index of their location
        f_idx,c_idx = np.nonzero(fragment_list[:,None] == f)
        
        # Populate feature vector only with neighborhoods observed in temp_df
        neighborhood_vector[f_idx] = c[c_idx]
        
        #### Finish with add_prop values
        prop_values = self._get_add_prop(struct)
        p_idx = np.nonzero(fragment_list[:,None] == self.add_prop)[0]
        neighborhood_vector[p_idx] = prop_values
        
        return neighborhood_vector
    
    
    def _get_add_prop(self, struct):
        """
        Small function for obtaining the add_prop values of the input 
        structure. 
        
        """
        prop_values = []
        for prop_name in self.add_prop:
            temp_prop_value = struct.get_property(prop_name)
            
            # Define behavior of not finding the property
            if not temp_prop_value:
                temp_prop_value = 0
                    
            prop_values.append(temp_prop_value)
        return prop_values
    
    
    def _get_target_df(self, struct_dict):
        """
        Returns dataframe of the target properties for the struct_dict
        
        """
        target_df = pd.DataFrame(
                     data=np.zeros((len(struct_dict), len(self.target_prop))),
                     columns=self.target_prop,
                     index=[x for x in struct_dict.keys()])
        
        for struct_id,struct in struct_dict.items():
            prop_values = []
            for prop_name in self.target_prop:
                temp_prop_value = struct.get_property(prop_name)
            
                # Define behavior of not finding the property
                if not temp_prop_value:
                    temp_prop_value = 0
                
                prop_values.append(temp_prop_value)
                
            target_df.loc[struct_id] = prop_values
        
        return target_df
    
    
    def _check_r_setting(self,r,method_name):
         # Check if self.neighborhood_df_r has been constructed
        if len(self.neighborhood_df_r) == 0:
            if len(self.neighborhood_df) == 0:
                raise Exception("The user called "+
                   "construct_bond_neighborhood_model.{} ".format(method_name) +
                   "before calling set_train. " +
                   "Please call set_train, and " +
                   "optionally regularize, before calling get_feature_vector.")
            elif r == True:
                raise Exception("The user called "+
                    "construct_bond_neighborhood_model.{} ".format(method_name) +
                    "with r=True before calling regularize. " +
                    "Either call regularize prior to get_neighborhood_vector "+
                    "or set r=False.")
        
        

if __name__ == "__main__":
#    from sklearn.linear_model import Ridge
#    regressor = Ridge()
#    
#    bn = BondNeighborhood()
#    constructor = construct_bond_neighborhood_model(bn,add_prop=["MC_volume"], 
#                                              target_prop=["molecule_volume"])
#    constructor.set_train(s_conquest)
#    constructor.regularize(tol=40)
#    constructor.set_test(s_paper)
    
    train_df,test_df = constructor.fit_model(regressor)
    train_err = np.abs(train_df["molecule_volume"] - train_df["predicted_molecule_volume"]) / train_df["predicted_molecule_volume"]
    test_err = np.abs(test_df["molecule_volume"] - test_df["predicted_molecule_volume"]) / test_df["predicted_molecule_volume"]