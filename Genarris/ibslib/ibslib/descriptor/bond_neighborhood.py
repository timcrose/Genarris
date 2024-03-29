# -*- coding: utf-8 -*-

import networkx as nx
import numpy as np 
import pandas as pd
from scipy.spatial.distance import pdist,squareform
from scipy.spatial.distance import __all__ as implemented_metrics

from ibslib.io import read,write
from ibslib.molecules.find_bonding import MoleculeBonding

import matplotlib.pyplot as plt


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
                       mb_kwargs={"natural_cutoff_mult": 1.3,
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
    Manages the construction of a regularized bond neighborhood model. Also
    has many functions which can be used to perform analysis of the learned
    model or on the dataset which was used, such as computing similarity of 
    molecules across the dataset. Models can be trained using any user defined 
    target value. Current use case is for the prediction of solid form volumes 
    of molecules from the Cambridge Structural Database. 
    
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
       
     
    Model Analysis
    --------------
    Functions which perform model analysis.
    
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
                        self.target_df.columns.tolist()]
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
        temp_df = pd.DataFrame(data=np.zeros((len(struct_dict),
                                              2+len(self.add_prop)),
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
    
    
    def model_complexity_analysis(self, regressor, choose_n_features=-1,
                                  cross_validation=0,
                                  neighborhood_df_train=pd.DataFrame(),
                                  target_train = pd.DataFrame(),
                                  neighborhood_df_test=pd.DataFrame(),
                                  target_test = pd.DataFrame()):
        """
        Generates results testing how the training and testing error
        change as a function of the number of features included in the model, 
        otherwise called the model complexity. Idea here is that as model
        complexity increases, the training error will decrease monoatomically,
        however, at some point, the testing error should increase indicating 
        overfitting of the model. Features are added to the model in a forward,
        greedy way. 
        
        Arguments
        ---------
        regressor: object
            Supervised learning estimator with a fit method that provides
            information about the feature importance either through a coef_ or
            feature_importances_ attribute. 
        cross_validation: int
            Number of cross validation calculations to be performed for each
            point of the model complexity testing. If set to 0 or 1, no
            cross validation will be used. If greater than 1, then a standard
            deviation of the predicted error can be computed. 
        choose_n_features: 
            Number of features which will be chosen greedily. 
        neighborhood_df_train/test: pd.DataFrame
            Allows users to pass in their desired dataframe. Otherwise, resorts
            to the default behavior, which is to use the dataframes stored in 
            self.neighborhood_df/self.neighborhood_df_t for train/test.
            
        """
        if len(neighborhood_df_train) == 0:
            if len(self.neighborhood_df) == 0:
                raise Exception("Called constructor.model_complexity_analysis "+
                        "without setting a training dataset.")
            neighborhood_df_train = self.neighborhood_df
            target_train = self.target_df.values
        else:
            pass
        if len(neighborhood_df_test) == 0:
            if len(self.neighborhood_df_t) == 0:
                raise Exception("Called constructor.model_complexity_analysis "+
                        "without setting a testing dataset." +
                        "Note that one could use the same training and testing "+
                        "dataset by calling set_train with the trainig set.")
            neighborhood_df_test = self.neighborhood_df_t
            target_test = self.target_df_t
        else:
            pass
        
        # Defines default behavior for choose_n_features
        if choose_n_features <= 0:
            choose_n_features = len(neighborhood_df_train)
        if choose_n_features > len(neighborhood_df_train):
            choose_n_features = len(neighborhood_df_train)
            
        # Build greedy feature model
        greedy_features = pd.DataFrame()
        greedy_model_err = np.zeros((choose_n_features,2))
        train = neighborhood_df_train
        target = target_train
        test = neighborhood_df_test
        target_t = target_test.values
        for n in range(0,choose_n_features):
            # Reset greedy values to begin process of adding a new feature
            model_err = []
            if n == 0:
                greedy_values = np.zeros((len(train.index),1))
            else:
                # Add column for new feature
                greedy_values = np.zeros((greedy_features.values.shape[0],
                                          greedy_features.values.shape[1]+1))
                greedy_values[:,0:-1] = greedy_features.values
                
            for j,feature in enumerate(train.columns):
                print(n, j,len(train.columns),feature)
                train_temp = train.iloc[:,j].values
                
                # Add new feature from training set to greedy values
                greedy_values[:,-1] = train_temp
                
                regressor.fit(greedy_values,target)
                pred = regressor.predict(greedy_values)
                err = np.mean(np.abs(target - pred) / pred)
                model_err.append(err)
            
            # Choose best new feature to add
            best_new_feature_idx = np.argmin(model_err)
            feature_name = train.columns[best_new_feature_idx]
            feature_values = train[feature_name].values
            
            # Remove this feature from future consideration in the train set
            train = train.drop(feature_name,axis=1)
            
            if n == 0:
                greedy_features = pd.DataFrame(data=feature_values,
                                               index=train.index,
                                               columns=[feature_name])
            else:
                greedy_features[feature_name] = feature_values
            
            # Now use testing set with the new greedy features
            regressor.fit(greedy_features.values,target)
            pred_t = regressor.predict(test.loc[:,greedy_features.columns])
            err_t = np.mean(np.abs(target_t - pred_t) / pred_t)
            
            # Store error
            greedy_model_err[n,0] = model_err[best_new_feature_idx]
            greedy_model_err[n,1] = err_t
            
            print(greedy_model_err[n,:])
            
        return greedy_features,greedy_model_err
        
    
    
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
    
    
    def plot_results(self,result_df):
        """
        Quickly plots the results from the constructed model for the results
        dataframe.
        """
        pass
    
    
    def plot_hist(self, neighborhood_df, 
                  exclude_add_prop=True, 
                  figname='',
                  regressor=None, 
                  most_important=-1, 
                  add_coef=False,
                  add_coef_text_kwargs =  {
                                            "whitespace_pad": 1.15,
                                            "labelpad_y": 10,
                                            "labelpad_x": -0.065,
                                            "fontsize": 16,
                                            "color": "tab:red",                                            
                                          },
                  add_obs=False,
                  add_obs_text_kwargs =   {
                                            "whitespace_pad": 1.15,
                                            "labelpad_y": 35,
                                            "labelpad_x": -0.075,
                                            "fontsize": 16,
                                            "color": "tab:purple"
                                          },
                  figsize=(12,8),
                  tick_params_kwargs_both={
                                             "axis": "both",
                                             "width": 3,
                                             "labelsize": 16,
                                             "length": 8,
                                           },
                  tick_params_kwargs_x =  {
                                             "axis": "x",
                                             "labelsize": 16,
                                             "labelrotation": 90,
                                          },
                  y_label_kwargs =        {
                                             "ylabel": "Number of Observations",
                                             "labelpad": 0,
                                             "fontsize": 18,
                                             "labelpad": 15,
                                          },
                  bar_kwargs =            {
                                            "edgecolor": "k",
                                            "linewidth": 2,
                                          },
                  ):
        """
        Plots a histogram of the frequency of bond neighborhoods identified
        in the input neighborhood dataframe
        
        Arguments
        ---------
        neighborhood_df: pandas.DataFrame
            Neighborhood dataframe from which to construct the histogram.
        exclude_add_prop: bool
            Whether to use columns from add prop in the histogram
        figname: str
            If a figname is provided, then the file is saved using this 
            string.
        Regressor: regressor object
            If a regressor is provided, then the histogram is constructed using
            only the N first most important bond neighborhoods.
        most_important: int
            If a values greater than 0 is provided, then the regressor must 
            also be provide. Uses on the N most_important features in 
            construction of historgram if greater than 0.
        add_coef: bool
            Adds the coefficient values from the regressor to the top of each
            bar plot.
        add_coef_text_kwargs: dict
            Keyword arguments for the ax.text call for the coeficient values.
        add_obs: bool
            Adds the number of observations for each plotted neighborhood
            to the top of each bar plot.
        add_obs_text_kwargs: dict
            Keyword arguments for the ax.text call for the observation values.
        figsize: (int,int)
            Figure size to pass into matplotlib.
        tick_params_kwargs_both: dict
            Dictionary of keyword arguments to pass to ax.tick_params
        tick_params_kwargs_both: dict
            Another dictionary of keyword arguments to pass to ax.tick_params.
            Idea here is that the *_both is used to format both axes while 
            this dictionary is used to format the x axis.
        y_label_kwargs: dict
            Dictionary of keyword arguments to control the label for the 
            y axis of the plot. 
        bar_kwargs: dict
            Dictionary of keyword arguments to control aesthetics of the bar
            plot such as color, opacity, and linewidth.
        
        """
        if exclude_add_prop:
            column_idx = np.arange(len(self.add_prop), 
                                   len(neighborhood_df.columns),
                                   1)
        else:
            column_idx = np.arange(0,
                                   len(neighborhood_df.columns),
                                   1)
        
        if most_important > 0:
            if not regressor:
                raise Exception("Argument most_important to plot_hist was {}."
                        .format(most_important) +
                        "However a regressor was not provided. Please provide "+
                        "A trained regressor object from sklearn.")
            
            if len(regressor.coef_.ravel()) != len(neighborhood_df.columns.values):
                raise Exception("The length of the regressor coeficients was {} "
                        .format(len(regressor.coef_)) +
                        "which does not match the number of columns of the "+
                        "the input neighborhood dataframe {}"
                        .format(len(neighborhood_df.columns.values)))
                
            coef = np.abs(regressor.coef_.ravel()[column_idx])
            important_idx = np.argsort(coef)[::-1]
            important_idx = important_idx[0:most_important]
            
            # Now reset column_idx to only be these most important indices
            column_idx = column_idx[important_idx]
                       
        column_names = neighborhood_df.columns.values[column_idx]
        sub_array = neighborhood_df.values[:,column_idx]
        column_totals = np.sum(sub_array, axis=0)
        
        sort_idx = np.argsort(column_totals)[::-1]
        column_idx = column_idx[sort_idx]
        column_totals = column_totals[sort_idx]
        column_names = column_names[sort_idx]
        x = np.arange(0,len(column_names),1)
        
        #### Making figure
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        ax.bar(x, column_totals, **bar_kwargs)
        ax.set_xticks(x)
        ax.set_xticklabels(column_names)
        ## Edit ticks
        ax.tick_params(**tick_params_kwargs_both)
        ax.tick_params(**tick_params_kwargs_x)
        
        # ylabel kwargs
        ax.set_ylabel(**y_label_kwargs)
        
        # Set outside of plot same size as tick_width
        outline_width = tick_params_kwargs_both["width"]
        ax.spines['top'].set_linewidth(outline_width)
        ax.spines['right'].set_linewidth(outline_width)
        ax.spines['bottom'].set_linewidth(outline_width)
        ax.spines['left'].set_linewidth(outline_width)
        
        if add_coef:
            # Increase the height of the whitespace in the plot
            whitespace_mult = add_coef_text_kwargs.pop("whitespace_pad")
            max_obs = np.max(column_totals)
            max_y = max_obs*whitespace_mult
            ax.set_ylim([0,max_y])
            
            # Get distance to add to each text y position
            labelpad_y = add_coef_text_kwargs.pop("labelpad_y")
            # Get distance for each text x position depending on the length
            # of the string. Because its a function of the length, this helps
            # to center the string perfectly.
            labelpad_x = add_coef_text_kwargs.pop("labelpad_x")
            for i,value in enumerate(regressor.coef_.ravel()[column_idx]):
                print_string = "{:.2f}".format(value)
                pos_x = x[i] + len(print_string)*labelpad_x
                pos_y = column_totals[i] + labelpad_y
                ax.text(pos_x,pos_y,print_string,
                        **add_coef_text_kwargs)
        
        # Now add coef if true
        if add_obs:
            # Increase the height of the whitespace in the plot
            whitespace_mult = add_obs_text_kwargs.pop("whitespace_pad")
            max_obs = np.max(column_totals)
            max_y = max_obs*whitespace_mult
            ax.set_ylim([0,max_y])
            
            # Get distance to add to each text y position
            labelpad_y = add_obs_text_kwargs.pop("labelpad_y")
            # Get distance for each text x position depending on the length
            # of the string. Because its a function of the length, this helps
            # to center the string perfectly.
            labelpad_x = add_obs_text_kwargs.pop("labelpad_x")
            for i,value in enumerate(column_totals):
                print_string = "{}".format(int(value))
                pos_x = x[i] + len(print_string)*labelpad_x
                pos_y = value + labelpad_y
                ax.text(pos_x,pos_y,print_string,
                        **add_obs_text_kwargs)
                
                
        plt.tight_layout()
        if len(figname) > 0:
            fig.savefig(figname)
        
    
    def similarity_matrix(self,neighborhood_df,metric="TC",
                          exclude_add_prop=True,
                          low_memory=False):
        """
        Constructs similarity matrix using the input neighborhood_df and 
        the input metric. 
        
        neighborhood_df: pandas.DataFrame
            Dataframe created by the constructor class. 
        metric: 
            Can be any scipy pdist metric or TC. 
        low_memory: bool
            Recommended to be true if using for a large number of molecules 
            (> 1000). If True, uses a low memory method should be used in the 
            computation of the TC metric. 
            
        """
        if "TC" not in implemented_metrics:   
            implemented_metrics.append("TC")
        if metric not in implemented_metrics:
            raise Exception("User used metric {} for similarity matrix. "
                    .format(metric) +
                    "Please use on of the implemented metrics: {}."
                    .format(implemented_metrics))
            
        if exclude_add_prop:
            column_idx = np.arange(len(self.add_prop), 
                                   len(neighborhood_df.columns),
                                   1)
        else:
            column_idx = np.arange(0,
                                   len(neighborhood_df.columns),
                                   1)
            
        values = neighborhood_df.values[:,column_idx]
            
        # Using any metric which is implemented for scipy pdist
        if metric != "TC":
            diff_matrix = squareform(pdist(values, metric=metric))
        
        # Using Tanimato Coefficient to compute similarity
        else: 
            diff_matrix = np.zeros((len(neighborhood_df.index),
                                    len(neighborhood_df.index)))
            if not low_memory:
                # First find NxN pairwise matrix of the number of fragments for 
                # each pair of molecules
                total_fragments = np.sum(values, axis=1)
                total_fragments_matrix = total_fragments + total_fragments[:,None]
                
                # Find where both have a nonzero value for the fragment
                nonzero_vectors = np.logical_and(values,values[:,None])
                # Find which has the minimum value for nonzero because that's  
                # the maximum number which can actually match 
                max_of_fragments = np.minimum(values,values[:,None])
                                
                # Mask maximum fragments by the nonzero mask
                c_array = np.ma.masked_array(max_of_fragments, 
                                             mask=nonzero_vectors, 
                                             fill_value=0).data
                
                # Now some the matching fragments for each pair of structure
                c_array = np.sum(c_array,axis=-1)
                
                # Compute the TC matrix as the diff_df
                diff_matrix = c_array / (total_fragments_matrix - c_array)
            
            # Low memory version follows the same steps as the above code
            # but loops over rows of the value matrix instead of computing
            # all similarities in memory at the same time.
            else:
                
                for i,row in enumerate(values):
                    print(i)
                    row = row[None,:]
                    total_fragments = np.sum(row)
                    # Get pairwise number of total fragments between the row
                    # and the rest of the molecules
                    total_fragments_matrix = total_fragments + np.sum(values,-1)
                    
                    # Find where the row and the other molecules have matching
                    # nonzero fragment values
                    nonzero_vectors = np.logical_and(row,values)
                                        
                    max_of_fragments = np.minimum(row,values)
                    c_array = np.ma.masked_array(max_of_fragments, 
                                             mask=nonzero_vectors, 
                                             fill_value=0).data
                                                 
                    c_array = np.sum(c_array,axis=-1)
                    if np.sum(c_array) == 0:
                        diff_matrix[i,:] = np.zeros(values.shape[0])
                    else:
                        diff_matrix[i,:] = c_array / (total_fragments_matrix - 
                                                      c_array)
                
        diff_df = pd.DataFrame(data=diff_matrix,
                               columns=neighborhood_df.index,
                               index=neighborhood_df.index)
        
        return diff_df
    
    
    def plot_similarity(self, diff_df, figname='',figsize=(10,10)):
        """
        Plot a similarity matrix as a color matrix. 
        """
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        im = ax.imshow(diff_df.values,interpolation='nearest')
        fig.colorbar(im)
        
        if len(figname) > 0:
            fig.savefig(figname)
            
        
        

if __name__ == "__main__":
    pass 
    
    ########## Initialization of model
#    from sklearn.linear_model import Ridge
#    regressor = Ridge()
#    
#    bn = BondNeighborhood()
#    constructor = construct_bond_neighborhood_model(bn,add_prop=["MC_volume"], 
#                                              target_prop=["molecule_volume"])
#    constructor.set_train(s_conquest)
#    constructor.regularize(tol=30)
#    constructor.set_test(s_paper)
#    
#    train_df,test_df = constructor.fit_model(regressor)
#    train_err = np.abs(train_df["molecule_volume"] - train_df["predicted_molecule_volume"]) / train_df["predicted_molecule_volume"]
#    test_err = np.abs(test_df["molecule_volume"] - test_df["predicted_molecule_volume"]) / test_df["predicted_molecule_volume"]
#    
#    constructor.regularize(tol=40)
#    constructor_test = construct_bond_neighborhood_model(bn,add_prop=["MC_volume"], 
#                                              target_prop=["molecule_volume"])
#    
#    
    ######### Plotting a histogrtam
#    constructor_test.plot_hist(constructor.neighborhood_df_r, 
#                              regressor=regressor,
#                              most_important=10,
#                              add_coef=True,
#                              add_obs=True)
    
    ########## Constructing a similarity matrix
#    diff_df = constructor_test.similarity_matrix(constructor.neighborhood_df_r,
#                                                 metric="TC",
#                                                 low_memory=True)
#    constructor_test.plot_similarity(diff_df)
    
#    from matplotlib.pyplot import cm
#    import matplotlib.colors as clr
#    import matplotlib.pyplot as plt
#    from matplotlib.colors import Normalize
#    from sklearn.decomposition import PCA
#    from sklearn.manifold import TSNE
#    pca = PCA(n_components=2)
##    pca = TSNE(n_components=2)
#    
#    reduced = pca.fit_transform(diff_df)
#    prop_values = constructor.target_df.values.ravel()
#    
#    # Set colormap to use
#    cmap = cm.hot
#    # Get min and max values
#    vmin = np.min(prop_values)
#    vmax = np.max(prop_values)
#    # Get normalized values for RGB conversion
#    norm = Normalize(vmin, vmax)
#    # Convert to RGB and use cmap in plotting
#    colors = norm(prop_values)
#    # Set scalar mappable to use when adding colorbar
#    SM = cm.ScalarMappable(norm, cmap)
#    # Set array for scalar mappable for colorbar on plot
#    SM.set_array(prop_values)
#    
#    
#    fig = plt.figure()
#    ax = fig.add_subplot(111)
#    ax.scatter(reduced[:,0],reduced[:,1], 
#                cmap=cmap, c=prop_values)
#    ax.set_title("Difference Matrix")
#    cbar_obj = plt.colorbar(SM)
#    cbar_obj.ax.set_ylabel('Molecule Volume, $\AA^3$', 
#                           rotation=270,
#                           labelpad=25)
    
    
    ################# Computing the complexity graph
#    from sklearn.linear_model import Ridge
#    regressor = Ridge(alpha=100)
#    
#    bn = BondNeighborhood()
#    constructor = construct_bond_neighborhood_model(bn,add_prop=["MC_volume"], 
#                                              target_prop=["molecule_volume"])
#    constructor.set_train(s_conquest)
#    constructor.regularize(10)
#    constructor.set_test(s_19,r=True)
#    
#    greedy_features,greedy_model_err = constructor.model_complexity_analysis(
#                                        regressor,
#                                        choose_n_features=500,
#                                        neighborhood_df_train=constructor.neighborhood_df_r,
#                                        target_train=constructor.target_df,
#                                        neighborhood_df_test=constructor.neighborhood_df_t,
#                                        target_test=constructor.target_df_t)
#    
#    fig = plt.figure()
#    ax = fig.add_subplot(111)
#    ax.plot(np.arange(0,greedy_model_err.shape[0]),greedy_model_err[:,0])
#    ax.plot(np.arange(0,greedy_model_err.shape[0]),greedy_model_err[:,1])
#    ax.set_xlim([00,450])
    
#    regressor.fit(greedy_features, constructor.target_df)
#    sort_idx = np.argsort(np.abs(regressor.coef_[0]))[::-1]
#    features_sorted = greedy_features.columns[sort_idx]
#    coef_sorted = regressor.coef_[0,sort_idx]         
    
#    regressor = Ridge(alpha=100)
#    temp_features = pd.DataFrame(data=greedy_features.values[:,0],
#                                 index=greedy_features.index,
#                                 columns=[greedy_features.columns[0]])
#    
#    test_df = constructor.neighborhood_df_t
#    target_t = constructor.target_df_t.values
#    
#    regressor.fit(temp_features,target)
#    pred = regressor.predict(temp_features)
#    pred_t = regressor.predict(test_df.loc[:,temp_features.columns])
#    err_temp = [np.mean(np.abs(target - pred) / pred)]
#    err_temp_test = [np.mean(np.abs(target_t - pred_t) / pred_t)]
#    total = len(greedy_features.columns)
##    
#    for column in greedy_features.columns[1:]:
#        temp_features[column] = greedy_features[column]
#        regressor.fit(temp_features,target)
#        pred = regressor.predict(temp_features)
#        pred_t = regressor.predict(test_df.loc[:,temp_features.columns])
#        err_temp.append(np.mean(np.abs(target - pred) / pred))
#        err_temp_test.append(np.mean(np.abs(target_t - pred_t) / pred_t))
#        total -= 1
#        print(total, regressor.coef_[0,0])
#
    
    
    
    ######### Making some nice plots
    from ibslib.analysis.pltlib import format_ticks
    fig = plt.figure()
    ax = fig.add_subplot(111)
    x=np.arange(0,len(greedy_features.columns),1)
    ax.plot(x,np.array(err_temp)*100, linewidth=3,
            label="Train")
    ax.plot(x,np.array(err_temp_test)*100, linewidth=3,
            label="Test")          
    format_ticks(ax,tick_size=14)
    ax.set_xlabel("Fragments in Model", fontsize=18)
    ax.set_ylabel("Percent MAE",fontsize=18)
    ax.axvline(64,c="tab:purple",linestyle="--",alpha=0.75,
               label="Model")
    ax.axvline(109,c="tab:red",linestyle="--",alpha=0.75,
               label="Overfit")
    plt.legend(bbox_to_anchor=(0.1,1))
    plt.tight_layout()
    fig.savefig("model_complexity.pdf")
    
    fig = plt.figure(figsize=(3,6))
    ax = fig.add_subplot(111)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
    format_ticks(ax)
    
    
    num_top_features = 10
    start=0.85
    stop=0
    spacing=(start-stop) / num_top_features
    xcol_1 = 0.05
    xcol_2 = 0.55
    
    fontsize=13
    ax.text(xcol_1,0.925,"Feature",fontsize=fontsize+2)
    ax.text(0.495,0.925,"Coefficient",fontsize=fontsize+2)
    for i,value in enumerate(greedy_features.columns[0:10]):
        xloc = xcol_1
        yloc = start - i*spacing
        print(yloc)
        if value == "MC_volume": value = "MC"
        ax.text(xloc,yloc,value,fontsize=fontsize)
        
        xloc = xcol_2
        ax.text(xloc,yloc,"{:.2f}".format(regressor.coef_[0,i]),
                fontsize=fontsize)
    ax.axvline(0.45, color='k', linestyle="--")
    ax.axhline(0.91, color='k', linestyle="--")
    plt.show()
        
    ######### Calc structure fragments and counts
#    struct = s_paper_r["GLYCIN_from_crystal"]
#    f,c = bn.calc(struct)
#    # Check which fragments in struct are in fragment list and return
#    # the index of their location
#    final_features = greedy_features.columns[0:]
#    f_idx,c_idx = np.nonzero(final_features.values[:,None] == f)
#    # Populate feature vector 
#    neighborhood_vector = np.zeros(final_features.shape[0])
#    neighborhood_vector[f_idx] = c[c_idx]
#    neighborhood_vector[0] = struct.properties["MC_volume"]
#
#        
#    # Calc structure fragments and counts
#    struct = s_paper_r["benzene_relaxed_geometry"]
#    f,c = bn.calc(struct)
#    # Check which fragments in struct are in fragment list and return
#    # the index of their location
#    final_features = greedy_features.columns[0:]
#    f_idx,c_idx = np.nonzero(final_features.values[:,None] == f)
#    # Populate feature vector 
#    neighborhood_vector = np.zeros(final_features.shape[0])
#    neighborhood_vector[f_idx] = c[c_idx]
#    neighborhood_vector[0] = struct.properties["MC_volume"]
#
#    
#    n=100
#    err_list = np.zeros(n)
#    g_err_list = np.zeros(n)
#    target = constructor.target_df.values
#    
#    for i in range(0,n):
#        regressor = Ridge(alpha=i,fit_intercept=False)
##        regressor.fit(test_df.loc[:,final_features],target_t)
#        regressor.fit(greedy_features.loc[:,final_features],target)
#        sort_idx = np.argsort(np.abs(regressor.coef_[0]))[::-1]
#        features_sorted = final_features[sort_idx]
#        coef_sorted = regressor.coef_[0,sort_idx]   
#        pred = regressor.predict(greedy_features.loc[:,final_features])
#        err_list[i] = np.std(np.abs(target - pred)/pred)
#        
#        glycine_pred = regressor.predict(neighborhood_vector[None,:64])
#        diff = np.abs(glycine_pred - 78) / 78
#        g_err_list[i] = diff
    
    ################# Computing greedy model with CV
#    from sklearn.linear_model import Ridge
#    from sklearn.model_selection import KFold
#    regressor = Ridge(alpha=100)
#    kf = KFold(n_splits=1,shuffle=False,random_state=1)
#    
#    train_dict = s_conquest
#    test_dict = s_19
#    
#    bn = BondNeighborhood()
#    constructor = construct_bond_neighborhood_model(bn,add_prop=["MC_volume"], 
#                                              target_prop=["molecule_volume"])
#    constructor.set_train(train_dict)
#    
#    constructor.regularize(tol=1)
#    temp_features = constructor.neighborhood_df_r
#    greedy_features_list = []
#    greedy_model_err_list = []
#    for train_index,val_index in kf.split(temp_features):
#        x_train = temp_features.iloc[train_index,:]
#        train_target = target[train_index,:]
#        
#        temp_dict = {}
#        for struct_id in x_train.index:
#            temp_dict[struct_id] = s_conquest[struct_id]
#        
#        constructor.set_train(temp_dict)
#        constructor.regularize(tol=1)
#        constructor.set_test(test_dict,r=True)
#        greedy_features,greedy_model_err = constructor.model_complexity_analysis(
#                                    regressor,
#                                    choose_n_features=500,
#                                    neighborhood_df_train=constructor.neighborhood_df_r,
#                                    neighborhood_df_test=constructor.neighborhood_df_t)
#        greedy_features_list.append(greedy_features)
#        greedy_model_err_list.append(greedy_model_err)
##    
#    sizes = [x.shape[0] for x in greedy_model_err_list]
#    same_size_train = [x[0:np.min(sizes),0][:,None] for x in greedy_model_err_list]
#    same_size_test = [x[0:np.min(sizes),1][:,None] for x in greedy_model_err_list]
#    greedy_model_err_array_train = np.hstack(same_size_train)
#    greedy_model_err_array_test = np.hstack(same_size_test)
#    
#    mean_train = np.zeros(np.min(sizes))
#    mean_test = np.zeros(np.min(sizes))
#    std_train = np.zeros(np.min(sizes))
#    std_test = np.zeros(np.min(sizes))
#    for i,err in enumerate(greedy_model_err_array_train[0:np.min(sizes)]):
#        mean_train[i] = np.mean(err)
#        std_train[i] = np.std(err)
#    for i,err in enumerate(greedy_model_err_array_test[0:np.min(sizes)]):
#        mean_test[i] = np.mean(err)
#        std_test[i] = np.std(err)
#        
#    fig = plt.figure()
#    ax = fig.add_subplot(111)
#    x = np.arange(0,np.min(sizes))
#    ax.plot(x,mean_train)
#    ax.fill_between(x,mean_train-std_train,mean_train+std_train,
#                    alpha=0.5)
#    ax.plot(x,mean_test)
#    ax.fill_between(x,mean_test-std_test,mean_test+std_test,
#                    alpha=0.5)
##    ax.set_xlim([0,200])
#    
#    fig = plt.figure()
#    ax = fig.add_subplot(111)
#    x = np.arange(0,np.min(sizes))
#    ax.plot(x,mean_train)
#    ax.fill_between(x,mean_train-std_train,mean_train+std_train,
#                    alpha=0.5)
#    ax.plot(x,mean_test)
#    ax.fill_between(x,mean_test-std_test,mean_test+std_test,
#                    alpha=0.5)
#    ax.set_xlim([-2.5,50])
#    
#    first_100_columns = [x.columns[0:50] for x in greedy_features_list]
#    shared_1_2 = np.intersect1d(first_100_columns[0], first_100_columns[1])
#    shared_2_3 = np.intersect1d(first_100_columns[1], first_100_columns[2])
#    shared_all = np.intersect1d(shared_1_2, shared_2_3)
#    
#    final_features = constructor.neighborhood_df.loc[:,shared_all]
#    
#    regressor.fit(final_features.values, constructor.target_df.values)
#    pred = regressor.predict(final_features)
#    final_err = np.abs(pred - constructor.target_df.values) / pred
#    
#    final_test = constructor.neighborhood_df_t.loc[:,shared_all]
#    pred_test = regressor.predict(final_test)
#    final_err_test = np.abs(pred_test - 
#                            constructor.target_df_t.values) / pred_test
#    
#    sort_idx = np.argsort(np.abs(regressor.coef_[0]))[::-1]
#    features_sorted = final_test.columns[sort_idx]
#    coef_sorted = regressor.coef_[0,sort_idx]                        
    
    
    
    ########### TEMP
    
    ### Will the greedy features constructed, I can make any graph I want
    ### with the regressor by iteratively adding columns to the model.
    ### However, I think the reason I'm not seeing overfitting is due to model
    ### Regularization. If the model is well regularized on the training set
    ### then it may not be overfit w.r.t. the testing set. Try again but 
    ### without the regularized model.
    
#    from sklearn.linear_model import LinearRegression,Lasso,Ridge
#    from sklearn.model_selection import KFold
#    regressor = Ridge(alpha=100)
#    
#    temp_features = pd.DataFrame(data=greedy_features.values[:,0],
#                                 index=greedy_features.index,
#                                 columns=[greedy_features.columns[0]])
#    
#    target = constructor.target_df.values
#    test_df = constructor.neighborhood_df_t
#    target_t = constructor.target_df_t.values
#    
#    regressor.fit(temp_features,target)
#    pred = regressor.predict(temp_features)
#    pred_t = regressor.predict(test_df.loc[:,temp_features.columns])
#    err_temp = [np.mean(np.abs(target - pred) / pred)]
#    err_temp_test = [np.mean(np.abs(target_t - pred_t) / pred_t)]
#    total = len(greedy_features.columns)
#    
#    n_splits=3
#    kf = KFold(n_splits=n_splits,shuffle=False,random_state=1)
#    err_train_std = []
#    err_train_mean = []
#    err_test_std = []
#    err_test_mean = []
#    err_val_std = []
#    err_val_mean = []
#    
#    
#    for column in greedy_features.columns:
#        temp_features[column] = greedy_features[column]
#        
#        #### Performing splitting of training set to construct std
#        err_train_for_std = []
#        err_test_for_std = []
#        err_val_for_std = []
#        for train_index,val_index in kf.split(temp_features):
#            x_train = temp_features.iloc[train_index,:]
#            train_target = target[train_index,:]
#            regressor.fit(x_train,train_target)
#            pred = regressor.predict(x_train)
#            train_err = np.mean(np.abs(train_target - pred) / pred)
#            err_train_for_std.append(train_err)
#            
#            pred_t = regressor.predict(test_df.loc[:,temp_features.columns])
#            test_err = np.mean(np.abs(target_t - pred_t) / pred_t)
#            err_test_for_std.append(test_err)
#            
#            x_val = temp_features.iloc[val_index,:]
#            val_target = target[train_index,:]
#            pred_v = regressor.predict(x_val)
#            val_err = np.mean(np.abs(train_target - pred) / pred)
#            err_val_for_std.append(val_err)
#            
#            
#        err_train_std.append(np.std(err_train_for_std))
#        err_train_mean.append(np.mean(err_train_for_std))
#        err_test_std.append(np.std(err_test_for_std))
#        err_test_mean.append(np.mean(err_test_for_std))
#        err_val_std.append(np.std(err_val_for_std))
#        err_val_mean.append(np.mean(err_val_for_std))
#        
#        
#        regressor.fit(temp_features,target)
#        pred = regressor.predict(temp_features)
#        pred_t = regressor.predict(test_df.loc[:,temp_features.columns])
#        err_temp.append(np.mean(np.abs(target - pred) / pred))
#        err_temp_test.append(np.mean(np.abs(target_t - pred_t) / pred_t))
#        total -= 1
#        print(total, regressor.coef_[0,0])
    
#    fig = plt.figure()
#    ax = fig.add_subplot(111)
#    ax.plot(np.arange(0,len(err_temp_test),1),err_temp)
#    ax.plot(np.arange(0,len(err_temp_test),1),err_temp_test)
#    
#    fig = plt.figure()
#    ax = fig.add_subplot(111)
#    err_train_mean = np.array(err_train_mean)
#    err_train_std = np.array(err_train_std)
#    err_test_mean = np.array(err_test_mean)
#    err_test_std = np.array(err_test_std)
#    x=np.arange(0,len(err_temp_test)-1,1)
#    ax.plot(x,err_train_mean)
#    ax.fill_between(x,err_train_mean-err_train_std,err_train_mean+err_train_std,
#                    alpha=0.5)
#    ax.plot(x,err_test_mean)
#    ax.fill_between(x,err_test_mean-err_test_std,err_test_mean+err_test_std,
#                    alpha=0.5)
#    ax.plot(x,err_val_mean)
#    ax.fill_between(x,err_val_mean-err_val_std,err_val_mean+err_val_std,
#                    alpha=0.5)
    