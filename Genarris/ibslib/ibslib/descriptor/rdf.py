

from ibslib.descriptor.R_descriptor import calc_R

import numpy as np
import torch
from torch import nn


class ACSF(nn.Module):
    """
    Initialized entire RDF layer.
    """
    def __init__(self, device, n_D_inter=4, init_scheme="centered", cutoff=12,
                   eta_range=[0.05,0.5], Rs_range=[0.1,12], learn_rep=False):
        """ 
        Representation layer: Current implementation is only for a 
            radial symmetry funciton. 

        Arguments
        ---------
        device: torch.device
        n_D_inter: int
            Number of symmetry functions to use for each interaction.
        init_scheme: str 
            One of "centered","shifted","random". Defines the initialization
              scheme which is used for the eta and Rs parameters. 
        cutoff: float
            Cutoff distance used in cutoff function.
        eta_range: list [start, end]
            Used in random initialization
        Rs_range: list [start, end]
            Used in random initialization
        learn_rep: bool
            Determines if gradients will be passed to the eta and Rs parameters

        """
        super(ACSF, self).__init__()
        self.device = device
        self.requires_grad = learn_rep
        self.cutoff = cutoff
        self.n_D_inter = n_D_inter
        self.eta_range = eta_range
        self.Rs_range = Rs_range

        self.avail_init_schemes = ["centered","shifted","random"]
        if init_scheme not in self.avail_init_schemes:
            raise Exception("Initializaiton scheme for representation layer "+
            "is not available. Available schemes are {}, however input was {}."
            .format(self.avail_init_schemes, init_scheme))
        
        # Performing parameter initialization
        eval("self.init_{}()".format(init_scheme))

        self.eta = nn.Parameter(self.eta)
        self.Rs = nn.Parameter(self.Rs)
    

    def init_centered(self):
        """
        Initialization scheme used for centered radial symmetry funciton as 
          defined in Gastegger et al. 2018. 
        """
        if self.n_D_inter > 1:
            del_R = (self.cutoff - 1.5) / (self.n_D_inter - 1)
            r = torch.range(start=1.0, end=(self.cutoff-0.5), step=del_R)
        else:
            del_R = self.cutoff/2
            r = torch.tensor([self.cutoff / 2])

        eta = 1 / (2*r.pow(2))
        self.eta = torch.tensor(eta, requires_grad=self.requires_grad,
                                device=self.device)
        
        Rs = torch.zeros(eta.size())
        self.Rs = torch.tensor(Rs, requires_grad=self.requires_grad,
                                device=self.device)

    
    def init_shifted(self):
        """
        Initialization scheme used for shifted radial symmetry funciton as 
          defined in Gastegger et al. 2018. 
        """
        if self.n_D_inter > 1:
            del_R = (self.cutoff - 1.5) / (self.n_D_inter - 1)
            r = torch.range(start=1.0, end=(self.cutoff-0.5), step=del_R)
        else:
            del_R = self.cutoff/2
            r = torch.tensor([self.cutoff / 2])

        self.Rs = torch.tensor(r, requires_grad=self.requires_grad,
                                device=self.device)
        eta = torch.zeros(self.Rs.size()) + 1 / (2*del_R**2)
        self.eta = torch.tensor(eta, requires_grad=self.requires_grad,
                                device=self.device)

    
    def init_random(self):
        """
        Random initialization of eta and Rs
        """
        self.eta = torch.FloatTensor(self.n_D_inter,).uniform_(
                                     self.eta_range[0],self.eta_range[1])
        self.eta = torch.tensor(self.eta, requires_grad=self.requires_grad, 
                                device=self.device)
        self.Rs = torch.FloatTensor(self.n_D_inter,).uniform_(
                                    self.Rs_range[0],self.Rs_range[1])
        self.Rs = torch.tensor(self.Rs, requires_grad=self.requires_grad, 
                               device=self.device)


    def cutoff_fn(self,R):
        """
        Computes the cutoff function for a batch of examples.

        cutoff = {
                    (1/2)[cos(pi*rij / cutoff) + 1] if rij <= rc
                    0
                 }

        Returns result from cutoff fn. This will be a torch tensor 
          with dimension [n_samples, n_atoms_per_system, R]
        """
        # Find Rij values less than the cutoff distance
        mask = R < self.cutoff
        result = (1/2)*(torch.cos(np.pi * R / self.cutoff) + 1)
        result = result * mask.float()
        return result
    
    
    def eval_struct(self, R):
        # Remove R[:,None] for single structure
        G = torch.exp(-self.eta[:,None,None]*(R - self.Rs[:,None,None]).pow(2))
        # Calculate cutoff function
        cutoff_result = self.cutoff_fn(R)
        # Broadcast cutoff function across all n_D_inter
        G = G * cutoff_result[None,:,:]
        # Now Sum
        G = torch.sum(G, dim=-1)
        # Sum here only performs sum across each atom, however, for a single 
        # evaluation of a structure, we are not conserned with the embedding of 
        # each individual atom. So, sum over th embedding of all atoms to give
        # index invariance.
        # However, make sure this works even if there's only a single atom of 
        # a specific type in the structure. 
        G = torch.sum(G, dim=-1)
        
        
        
        return G
    
    
    def forward(self, R):
        """
        Arguments:
        R: torch.tensor
            Tensor with dimension [n_samples, n_atoms_per_system, R]
        
        Returns final representation vector for each atom in the systems
            with dimensions [n_samples, n_D_inter, n_atom]
        """

        # Verified for minibatch
        G = torch.exp(-self.eta[None,:,None,None]*(R[:,None] - self.Rs[:,None,None]).pow(2))
        self.eta.retain_grad()
        self.Rs.retain_grad() 

        # Calculate cutoff function
        cutoff_result = self.cutoff_fn(R)
        # Broadcast cutoff function across all n_D_inter
        G = G * cutoff_result[:,None,:,:]
        # Now Sum
        G = torch.sum(G, dim=-1)
        
        # Batch Norm
        mean = torch.mean(G, dim=1, keepdim=True)
        std = torch.std(G, dim=1, keepdim=True)

        # Atom batch norm
        mean = torch.mean(mean,dim=-1, keepdim=True)
        std = torch.mean(std, dim=-1, keepdim=True)

        G = (G - mean) / std

        return G
    
#    
#    def __init__(self, nodes=1, init_type='shifted',
#                 eta=torch.tensor([]), Rs=torch.tensor([])):
#        """
#        Function initialized with eta and Rs values
#        """
#        if eta.nelement() == 0:
#            self.eta = 0
#        self.eta = eta
#        self.Rs = Rs
#        
#    
#    def _initialize(self):
#        pass
#    
#    
#    def eval_struct(self, struct):
#        """
#        For the purpose of evaluation descriptor of single structure
#        """
#        pass
#    
#    
#    def forward(self, R):
#        """
#        Passes through entire batch
#        """
#        pass
        
if __name__ == "__main__":
    pass
    