

import pandas as pd
import numpy as np

from ibslib.io import read,write
from ibslib.structures.volume_estimator import MoleculeVolumeEstimator
from ibslib.motif.identify_single_molecule import UniqueMolecules

from scipy.stats import linregress

#try: struct_dict
#except: struct_dict = read("../test_cif")

#struct_dict = read("../test_cif")
#
#columns = ["Experimental", "vdW Spheres", "MC Volume"]
#results = pd.DataFrame(0, index=struct_dict.keys(), columns=columns)
#tracking = []
#i=0
#for struct_id,struct in struct_dict.items():
#    unique = UniqueMolecules(struct)
#    if len(unique.unique_molecule_list) > 1:
#        print("THIS STRUCTURE HAS MOLECULE ISSUES {}".format(struct_id))
#    mve = MoleculeVolumeEstimator(iterations=1e7,batch=1000)
#    mve.MC(unique.unique_molecule_list[0])
#    
#    vdW_spheres = np.sum(np.power(mve.radii,3)*np.pi*(4.0/3.0))
#    
#    experimental = struct.get_unit_cell_volume() / len(unique.molecule_struct_list)
#    print("{} experimental; MC; vdW :  {}; {}; {}"
#          .format(struct_id, experimental, mve.molecule_volume, vdW_spheres))
#    
#    tracking.append(mve.volume_tracking)
#    
#    results["Experimental"][struct_id] = experimental
#    results["vdW Spheres"][struct_id] = vdW_spheres
#    results["MC Volume"][struct_id] = mve.molecule_volume
#    
#    i+=1

MC_results = linregress(results["Experimental"], results["MC Volume"])
sphere_results = linregress(results["Experimental"], results["vdW Spheres"])


from ibslib.analysis.pltlib import format_ticks
import matplotlib.pyplot as plt
from scipy.stats import linregress
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(results["Experimental"], results["MC Volume"], label="MC; r$^2$={:.3f}".format(MC_results[2]))
ax.scatter(results["Experimental"], results["vdW Spheres"], label="Spheres; r$^2$={:.3f}".format(sphere_results[2]))

ax.set_xlabel("Experimental Volume per Molecule", fontsize=14)
ax.set_ylabel("Predicted Volume per Molecule", fontsize=14)

plt.legend()
#ax.set_xlim([0,1500])
#ax.set_xticks(np.arange(0,1500,200))
#ax.set_ylim([0,1500])
#ax.set_yticks(np.arange(0,1500,200))

format_ticks(ax,grid=True)


plt.tight_layout()


line_x = np.arange(0,max(results["Experimental"])+200,100)
line_y = line_x*MC_results[0]+MC_results[1]
line_vdw = line_x*sphere_results[0] + sphere_results[1]
ax.plot(line_x,line_y,c="tab:blue")
ax.plot(line_x,line_vdw,c="tab:orange")

fig.savefig("PAHS_results.pdf")  

#fig = plt.figure()
#ax = fig.add_subplot(111)
#format_ticks(ax,grid=True)
#
#for tracked in tracking:
#    iterations = np.arange(0,1e7,(1e7/len(tracked)))
#    ax.plot(iterations,tracked)
#    
#ax.set_xlabel("MC Iteration", fontsize=14)
#ax.set_ylabel("Volume", fontsize=14)
#
#plt.tight_layout()
#
#fig.savefig("MC_Convergence.pdf")
