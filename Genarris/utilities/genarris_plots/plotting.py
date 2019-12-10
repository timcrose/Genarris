from ibslib.analysis.plotting import plot_IP_hist
from ibslib.io import read,write
from ibslib.analysis import get
from ibslib.analysis.plotting import plot_IP_hist, plot_IP_hist_dict, plot_IP_hist_dir
import numpy as np
import sys, os
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, IndexLocator, FixedLocator, FixedFormatter


from ibslib.io import import_structures,output_struct_dict

try: plot_GA_prop
except: from ibslib.analysis.plotting import plot_GA_prop,motif_prop_plot

#try: plot_energy_property_motif
#except: from ibslib.analysis.plotting import plot_energy_property_motif

from ibslib.analysis.plotting import correct_energy,prepare_prop_from_dict
from ibslib.descriptor.clustering import run_KMeans

from sklearn.cluster import KMeans

from ase.data import atomic_masses, atomic_numbers

import matplotlib
import matplotlib.colors as mcolors
import matplotlib.colors as clr
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import axes3d, Axes3D
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
import matplotlib.patheffects as PathEffects

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=10000):
    new_cmap = clr.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def lattice_parameter_plot(structure_dir, expt_fpath=None, energy_name='energy',
              xlim=-1, ylim=-1, zlim=-1,
              xlabel='a ($\AA$)', ylabel='b ($\AA$)', zlabel='c ($\AA$)',
              cbar_title='Relative Energy (kJ/mol/molecule)',
              fig_size=(16,10), global_min=None, global_max=None,
              experimental_energy_name='energy',s=200, alpha=0.6,
              ax_x_labelpad=12, ax_y_labelpad=12, ax_z_labelpad=12, 
              ax_x_tickpad=0, ax_y_tickpad=0, ax_z_tickpad=0,
              cbar_labelpad=30,cmap_truncation=[0.1,0.8],
              cbar_rotation=270,
              ax_fontdict={ 'family' : 'DejaVu Sans',
                            'weight' : 'normal',
                            'size'   : 22},
                expt_fontdict={ 'family': 'calibri',
                                'color':  'green',
                                'weight': 100,
                                'size': 28},
                figname='lattice_parameter_plot.png', dpi=1200, show_x=False,
                plt_show=False, colormap=cm.hot, show_cbar=False, cbar_orientation='horizontal'):

    plt.close()
    plt.cla()
    plt.clf()
    '''
    SMALL_SIZE = 8
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    '''
    plt.rc('font', **ax_fontdict)
    struct_dict = read(structure_dir)
    
    aval, bval, cval, energyval = ([] for i in range(4))
    directions = ['a', 'b', 'c']

    for struct_id,struct in struct_dict.items():
        Z = struct.get_property('Z')
        aval.append(struct.get_property('a'))
        bval.append(struct.get_property('b'))
        cval.append(struct.get_property('c'))
        if energy_name in struct.properties and struct.properties[energy_name] != 'none':
            energyval.append(struct.get_property(energy_name))
    
    if expt_fpath is not None:
        expt_struct = read(expt_fpath)
        #struct_dict['expt_struct'] = expt_struct
        if experimental_energy_name in expt_struct.properties:
            experimental_eV = expt_struct.get_property(experimental_energy_name)
            #energyval.append(experimental_eV)
        experimental_cell = [expt_struct.get_property(d) for d in directions]
        #aval.append(experimental_cell[0])
        #bval.append(experimental_cell[1])
        #cval.append(experimental_cell[2])

    if len(energyval) <= 1:
        energyval = None
        colormap = None
    else:
        energyval = correct_energy(energyval, nmpc=Z, global_min=global_min)
        colormap = truncate_colormap(colormap, cmap_truncation[0], cmap_truncation[1])
        aval, bval, cval, energyval = ([] for i in range(4))
        for struct_id,struct in struct_dict.items():
            if energy_name in struct.properties and struct.properties[energy_name] != 'none':
                aval.append(struct.get_property('a'))
                bval.append(struct.get_property('b'))
                cval.append(struct.get_property('c'))
                energyval.append(struct.get_property(energy_name))

    if (global_max is not None) and (energyval is not None):
        aval.append(-1)
        bval.append(-1)
        cval.append(-1)
        energyval.append(global_max)

    if energyval is not None:
        energyval = correct_energy(energyval, nmpc=Z, global_min=global_min)
    
    fig = plt.figure(figsize=fig_size)
    ax = fig.add_subplot(111, projection='3d')
    
    if expt_fpath is not None and show_x:
        txt = ax.text(experimental_cell[0], 
                experimental_cell[1], 
                experimental_cell[2],  
                'X', fontdict=expt_fontdict,
                horizontalalignment='center',
                verticalalignment='center')
            
        txt.set_path_effects([PathEffects.withStroke(linewidth=3, foreground='k')])
    
    p = ax.scatter(aval, bval, cval, c=energyval, 
                    cmap=colormap, s=s, alpha=alpha) 
    
    if colormap is not None and show_cbar:
        cbar = plt.colorbar(p, orientation=cbar_orientation)
        cbar.ax.set_ylabel(cbar_title, rotation=cbar_rotation, labelpad=cbar_labelpad)
    if xlim != -1:
        ax.set_xlim(xlim[0], xlim[1])
    if ylim != -1:
        ax.set_ylim(ylim[0], ylim[1])
    if zlim != -1:
        ax.set_zlim(zlim[0], zlim[1])

    ax.set_xlabel(xlabel, labelpad=ax_x_labelpad)
    ax.set_ylabel(ylabel, labelpad=ax_y_labelpad)
    ax.set_zlabel(zlabel, labelpad=ax_z_labelpad)

    ax.tick_params(axis='x', pad=ax_x_tickpad)
    ax.tick_params(axis='y', pad=ax_y_tickpad)
    ax.tick_params(axis='z', pad=ax_z_tickpad)
    plt.tight_layout()
    fig.savefig(figname, dpi=dpi)
    if plt_show:
        plt.show()
    
'''
def cluster_plot(struct_dict, value='data'):
    aval, bval, cval, energyval = ([] for i in range(4))

    for struct_id,struct in struct_dict.items():
        aval.append(struct.get_property('a'))
        bval.append(struct.get_property('b'))
        cval.append(struct.get_property('c'))
    
    feature_matrix = np.array([aval,bval,cval])
    feature_matrix = feature_matrix.T
    
    
    #k_model = run_KMeans(feature_matrix, **{'n_clusters':12})
    k_model = KMeans(n_clusters=24,random_state=2)
    k_model.fit(feature_matrix)
    
    fig = plt.figure(figsize=(6,5))
    ax = fig.add_subplot(111, projection='3d')
    
    colormap = cm.viridis
    colormap = truncate_colormap(colormap, 0.1, 0.9)
    
    if value == 'data':
        plt1 = ax.scatter(feature_matrix[:,0],
                          feature_matrix[:,1],
                          feature_matrix[:,2],
                          c=k_model.labels_,
                          s=100)
        ax.set_zlim([10,27.5])
        ax.set_ylim([5,15.5])
        ax.set_yticks(np.arange(5,15.5,2))
    elif value =='center':
        plt1 = ax.scatter(k_model.cluster_centers_[:,0],
                          k_model.cluster_centers_[:,1],
                          k_model.cluster_centers_[:,2],
                          c=range(0,k_model.cluster_centers_.shape[0]),
                          s=100)
        ax.set_zlim([10,27.5])
        ax.set_ylim([5,15.5])
        ax.set_yticks(np.arange(5,15.5,2))
        
    
    ax.tick_params(axis='both', width=10, length=7, labelsize=16, pad=0)
    ax.zaxis.set_major_formatter(FormatStrFormatter('%i'))
    ax.tick_params(axis='z', pad=1, )
    ax.set_xlabel(' Length a, $\AA$', labelpad=5, fontsize=18)
    ax.set_ylabel('Length b, $\AA$', labelpad=5, fontsize=18)
    ax.set_zlabel('Length c, $\AA$', labelpad=5, fontsize=18)
    
    plt.tight_layout()
    plt.show()
    
    fig.savefig('Clustering_Plot_Center_20190516.pdf')
    
    return k_model
    

    #example usage
    #struct_dir = "../density_scripts/FUQJIK01/ensemble"
    #output_dir = "FUQJIK01/Energy_Cross075"
    #struct_dict = import_structures(struct_dir)

    k_model = cluster_plot(struct_dict, value='center')
'''


def plot_property(struct_dir, prop='unit_cell_volume', nmpc=2, figname='property_histogram.pdf', 
                    xlabel='Structure Volume, $\AA^3$', ylabel='Number of structures', figure_size=(12,8),
                        label_size=24, tick_size=18, tick_width=3, GAtor_IP=False, expt_fpath=None,
                        bins='sqrt', predicted_volume=None, hist_range=None, y_axis_side='left', dpi=1200,
                        xtick_label_frequency='all', xtick_labelsize=26, xtick_labelrotation=0):
    # Can change any of the following key word arguments to adjust plotting
    kwargs = {
          # Property to extract from structures
          'prop': prop,
          # nmpc is necessary for some property calculations such as
          #   for space groups
          'nmpc': nmpc,
          # Change figname to save the figure that's produced.
          #   It's recommended to save figures in a .pdf format
          'figname': figname,
          # Other graph plotting parameters
          'xlabel': xlabel,
          'ylabel': ylabel,
          'figure_size': figure_size,
          # fed into matplotlib plt.hist (so could be int)
          'hist_kwargs':{"density": False,
                        "bins": bins,
                        "edgecolor": "k",
                        'range': hist_range,
                        },
          #'density': False,
          'label_size': label_size,
          'tick_size': tick_size,
          'tick_width': tick_width,
          # GAtor_IP controls if the histogram for the initial pool should be plotted 
          #   alongside the entire GA pool of structures for GA plotting. 
          'GAtor_IP': GAtor_IP
          }
    # Calling the function with just the directory of structure files
    plt.close()
    ax1 = plot_IP_hist_dir(struct_dir, **kwargs)
    if expt_fpath is not None:
        expt_struct = read(expt_fpath)
        expt_vol = expt_struct.properties['volume_before_relaxation']
        ax1.axvline(x=expt_vol, color='lightgreen', linestyle='-', linewidth=3)
        if predicted_volume is not None:
            ax1.axvline(x=predicted_volume, color='orange', linestyle='--', linewidth=3)
    if y_axis_side == 'right':
        ax1.yaxis.set_label_position("right")
        ax1.yaxis.tick_right()

    ax1.tick_params(axis='x', labelsize=xtick_labelsize, labelrotation=xtick_labelrotation)

    if xtick_label_frequency != 'all':
        ax1.xaxis.set_minor_locator(AutoMinorLocator())
        ax1.tick_params(which='minor', length=7)
        ax1.tick_params(which='major', length=10)
        # Keeps every xtick_label_frequency-th label
        [l.set_visible(False) for (i,l) in enumerate(ax1.xaxis.get_ticklabels()) if i % xtick_label_frequency != 0]

    #Round to nearest 10 upwards
    y_lim = plt.gca().get_ylim()
    print('y_lim', y_lim)
    max_count = int(np.ceil(y_lim[1]))
    max_tick = 10 - (max_count % 10) + max_count
    ytick_interval = max_tick / 10
    ax1.set_yticks(np.arange(0, max_tick + ytick_interval, ytick_interval))
    plt.tight_layout()
    plt.savefig(figname, dpi=dpi)


def parse_spg_outfile(pygenarris_outfile='outfile'):
    '''
    outfile: str
        file that pygenarris writes to when generating structures (not geometry output file)
    '''
    with open(pygenarris_outfile) as f:
        lines = f.readlines()
    all_allowed_spgs = []
    for line in lines:
        split_line = line.split()
        if len(split_line) == 0:
            continue
        if 'spg' in split_line[0] and len(split_line) == 11:
            spg = split_line[1]
        else:
            continue
        if spg in all_allowed_spgs:
            continue
        all_allowed_spgs.append(spg)
        
    return all_allowed_spgs


def plot_spg_bar_chart(struct_dir, pygenarris_outfile='outfile', spg_bar_chart_fname='spg_bar_chart.pdf',
                        width=0.5, ylabel='Number of structures', xlabel='Allowed space groups',
                        title='spg_histogram', dpi=1200, legend_fontsize=14,
                        show_legend=True, fontdict={'family' : 'DejaVu Sans',
                                                    'weight' : 'normal',
                                                    'size'   : 14}, expt_fpath=None, fig_size=None, y_axis_side='left',
                        xtick_label_frequency='all', xtick_labelsize=26, xtick_labelrotation=0,
                        ax_x_labelpad=1, ax_y_labelpad=2, spg_arrow_width=0.09, ytick_labelsize=24, spg_list=[],
                        spg_show_arrow=True, one_color_for_spg_plot=False):
    all_allowed_spgs = parse_spg_outfile(pygenarris_outfile=pygenarris_outfile)
    # Directory of json files is struct_dir
    # Run using the struct_dict or by calling the function with a directory
    struct_dict = read(struct_dir)
    results_df = get(struct_dict, 'prop', ['spg'])
    spgs_gotten = list(results_df['spg'].values)
    unique_spgs_allowed = list(map(int, set(all_allowed_spgs)))
    unique_spgs_allowed += spg_list
    unique_spgs_allowed = list(set(unique_spgs_allowed))
    unique_spgs_allowed.sort()

    N = len(unique_spgs_allowed)
    #spgs_gotten_counts = [spgs_gotten.count(spg) for spg in unique_spgs_allowed]
    spgs_general_pos_counts = [0] * N
    spgs_special_pos_counts = [0] * N
    spgs_all_pos_counts = [0] * N
    for i,spg in enumerate(unique_spgs_allowed):
        for struct_id in struct_dict:
            if int(struct_dict[struct_id].properties['spg']) == spg:
                site_symm = struct_dict[struct_id].properties['site_symmetry_group']
                try:
                    site_symm = int(site_symm)
                except:
                    pass
                # Right now the site_symm might be incorrect if Z' > 1 so manually make spg 1 be classified as general position and also glycine
                if site_symm == 1 or struct_dict[struct_id].properties['spg'] == 1 or 'glycine' in struct_dir or 'glycine' in struct_id or (expt_fpath is not None and 'glycine' in expt_fpath):
                    # general position
                    spgs_general_pos_counts[i] = spgs_general_pos_counts[i] + 1
                else:
                    # special position
                    spgs_special_pos_counts[i] = spgs_special_pos_counts[i] + 1
                spgs_all_pos_counts[i] = spgs_all_pos_counts[i] + 1

    ind = np.arange(N)    # the x locations for the groups
    width = width       # the width of the bars: can also be len(x) sequence
    plt.close()
    plt.cla()
    plt.clf()
    plt.close()

    plt.rc('font', **fontdict)
    plt.rcParams['figure.figsize'] = fig_size
    _, ax = plt.subplots()
    if one_color_for_spg_plot:
        print('plotting purple bars', flush=True)
        p0 = plt.bar(ind,spgs_all_pos_counts, width, color='purple')
        if show_legend:
            p1 = plt.bar(ind, np.zeros(len(spgs_general_pos_counts)), width)
            p2 = plt.bar(ind, np.zeros(len(spgs_special_pos_counts)), width, bottom=spgs_general_pos_counts)
            plt.rc('legend', fontsize=legend_fontsize)    # legend fontsize
            plt.legend((p1[0], p2[0], p0[0]), ('General positions', 'Special positions', 'Either position'))
    else:
        p1 = plt.bar(ind, spgs_general_pos_counts, width)
        p2 = plt.bar(ind, spgs_special_pos_counts, width, bottom=spgs_general_pos_counts)
        if show_legend:
            plt.rc('legend', fontsize=legend_fontsize)    # legend fontsize
            plt.legend((p1[0], p2[0]), ('General positions', 'Special positions'))

    plt.ylabel(ylabel, labelpad=ax_y_labelpad)
    plt.xlabel(xlabel, labelpad=ax_x_labelpad)
    if title is not None:
        plt.title(title)
    ax.set_xticks(ind)
    max_count = max([max(spgs_general_pos_counts), max(spgs_special_pos_counts)])
    #Round to nearest 10 upwards
    max_tick = 10 - (max_count % 10) + max_count
    ytick_interval = max_tick / 10
    ax.set_yticks(np.arange(0, max_tick + ytick_interval, ytick_interval))
    labels = list(map(str, unique_spgs_allowed))
    ax.set_xticklabels(labels)

    y_lim = plt.gca().get_ylim()
    ax.set_ylim([y_lim[0], y_lim[1] * 1.1])
    y_lim = plt.gca().get_ylim()

    if expt_fpath is not None:
        expt_struct = read(expt_fpath)
        expt_spg = expt_struct.properties['spg_before_relaxation']
        
        x = unique_spgs_allowed.index(expt_spg)
        if spg_show_arrow:
            plt.arrow(x, y_lim[1], 0, max_count - y_lim[1], color='lightgreen', width=spg_arrow_width, head_length=y_lim[1] / 15.0, alpha=0.75)


    if y_axis_side == 'right':
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()

    ax.tick_params(axis='x', labelsize=xtick_labelsize, labelrotation=xtick_labelrotation)
    ax.tick_params(axis='y', labelsize=ytick_labelsize)

    if xtick_label_frequency != 'all':
        ax.set_xticklabels([lab for i,lab in enumerate(labels) if i % xtick_label_frequency == 0])
        ax.xaxis.set_major_locator(FixedLocator([i for i in ind if i % xtick_label_frequency == 0]))
        ax.xaxis.set_minor_locator(FixedLocator([i for i in ind if i % xtick_label_frequency != 0]))
        ax.tick_params(which='major', length=7)
        ax.tick_params(which='minor', length=4)

    plt.tight_layout()
    plt.savefig(spg_bar_chart_fname, dpi=dpi)
