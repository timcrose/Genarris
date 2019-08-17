from ibslib.analysis.plotting import plot_IP_hist
from ibslib.io import read,write
from ibslib.analysis import get
from ibslib.analysis.plotting import plot_IP_hist, plot_IP_hist_dict, plot_IP_hist_dir
import numpy as np
import sys, os
import matplotlib.pyplot as plt

def plot_property(struct_dir, prop='unit_cell_volume', nmpc=2, figname='property_histogram.pdf', 
                    xlabel='Structure Volume, $\AA^3$', ylabel='Counts', figure_size=(12,8),
                        label_size=24, tick_size=18, tick_width=3, GAtor_IP=False):
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
          #'bins': 'auto',
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
    plot_IP_hist_dir(struct_dir, **kwargs)


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
        if spg in all_allowed_spgs:
            continue
        all_allowed_spgs.append(spg)
        
    return all_allowed_spgs


def plot_spg_bar_chart(struct_dir, pygenarris_outfile='outfile', spg_bar_chart_fname='spg_bar_chart.pdf',
                        width=0.5, ylabel='Count', xlabel='Allowed space groups',
                        title='spg_histogram', tick_rotation='vertical'):
    all_allowed_spgs = parse_spg_outfile(pygenarris_outfile=pygenarris_outfile)
    # Directory of json files is struct_dir
    # Run using the struct_dict or by calling the function with a directory
    struct_dict = read(struct_dir)
    results_df = get(struct_dict, 'prop', ['spg'])
    spgs_gotten = list(results_df['spg'].values)
    unique_spgs_allowed = list(map(int, set(all_allowed_spgs)))
    unique_spgs_allowed.sort()

    N = len(unique_spgs_allowed)
    spgs_gotten_counts = [spgs_gotten.count(spg) for spg in unique_spgs_allowed]

    ind = np.arange(N)    # the x locations for the groups
    width = width       # the width of the bars: can also be len(x) sequence
    plt.close()
    p1 = plt.bar(ind, spgs_gotten_counts, width)

    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    if title is not None:
        plt.title(title)
    plt.xticks(ind, unique_spgs_allowed, rotation=tick_rotation)

    plt.savefig(spg_bar_chart_fname)