from ibslib.io import read,write
from ibslib.analysis.plotting import plot_IP_hist, plot_IP_hist_dict

# Directory of json files
struct_dir = 'delta_glycine_raw_pool'

# Keyword arguments to be passed to matplotlib.pyplot.hist
hist_kwargs = \
    {
      "density": False,
      "bins": "auto",
      "edgecolor": "k",
    }

# Can change any of the following key word arguments to adjust plotting
kwargs = {
          "hist_kwargs": hist_kwargs,
          # Property to extract from structures
          'prop': 'cell_vol',
          # nmpc is necessary for some property calculations such as
          #   for space groups
          'nmpc': 4,
          # Change figname to save the figure that's produced.
          #   It's recommended to save figures in a .pdf format
          'figname': None,
          # Other graph plotting parameters
          'xlabel': 'Structure Volume',
          'ylabel': 'Observed Distribution',
          'figure_size': (12,8),
          'label_size': 24,
          'tick_size': 18,
          'tick_width': 3,
          # GAtor_IP controls if the histogram for the initial pool should be plotted 
          #   alongside the entire GA pool of structures for GA plotting. 
          'GAtor_IP': False
          }

# Run using the struct_dict or by calling the function with a directory
struct_dict = read(struct_dir)
plot_IP_hist_dict(struct_dict, **kwargs)

# Calling the function with just the directory of structure files
#plot_IP_hist(struct_dir, **kwargs)
