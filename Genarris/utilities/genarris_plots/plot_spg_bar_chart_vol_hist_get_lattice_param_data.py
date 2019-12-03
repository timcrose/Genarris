from Genarris.utilities import  file_utils
from Genarris.core import instruct
from ibslib.io import read, write
import sys, os
import numpy as np
import plotting

def insert_lattice_lengths(struct_dict, struct_dir):
    directions = ['a', 'b', 'c']
    changed = False
    for struct_id in struct_dict:
        for d in directions:
            if d not in struct_dict[struct_id].properties:
                changed = True
                struct_dict[struct_id].properties[d] = np.linalg.norm(
                        np.array(struct_dict[struct_id].properties['lattice_vector_' + d]))

    if changed:
        write(struct_dir, struct_dict, overwrite=True)
    return struct_dict

def read_dir_into_dict(struct_dir):
    if os.path.exists(struct_dir):
        struct_dict = read(struct_dir)
    else:
        struct_dict = {}
    return struct_dict


def plot_lattice_parameter_space(inst, plot_list):
    for sname in plot_list:
        structure_dir = inst.get(sname, 'structure_dir')
        raw_pool_dir = inst.get('plotting', 'raw_pool_dir')
        energy_dir = inst.get_or_none('plotting', 'energy_dir')
        relaxed_pool_dir = inst.get_or_none('plotting', 'relaxed_pool_dir')
        expt_fpath = inst.get_or_none('plotting', 'expt_fpath')
        xlim = inst.get_with_default(sname, 'lp_xlim', 'automatic')
        ylim = inst.get_with_default(sname, 'lp_ylim', 'automatic')
        zlim = inst.get_with_default(sname, 'lp_zlim', 'automatic')
        xlabel = inst.get_with_default(sname, 'lp_xlabel', 'a ($\AA$)')
        ylabel = inst.get_with_default(sname, 'lp_ylabel', 'b ($\AA$)')
        zlabel = inst.get_with_default(sname, 'lp_zlabel', 'c ($\AA$)')
        cbar_title = inst.get_with_default(sname, 'lp_cbar_title', 'Relative Energy (kJ/mol/molecule)')
        fig_size = inst.get_with_default(sname, 'lp_fig_size', (16,10), eval=True)
        global_min = inst.get_with_default(sname, 'lp_global_min', 'automatic')
        global_max = inst.get_with_default(sname, 'lp_global_max', 'automatic')
        energy_name = inst.get_with_default(sname, 'lp_energy_name', 'energy')
        experimental_energy_name = inst.get_with_default(sname, 'lp_experimental_energy_name', 'energy')
        s = inst.get_with_default(sname, 'lp_s', 200, eval=True)
        alpha = inst.get_with_default(sname, 'lp_alpha', 0.6, eval=True)
        ax_x_labelpad = inst.get_with_default(sname, 'lp_ax_x_labelpad', 35, eval=True)
        ax_y_labelpad = inst.get_with_default(sname, 'lp_ax_y_labelpad', 28, eval=True)
        ax_z_labelpad = inst.get_with_default(sname, 'lp_ax_z_labelpad', 47, eval=True)
        ax_x_tickpad = inst.get_with_default(sname, 'lp_ax_x_tickpad', 0, eval=True)
        ax_y_tickpad = inst.get_with_default(sname, 'lp_ax_y_tickpad', 0, eval=True)
        ax_z_tickpad = inst.get_with_default(sname, 'lp_ax_z_tickpad', 19, eval=True)
        cbar_labelpad = inst.get_with_default(sname, 'lp_cbar_labelpad', 30, eval=True)
        cmap_truncation = inst.get_with_default(sname, 'lp_cmap_truncation', [0.1,0.8], eval=True)
        cbar_rotation = inst.get_with_default(sname, 'lp_cbar_rotation', 0, eval=True)
        figname = inst.get_with_default(sname, 'lp_figname', 'lattice_parameter_plot.png')
        dpi = inst.get_with_default(sname, 'dpi', 1200, eval=True)
        dpi = None
        ax_font_size = inst.get_with_default(sname, 'lp_ax_font_size', 35, eval=True)
        show_x = inst.get_boolean(sname, 'lp_show_x')
        additional_structure_dir = inst.get_or_none('plotting', 'additional_structure_dir')
        show_cbar = inst.get_boolean(sname, 'show_cbar')
        cbar_orientation = inst.get_with_default(sname, 'cbar_orientation', 'horizontal')

        print('sname', sname)
        expt_font_size = inst.get_with_default(sname, 'lp_expt_font_size', 24, eval=True)

        ax_fontdict={   'family' : 'DejaVu Sans',
                        'weight' : 'normal',
                        'size'   : ax_font_size}

        expt_fontdict={ 'family': 'DejaVu Sans',
                                'color':  'lightgreen',
                                'weight': 1000,
                                'size': expt_font_size}

        struct_dict = read(structure_dir)
        raw_pool_dict = read(raw_pool_dir)
        if energy_dir is not None:
            energy_pool_dict = read(energy_dir)
        if relaxed_pool_dir is not None:
            relaxed_pool_dict = read(relaxed_pool_dir)

        if additional_structure_dir is not None:
            #add additional structure directory to struct_dict to get max and min energies and lattice parameters for consistent scales on plot
            additional_raw_pool_dir = os.path.join(additional_structure_dir, os.path.basename(raw_pool_dir))
            additional_energy_dir = os.path.join(additional_structure_dir, os.path.basename(energy_dir))
            additional_relaxed_pool_dir = os.path.join(additional_structure_dir, os.path.basename(relaxed_pool_dir))

            additional_raw_pool_dict = read_dir_into_dict(additional_raw_pool_dir)
            additional_energy_dict = read_dir_into_dict(additional_energy_dir)
            additional_relaxed_pool_dict = read_dir_into_dict(additional_relaxed_pool_dir)

            additional_raw_pool_dict = insert_lattice_lengths(additional_raw_pool_dict, additional_raw_pool_dir)
            additional_energy_dict = insert_lattice_lengths(additional_energy_dict, additional_energy_dir)
            additional_relaxed_pool_dict = insert_lattice_lengths(additional_relaxed_pool_dict, additional_relaxed_pool_dir)

        struct_dict = insert_lattice_lengths(struct_dict, structure_dir)
        raw_pool_dict = insert_lattice_lengths(raw_pool_dict, raw_pool_dir)
        if energy_dir is not None:
            energy_pool_dict = insert_lattice_lengths(energy_pool_dict, energy_dir)
        if relaxed_pool_dir is not None:
            relaxed_pool_dict = insert_lattice_lengths(relaxed_pool_dict, relaxed_pool_dir)

        print('additional_structure_dir', additional_structure_dir, flush=True)
    
        if global_min == 'automatic' or global_max == 'automatic':
            if energy_dir is not None and relaxed_pool_dir is not None:
                energies = [energy_pool_dict[struct_id].properties[energy_name] for struct_id in energy_pool_dict if energy_name in energy_pool_dict[struct_id].properties and energy_pool_dict[struct_id].properties[energy_name] != 'none'] + \
                        [relaxed_pool_dict[struct_id].properties[energy_name] for struct_id in relaxed_pool_dict if energy_name in relaxed_pool_dict[struct_id].properties and relaxed_pool_dict[struct_id].properties[energy_name] != 'none']
            if additional_structure_dir is not None:
                energies += \
                    [additional_energy_dict[struct_id].properties[energy_name] for  struct_id in additional_energy_dict if energy_name in additional_energy_dict[struct_id].properties and additional_energy_dict[struct_id].properties[energy_name] != 'none'] + \
                    [additional_relaxed_pool_dict[struct_id].properties[energy_name] for  struct_id in additional_relaxed_pool_dict if energy_name in additional_relaxed_pool_dict[struct_id].properties and additional_relaxed_pool_dict[struct_id].properties[energy_name] != 'none']
            if expt_fpath is not None:
                expt_struct = read(expt_fpath)
                if experimental_energy_name in expt_struct.properties:
                    expt_energy = expt_struct.properties[experimental_energy_name]
                    if expt_energy != 'none':
                        energies.append(expt_energy)
            if energy_dir is not None and relaxed_pool_dir is not None:            
                if global_min == 'automatic':
                    global_min = min(energies)
                    print('global_min energy', global_min)
                if global_max == 'automatic':
                    global_max = max(energies)
                    print('global_max energy', global_max)

        if xlim == 'automatic':
            if energy_dir is not None and relaxed_pool_dir is not None:
                a = [raw_pool_dict[struct_id].properties['a'] for struct_id in raw_pool_dict] + \
                    [relaxed_pool_dict[struct_id].properties['a'] for struct_id in relaxed_pool_dict]
            else:
                a = [raw_pool_dict[struct_id].properties['a'] for struct_id in raw_pool_dict]
            if additional_structure_dir is not None:
                a += \
                    [additional_raw_pool_dict[struct_id].properties['a'] for struct_id in additional_raw_pool_dict] + \
                    [additional_relaxed_pool_dict[struct_id].properties['a'] for struct_id in additional_relaxed_pool_dict]
            xlim = [int(min(a)), int(np.ceil(max(a)))]
            print('xlim (a)', xlim)
        if ylim == 'automatic':
            if energy_dir is not None and relaxed_pool_dir is not None:
                b = [raw_pool_dict[struct_id].properties['b'] for struct_id in raw_pool_dict] + \
                    [relaxed_pool_dict[struct_id].properties['b'] for struct_id in relaxed_pool_dict]
            else:
                b = [raw_pool_dict[struct_id].properties['b'] for struct_id in raw_pool_dict]
            if additional_structure_dir is not None:
                b += \
                [additional_raw_pool_dict[struct_id].properties['b'] for struct_id in additional_raw_pool_dict] + \
                [additional_relaxed_pool_dict[struct_id].properties['b'] for struct_id in additional_relaxed_pool_dict]
            ylim = [int(min(b)), int(np.ceil(max(b)))]
            print('ylim (b)', ylim)
        if zlim == 'automatic':
            if energy_dir is not None and relaxed_pool_dir is not None:
                c = [raw_pool_dict[struct_id].properties['c'] for struct_id in raw_pool_dict] + \
                    [relaxed_pool_dict[struct_id].properties['c'] for struct_id in relaxed_pool_dict]
            else:
                c = [raw_pool_dict[struct_id].properties['c'] for struct_id in raw_pool_dict]
            if additional_structure_dir is not None:
                c += \
                [additional_raw_pool_dict[struct_id].properties['c'] for struct_id in additional_raw_pool_dict] + \
                [additional_relaxed_pool_dict[struct_id].properties['c'] for struct_id in additional_relaxed_pool_dict]
            zlim = [int(min(c)), int(np.ceil(max(c)))]
            print('zlim (c)', zlim)
        
            
    
        plotting.lattice_parameter_plot(structure_dir, expt_fpath=expt_fpath, energy_name=energy_name,
              xlim=xlim, ylim=ylim, zlim=zlim,
              xlabel=xlabel, ylabel=ylabel, zlabel=zlabel,
              cbar_title=cbar_title,
              fig_size=fig_size, global_min=global_min, global_max=global_max,
              experimental_energy_name=experimental_energy_name,s=s, alpha=alpha,
              ax_x_labelpad=ax_x_labelpad, ax_y_labelpad=ax_y_labelpad, ax_z_labelpad=ax_z_labelpad, 
              ax_x_tickpad=ax_x_tickpad, ax_y_tickpad=ax_y_tickpad, ax_z_tickpad=ax_z_tickpad,
              cbar_labelpad=cbar_labelpad,cmap_truncation=cmap_truncation,
              cbar_rotation=cbar_rotation,
              show_cbar=show_cbar,
              ax_fontdict=ax_fontdict,
                expt_fontdict=expt_fontdict,
                figname=figname, dpi=dpi, show_x=show_x, cbar_orientation=cbar_orientation)



def plot_property(inst, plot_list):
    for sname in plot_list:
        structure_dir = inst.get(sname, 'structure_dir')
        prop = inst.get_with_default(sname, 'prop', 'unit_cell_volume')
        prop_figname = inst.get_with_default(sname, 'prop_figname', 'raw_pool_volume_histogram.pdf')
        prop_xlabel = inst.get_with_default(sname, 'prop_xlabel', 'Structure Volume, $\AA^3$')
        prop_ylabel = inst.get_with_default(sname, 'prop_ylabel', 'Number of structures')
        prop_figure_size = inst.get_with_default(sname, 'prop_figure_size', (12,8), eval=True)
        prop_label_size = inst.get_with_default(sname, 'prop_label_size', 30, eval=True)
        prop_tick_size = inst.get_with_default(sname, 'prop_tick_size', 26, eval=True)
        prop_tick_width = inst.get_with_default(sname, 'prop_tick_width', 3, eval=True)
        prop_GAtor_IP = inst.get_boolean(sname, 'prop_GAtor_IP')
        expt_fpath = inst.get_or_none('plotting', 'expt_fpath')
        bins = inst.get_with_default(sname, 'prop_bins', 'sqrt')
        predicted_volume = inst.get_or_none('plotting', 'predicted_volume', eval=True)
        raw_pool_dir = inst.get('plotting', 'raw_pool_dir')
        relaxed_pool_dir = inst.get_or_none('plotting', 'relaxed_pool_dir')
        expt_fpath = inst.get_or_none('plotting', 'expt_fpath')
        additional_structure_dir = inst.get_or_none('plotting', 'additional_structure_dir')
        y_axis_side = inst.get_with_default(sname, 'prop_y_axis_side', 'left')
        dpi = inst.get_with_default(sname, 'dpi', 1200, eval=True)
        dpi = None
        xtick_label_frequency = inst.get_with_default(sname, 'prop_xtick_label_frequency', 'all')
        xtick_labelsize = inst.get_with_default(sname, 'prop_xtick_labelsize', 24, eval=True)
        xtick_labelrotation = inst.get_with_default(sname, 'prop_xtick_labelrotation', 0, eval=True)
        if xtick_label_frequency != 'all':
            xtick_label_frequency = eval(xtick_label_frequency)

        if bins.isnumeric():
            bins = int(bins)

        struct_dict = read_dir_into_dict(structure_dir)
        struct_dict = insert_lattice_lengths(struct_dict, structure_dir)
        vol_list = [struct_dict[struct_id].properties['unit_cell_volume'] for struct_id in struct_dict]

        if os.path.exists(raw_pool_dir):
            raw_pool_dict = read(raw_pool_dir)
            raw_pool_dict = insert_lattice_lengths(raw_pool_dict, raw_pool_dir)
            vol_list += [raw_pool_dict[struct_id].properties['unit_cell_volume'] for struct_id in raw_pool_dict]

        if relaxed_pool_dir is not None:
            relaxed_pool_dict = read(relaxed_pool_dir)
            relaxed_pool_dict = insert_lattice_lengths(relaxed_pool_dict, relaxed_pool_dir)
            vol_list += [relaxed_pool_dict[struct_id].properties['unit_cell_volume'] for struct_id in relaxed_pool_dict]

        if additional_structure_dir is not None:
            #add additional structure directory to struct_dict to get max and min volume range for consistent scales on plot
            additional_raw_pool_dir = os.path.join(additional_structure_dir, os.path.basename(raw_pool_dir))
            additional_relaxed_pool_dir = os.path.join(additional_structure_dir, os.path.basename(relaxed_pool_dir))
            if os.path.exists(additional_raw_pool_dir):
                additional_raw_pool_dict = read_dir_into_dict(additional_raw_pool_dir)
                additional_raw_pool_dict = insert_lattice_lengths(additional_raw_pool_dict, additional_raw_pool_dir)
                vol_list += [additional_raw_pool_dict[struct_id].properties['unit_cell_volume'] for struct_id in additional_raw_pool_dict]

            if os.path.exists(additional_relaxed_pool_dir):
                additional_relaxed_pool_dict = read_dir_into_dict(additional_relaxed_pool_dir)
                additional_relaxed_pool_dict = insert_lattice_lengths(additional_relaxed_pool_dict, additional_relaxed_pool_dir)
                vol_list += [additional_relaxed_pool_dict[struct_id].properties['unit_cell_volume'] for struct_id in additional_relaxed_pool_dict]

        if expt_fpath is not None:
            expt_struct = read(expt_fpath)
            vol_list.append(expt_struct.properties['volume_before_relaxation'])

        hist_range = (min(vol_list), max(vol_list))

        plotting.plot_property(structure_dir, prop=prop, nmpc=-1, figname=prop_figname, 
                            xlabel=prop_xlabel, ylabel=prop_ylabel, figure_size=prop_figure_size,
                            label_size=prop_label_size, tick_size=prop_tick_size, tick_width=prop_tick_width, GAtor_IP=prop_GAtor_IP,
                            expt_fpath=expt_fpath, bins=bins, predicted_volume=predicted_volume, hist_range=hist_range, y_axis_side=y_axis_side,
                            xtick_label_frequency=xtick_label_frequency, xtick_labelsize=xtick_labelsize, xtick_labelrotation=xtick_labelrotation)


def plot_spg_bar_chart(inst, plot_list):
    for sname in plot_list:
        structure_dir = inst.get(sname, 'structure_dir')
        pygenarris_outfile = inst.get_with_default(sname, 'pygenarris_outfile', 'outfile')
        spg_bar_chart_fname = inst.get_with_default(sname, 'spg_bar_chart_fname', 'raw_pool_spg_bar_chart.pdf')
        spg_bar_width = inst.get_with_default(sname, 'spg_bar_width', 0.5, eval=True)
        spg_bar_xlabel = inst.get_with_default(sname, 'spg_bar_xlabel', 'Allowed space groups')
        spg_bar_ylabel = inst.get_with_default(sname, 'spg_bar_ylabel', 'Number of structures')
        spg_bar_title = inst.get_or_none(sname, 'spg_bar_title')
        dpi = inst.get_with_default(sname, 'dpi', 1200, eval=True)
        dpi = None
        legend_fontsize = inst.get_with_default(sname, 'spg_legend_fontsize', 14, eval=True)
        show_legend = inst.get_boolean(sname, 'spg_show_legend')
        fontsize = inst.get_with_default(sname, 'spg_fontsize', 30, eval=True)
        raw_pool_dir = inst.get('plotting', 'raw_pool_dir')
        relaxed_pool_dir = inst.get_or_none('plotting', 'relaxed_pool_dir')
        expt_fpath = inst.get_or_none('plotting', 'expt_fpath')
        additional_structure_dir = inst.get_or_none('plotting', 'additional_structure_dir')
        fig_size = inst.get_with_default(sname, 'spg_fig_size', (12,8), eval=True)
        y_axis_side = inst.get_with_default(sname, 'spg_y_axis_side', 'left')
        ytick_labelsize = inst.get_with_default(sname, 'spg_ytick_labelsize', 24, eval=True)
        xtick_label_frequency = inst.get_with_default(sname, 'spg_xtick_label_frequency', 'all')
        xtick_labelsize = inst.get_with_default(sname, 'spg_xtick_labelsize', 24, eval=True)
        xtick_labelrotation = inst.get_with_default(sname, 'spg_xtick_labelrotation', 0, eval=True)
        ax_y_labelpad = inst.get_with_default(sname, 'spg_y_labelpad', 10, eval=True)
        ax_x_labelpad = inst.get_with_default(sname, 'spg_x_labelpad', 3, eval=True)
        spg_arrow_width = inst.get_with_default(sname, 'spg_arrow_width', 0.09, eval=True)
        one_color_for_spg_plot = inst.get_boolean(sname, 'one_color_for_spg_plot')
        use_additional_dir = inst.get_boolean(sname, 'spg_use_additional_dir')
        use_relaxed_dir = inst.get_boolean(sname, 'spg_use_relaxed_dir')
        #spg_show_arrow = inst.get_boolean(sname, 'spg_show_arrow')


        if xtick_label_frequency != 'all':
            xtick_label_frequency = eval(xtick_label_frequency)

        spg_list = []
        if os.path.exists(raw_pool_dir):
            raw_pool_dict = read(raw_pool_dir)
            raw_pool_dict = insert_lattice_lengths(raw_pool_dict, raw_pool_dir)
            spg_list += [raw_pool_dict[struct_id].properties['spg'] for struct_id in raw_pool_dict]

        if relaxed_pool_dir is not None and use_relaxed_dir:
            relaxed_pool_dict = read(relaxed_pool_dir)
            relaxed_pool_dict = insert_lattice_lengths(relaxed_pool_dict, relaxed_pool_dir)
            spg_list += [relaxed_pool_dict[struct_id].properties['spg'] for struct_id in relaxed_pool_dict]

        if additional_structure_dir is not None and use_additional_dir:
            #add additional structure directory to struct_dict to get max and min volume range for consistent scales on plot
            additional_raw_pool_dir = os.path.join(additional_structure_dir, os.path.basename(raw_pool_dir))
            additional_relaxed_pool_dir = os.path.join(additional_structure_dir, os.path.basename(relaxed_pool_dir))
            if os.path.exists(additional_raw_pool_dir):
                additional_raw_pool_dict = read_dir_into_dict(additional_raw_pool_dir)
                additional_raw_pool_dict = insert_lattice_lengths(additional_raw_pool_dict, additional_raw_pool_dir)
                spg_list += [additional_raw_pool_dict[struct_id].properties['spg'] for struct_id in additional_raw_pool_dict]

            if os.path.exists(additional_relaxed_pool_dir):
                additional_relaxed_pool_dict = read_dir_into_dict(additional_relaxed_pool_dir)
                additional_relaxed_pool_dict = insert_lattice_lengths(additional_relaxed_pool_dict, additional_relaxed_pool_dir)
                spg_list += [additional_relaxed_pool_dict[struct_id].properties['spg'] for struct_id in additional_relaxed_pool_dict]

        if expt_fpath is not None and use_additional_dir:
            expt_struct = read(expt_fpath)
            spg_list.append(expt_struct.properties['spg'])

        spg_list = list(set(spg_list))
        spg_list.sort()


        plotting.plot_spg_bar_chart(structure_dir, pygenarris_outfile=pygenarris_outfile, spg_bar_chart_fname=spg_bar_chart_fname,
                            width=spg_bar_width, ylabel=spg_bar_ylabel, xlabel=spg_bar_xlabel,
                            title=spg_bar_title, dpi=dpi, legend_fontsize=legend_fontsize, show_legend=show_legend,
                            fontdict={'family' : 'DejaVu Sans',
                                                    'weight' : 'normal',
                                                    'size'   : fontsize}, expt_fpath=expt_fpath, fig_size=fig_size, y_axis_side=y_axis_side,
                            xtick_label_frequency=xtick_label_frequency, xtick_labelsize=xtick_labelsize, xtick_labelrotation=xtick_labelrotation,
                            ax_x_labelpad=ax_x_labelpad, ax_y_labelpad=ax_y_labelpad, spg_arrow_width=spg_arrow_width, ytick_labelsize=ytick_labelsize,
                            spg_list=spg_list, one_color_for_spg_plot=one_color_for_spg_plot)


def convert_ev_to_kj_per_mol_per_molecule(energies, Z):
    #1 eV = 96.49 kJ/mol
    return energies * 96.49 / float(Z)

def get_lattice_param_data(inst, plot_list):
    sections_with_energies = ['AP_2', 'relaxed_pool']
    Z = inst.get_eval('plotting', 'Z')
    for sname in plot_list:
        energy_name = inst.get_with_default(sname, 'energy_name', 'energy')
        structure_dir = inst.get(sname, 'structure_dir')
        
        lattice_param_data_fname = inst.get(sname, 'lattice_param_data_fname')
        struct_dict = read(structure_dir, file_format='json')
        struct_dict = insert_lattice_lengths(struct_dict, structure_dir)
        
        lattice_param_list = []
        energies = []
        for struct_id in struct_dict:
            if sname in sections_with_energies:
                if energy_name in struct_dict[struct_id].properties:
                    energy = struct_dict[struct_id].properties[energy_name]
                    if energy != 'none':
                        energies.append(energy)
                        a = struct_dict[struct_id].properties['a']
                        b = struct_dict[struct_id].properties['b']
                        c = struct_dict[struct_id].properties['c']
                        lattice_param_list.append([a, b, c])
            else:
                a = struct_dict[struct_id].properties['a']
                b = struct_dict[struct_id].properties['b']
                c = struct_dict[struct_id].properties['c']
                lattice_param_list.append([a, b, c])

        energies = np.array(energies)
        if sname in sections_with_energies:
            kj_per_mol_per_molecule = convert_ev_to_kj_per_mol_per_molecule(energies, Z)    
            # This makes the maximum map to 0, minimum map to 1, and everything else is proportionally in between
            #print(np.max(kj_per_mol_per_molecule), kj_per_mol_per_molecule, np.min(kj_per_mol_per_molecule))
            
            relative_energies = (np.max(kj_per_mol_per_molecule) - kj_per_mol_per_molecule) / (np.max(kj_per_mol_per_molecule) - np.min(kj_per_mol_per_molecule))
        else:
            kj_per_mol_per_molecule = np.zeros(len(lattice_param_list))
            relative_energies = np.zeros(len(lattice_param_list))

        
        data = [lattice_params + [relative_energies[i]] + [kj_per_mol_per_molecule[i]] for i,lattice_params in enumerate(lattice_param_list)]
        file_utils.write_rows_to_csv(lattice_param_data_fname, data)


def main():
    data_dir = sys.argv[2]
    print('data_dir', data_dir)
    os.chdir(data_dir)
    inst_path = sys.argv[1]
    inst = instruct.Instruct()
    #enable getting values from the conf file. In other words, this function
    # will read the config file in a format that SafeConfigParser can
    # parse.
    inst.load_instruction_from_file(inst_path)

    plot_list = inst.get_eval('plotting', 'plot_list')

    plot_property(inst, plot_list)
    plot_spg_bar_chart(inst, plot_list)
    get_lattice_param_data(inst, plot_list)
    plot_lattice_parameter_space(inst, plot_list)

    
main()
