import os, sys
import get_close_structs_to_expt, niggli_reduce, update_json_with_relaxed_geo
from update_spg import update_jsons_with_spg

def main():
    pwd = os.path.abspath(os.getcwd())
    input_params = [
        ['glycine',
        '2',
        '/home/trose/genarris_run_calcs/genarris_v2/glycine/2mpc/k333_full_relaxed_rdf',
        'plotting_glycine_2mpc_rdf.conf',
        '/home/trose/genarris_run_calcs/genarris_v2/expts/relax_with_k333/glycine/2mpc/glycine_2mpc_expt_relaxed_niggli_reduced/glycine_2mpc_expt.json'
        ],
        ['glycine',
        '3',
        '/home/trose/genarris_run_calcs/genarris_v2/glycine/3mpc/k333_full_relaxed_rdf',
        'plotting_glycine_3mpc_rdf.conf',
        '/home/trose/genarris_run_calcs/genarris_v2/expts/relax_with_k333/glycine/3mpc/glycine_3mpc_expt_relaxed_niggli_reduced/glycine_3mpc_expt.json'
        ],
        ['glycine',
        '4',
        '/home/trose/genarris_run_calcs/genarris_v2/glycine/4mpc/k333_full_relaxed_rdf',
        'plotting_glycine_4mpc_rdf.conf',
        '/home/trose/genarris_run_calcs/genarris_v2/expts/relax_with_k333/glycine/4mpc/glycine_4mpc_expt_relaxed_niggli_reduced/glycine_4mpc_expt.json'
        ],
        ['benzene',
        '2',
        '/home/trose/genarris_run_calcs/genarris_v2/benzene/2mpc/k333_full_relaxed_pressure_rdf',
        'plotting_benzene_2mpc_rdf.conf',
        '/home/trose/genarris_run_calcs/genarris_v2/expts/relax_with_k333/benzene/2mpc/benzene_2mpc_expt_relaxed_niggli_reduced/benzene_2mpc_expt.json'
        ],
        ['benzene',
        '4',
        '/home/trose/genarris_run_calcs/genarris_v2/benzene/4mpc/k333_full_relaxed_rdf',
        'plotting_benzene_4mpc_rdf.conf',
        '/home/trose/genarris_run_calcs/genarris_v2/expts/relax_with_k333/benzene/4mpc/benzene_4mpc_expt_relaxed_niggli_reduced/benzene_4mpc_expt.json'
        ],
        ['glycine',
        '2',
        '/home/trose/genarris_run_calcs/genarris_v2/glycine/2mpc/k333_full_relaxed_rcd',
        'plotting_glycine_2mpc_rcd.conf',
        '/home/trose/genarris_run_calcs/genarris_v2/expts/relax_with_k333/glycine/2mpc/glycine_2mpc_expt_relaxed_niggli_reduced/glycine_2mpc_expt.json'
        ],
        ['glycine',
        '3',
        '/home/trose/genarris_run_calcs/genarris_v2/glycine/3mpc/k333_full_relaxed_rcd',
        'plotting_glycine_3mpc_rcd.conf',
        '/home/trose/genarris_run_calcs/genarris_v2/expts/relax_with_k333/glycine/3mpc/glycine_3mpc_expt_relaxed_niggli_reduced/glycine_3mpc_expt.json'
        ],
        ['glycine',
        '4',
        '/home/trose/genarris_run_calcs/genarris_v2/glycine/4mpc/k333_full_relaxed_rcd',
        'plotting_glycine_4mpc_rcd.conf',
        '/home/trose/genarris_run_calcs/genarris_v2/expts/relax_with_k333/glycine/4mpc/glycine_4mpc_expt_relaxed_niggli_reduced/glycine_4mpc_expt.json'
        ],
        ['benzene',
        '2',
        '/home/trose/genarris_run_calcs/genarris_v2/benzene/2mpc/k333_full_relaxed_pressure_rcd',
        'plotting_benzene_2mpc_rcd.conf',
        '/home/trose/genarris_run_calcs/genarris_v2/expts/relax_with_k333/benzene/2mpc/benzene_2mpc_expt_relaxed_niggli_reduced/benzene_2mpc_expt.json'
        ],
        ['benzene',
        '4',
        '/home/trose/genarris_run_calcs/genarris_v2/benzene/4mpc/k333_full_relaxed_rcd',
        'plotting_benzene_4mpc_rcd.conf',
        '/home/trose/genarris_run_calcs/genarris_v2/expts/relax_with_k333/benzene/4mpc/benzene_4mpc_expt_relaxed_niggli_reduced/benzene_4mpc_expt.json'
        ],
        ['glycine',
        '4',
        '/home/trose/genarris_run_calcs/genarris_v2/glycine/4mpc/k333_full_relaxed_weakH_rdf',
        'plotting_glycine_4mpc_weakH_rdf.conf',
        '/home/trose/genarris_run_calcs/genarris_v2/expts/relax_with_k333/glycine/4mpc/glycine_4mpc_expt_relaxed_niggli_reduced/glycine_4mpc_expt.json'
        ],
        ['glycine',
        '4',
        '/home/trose/genarris_run_calcs/genarris_v2/glycine/4mpc/k333_full_relaxed_weakH_rcd',
        'plotting_glycine_4mpc_weakH_rcd.conf',
        '/home/trose/genarris_run_calcs/genarris_v2/expts/relax_with_k333/glycine/4mpc/glycine_4mpc_expt_relaxed_niggli_reduced/glycine_4mpc_expt.json'
        ]
    ]

    for input_param in input_params:
        molecule_name, Z, data_dir, plot_conf_file, expt_fpath = input_param
        if plot_conf_file != 'plotting_glycine_4mpc_weakH_rcd.conf' and plot_conf_file != 'plotting_glycine_4mpc_weakH_rdf.conf':
            continue
        #if 'glycine' != molecule_name:
        #    continue
        print(molecule_name, Z)
        if 'benzene' in data_dir:
            napm = 12
        elif 'glycine' in data_dir:
            napm = 10
        else:
            raise Exception('only benzene and glycine supported')
        output_dir = molecule_name + '_' + Z + 'mpc_relaxed_niggli_reduced'
        aims_output_dir = 'aims_' + molecule_name + '_' + Z + 'mpc_relaxed'
        if not os.path.exists(os.path.join(data_dir, aims_output_dir)):
            continue
        
        plot_conf_file = os.path.abspath(plot_conf_file)
        os.chdir(data_dir)
        '''
        niggli_reduce.main(molecule_name, Z, napm)
        update_json_with_relaxed_geo.main(Z, aims_output_dir, output_dir)
        '''
        os.system('python ' + os.path.join(pwd, 'plot_spg_bar_chart_vol_hist_get_lattice_param_data.py') + ' '+ plot_conf_file + ' ' + data_dir)
        os.chdir(data_dir)
        '''
        get_close_structs_to_expt.main(molecule_name, Z, expt_fpath)
        raw_pool_dir = molecule_name + '_' + Z + 'mpc_raw_jsons_niggli_reduced'
        AP1_dir = 'sample_structures_exemplars_niggli_reduced'
        AP2_dir = 'sample_structures_exemplars_2_niggli_reduced'
        relaxed_dir = molecule_name + '_' + Z + 'mpc_relaxed_niggli_reduced'
        structure_dirs = [raw_pool_dir, AP1_dir, AP2_dir, relaxed_dir]
        symprecs = [0.1, 0.01, 0.001, 0.0001]
        for structure_dir in structure_dirs:
            for symprec in symprecs:
                update_jsons_with_spg(structure_dir, int(Z), symprec=symprec, rm_geo=True)
        '''
        os.chdir(pwd)

main()