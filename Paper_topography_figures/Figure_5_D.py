### Figure 5D - Obenhaus et al. 
# k-NN (k nearest neighbours) over various NN group sizes

import sys, os
import os.path
import numpy as np
import pandas as pd 
import datajoint as dj
import cmasher as cmr

# Make plots pretty 

import seaborn as sns
sns.set(style='white')

# Prevent bug in figure export as pdf: 
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

##### IMPORTS ###########################################################################
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
from helpers_topography.notebooks.pairw_distances import norm_nn_df, plot_mean_nn_over_nn
from general import print_mannwhitneyu, print_wilcoxon

##### LOAD SCHEMA COMPONENTS #############################################################
from dj_schemas.dj_conn import *

##### EXPORT LOCATION ####################################################################
figure_dir = 'YOUR_EXPORT_DIRECTORY/'


def knn_dist(animals,
            param_hash_session='cf83e1357eefb8bd',
            param_hash_id_cell = 'standard',
            region='MEC',
            pairwise_dist_param='A',
            cutoff_n_starters=0, 
            ):


    # Brain region filter 
    assert region in ['MEC','PAS'], f'Region "{region}" not understood. Choose "MEC" or "PAS"'

    all_sessions = (Session.proj('animal_name') * FilteredSessions 
                    & [f'animal_name = "{animal}"' for animal in animals]
                    & f'param_hash_session = "{param_hash_session}"'
                    )
    # Print pairw dist. parameter 
    score_, score_cutoff_ = (PairwDistParams & f'pairwise_dist_param = "{pairwise_dist_param}"').fetch1('score','score_cutoff')
    print(f'Filtering pairwise distances by {score_} > {score_cutoff_}')

    nns = (Session.proj('animal_name') * PairwDist.Cells * PairwDist.NN 
                    & all_sessions.proj() 
                    & f'param_hash_id_cell = "{param_hash_id_cell}"' 
                    & f'pairwise_dist_param = "{pairwise_dist_param}"' 
                    & f'region = "{region}"'
                    & f'n_startr > {cutoff_n_starters}')

    nns_df   = pd.DataFrame(nns.fetch(as_dict=True))
    # Print statistics 
    print(f'Data over {len(nns_df.session_name)} datasets (Careful! Can be multiplane!) ({len(set(nns_df.animal_name))} animals)')
    print(f'{set(nns_df.animal_name)}')


    # Normalize
    mean_nns, mean_nn_shuff_refs = norm_nn_df(nns_df, norm_to='mean_nn_shuff_all')
    # Plot 
    plot_mean_nn_over_nn(mean_nns, mean_nn_shuff_refs, save_path=figure_dir)

    # Stats
    print('Statistics NN:')
    for col in np.arange(mean_nns.shape[1]):
        print_mannwhitneyu(mean_nns[:,col], mean_nn_shuff_refs[:,col], label_A=pairwise_dist_param, label_B='Ref')
    for col in np.arange(mean_nns.shape[1]):
        try:
            print_wilcoxon(mean_nns[:,col] - 1., label=pairwise_dist_param)
        except ValueError:
            print(f'Skipping column {col} (all zero?)')
        
    return





if __name__ == "__main__":
  
  
    grid_mice = [ 
                '82913','88592', '87244', '60480',
                '97046','89841' 
                ]
    ov_mice   = [ 
                '87187','88106','87245','90222',
                '94557','89622'
                ]

    all_animals = [
                '90222','90218','90647',
                '82913','88592','89622',
                '87244','89841','60480',
                '87245','87187','88106',
                '94557','97045','97046',
                ] 

    animals = grid_mice 
    pairwise_dist_param = "A"
    region = 'MEC'
    param_hash_id_cell = 'standard'

    print(f'Creating k-NN distance figure for {len(animals)} animal(s)')
    print(animals)
    print('\n')

    # Cutoff number of cells 
    cutoff_n_starters = 15. 

    knn_dist(animals,
            param_hash_session='cf83e1357eefb8bd',
            param_hash_id_cell=param_hash_id_cell,
            region=region,
            pairwise_dist_param=pairwise_dist_param,
            cutoff_n_starters=cutoff_n_starters,
            )


    print(figure_dir)
    print('Success.')