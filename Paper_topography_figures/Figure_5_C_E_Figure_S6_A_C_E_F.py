### Figure 5 C and E - Obenhaus et al. 
# Figure S6 A, C, E and F - Obenhaus et al. 
# 
# NN distance analysis 
# Pairwise distance analysis 
#  
import sys, os
import os.path
import numpy as np
import pandas as pd 
import datajoint as dj
import cmasher as cmr
from tabulate import tabulate
import itertools

# Make plots pretty 
import seaborn as sns
sns.set(style='white')

# Prevent bug in figure export as pdf: 
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

##### IMPORTS ###########################################################################
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
from general import print_mannwhitneyu, print_wilcoxon
from dj_plotter.helpers.plotting_helpers import make_linear_colormap
from helpers_topography.notebooks.pairw_distances import norm_pairw_nn_df, plot_pairw_nn_summary

##### LOAD SCHEMA COMPONENTS #############################################################
from dj_schemas.dj_conn import *

##### EXPORT LOCATION ####################################################################
figure_dir = 'YOUR_EXPORT_DIRECTORY/'



def pairw_dist(animals,
               col_dict,
               param_hash_session='cf83e1357eefb8bd',
               param_hash_id_cell='standard',
               region='MEC',
               pairwise_dist_param='A',
               cutoff_n_starters=0, 
               plot=True
               ):
    # Print col_dict 
    print(f'\nReceived the following column dictionary \n{col_dict}\n')

    # Brain region filter 
    assert region in ['MEC','PAS'], f'Region "{region}" not understood. Choose "MEC" or "PAS"'

    all_sessions = (Session.proj('animal_name') * FilteredSessions 
                    & [f'animal_name = "{animal}"' for animal in animals]
                    & f'param_hash_session = "{param_hash_session}"'
                    )
    # Print pairw dist. parameter 
    score_, score_cutoff_ = (PairwDistParams & f'pairwise_dist_param = "{pairwise_dist_param}"').fetch1('score','score_cutoff')
    print(f'Filtering pairwise distances by {score_} > {score_cutoff_}')

    pairw = (Session.proj('animal_name') * PairwDist.Cells * PairwDist.PairwD 
                    & all_sessions.proj() \
                    & f'param_hash_id_cell = "{param_hash_id_cell}"' 
                    & f'pairwise_dist_param = "{pairwise_dist_param}"' 
                    & f'region = "{region}"'
                    & f'n_startr > {cutoff_n_starters}')

    pairw_df = pd.DataFrame(pairw.fetch(as_dict=True))
    pairw_df.dropna(inplace=True)
    
    colors   = make_linear_colormap(pairw_df.animal_name, categorical=True, cmap='cmr.guppy')


    ### COLS TO NORMALIZE #################################################################################
    cols_to_norm        = col_dict['cols_to_norm'] # ['mean_pairw_dist_shuffref', 'mean_pairw_dist']
    cols_to_norm_label  = col_dict['cols_to_norm_label'] # ['Ref', 'Data']
    norm_to             = col_dict['norm_to'] # mean_pairw_dist_shuffall'   
    cols                = col_dict['cols'] # 'animal_name'

    # Normalize
    pairw_df_norm = norm_pairw_nn_df(pairw_df, cols_to_norm, cols, norm_to)
    pairw_df_norm.reset_index(drop=True, inplace=True)

    # Plot 
    if plot:
        plot_pairw_nn_summary(pairw_df_norm, 
                              cols_to_norm, 
                              colors=colors, 
                              xlabels=cols_to_norm_label, 
                              save_path=figure_dir,
                              label='PairwD')

    # Print statistics 
    print(f'Data over {len(pairw_df.session_name)} datasets (careful! Can be multiplane!) ({len(set(pairw_df.animal_name))} animals)')
    print(f'{set(pairw_df.animal_name)}')


    # Calculate p values MannWhithney and 1 sample Wilcoxon rank
    pairw_df_norm_ = pairw_df_norm[cols_to_norm]
    results = pd.DataFrame(columns = pairw_df_norm_.columns, 
                           index = pairw_df_norm_.columns)
    for (label1, column1), (label2, column2) in itertools.combinations(pairw_df_norm_.items(), 2):
        _ ,results.loc[label1, label2] = _ ,results.loc[label2, label1] = print_mannwhitneyu(column1, column2, label_A=label1, label_B=label2) 
    #print(tabulate(results, headers='keys', tablefmt='psql'))
    print('\nWilcoxon signed rank test (against 1.):')
    for col in cols_to_norm:
        try:
            print_wilcoxon(pairw_df_norm[col] - 1., label=col)
        except ValueError:
            print(f'Skipping column {col} (all zero?)')

    # Print some more stats
    print('Mean and SEM for PairwDist results')
    for col in cols_to_norm:
        mean_col, sem_col = np.nanmean(pairw_df_norm[col]), np.std(pairw_df_norm[col]) / np.sqrt(len(pairw_df_norm[col]))
        print(f'{col:<30} | Mean ± SEM: {mean_col:.2f} ± {sem_col:.2f}')

    return pairw_df_norm, len(set(pairw_df.animal_name)), len(pairw_df.session_name)


def group_nn_dist(animals,
                  col_dict,
                  param_hash_session='cf83e1357eefb8bd',
                  param_hash_id_cell = 'standard',
                  region='MEC',
                  pairwise_dist_param='A',
                  cutoff_n_starters=0, 
                  nn_group_number=5,
                  plot=True
                  ):

    ''' 
    Like pairw_dist() but for PairwDist.NN instead of PairwDist.PairwD, i.e. grouped NN results
    
    nn_group_number : default 5 : Number of NN to consider (group size). 
                            Careful: Zero indexed! 0 = first nearest neighbour
    '''    
    # Print col_dict 
    print(f'\nReceived the following column dictionary \n{col_dict}\n')

    # Brain region filter 
    assert region in ['MEC','PAS'], f'Region "{region}" not understood. Choose "MEC" or "PAS"'

    all_sessions = (Session.proj('animal_name') * FilteredSessions 
                    & [f'animal_name = "{animal}"' for animal in animals]
                    & f'param_hash_session = "{param_hash_session}"'
                    )
    # Print pairw dist. parameter 
    score_, score_cutoff_ = (PairwDistParams & f'pairwise_dist_param = "{pairwise_dist_param}"').fetch1('score','score_cutoff')
    print(f'Filtering pairwise distances by {score_} > {score_cutoff_}')

    nn = (Session.proj('animal_name') * PairwDist.Cells * PairwDist.NN 
                    & all_sessions.proj() 
                    & f'param_hash_id_cell = "{param_hash_id_cell}"' 
                    & f'pairwise_dist_param = "{pairwise_dist_param}"' 
                    & f'region = "{region}"'
                    & f'n_startr > {cutoff_n_starters}')

    nn_df    = pd.DataFrame(nn.fetch(as_dict=True))
    nn_df.dropna(inplace=True) # Important here because apparently some of the stuff can be None

    colors   = make_linear_colormap(nn_df.animal_name, categorical=True, cmap='cmr.guppy')

    # Subselect a specific nn_number = number of NN in result (group size)
    data_cols_pairwDist_NN = ['mean_nn','mean_nn_shuff_all',
                              'mean_nn_shuff_ref','mean_nn_csr'] # All data columns in table
    for col in data_cols_pairwDist_NN:
        nn_df[col] = [res[nn_group_number] for res in nn_df[col]]
    
    ### COLS TO NORMALIZE #################################################################################
    cols_to_norm        = col_dict['cols_to_norm'] 
    cols_to_norm_label  = col_dict['cols_to_norm_label'] 
    norm_to             = col_dict['norm_to'] 
    cols                = col_dict['cols'] 

    # Normalize
    nn_df_norm = norm_pairw_nn_df(nn_df, cols_to_norm, cols, norm_to)
    nn_df_norm.reset_index(drop=True, inplace=True)


    # Plot 
    if plot:
        plot_pairw_nn_summary(nn_df_norm, 
                              cols_to_norm, 
                              colors=colors, 
                              xlabels=cols_to_norm_label, 
                              save_path=figure_dir,
                              label='NN')

    # Print statistics 
    print(f'Data over {len(nn_df.session_name)} datasets (careful! Can be multiplane!) ({len(set(nn_df.animal_name))} animals)')
    print(f'{set(nn_df.animal_name)}')


    # Calculate p values MannWhithney and 1 sample Wilcoxon rank
    nn_df_norm_ = nn_df_norm[cols_to_norm]
    results = pd.DataFrame(columns = nn_df_norm_.columns, 
                           index = nn_df_norm_.columns)
    for (label1, column1), (label2, column2) in itertools.combinations(nn_df_norm_.items(), 2):
        _ ,results.loc[label1, label2] = _ ,results.loc[label2, label1] = print_mannwhitneyu(column1, column2, label_A=label1, label_B=label2) 
    #print(tabulate(results, headers='keys', tablefmt='psql'))
    print('\nWilcoxon signed rank test (against 1.):')
    for col in cols_to_norm:
        #_, onesample_p_data = ttest_1samp(pairw_df_norm[col], 1.)
        try:
            print_wilcoxon(nn_df_norm[col] - 1., label=col)
        except ValueError:
            print(f'Skipping column {col} (all zero?)')
        
    # Print some more stats
    print('Mean and SEM for NN results')
    for col in cols_to_norm:
        mean_col, sem_col = np.nanmean(nn_df_norm[col]), np.std(nn_df_norm[col]) / np.sqrt(len(nn_df_norm[col]))
        print(f'{col:<30} | Mean ± SEM: {mean_col:.2f} ± {sem_col:.2f}')

    return nn_df_norm, len(set(nn_df.animal_name)), len(set(nn_df.session_name))



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
    param_hash_id_cell = 'standard'
    region = 'MEC'

    # Cutoff number of cells 
    cutoff_n_starters = 15. 

    # For NN 
    nn_group_number = 5

    ###### PAIRWISE DISTANCES ####################################################################################

    # print(f'Creating pairwise distance figure for {len(animals)} animal(s)')
    # print(animals)
    # print('\n')

    # # Create column dictionary 
    # col_dict = {}
    # col_dict['cols_to_norm']       = ['mean_pairw_dist_shuffall', 'mean_pairw_dist_shuffref', 'mean_pairw_dist'] 
    #                                   #mean_pairw_dist_shuffref, mean_pairw_dist_shuffall
    # col_dict['cols_to_norm_label'] = ['All', 'Ref', 'Data']
    # col_dict['norm_to']            = 'mean_pairw_dist_shuffall'
    # col_dict['cols']               = 'animal_name'

    # pairw_dist(animals,
    #            col_dict,
    #            param_hash_session='cf83e1357eefb8bd',
    #            param_hash_id_cell=param_hash_id_cell,
    #            region=region,
    #            pairwise_dist_param=pairwise_dist_param,
    #            cutoff_n_starters=cutoff_n_starters,
    #            )

    ####### NN DISTANCES ###########################################################################################
    print('\n########################################################################################################')
    print(f'\nCreating NN distance figure for {len(animals)} animal(s)')
    print(animals)
    print('\n')

    # Create column dictionary 
    col_dict = {}
    col_dict['cols_to_norm']       = ['mean_nn_shuff_all', 'mean_nn_shuff_ref', 'mean_nn'] 
    col_dict['cols_to_norm_label'] = ['All', 'Ref', 'Data']
    col_dict['norm_to']            = 'mean_nn_shuff_all'
    col_dict['cols']               = 'animal_name'
    group_nn_dist(animals,
                  col_dict,
                  param_hash_session='cf83e1357eefb8bd',
                  param_hash_id_cell=param_hash_id_cell,
                  region=region,
                  pairwise_dist_param=pairwise_dist_param,
                  cutoff_n_starters=cutoff_n_starters, 
                  nn_group_number=nn_group_number,
                  plot=True)

    print(figure_dir)
    print('Success.')