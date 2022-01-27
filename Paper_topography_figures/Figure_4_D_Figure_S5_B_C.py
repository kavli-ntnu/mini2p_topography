### Figure 4D and Figure S5 B and C - Obenhaus et al. 
# NN distance box pots and heatmaps 

from general import print_wilcoxon
import sys, os
import os.path
import pandas as pd 
import datajoint as dj
import cmasher as cmr
from pathlib import Path

# Make plots pretty 
import seaborn as sns
sns.set(style='white')

# Prevent bug in figure export as pdf: 
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

##### IMPORTS ###########################################################################
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
import numpy as np
from general import print_mannwhitneyu, print_ttest_1samp, print_kruskal
from itertools import combinations
from matplotlib import pyplot as plt
from tabulate import tabulate

##### LOAD SCHEMA COMPONENTS #############################################################
from dj_schemas.dj_conn import *

##### GLOBAL PLOTTING OPTIONS ############################################################

##### EXPORT LOCATION ####################################################################
figure_dir = 'YOUR_EXPORT_DIRECTORY/'
figure_dir = Path(figure_dir)




def get_norm_nn_dict(animals,
                     params,
                     param_hash_id_cell='standard',
                     param_hash_session='cf83e1357eefb8bd',
                     norm_col='shuffall',
                     region='MEC',
                     cutoff_n_cells=5,
                     ):
    '''
    NEAREST NEIGHBOUR (NN)
    Loop over parameter combinations and extract normalized NN metrics

    '''

    norm_nns_dict   = {}
    norm_nns_matrix = np.zeros((len(params),len(params))) # Save the average into a matrix to plot heatmap later

    for param_set in list(combinations(params.keys(), 2)):
        sessions = _get_nn_sessions(animals, 
                                    param_set[0], 
                                    param_set[1], 
                                    param_hash_id_cell=param_hash_id_cell,
                                    param_hash_session=param_hash_session,
                                    norm_col=norm_col,
                                    region=region,
                                    cutoff_n_cells=cutoff_n_cells
                                    )
        if not sessions: # Just reverse and try again 
            sessions = _get_nn_sessions(animals, 
                                        param_set[1], 
                                        param_set[0], 
                                        param_hash_id_cell=param_hash_id_cell,
                                        param_hash_session=param_hash_session,
                                        norm_col=norm_col,
                                        region=region,
                                        cutoff_n_cells=cutoff_n_cells
                                        )
            
        # Get average norm NN (ab and ba)
        norm_nns_ab_, norm_nns_ba_ = sessions.fetch('norm_nns_ab','norm_nns_ba')
        norm_nns_dict[f'{labels[param_set[0]]} - {labels[param_set[1]]}'] = np.mean(np.stack([norm_nns_ab_, norm_nns_ba_]), axis=0)
        norm_nns_matrix[params[param_set[0]], params[param_set[1]]] = np.mean(np.mean(np.stack([norm_nns_ab_, norm_nns_ba_]), axis=0))
    
    norm_nns_matrix = norm_nns_matrix.T

    return norm_nns_dict, norm_nns_matrix


def get_interintra_dict(animals,
                        params,
                        param_hash_id_cell='standard',
                        param_hash_session='cf83e1357eefb8bd',
                        region='MEC',
                        cutoff_n_cells=5,
                        ):
    '''
    INTER INTRA (RATIO)
    Loop over parameter combinations and extract inter / intra ratio (NN) metrics

    '''

    ratio_dict   = {}
    ratio_matrix = np.zeros((len(params),len(params))) # Save the average into a matrix to plot heatmap later

    for param_set in list(combinations(params.keys(), 2)):
        sessions = _get_interintra_sessions(animals, 
                                    param_set[0], 
                                    param_set[1], 
                                    param_hash_id_cell=param_hash_id_cell,
                                    param_hash_session=param_hash_session,
                                    region=region,
                                    cutoff_n_cells=cutoff_n_cells
                                    )
        if not sessions: # Just reverse and try again 
            sessions = _get_interintra_sessions(animals, 
                                        param_set[1], 
                                        param_set[0], 
                                        param_hash_id_cell=param_hash_id_cell,
                                        param_hash_session=param_hash_session,
                                        region=region,
                                        cutoff_n_cells=cutoff_n_cells
                                        )
            
        # Get average norm NN (ab and ba)
        ratio_ab_, ratio_ba_ = sessions.fetch('ratio_ab','ratio_ba')
        ratio_dict[f'{labels[param_set[0]]} - {labels[param_set[1]]}'] = np.mean(np.stack([ratio_ab_, ratio_ba_]), axis=0)
        ratio_matrix[params[param_set[0]], params[param_set[1]]]       = np.mean(np.mean(np.stack([ratio_ab_, ratio_ba_]), axis=0))
    
    ratio_matrix = ratio_matrix.T

    return ratio_dict, ratio_matrix



def _get_interintra_sessions(animals, 
                             param_A, 
                             param_B, 
                             param_hash_id_cell,
                             param_hash_session,
                             region,
                             cutoff_n_cells,
                             ):
    '''
    Helper for finding sessions and projecting out inter / intra ratios (NN)

    '''

    sessions = (Session.proj('animal_name') * NNeighbourInterIntra.RatioSub * NNeighbourInterIntra.Cells  
                            & f'param_hash_id_cell = "{param_hash_id_cell}"'
                            & f'param_hash_session = "{param_hash_session}"'
                            & f'pairwise_dist_param_A = "{param_A}"' 
                            & f'pairwise_dist_param_B = "{param_B}"'
                            & f'region="{region}"' 
                            & f'n_startr_a > {cutoff_n_cells}'
                            & f'n_startr_b > {cutoff_n_cells}' 
                            & [f'animal_name = "{animal}"' for animal in animals])

    
    unique_animals = set(sessions.fetch('animal_name'))
    print(f'{labels[param_A]:>10} vs {labels[param_B]:<10} | Found {len(sessions):>3} sessions over {len(unique_animals):>3} animals  |  {list(unique_animals)}')
    
    # NO NORMALIZATION !
    sessions_proj = sessions.proj('ratio_ab','ratio_ba')    

    return sessions_proj


def _get_nn_sessions(animals, 
                     param_A, 
                     param_B, 
                     param_hash_id_cell,
                     param_hash_session,
                     norm_col,
                     region,
                     cutoff_n_cells,
                     ):
    '''
    Helper for finding sessions and projecting out a normalized distance metric (NN)

    '''
    #filtered_sessions = get_filtered_sessions(animals, verbose=False).proj()
    sessions = (Session.proj('animal_name') * NNeighbourInterIntra.NNSub * NNeighbourInterIntra.Cells  
                                & f'param_hash_id_cell = "{param_hash_id_cell}"'
                                & f'param_hash_session = "{param_hash_session}"'
                                & f'pairwise_dist_param_A = "{param_A}"' 
                                & f'pairwise_dist_param_B = "{param_B}"'
                                & f'region="{region}"' 
                                & f'n_startr_a > {cutoff_n_cells}'
                                & f'n_startr_b > {cutoff_n_cells}' 
                                & [f'animal_name = "{animal}"' for animal in animals])

    
    unique_animals = set(sessions.fetch('animal_name'))
    print(f'{param_A:>3} vs {param_B:<3} | Found {len(sessions):>3} sessions over {len(unique_animals):>3} animals  |  {list(unique_animals)}')
    
    # Project 
    ab_norm_col = f'nns_ab_{norm_col}' 
    ba_norm_col = f'nns_ba_{norm_col}'
    
    sessions_proj = sessions.proj(norm_nns_ab=f'nns_ab/{ab_norm_col}', 
                                  norm_nns_ba=f'nns_ba/{ba_norm_col}',
                                  )   
    return sessions_proj





def plot_boxplot_summary(results,
                         save_name = 'nn boxplot',
                         population_mean = 1.,
                         ylim = None,
                         figure_dir = figure_dir
                         ):
    '''
    Draw box / scatter plot of results

    Parameter
    ---------
    results : dictionary

    '''

    ## SCATTER PLOT 
    sns.set(style='white', font_scale=2)

    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['xtick.bottom'] = True
    plt.rcParams['ytick.left'] = True

    # Box plot properties 
    boxprops    = dict(color='k', linewidth=1.5)
    medianprops = dict(color='orange', linewidth=0, ls='-', solid_capstyle='butt') # HIDE MEDIAN LINE! 
    meanprops   = dict(color='k', linewidth=2, ls='-', solid_capstyle='butt')


    figure = plt.figure(figsize=(5,5))
    ax = figure.add_subplot(111)

    xlabels = []
    values  = []

    counter=0
    for key, v in results.items():
        
        if len(v) and len(v)>1:
            # random jitter for display
            rnds = np.random.rand(len(v))-.5
            rnds /= 4.5
            ax.scatter(np.zeros_like(v)+counter+.5+rnds,v, color='#333', alpha=.2, s=60, lw=0)
            ax.boxplot(v,
                       positions=[counter+.85], 
                       widths=.25,
                       showfliers=False,
                       whis=[1,99],
                       showmeans=True,
                       meanline=True,
                       boxprops=boxprops,
                       medianprops=medianprops,
                       meanprops=meanprops,
                       )

            xlabels.append(key)
            values.append(v)

            # 1 Sample ttest to theoretical mean of 1.
            print_ttest_1samp(v,population_mean,label=key)
            print_wilcoxon(np.array(v)-population_mean,label=key)
            v = np.array(v).astype(float)
            print(f'Mean ± SD: {np.nanmean(v):.2f}±{np.nanstd(v):.2f}')
            counter+=1

    x_pos = np.arange(len(xlabels))+.75
    ax.set_xticks(x_pos)
    ax.set_xticklabels(xlabels, rotation=40, ha='right')
    ax.axhline(y=population_mean, color='k', ls=':') 
    ax.set_xlim(0, len(x_pos)+.5)
    if ylim is not None: 
        ax.set_ylim(ylim)

    sns.despine(left=True,bottom=True)

    figure.savefig(figure_dir / f'{save_name}.pdf', dpi=300, bbox_inches='tight')
    return values, xlabels



def plot_heatmap(matrix, 
                 labels, 
                 params,
                 save_name='NN heatmap'):

    ## HEATMAP 

    sns.set(style='white', font_scale=2.0)
    # Generate a mask for the upper triangle
    mask = np.triu(np.ones_like(matrix, dtype=bool))
    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(250, 24, s=100, l=70, as_cmap=True)

    figure = plt.figure(figsize=(5,5))
    ax = figure.add_subplot(111)


    # Draw the heatmap with the mask and correct aspect ratio
    sns.heatmap(matrix, 
                mask=mask, 
                cmap=cmap, 
                vmin=0.9, 
                vmax=1.1, 
                center=1,
                square=True, 
                linewidths=1, 
                cbar_kws={"shrink": .5}, 
                ax=ax
                )

    ax.set_xticks(np.arange(4)+.5)
    ax.set_xticklabels([labels[label] for label in params.keys()], rotation=40)
    ax.set_yticks(np.arange(4)+.5)
    ax.set_yticklabels([labels[label] for label in params.keys()], rotation=0)

    figure.savefig(figure_dir / f'{save_name}.pdf', bbox_inches='tight')

    return 


def _quick_stats(values, labels):

    ####### Kruskal-Wallis H-test for independent samples
    print_kruskal(values)

    ####### Pairwise stats (results from 'plot_boxplot_summary')
    print('Pairwise stats (Mann Whitney) of NN')    

    pairw_stats_p = np.zeros((len(values),len(values)))
    pairw_stats_U = np.zeros((len(values),len(values)))

    for (row, value_row), label_A in zip(enumerate(values), labels):
        for (col, value_col), label_B in zip(enumerate(values), labels):
            u, p =  print_mannwhitneyu(value_row, value_col, label_A=label_A, label_B=label_B)
            pairw_stats_p[row, col] = p
            pairw_stats_U[row, col] = u

    pairw_stats_p = pd.DataFrame(pairw_stats_p)
    pairw_stats_U = pd.DataFrame(pairw_stats_U)

    pairw_stats_p.columns   = labels
    pairw_stats_p['labels'] = labels
    pairw_stats_p.set_index('labels', inplace=True)
    print('Mann-Whitney p-values')
    print(tabulate(pairw_stats_p, headers='keys', tablefmt='psql'))
    
    pairw_stats_U.columns   = labels
    pairw_stats_U['labels'] = labels
    pairw_stats_U.set_index('labels', inplace=True)
    print('Mann-Whitney U values')
    print(tabulate(pairw_stats_U, headers='keys', tablefmt='psql'))

    return






if __name__ == "__main__":
  
    animals = ['88592','87244','82913','60480',
               '87187','88106','87245','90222',
               '90218','90647','89841','89622',
               '94557','97045','97046']
    region = 'MEC'

    # For figure legends ... 
    labels = {
                'A' : 'Grid95',
                'B' : 'Grid99',
                'C' : 'HD95',
                'D' : 'HD99',
                'E' : 'OV',
                'F' : 'Border95',
                'G' : 'BVS95',
                'H' : 'Info95',  
                'I' : 'Border99',
                'J' : 'BVS99',
            }

    cutoff_n_cells = 5
    region = 'MEC'
    nn_norm_column = 'shuffall'



    ###############################################################################################################
    ###############################################################################################################
    ####### GET NN METRICS  #######################################################################################
    # Key table: NNeighbourInterIntra.NNSub
    print(f'{"#"*200}\nNN (Nearest neighbour)')
    params = {
              'A' : 0,
              'C' : 1,
              'E' : 2,
              'F' : 3,
              }

    # Get dictionary of session averages for each score combination (norm_nns_dict)
    # and a matrix of averages over all sessions (norm_nns_matrix)
    norm_nns_dict, norm_nns_matrix = get_norm_nn_dict(animals,
                                                      params,
                                                      param_hash_id_cell='standard',
                                                      param_hash_session='cf83e1357eefb8bd',
                                                      norm_col=nn_norm_column,
                                                      region=region,
                                                      cutoff_n_cells=cutoff_n_cells,
                                                      )

    ####### PLOT BOX PLOTS ########################################################################################
    print('Generating NN summary plot')
    nns_, labels_ = plot_boxplot_summary(norm_nns_dict, 
                                         save_name= f'Norm NN stat {region}'
                                         )
    print('Stats')
    _quick_stats(nns_, labels_)

    ####### PLOT HEATMAP ##########################################################################################
    print('Generating heatmap')
    plot_heatmap(norm_nns_matrix, labels, params, save_name= f'Norm NN heatmap {region}')



    # ##############################################################################################################
    # ##############################################################################################################
    # ####### GET RATIO (INTER/INTRA) METRICS  #####################################################################
    # # Key table: NNeighbourInterIntra.RatioSub
    # print(f'{"#"*200}\nINTER INTRA (inter nearest neighbour / intra nearest neighbour')

    # # Get dictionary of session averages for each score combination (ratio_dict)
    # # and a matrix of averages over all sessions (ratio_matrix)
    # ratio_dict, ratio_matrix       = get_interintra_dict(animals,
    #                                                      params,
    #                                                      param_hash_id_cell='standard',
    #                                                      param_hash_session='cf83e1357eefb8bd',
    #                                                      region=region,
    #                                                      cutoff_n_cells=cutoff_n_cells,
    #                                                      )

    # ####### PLOT BOX PLOTS ########################################################################################
    # print('Generating Inter/Intra ratio summary plot')
    # ratios_, labels_ = plot_boxplot_summary(ratio_dict, 
    #                                         save_name= f'InterIntra stat {region}',
    #                                         )
    # print('Stats')
    # _quick_stats(ratios_, labels_)

    # ####### PLOT HEATMAP ##########################################################################################
    # print('Generating heatmap')
    # plot_heatmap(ratio_matrix, labels, params, save_name= f'InterIntra ratio heatmap {region}')



    print(figure_dir)
    print('Success.')