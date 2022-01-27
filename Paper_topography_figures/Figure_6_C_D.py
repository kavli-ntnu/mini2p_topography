### Figure 6 C and D - Obenhaus et al. 
# Scoremap correlation results
# in boxplot and matrix form 
 
import sys, os
import os.path
from pathlib import Path
import numpy as np 
import datajoint as dj
import cmasher as cmr

# Make plots pretty 
import seaborn as sns
sns.set(style='white')

# Prevent bug in figure export as pdf: 
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

##### IMPORTS #########################################################################################
from itertools import combinations
from matplotlib import pyplot as plt
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
from Figure_4_D_Figure_S5_B_C import plot_boxplot_summary
from mpl_toolkits.axes_grid1 import make_axes_locatable


##### LOAD SCHEMA COMPONENTS ########################################################################## 
from dj_schemas.dj_conn import *


##### EXPORT LOCATION #################################################################################
figure_dir = 'YOUR_EXPORT_DIRECTORY/'
figure_dir = Path(figure_dir)


def _plot_matrix_xcorr(corr_matrix, binning_parameters ):
    '''
    Plot color coded matrix of 

    '''
       
    sns.set(style='white', font_scale=1.8)
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['xtick.bottom'] = True
    plt.rcParams['ytick.left'] = True

    figure = plt.figure(figsize=(8,4))

    ax = figure.add_subplot(111)
    g = ax.imshow(corr_matrix, cmap='bwr', vmin=-.65, vmax=.65, extent=[-.5,3.5,-.5,3.5], aspect='equal')

    ax.set_xticks([0,1,2,3]); ax.set_yticks([0,1,2,3])
    ax.set_xticklabels(list(binning_parameters.keys()),       rotation=50, ha='right', va='top') # X 
    ax.set_yticklabels(list(binning_parameters.keys())[::-1], rotation=0,  ha='right', va='center') # Y

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4%", pad=0.1)
    plt.colorbar(g, cax=cax)

    sns.despine(left=True,bottom=True)
    figure.savefig(figure_dir / f'square matrix correlation.pdf', 
                            dpi=600, bbox_inches='tight')





def _plot_barchart_xcorr_sign(no_session_dict):
    '''
    Plot horizontal bar chart showing the amount of 
    xcorr results 
    > threshold high 
    < threshold low 

    Sizes are saved as [len(pos), len(neg), len(rest)]

    ''' 

    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['xtick.bottom'] = True
    plt.rcParams['ytick.left'] = True

    figure = plt.figure(figsize=(4,5))

    N = len(no_session_dict)
    pos  = np.array([v[0] for v in no_session_dict.values()])
    neg  = np.array([v[1] for v in no_session_dict.values()])
    neut = np.array([v[2] for v in no_session_dict.values()])

    ind = np.arange(N)    # the x locations for the groups
    width = 0.9       # the width of the bars: can also be len(x) sequence

    plt.barh(ind, pos,  width, color='r')
    plt.barh(ind, neut, width, left=pos, color='#ccc')
    plt.barh(ind, neg,  width, left=pos+neut, color='cornflowerblue')

    plt.yticks(ind, no_session_dict.keys(), rotation=0, ha='right')
    plt.ylim=(0,.5)
    plt.xlabel('# animals')
    sns.despine(left=True)

    print('Saving figure under {}'.format(str(figure_dir)))
    figure.savefig(figure_dir / f'xcorr no sessions.pdf', dpi=300, bbox_inches='tight')








def get_xcorr_results(animals,
                      proj_param_id,
                      binning_parameters,
                      cutoff_low,
                      cutoff_high,
                      param_hash_session='cf83e1357eefb8bd',
                      param_hash_id_cell='standard',
                      ):

    '''
    Get pairwise FOV scoremap XCORR results for 
    'grid','border_bvs','mvl','ovc'

    '''



    filtered_sessions  = (Session.proj('animal_name') * FilteredSessions 
                            & [f'animal_name = "{animal}"' for animal in animals]
                            & f'param_hash_session = "{param_hash_session}"'
                            )
    corr_entries       = (ScoremapCorr 
                            * ScoremapFOVParams.proj(binning_param_A='binning_param_short', binning_param_set_A='binning_param_set') 
                            * ScoremapFOVParams.proj(binning_param_B='binning_param_short', binning_param_set_B='binning_param_set')
                            & filtered_sessions.proj(metasession_ref='metasession_name') 
                            & f'param_hash_id_cell = "{param_hash_id_cell}"'
                            & f'proj_param_id = "{proj_param_id}"'
                            )


    # Loop over parameter set combinations
    norm_xcorr_dict = {}
    no_session_dict = {}
    corr_matrix = np.zeros((4,4))

    for param_set in list(combinations(binning_parameters.keys(), 2)):
    
        scoremap_corrs = (corr_entries  
                            & 'xcorr is not null' 
                            & f'binning_param_set_A = "{param_set[0]}"' 
                            & f'binning_param_set_B = "{param_set[1]}"')
        if not scoremap_corrs:
            # Turn around and try again 
            scoremap_corrs = (corr_entries  
                            & 'xcorr is not null' 
                            & f'binning_param_set_A = "{param_set[1]}"' 
                            & f'binning_param_set_B = "{param_set[0]}"')
    
        # Print some stats and save number of animals in dictionary: 
        animals_in_result = set((Session.proj(..., metasession_ref='metasession_name') & scoremap_corrs).fetch('animal_name'))
        
        


        # Loop over xcorr results and subtract shuffled mean from xcorr result
        norm_xcorrs = []
        for result in scoremap_corrs:
            # Be aware of positive and negative results
            norm_xcorr = result['xcorr'] - np.nanmedian(result['xcorr_shuffles'])    
            norm_xcorrs.append(norm_xcorr)
        norm_xcorr_dict[f'{param_set[0]}__{param_set[1]}'] = norm_xcorrs
        
        # Above / below shuffling cutoff? 
        pos     = scoremap_corrs & f'xcorr > {cutoff_high}'
        neg     = scoremap_corrs & f'xcorr < {cutoff_low}'
        rest    = scoremap_corrs - pos.proj() - neg.proj()
        sizes   = [len(pos), len(neg), len(rest)]

        no_session_dict[f'{param_set[0]}__{param_set[1]}'] = sizes
        
        # Write to matrix 
        discrim = (len(pos)-len(neg)) / len(scoremap_corrs)
        corr_matrix[binning_parameters[param_set[1]], binning_parameters[param_set[0]]] = discrim

        print(f'{param_set[0]} X {param_set[1]}: Found {len(scoremap_corrs)} results across {len(animals_in_result)} animals')
        print(f'... Discrimination ratio: {discrim:.3f}')

    return norm_xcorr_dict, no_session_dict, corr_matrix





if __name__ == "__main__":

    all_animals = [
               '90222','90218','90647',
               '82913','88592','89622',
               '87244','89841','60480',
               '87245','87187','88106',
               '94557','97045','97046',
               ] 
    print(f'Including {len(all_animals)} animals')
    animals = all_animals
    
    proj_param_id = 'A' #ScoremapFOVProjParams() - printed above
    param_hash_id_cell ='standard'
    param_hash_session ='cf83e1357eefb8bd'

    cutoff_high = 'xcorr_95'
    cutoff_low = 'xcorr_5'

    binning_parameters = {
                          'grid': 0,
                          'border_bvs': 1,
                          'mvl': 2,
                          'ovc': 3,
                          } 



    norm_xcorr_dict, no_session_dict, corr_matrix = get_xcorr_results(animals,
                                                                      proj_param_id,
                                                                      binning_parameters,
                                                                      cutoff_low,
                                                                      cutoff_high,
                                                                      param_hash_session='cf83e1357eefb8bd',
                                                                      param_hash_id_cell='standard',
                                                                      )                 


    
    print('\nPlot norm xcorr results')
    xcorrs_, labels_ = plot_boxplot_summary(norm_xcorr_dict, 
                                            save_name = 'norm xcorr boxplot', 
                                            population_mean = 0.,
                                            figure_dir = figure_dir
                                            )

    #print('Stats')
    #_quick_stats(xcorrs_, labels_)

    # print('Plotting horizontal bar chart of xcorr significances')
    # _plot_barchart_xcorr_sign(no_session_dict)


    print('Plotting color coded matrix of xcorr results')
    _plot_matrix_xcorr(corr_matrix, binning_parameters)


    print(figure_dir)
    print('Success.')