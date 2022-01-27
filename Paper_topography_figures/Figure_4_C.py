### Figure 4 - Panel C - Obenhaus et al.
# Cell maps and distribution plots

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
from scipy.stats import mannwhitneyu, ks_2samp
from helpers_topography.nn_dists import get_nn_pwd_dists, plot_blob_set
from helpers_topography.plotting_helpers import plot_cumulative_horizontal

##### LOAD SCHEMA COMPONENTS #############################################################
from dj_schemas.dj_conn import *

##### EXPORT LOCATION ####################################################################
figure_dir = 'YOUR_EXPORT_DIRECTORY/'
figure_dir = Path(figure_dir)



def get_cell_centers(key, 
                     params_A, 
                     params_B, 
                     param_hash_id_cell='standard',
                     region='MEC',
                     cutoff_n_cells=20
                     ):


    '''
    Return cell centers for populations A and B. 
    These are filtered cells which match a certain score cutoff criterion. 


    '''

    # Get cell parameter dictionary 
    cell_parameter_dict = (FilteredCellsParams & f'param_hash_id_cell = "{param_hash_id_cell}"').fetch1('parameter_dict_cell')


    # Brain region filter 
    assert region in ['MEC','PAS'], f'Region "{region}" not understood. Choose "MEC" or "PAS"'
    if region == 'MEC':
        print('Choosing brain region filter MEC')
        brainregion = RoisCorrBrainLoc.MEC 
    elif region == 'PAS':
        print('Choosing brain region filter PAS')
        brainregion = RoisCorrBrainLoc.PAS



    # Get filtered cells 
    cells_session_part = (eval(params_A['scoretables']) * eval(params_B['scoretables']) * RoisCorr 
                                    & FilteredCells 
                                    & key
                                    & cell_parameter_dict 
                                    & brainregion)

 
    # Define starter cells
    starter_cells_A = cells_session_part \
                            & f'{params_A["score"]} >  {params_A["score_cutoff"]}' \
                            & f'{params_B["score"]} <= {params_B["score_cutoff"]}'
    starter_cells_B = cells_session_part \
                            & f'{params_B["score"]} >  {params_B["score_cutoff"]}' \
                            & f'{params_A["score"]} <= {params_A["score_cutoff"]}'

    if (len(starter_cells_A) < cutoff_n_cells) \
        or (len(starter_cells_B) < cutoff_n_cells):
        print('Too few cells found')
        return
        
    # Check exclusivity of sets
    cell_ids_A = set(starter_cells_A.fetch('cell_id'))
    cell_ids_B = starter_cells_B.fetch('cell_id')
    if len(cell_ids_A.intersection(cell_ids_B)):
        raise ValueError('Cell IDs of set A appear in set B')


    cell_centers_starters_A  = np.stack(starter_cells_A.fetch('center_x_corr','center_y_corr')).T
    cell_centers_starters_B  = np.stack(starter_cells_B.fetch('center_x_corr','center_y_corr')).T

    return cell_centers_starters_A, cell_centers_starters_B





def plot_center_blobs(centers_A, 
                      centers_B, 
                      title,
                      label_A,
                      label_B,
                      min_dist_thresh=10):
    '''
    Blob centers of (subsampled) population centers of A vs. B 

    '''

    centersAB = np.concatenate((centers_A, centers_B))
    labelsAB  = np.concatenate((np.zeros_like(centers_A[:,0]), 
                                np.ones_like( centers_B[:,0]))).astype(int)

    chunk_length = np.min([len(centers_A), len(centers_B)])
    print(f'Extracted chunk length: {chunk_length}')


    # ... this here is "only" needed to get a representative / plottable sample back
    _, _, pop_A_sub_sample, pop_B_sub_sample = get_nn_pwd_dists(
                                                                centers_A, 
                                                                centers_B, 
                                                                min_dist_thresh=min_dist_thresh, 
                                                                chunk_length=chunk_length,
                                                                permutations=1
                                                                )

    # Create figure                                                                 
    blobfig = plot_blob_set(
                            centersAB, 
                            labelsAB,
                            pop_A_sub=pop_A_sub_sample, 
                            pop_B_sub=pop_B_sub_sample,
                            legend=True,
                            label_A=label_A,
                            label_B=label_B,
                            )

    blobfig.savefig(figure_dir / f'{title}.pdf', bbox_inches='tight')


    return


def plot_summary_nn(key,
                    params_A, 
                    params_B,
                    title=None,
                    param_hash_id_cell='standard',
                    region='MEC'):


    # Retrieve from SUBSAMPLED results
    # Populations / Distributions
    dist_sub_results = (NNeighbourInterIntra.DistSub 
                        & key 
                        & 'pairwise_dist_param_A = "{}"'.format(params_A['pairwise_dist_param'])
                        & 'pairwise_dist_param_B = "{}"'.format(params_B['pairwise_dist_param'])
                        & f'param_hash_id_cell = "{param_hash_id_cell}"'
                        & f'region="{region}"').fetch1('nns_ab','nns_ba', 'nns_ab_shuffall','nns_ba_shuffall')
    nns_ab_pop, nns_ba_pop, nns_ab_shuff_pop, nns_ba_shuff_pop = dist_sub_results

    # NN means
    # This is only needed to decide which population result to plot:
    # Take the one with largest mean NN 
    nn_sub_results  = (NNeighbourInterIntra.NNSub 
                        & key 
                        & 'pairwise_dist_param_A = "{}"'.format(params_A['pairwise_dist_param'])
                        & 'pairwise_dist_param_B = "{}"'.format(params_B['pairwise_dist_param'])
                        & f'param_hash_id_cell = "{param_hash_id_cell}"'
                        & f'region="{region}"').fetch1('nns_ab','nns_ba')
    nns_ab, nns_ba = nn_sub_results

    if nns_ba > nns_ab:
        pop = nns_ba_pop
        pop_shuff = nns_ba_shuff_pop
    else:
        pop = nns_ab_pop
        pop_shuff = nns_ab_shuff_pop


    bin_width = 5 # microns
    print(f'Choosing a bin width of {bin_width} microns')
    bins = np.arange(0., np.max(pop_shuff), bin_width)

    plot_cumulative_horizontal(values_A=pop_shuff,
                               values_B=pop,
                               bins=bins, 
                               Alabel='Shuff',
                               Blabel='Data', 
                               xlabel=u'NN dist [\u03bcm]', 
                               median_mean='mean',
                               save_path=figure_dir,
                               save_title=title)

    # Quick stats
    _, mannwhitneyp = mannwhitneyu(pop, pop_shuff)
    _, kolmogorovp  = ks_2samp(pop, pop_shuff)
    print(f'\np Mann Whitney: {mannwhitneyp:.6f}')
    print(f'p Kolmogorov  : {kolmogorovp:.6f}')
    return







if __name__ == "__main__":
  
    # A->C: 1f20835f09e28706
    # A->F: 1f20835f09e28706
    # A->E: 0ed19ecd643fdafa
    # E->F: 79cf52ce2fd423a6
    # C->F: 7bc761d3dc5dcae1
    # C->E: bdad271945828f20

    pairwise_dist_param_A = 'A'
    pairwise_dist_param_B = 'C'
    session_name          = '1f20835f09e28706'

    param_hash_id_cell = 'standard'
    cutoff_n_cells = 5
    region = 'MEC'
    min_dist_thresh = 10.

    ###############################################################################################################
    session = Session & f'session_name = "{session_name}"'
    params_A = (PairwDistParams & f'pairwise_dist_param = "{pairwise_dist_param_A}"')
    params_B = (PairwDistParams & f'pairwise_dist_param = "{pairwise_dist_param_B}"')

    # Fetch keys and parameters A and B
    key = session.fetch1('KEY')
    params_A = params_A.fetch1()
    params_B = params_B.fetch1()

    # ... print some info 
    print('Creating cell center figure and NN distribution plots for')
    print(f'Session {key["session_name"]}')
    print('Parameter Set A')
    print(f'Score: {params_A["score"]} | Cutoff: {params_A["score_cutoff"]}')
    print('Parameter Set B')
    print(f'Score: {params_B["score"]} | Cutoff: {params_B["score_cutoff"]}')

    print('\n')

    centers_A, centers_B = get_cell_centers(key, 
                                            params_A, 
                                            params_B,
                                            param_hash_id_cell=param_hash_id_cell,
                                            region=region,
                                            cutoff_n_cells=cutoff_n_cells,
                                            )


    print(f'Found \n{len(centers_A)} centers for population A and \n{len(centers_B)} centers for population B\n')

    #### PLOT BLOB FIGURE #######################################################################################
    title = '_'.join([key['session_name'], params_A['score'], params_B['score']])
    plot_center_blobs(centers_A, 
                      centers_B, 
                      title=title,
                      label_A= params_A['score'],
                      label_B= params_B['score'],
                      min_dist_thresh=min_dist_thresh)


    #### PLOT NN DIST FIGURE ####################################################################################
    plot_summary_nn(key,
                    params_A, 
                    params_B,
                    title=title,
                    param_hash_id_cell=param_hash_id_cell,
                    region=region)



    print(figure_dir)
    print('Success.')