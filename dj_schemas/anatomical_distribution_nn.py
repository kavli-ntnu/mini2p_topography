### Scores, etc. over anatomical space - second 
import sys, os
from itertools import combinations
from datetime import datetime, timedelta
import warnings # Disable np.nanmean Runtime warning

import numpy as np
import datajoint as dj
from tqdm.auto import tqdm
from pointpats import PointPattern, PoissonPointProcess
from sklearn.metrics import pairwise_distances
from scipy.stats import ks_2samp, wasserstein_distance

import random
from helpers_topography.utils import chunks

#### LOAD DATABASE #########################################

# Load base schema
from .dj_conn import * 
imhotte = dj.schema(horst_imaging_db)


from .anatomical_distribution4 import PairwDistParams
from helpers_topography.nn_dists import get_nn_pwd_dists, shuffle_nn_pwd_dists, csr_nn_pwd_dists

from dj_schemas.constants import NN_CUTOFF_NO_CELLS, NN_CUTOFF_NO_STARTER_CELLS, MIN_NN_DIST_THRESH



@imhotte 
class NNeighbourInterIntra(dj.Computed):
    definition = """
    # NN scores inter vs. intra score
    -> FilteredSessions
    -> FilteredCellsParams
    -> PairwDistParams.proj(pairwise_dist_param_A = "pairwise_dist_param") 
    -> PairwDistParams.proj(pairwise_dist_param_B = "pairwise_dist_param") 
    -> ProjectionCorr
    ---    

    """
    class Cells(dj.Part):
        definition = """
        # Cell numbers
        -> master
        region                                  :  char(3)         # Brain region (3 letter abbreviation)
        --- 
        n_all                                   :  smallint        # How many cells in total were considered?
        n_startr_a                              :  smallint        # How many starter cells was the average calculated over (population A)? 
        n_startr_b                              :  smallint        # How many starter cells was the average calculated over (population B)?
        """


    #### DISTANCES (POPULATION RESULTS) ####################################################################
    class DistAll(dj.Part):
        definition = """
        # NN distance results per region all cells - raw results
        -> master
        region                                  :  char(3)       
        --- 
        nns_ab                                  :  blob@hottestore 
        nns_ba                                  :  blob@hottestore 
        nns_ab_shuffab                          :  blob@hottestore 
        nns_ba_shuffab                          :  blob@hottestore 
        nns_ab_shuffall                         :  blob@hottestore 
        nns_ba_shuffall                         :  blob@hottestore 
        nns_ab_csr                              :  blob@hottestore 
        nns_ba_csr                              :  blob@hottestore 
        """
    class DistSub(dj.Part):
        definition = """
        # NN distance results per region subsampled populations (size matched) - raw results
        -> master
        region                                  :  char(3)       
        --- 
        nns_ab                                  :  blob@hottestore 
        nns_ba                                  :  blob@hottestore 
        nns_ab_shuffab                          :  blob@hottestore 
        nns_ba_shuffab                          :  blob@hottestore 
        nns_ab_shuffall                         :  blob@hottestore 
        nns_ba_shuffall                         :  blob@hottestore 
        nns_ab_csr                              :  blob@hottestore 
        nns_ba_csr                              :  blob@hottestore 
        """

    #### Nearest neighbour means ###########################################################################
    class NNAll(dj.Part):
        definition = """
        # Mean NN distance results per region all cells
        -> master
        region                                  :  char(3)         
        --- 
        nns_ab                                  :  double
        nns_ba                                  :  double 
        nns_ab_shuffab                          :  double
        nns_ba_shuffab                          :  double
        nns_ab_shuffall                         :  double 
        nns_ba_shuffall                         :  double
        nns_ab_csr                              :  double
        nns_ba_csr                              :  double 
        """
    class NNSub(dj.Part):
        definition = """
        # Mean NN distance results per region subsampled populations (size matched) 
        -> master
        region                                  :  char(3)         
        --- 
        nns_ab                                  :  double 
        nns_ba                                  :  double 
        nns_ab_shuffab                          :  double 
        nns_ba_shuffab                          :  double 
        nns_ab_shuffall                         :  double 
        nns_ba_shuffall                         :  double 
        nns_ab_csr                              :  double 
        nns_ba_csr                              :  double 
        """

    #### Inter to intra ratio #############################################################################
    class RatioAll(dj.Part):
        definition = """
        # Inter to intra NN distances per region all cells
        -> master
        region                                  :  char(3)      
        ---  
        ratio_ab                                :  double
        ratio_ba                                :  double
        ratio_ab_shuffab                        :  double
        ratio_ba_shuffab                        :  double
        ratio_ab_shuffall                       :  double
        ratio_ba_shuffall                       :  double
        ratio_ab_csr                            :  double
        ratio_ba_csr                            :  double
        """
    class RatioSub(dj.Part):
        definition = """
        # Inter to intra NN distances per region subsampled populations (size matched) 
        -> master
        region                                  :  char(3)      
        ---  
        ratio_ab                                :  double
        ratio_ba                                :  double
        ratio_ab_shuffab                        :  double
        ratio_ba_shuffab                        :  double
        ratio_ab_shuffall                       :  double
        ratio_ba_shuffall                       :  double
        ratio_ab_csr                            :  double
        ratio_ba_csr                            :  double
        """


    class Unprocessed(dj.Part):
        definition = """
        # Inter to intra NN distances per region (and cell numbers)
        -> master
        --- 
        """


    @property
    def key_source(self):

        keys = super().key_source.proj() - (super().key_source.proj() & 'pairwise_dist_param_A = pairwise_dist_param_B')

        # Exclude all "reverse" pairs. Unelegant but it works. 
        available_dist_params = PairwDistParams.fetch('pairwise_dist_param')
        available_combinations = list(combinations(available_dist_params, 2))
        key2exclude = []
        for comb in available_combinations:
            exc = (keys & f'pairwise_dist_param_A = "{comb[1]}"' & f'pairwise_dist_param_B = "{comb[0]}"').fetch('KEY')
            key2exclude.append(exc)
        key2exclude = [item for sublist in key2exclude for item in sublist]
        keys = keys - key2exclude - self.Unprocessed

        return keys

    def make(self, key): 
        
        # Retrieve parameters
        params_A = (PairwDistParams.proj(pairwise_dist_param_A='pairwise_dist_param',
                                         score='score',
                                         score_cutoff='score_cutoff',
                                         scoretables='scoretables')\
                                         & key).fetch1()

        params_B = (PairwDistParams.proj(pairwise_dist_param_B='pairwise_dist_param',
                                         score='score',
                                         score_cutoff='score_cutoff',
                                         scoretables='scoretables')\
                                         & key).fetch1()
        
        # There are parameter sets with the same score but different cutoffs. Exclude those
        if params_A['score'] == params_B['score']:
            self.insert1(key, ignore_extra_fields=True)
            self.Unprocessed.insert1(key, ignore_extra_fields=True)
            return

        cell_parameter_dict = (FilteredCellsParams & key).fetch1('parameter_dict_cell')
     

        filtered_cells_AB = eval(params_A['scoretables']) * eval(params_B['scoretables']) \
                         * RoisCorr & FilteredCells & key & cell_parameter_dict

        # GET MEC, PAS, ... ANATOMICAL FILTERS
        mec_filter = RoisCorrBrainLoc.MEC & key
        pas_filter = RoisCorrBrainLoc.PAS & key

        # Define entry dictionaries to loop over for part tables
        entry_dicts = [
                    {
                'cell_filter' : mec_filter,
                'name'        : 'MEC',
                    },
                    {
                'cell_filter' : pas_filter,
                'name'        : 'PAS',
                    },
                ]

        # Insert into master
        self.insert1(key, ignore_extra_fields=True)

        ### LOOP OVER BRAIN REGIONS 
        for entry_dict in entry_dicts:

            cells_session_part = filtered_cells_AB & entry_dict['cell_filter'] # i.e. MEC cells, ...
            n_cells = len(cells_session_part)
            
            if (n_cells < NN_CUTOFF_NO_CELLS): 
                continue

            # Define starter cells
            starter_cells_A = cells_session_part \
                                    & f'{params_A["score"]} > {params_A["score_cutoff"]}' \
                                    & f'{params_B["score"]} <= {params_B["score_cutoff"]}'
            starter_cells_B = cells_session_part \
                                    & f'{params_B["score"]} > {params_B["score_cutoff"]}' \
                                    & f'{params_A["score"]} <= {params_A["score_cutoff"]}'

            if (len(starter_cells_A) < NN_CUTOFF_NO_STARTER_CELLS) \
                or (len(starter_cells_B) < NN_CUTOFF_NO_STARTER_CELLS):
                continue
                
            # Check exclusivity of sets
            cell_ids_A = set(starter_cells_A.fetch('cell_id'))
            cell_ids_B = starter_cells_B.fetch('cell_id')
            if len(cell_ids_A.intersection(cell_ids_B)):
                raise ValueError('Cell IDs of set A appear in set B')
            
            # Define the three populations: 
            # - All centers
            # - A 
            # - B 

            cell_centers_all         = np.stack(cells_session_part.fetch('center_x_corr','center_y_corr')).T
            cell_centers_starters_A  = np.stack(starter_cells_A.fetch('center_x_corr','center_y_corr')).T
            cell_centers_starters_B  = np.stack(starter_cells_B.fetch('center_x_corr','center_y_corr')).T
            cell_centers_starters_AB = np.concatenate((cell_centers_starters_A, cell_centers_starters_B))
               

            # DATA ... 

            # Loop over permuted data 
            # Rational: 
            # Cutting out chunks introduces biased sampling of the underlying population 
            # i.e. when pop A contains 15 samples and pop B contains 40 samples, 10 samples will be missed. 
            # Solution: 
            # Permute the order of centers of population A and B, and calculate population and median from that.

            chunk_length = np.min([len(cell_centers_starters_A), len(cell_centers_starters_B)])

            # FIX GLOBAL # OF ITERATIONS (SHUFFLES)
            n_shuffles = 1000


            # 1. For all data 
            nns_all, ratios_all, _, _ = get_nn_pwd_dists(cell_centers_starters_A, 
                                                         cell_centers_starters_B, 
                                                         min_dist_thresh=MIN_NN_DIST_THRESH
                                                         )

            # 2. For subsampled data
            nns_sub, ratios_sub, _, _ = get_nn_pwd_dists(cell_centers_starters_A, 
                                                         cell_centers_starters_B, 
                                                         min_dist_thresh=MIN_NN_DIST_THRESH,
                                                         chunk_length=chunk_length,
                                                         permutations=n_shuffles
                                                         )
            

            # SHUFFLE ... 
            # SHUFFLE WITH POP A AND POP B AS START 
            # 1. For all data 
            # i.e. asymmetric population sizes 
            nns_shuff_AB, ratios_shuff_AB = shuffle_nn_pwd_dists(cell_centers_starters_AB, 
                                                                 min_dist_thresh=MIN_NN_DIST_THRESH,
                                                                 length_A=len(cell_centers_starters_A), 
                                                                 length_B=len(cell_centers_starters_B), 
                                                                 shuffle_iter=1000
                                                                 )

            # 2. For subsampled data
            nns_shuff_AB_sub, ratios_shuff_AB_sub = shuffle_nn_pwd_dists(cell_centers_starters_AB, 
                                                                         min_dist_thresh=MIN_NN_DIST_THRESH,
                                                                         chunk_length=chunk_length,
                                                                         shuffle_iter=1000
                                                                         )

            # SHUFFLE WITH ALL CELLS AS START 
            # 1. For all data 
            # i.e. asymmetric population sizes 
            nns_shuff_all, ratios_shuff_all = shuffle_nn_pwd_dists(cell_centers_all, 
                                                                   min_dist_thresh=MIN_NN_DIST_THRESH,
                                                                   length_A=len(cell_centers_starters_A), 
                                                                   length_B=len(cell_centers_starters_B), 
                                                                   shuffle_iter=1000
                                                                   )

            # 2. For subsampled data
            nns_shuff_all_sub, ratios_shuff_all_sub = shuffle_nn_pwd_dists(cell_centers_all, 
                                                                           min_dist_thresh=MIN_NN_DIST_THRESH,
                                                                           chunk_length=chunk_length,
                                                                           shuffle_iter=1000
                                                                           )


            # SHUFFLE CSR ... 
            # Generate a CSR (complete spatial randomness) pattern
            # 1. For all data 
            nns_csr_all, ratios_csr_all = csr_nn_pwd_dists(cell_centers_starters_A, 
                                                           cell_centers_starters_B, 
                                                           min_dist_thresh=MIN_NN_DIST_THRESH,
                                                           length_A=len(cell_centers_starters_A), 
                                                           length_B=len(cell_centers_starters_B),
                                                           shuffle_iter=1000
                                                           )
            # 2. For subsampled data
            nns_csr_sub, ratios_csr_sub = csr_nn_pwd_dists(cell_centers_starters_A, 
                                                           cell_centers_starters_B, 
                                                           min_dist_thresh=MIN_NN_DIST_THRESH,
                                                           chunk_length=chunk_length,
                                                           shuffle_iter=1000
                                                           )

            # SAVE ...

            # Cell numbers 
            entry_dict_celln = {
                 'region'                  : entry_dict['name'],
                 'n_all'                   : n_cells,
                 'n_startr_a'              : len(cell_centers_starters_A),
                 'n_startr_b'              : len(cell_centers_starters_B),
            }

            # Populations
            entry_dict_dist_all = {
                 'region'                  : entry_dict['name'],
                 'nns_ab'                  : nns_all['AB'],
                 'nns_ba'                  : nns_all['BA'],
                 'nns_ab_shuffab'          : nns_shuff_AB['AB'],
                 'nns_ba_shuffab'          : nns_shuff_AB['BA'],
                 'nns_ab_shuffall'         : nns_shuff_all['AB'],
                 'nns_ba_shuffall'         : nns_shuff_all['BA'],
                 'nns_ab_csr'              : nns_csr_all['AB'],
                 'nns_ba_csr'              : nns_csr_all['BA'],
            }

            entry_dict_dist_sub = {
                 'region'                  : entry_dict['name'],
                 'nns_ab'                  : nns_sub['AB'],
                 'nns_ba'                  : nns_sub['BA'],
                 'nns_ab_shuffab'          : nns_shuff_AB_sub['AB'],
                 'nns_ba_shuffab'          : nns_shuff_AB_sub['BA'],
                 'nns_ab_shuffall'         : nns_shuff_all_sub['AB'],
                 'nns_ba_shuffall'         : nns_shuff_all_sub['BA'],
                 'nns_ab_csr'              : nns_csr_sub['AB'],
                 'nns_ba_csr'              : nns_csr_sub['BA'],
            }

            # Mean NN 
            entry_dict_nn_all = {
                'region'                   : entry_dict['name'],
                'nns_ab'                   : np.nanmean(nns_all['AB']),
                'nns_ba'                   : np.nanmean(nns_all['BA']),
                'nns_ab_shuffab'           : np.nanmean(nns_shuff_AB['AB']),
                'nns_ba_shuffab'           : np.nanmean(nns_shuff_AB['BA']),
                'nns_ab_shuffall'          : np.nanmean(nns_shuff_all['AB']),
                'nns_ba_shuffall'          : np.nanmean(nns_shuff_all['BA']),
                'nns_ab_csr'               : np.nanmean(nns_csr_all['AB']),
                'nns_ba_csr'               : np.nanmean(nns_csr_all['BA']),
            }

            entry_dict_nn_sub = {
                'region'                   : entry_dict['name'],
                'nns_ab'                   : np.nanmean(nns_sub['AB']),
                'nns_ba'                   : np.nanmean(nns_sub['BA']),
                'nns_ab_shuffab'           : np.nanmean(nns_shuff_AB_sub['AB']),
                'nns_ba_shuffab'           : np.nanmean(nns_shuff_AB_sub['BA']),
                'nns_ab_shuffall'          : np.nanmean(nns_shuff_all_sub['AB']),
                'nns_ba_shuffall'          : np.nanmean(nns_shuff_all_sub['BA']),
                'nns_ab_csr'               : np.nanmean(nns_csr_sub['AB']),
                'nns_ba_csr'               : np.nanmean(nns_csr_sub['BA']),
            }

            # Inter to intra ratios
            entry_dict_ratios_all = {
                'region'                  : entry_dict['name'],
                'ratio_ab'                : ratios_all['AB'],
                'ratio_ba'                : ratios_all['BA'],
                'ratio_ab_shuffab'        : ratios_shuff_AB['AB'],
                'ratio_ba_shuffab'        : ratios_shuff_AB['BA'],
                'ratio_ab_shuffall'       : ratios_shuff_all['AB'],
                'ratio_ba_shuffall'       : ratios_shuff_all['BA'],
                'ratio_ab_csr'            : ratios_csr_all['AB'],
                'ratio_ba_csr'            : ratios_csr_all['BA'],
            }

            entry_dict_ratios_sub = {
                'region'                  : entry_dict['name'],
                'ratio_ab'                : ratios_sub['AB'],
                'ratio_ba'                : ratios_sub['BA'],
                'ratio_ab_shuffab'        : ratios_shuff_AB_sub['AB'],
                'ratio_ba_shuffab'        : ratios_shuff_AB_sub['BA'],
                'ratio_ab_shuffall'       : ratios_shuff_all_sub['AB'],
                'ratio_ba_shuffall'       : ratios_shuff_all_sub['BA'],
                'ratio_ab_csr'            : ratios_csr_sub['AB'],
                'ratio_ba_csr'            : ratios_csr_sub['BA'],
            }

            # Write into table 
            self.Cells.insert1({**key, **entry_dict_celln})

            self.DistAll.insert1({**key, **entry_dict_dist_all})
            self.DistSub.insert1({**key, **entry_dict_dist_sub})

            self.NNAll.insert1({**key, **entry_dict_nn_all}) 
            self.NNSub.insert1({**key, **entry_dict_nn_sub})

            self.RatioAll.insert1({**key, **entry_dict_ratios_all})
            self.RatioSub.insert1({**key, **entry_dict_ratios_sub})