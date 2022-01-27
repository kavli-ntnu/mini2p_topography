### Figure S5 - D and E - Obenhaus et al. 
# Shuffling of NN graphs

import sys, os
import os.path
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

#from fig_4_C_E import _get_nn_sessions
from Figure_4_E_and_Figure_S5_D import get_nn_weights, _make_nn_graph, _spring_loaded_model, _get_distance_metric
import numpy as np
from tqdm.auto import trange
import pickle

##### LOAD SCHEMA COMPONENTS #############################################################
from dj_schemas.dj_conn import *


##### GLOBAL PLOTTING OPTIONS ############################################################



def _return_simulation_points(weights,
                              params_A,
                              params_B,
                              iterations,
                              threshold,
                              ):
    '''
    Helper for building graph and running spring loaded simulation.
    Returns node positions after simulation has run through.
    '''

    # Build graph
    nn_graph = _make_nn_graph(weights, params_A, params_B)
    # ... and run simulation
    model_node_pos = _spring_loaded_model(nn_graph, iterations = iterations, threshold = threshold)

    return model_node_pos
 





if __name__ == "__main__":
  
    animals = ['88592','87244','82913','60480',
               '87187','88106','87245','90222',
               '90218','90647','89841','89622',
               '94557','97045','97046']
    region = 'MEC'

    # For figure legends ... 
    # SIMPLIFIED VERSION
    labels = {
                'A' : 'Grid',
                'B' : 'Grid',
                'C' : 'HD',
                'D' : 'HD',
                'E' : 'OV',
                'F' : 'Border',
                'G' : 'Border',
                'H' : 'Info',  
                'I' : 'Border',
                'J' : 'Border',
            }

    cutoff_n_cells = 5
    region = 'MEC'
    nn_norm_column = 'shuffall'



    ##############################################################################################################
    ##############################################################################################################
    ####### GET NN GRAPH  ########################################################################################
    # Key table: NNeighbourInterIntra.NNSub
    print(f'{"#"*200}\nNN (Nearest neighbour) Graph')
    params = {
              'A' : 0,
              'C' : 1,
              'E' : 2,
              'F' : 3,
              }
    grid_key  = 'A'
    other_keys= ('C','E','F')

    export_folder = 'data_exports/NN_graph_shuffling/'

    ###### GET WEIGHTS ##########################################################################################
    weights, params_A, params_B = get_nn_weights(animals,
                                                 params,
                                                 param_hash_id_cell='standard',
                                                 param_hash_session='cf83e1357eefb8bd',
                                                 norm_col=nn_norm_column,
                                                 region=region,
                                                 cutoff_n_cells=cutoff_n_cells,
                                                 )


    ###### GET METRICS AFTER SIMULATION ##########################################################################
    # Loop over simulation lengths (iterations)
    for iter in [250, 500, 1000, 5000, 10000]:   #np.arange(0, 10000, 250):
        print(f'Shuffling for {iter} simulation steps')
        shuffle_iter = 1000

        data_dists = [] 
        for i in trange(shuffle_iter, desc='Data'):
            # Helper for getting graph / simulation
            model_node_pos = _return_simulation_points(weights,
                                                       params_A,
                                                       params_B,
                                                       iterations=iter,
                                                       threshold=0.,
                                                       )

            ratio_inter_intra = _get_distance_metric(model_node_pos=model_node_pos,
                                                     grid_key=grid_key,
                                                     other_keys=other_keys,
                                                     )
            data_dists.append(ratio_inter_intra)

        # Save results
        with open(f'{export_folder}/fig_4_E_data_dists_iter_{iter}.pkl', 'wb') as f:
            pickle.dump(data_dists, f)

        ####### SHUFFLING #########################################################################################
        shuffled_dists = [] 
        for i in trange(shuffle_iter, desc='Shuffling'):
            # SHUFFLE WEIGHTS
            weights_ = weights.copy()
            weights_ = np.random.permutation(weights_)
            # Helper for getting graph / simulation
            model_node_pos = _return_simulation_points(weights_,
                                                       params_A,
                                                       params_B,
                                                       iterations=iter,
                                                       threshold=0.,
                                                       )

            ratio_inter_intra = _get_distance_metric(model_node_pos=model_node_pos,
                                                     grid_key=grid_key,
                                                     other_keys=other_keys,
                                                     )
            shuffled_dists.append(ratio_inter_intra)

        # Save results
        with open(f'{export_folder}/fig_4_E_shuff_dists_iter_{iter}.pkl', 'wb') as f:
            pickle.dump(shuffled_dists, f)
        
    

    print('Success.')