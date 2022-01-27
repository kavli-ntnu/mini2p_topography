### Figure 4 E and S5D: NN graph and supplementary - Obenhaus et al. 

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

from Figure_4_D_Figure_S5_B_C import _get_nn_sessions

from itertools import combinations
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle
from sklearn.metrics import pairwise_distances
import networkx as nx
from dj_plotter.helpers.plotting_helpers import make_linear_colormap


##### LOAD SCHEMA COMPONENTS #############################################################
from dj_schemas.dj_conn import * 

##### EXPORT LOCATION ####################################################################
figure_dir = 'YOUR_EXPORT_DIRECTORY/'
figure_dir = Path(figure_dir)



def get_nn_weights(animals,
                   params,
                   param_hash_id_cell='standard',
                   param_hash_session='cf83e1357eefb8bd',
                   norm_col='shuffall',
                   region='MEC',
                   cutoff_n_cells=5,
                   ):
    '''
    NEAREST NEIGHBOUR (NN) GRAPH BUILDING
    Loop over parameter combinations and extract normalized NN metrics and build network

    Returns 
    -------
    nn_graph : NetworkX Multigraph instance
    weights  : Weights between nodes in the graph 

    '''
    ###### COLLECT WEIGHTS ##########################################################################

    weights  = []
    params_A = []
    params_B = []
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
            
        ### LOOP OVER SESSIONS AND ADD EDGES 
        for session in sessions:          
            weight = (1 / np.mean([session['norm_nns_ab'], session['norm_nns_ba']])) - 1
            weights.append(weight)
            params_A.append(session['pairwise_dist_param_A'])
            params_B.append(session['pairwise_dist_param_B'])
    
    return weights, params_A, params_B



##### GRAPH RELATED ###############################################################################################



def _make_nn_graph(weights, params_A, params_B):
    '''
    Taking lists of weights and corresponding nodes (param A vs. B),
    build a networkx MultiGraph()

    Returns
    -------
    graph   : nx.MultiGraph()
              Weight attribute is called 'weight'

    '''

    ###### BUILD GRAPH ##############################################################################
    graph = nx.MultiGraph() # Implementation of undirected graph that can hold multiple nodes
    for w, a, b in zip(weights, params_A, params_B):
        graph.add_edge(a, b, weight = w)

    return graph


def _generate_initialpos(graph):
    '''
    Randomly initiate graph node positions in a square centered on zero
    '''

    total = len(graph.nodes)
    assert total == 4, f'{total} graph nodes found. This is not supported currently (needs to be 4 exactly)'

    # Initialize 4 points 
    rnd = np.random.random
    start_coords = np.array([(-1*rnd(),-1*rnd()),
                             (-1*rnd(),rnd()),
                             (rnd(),-1*rnd()),
                             (rnd(),rnd())]
                             )
    initial_pos = {}
    for node_n, coord in zip(list(graph.nodes), np.random.permutation(start_coords)):
        initial_pos[node_n] = coord.tolist()
        
    return initial_pos


def _spring_loaded_model(graph, iterations, threshold):
    '''  
    Simulate "spring loaded model" 
    https://networkx.org/documentation/stable/reference/generated/networkx.drawing.layout.spring_layout.html
    
    Initialize positions around zero. The initialization itself is random (permutation of possible positions)

    Parameters
    ----------
    graph     : networkx graph object, where weight attribute is called 'weight'
    iterations: int
                Number of iterations to run for spring loaded model simulation
    threshold : float
                See doc of layout.spring_layout 


    Returns
    -------
    model_nodel_pos : networkx dictionary of keys:positions after simulation has stopped
    '''

    # Randomly initialize positions in a square centered on zero
    initialpos = _generate_initialpos(graph)

    model_node_pos = nx.spring_layout(graph, 
                                      weight='weight', 
                                      pos=initialpos, 
                                      iterations=iterations, 
                                      threshold=threshold,
                                      center=[0,0])

    return model_node_pos



def plot_graph(graph, weights, model_node_pos, iterations, ratio): 
    '''
    Plot outcome of spring loaded model simulation

    Create colormap over all weights, 
    then loop over edges and draw node positions / vertices

    '''

    ####### SIMULATE + PLOT #########################################################################

    cmap = 'coolwarm' # #RdBu_r
    colormap_boundaries = np.array([-np.max(weights), np.max(weights)])

    weight_colors = make_linear_colormap(weights, 
                                         cmap=cmap, 
                                         percentile=100.,
                                         reference_numbers=colormap_boundaries)

    sns.set(style='white', font_scale=1.6)
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.size'] = 10
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['xtick.bottom'] = True
    plt.rcParams['ytick.left'] = True

    figure = plt.figure(constrained_layout=False, figsize=(6,7))
    gs = figure.add_gridspec(nrows=6, ncols=1, hspace=1.5)

    ax_graph = figure.add_subplot(gs[:-2, :])   
 
    # Loop over graph edges and draw them, color coding by weight
    for w, e, color in zip(weights, graph.edges, weight_colors):
        ax_graph.annotate("",
                    xy = model_node_pos[e[0]], xycoords='data',
                    xytext = model_node_pos[e[1]], textcoords='data',
                    arrowprops = dict(arrowstyle="-", color=color,
                                      alpha=.1*(np.abs(w)+1),
                                      lw=(w+1)*10,
                                      connectionstyle="arc3,rad=rrr".replace('rrr',str(0.0020*e[2])
                                      ),
                                      ),
                    )
    # Add nodes and labels
    nx.draw_networkx_nodes(graph,  model_node_pos, node_color='#ccc', node_size=2000, alpha=1,)
    nx.draw_networkx_labels(graph, model_node_pos)

    # ... take care of colorbar
    norm = mpl.colors.Normalize(vmin=colormap_boundaries[0], vmax=colormap_boundaries[1])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    divider = make_axes_locatable(ax_graph)
    cax = divider.append_axes("right", size="3%", pad=0.1)
    cbar = plt.colorbar(sm, 
                        cax=cax, 
                        ticks=np.linspace(colormap_boundaries[0], colormap_boundaries[1], 5), 
                        )
    cbar.outline.set_visible(False)
    

    xlim = ax_graph.get_xlim()
    x_stretch = .1* np.abs([xlim[1]-xlim[0]])
    ylim = ax_graph.get_ylim()
    y_stretch = .1* np.abs([ylim[1]-ylim[0]])

    ax_graph.set_xlim(xlim[0]-x_stretch,xlim[1]+x_stretch)
    ax_graph.set_ylim(ylim[0]-y_stretch,ylim[1]+y_stretch)

    # Display ratio 
    ax_graph.text(.8*xlim[0],1.1*ylim[0], f'Ratio: {ratio:.2f}')


    ##### Load in shuffled results and data from previous run of simulations
    # This is the output of fig_4_D_shuffling.py 

    export_folder = 'data_exports/NN_graph_shuffling/'
    print(f'Opening files for {iterations} iterations (=simulation steps)')
    with open(f'{export_folder}/fig_4_E_data_dists_iter_{iterations}.pkl', 'rb') as f:
        data_dists  = pickle.load(f)
    with open(f'{export_folder}/fig_4_E_shuff_dists_iter_{iterations}.pkl', 'rb') as f:
        shuff_dists = pickle.load(f)

    # ... add new axis
    ax_hist_data_shuff = figure.add_subplot(gs[-2:, :])   

    all_results = data_dists + shuff_dists
    median_shuff = np.median(shuff_dists)
    median_data  = np.median(data_dists)
    min = .8  * np.min(all_results)
    max = 1.2 * np.max(all_results)
    bins = np.linspace(min, max, 200)
    # Histogram
    ax_hist_data_shuff.hist(shuff_dists, bins=bins, lw=0, color='#333', alpha=.6, label='Shuffled')
    ax_hist_data_shuff.hist(data_dists, bins=bins, lw=0, color='k', alpha=.8, label='Data')
    ax_hist_data_shuff.axvline(x=median_shuff, ls=':', color='k',alpha=.9, label=f'Med. shuff: {median_shuff:.1f}')
    ax_hist_data_shuff.axvline(x=median_data, color='k', label=f'Med. data: {median_data:.1f}')
    ax_hist_data_shuff.legend()
    ax_hist_data_shuff.set_xlim(0, 15.)
    ax_hist_data_shuff.set_xlabel('Dist. ratio grid (inter) vs. rest (intra)')

    # # TEMPORARY: ADD VERTICAL LINES 
    # ax_hist_data_shuff.axvline(x=1.68, ls=':', color='r', alpha=.9)
    # ax_hist_data_shuff.axvline(x=4.78, ls=':', color='r', alpha=.9)
    # ax_hist_data_shuff.axvline(x=8.42, ls=':', color='r', alpha=.9)

    sns.despine(left=True,bottom=True)

    # # Draw weight histogram
    # figure = plt.figure(figsize=(5,2))
    # ax = figure.add_subplot(111)
    # ax.hist(weights, lw=0, color='#333', bins=20);
    # ax.set_xlabel('Weights');
    plt.savefig(
            f'{figure_dir}/NN graph ratio {int(np.round(ratio))}.pdf', 
            dpi=300, 
            facecolor='w', 
            edgecolor='w',
            orientation='portrait',
            format=None,
            transparent=False, 
            bbox_inches='tight', 
            pad_inches=0,
            )

    #plt.show()

    return



def _get_distance_metric(model_node_pos,
                         grid_key='A',
                         other_keys=('C','E','F'),
                         ):
    '''
    Get distance metric

    Parameters
    ----------
    model_node_pos : Positions after simulation


    '''

    # Work with position output from simulation
    pairw_ov_border_hd  = pairwise_distances([model_node_pos[other_keys[0]],
                                              model_node_pos[other_keys[1]],
                                              model_node_pos[other_keys[2]]],
                                              metric='euclidean')
    pairw_ov_border_hd = np.mean(pairw_ov_border_hd[np.triu_indices_from(pairw_ov_border_hd, k=1)])

    # Center point of 
    center_ov_border_hd = np.mean([model_node_pos[other_keys[0]],
                                   model_node_pos[other_keys[1]],
                                   model_node_pos[other_keys[2]]])

    dist_grid_to_rest   = np.linalg.norm(center_ov_border_hd - model_node_pos[grid_key])
    # "Inter to intra"
    ratio_inter_intra   = dist_grid_to_rest / pairw_ov_border_hd

    ##############################################################################################################
    #print(f'Distance Grid vs. Rest: {ratio_inter_intra}\n\n')

    return ratio_inter_intra









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

    ### SHUFFLE WEIGHTS? 
    shuffle = False

    weights, params_A, params_B = get_nn_weights(animals,
                                                 params,
                                                 param_hash_id_cell='standard',
                                                 param_hash_session='cf83e1357eefb8bd',
                                                 norm_col=nn_norm_column,
                                                 region=region,
                                                 cutoff_n_cells=cutoff_n_cells,
                                                 )
    
    if shuffle: 
        print('YES ! SHUFFLING WEIGHTS')
        weights = np.random.permutation(weights)
    else:
        print('NOT ! SHUFFLING WEIGHTS')
        weights = np.array(weights)


    ####### BUILD GRAPH + SIMULATE  #############################################################################
    iterations = 1000

    nn_graph = _make_nn_graph(weights, params_A, params_B)
    model_node_pos = _spring_loaded_model(nn_graph, iterations = iterations, threshold = 0.)
    # Stats
    ratio_inter_intra = _get_distance_metric(model_node_pos,
                                             grid_key='A',
                                             other_keys=('C','E','F'),
                                             )

    #### PLOT  ##################################################################################################
    plot_graph(nn_graph, weights, model_node_pos, iterations = iterations, ratio=ratio_inter_intra)


    print(figure_dir)
    print('Success.')