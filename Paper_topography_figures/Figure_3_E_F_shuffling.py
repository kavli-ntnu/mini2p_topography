### Figure 3E/F Shuffling - Obenhaus et al. 
# Shuffles data and exports to disk

import sys, os
import os.path
import numpy as np 
import datajoint as dj
import pickle
import cmasher as cmr
from tqdm.auto import tqdm

# Make plots pretty 
from matplotlib import pyplot as plt
import seaborn as sns
sns.set(style='white')

# Prevent bug in figure export as pdf: 
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

##### IMPORTS ###########################################################################
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

##### LOAD SCHEMA COMPONENTS #############################################################
from dj_schemas.dj_conn import * 

##### GLOBAL PLOTTING OPTIONS ############################################################


def _shuffle_fractions(n_all, pop_size, n_functional, verbose=False):
    '''
    GRID vs. OV ratio shuffling

    Given a total population size (n_all) and the number of possibly functional
    cell types (n_functional), create a 50/50 distribution of cell types (OV vs. Grid), 
    and subselect n=pop_size cells (random picks after permutation)    
    
    '''
    if np.mod(n_functional, 2): 
        n_functional+=1
        if verbose: print(f'n_functional = {n_functional}')
        
    # Choose a 50 / 50 distribution 
    grid_ = np.zeros(n_all)
    grid_[:int(n_functional*.5)] = 1
    ov_ = np.zeros(n_all)
    ov_[:int(n_functional*.5)] = 1

    if verbose : print(len(ov_[ov_>0]), len(grid_[grid_>0]))
    
    # Shuffle vectors to prepare subselection ... 
    random_grid_ids = np.random.permutation(grid_)
    random_ov_ids   = np.random.permutation(ov_) 
    # Subselect n = pop_size
    chunk_grid = random_grid_ids[:pop_size]
    chunk_ov   = random_ov_ids[:pop_size]

    # Which are the grid and the ov cell ids? 
    grid_ids = np.argwhere(chunk_grid>0).squeeze()
    ov_ids   = np.argwhere(chunk_ov>0).squeeze()

    # Clean up those mixed cells ... 
    ov_pure_   = np.array(list(set(ov_ids)   - set(grid_ids)))
    grid_pure_ = np.array(list(set(grid_ids) - set(ov_ids)))
    ov_ids = ov_pure_
    grid_ids = grid_pure_
    
    fraction_ov   = len(ov_ids)/pop_size
    fraction_grid = len(grid_ids)/pop_size
    
    # Discrimination index
    disc_grid_ov = (len(grid_ids) - len(ov_ids)) / (len(grid_ids) + len(ov_ids))

    if verbose : print(fraction_ov, fraction_grid, disc_grid_ov)
    return fraction_ov, fraction_grid, disc_grid_ov



def _extract_info_grid_ov(animals, 
                          pure_only=True,
                          param_hash_id_cell='standard',
                          param_hash_session='cf83e1357eefb8bd',
                          ov_cutoff_id = 'A',
                          region='MEC',
                          cutoff_n_cells=20):


    cell_parameter_dict = (FilteredCellsParams & f'param_hash_id_cell = "{param_hash_id_cell}"').fetch1('parameter_dict_cell')
    
    # Get rid of ov_cutoff_id -> set this custom 
    ov_cutoff_id_previous = cell_parameter_dict.pop('ov_cutoff_id')
    cell_parameter_dict['ov_cutoff_id'] = ov_cutoff_id
    print(f'Previous cutoff OV cells (ov_cutoff_id): {ov_cutoff_id_previous} | Set to: {ov_cutoff_id}\n')
    
    # Brain region filter 
    assert region in ['MEC','PAS'], f'Region "{region}" not understood. Choose "MEC" or "PAS"'
    if region == 'MEC':
        print('Choosing brain region filter MEC')
        brainregion = RoisCorrBrainLoc.MEC 
    elif region == 'PAS':
        print('Choosing brain region filter PAS')
        brainregion = RoisCorrBrainLoc.PAS

    # Retrieve sessions
    all_sessions = (Session.proj('animal_name') * FilteredSessions 
                    & [f'animal_name = "{animal}"' for animal in animals]
                    & f'param_hash_session = "{param_hash_session}"'
                    )

    # All OVC entries
    ov_sessions = all_sessions & OVC.proj(session_name='base_session')

    ############################################################################################################################################################
    ### STEP 1: Collect all valid sessions (based on cell number cutoff)
    animal_names  = []
    session_names = []
    for sess in ov_sessions.proj('animal_name'):

        filtered_cells = FilteredCells.proj() & sess & brainregion & f'param_hash_id_cell = "{param_hash_id_cell}"'
        if len(filtered_cells) < cutoff_n_cells:
            # Skip this session
            continue 
        else: 
            animal_names.append(sess['animal_name'])
            session_names.append(sess['session_name'])

    # Quick stats for display: 
    animals_considered  = set(animal_names)
    sessions_considered = set(session_names)
    print(f'\n{len(animals_considered)} animals across {len(sessions_considered)} sessions\n')
    print(f'Animals: {list(animals_considered)}')
    
    ############################################################################################################################################################
    ### STEP 2: Collect all collected sessions together and extract fraction grid / OV etc. 
    filtered_cells = FilteredCells.proj() & brainregion & f'param_hash_id_cell = "{param_hash_id_cell}"' & [f'session_name = "{sess}"' for sess in session_names]
    # Total number of cells
    n_all = len(filtered_cells)
    # Number of cells per session
    pop_size = int(np.round(n_all/len(session_names)))

    # Filter further for Grid and OV
    ov_cells_all_sessions  = OVC.proj(session_name='base_session', is_ovc='is_ovc') & filtered_cells & cell_parameter_dict & 'is_ovc = 1'
    gridcells_all_sessions = GridScore * CutoffsGridscore & filtered_cells & cell_parameter_dict & 'gridscore > gridscore_95'

    if pure_only:
        # Clean up those mixed cells ... 
        ov_cells_session_   = ov_cells_all_sessions   - gridcells_all_sessions.proj()
        gridcells_session_  = gridcells_all_sessions  - ov_cells_all_sessions.proj()
        ov_cells_all_sessions  = ov_cells_session_
        gridcells_all_sessions = gridcells_session_

    # What is the total fraction of Grid and OV cells combined in the population?
    n_functional = len(gridcells_all_sessions) + len(ov_cells_all_sessions)
    fraction_functional = n_functional / n_all
    
    # Print quick results 
    print(f'n_all: {n_all}')
    print(f'fraction_functional: {fraction_functional}')
    print(f'n_functional: {n_functional}\n')
    print(f'pop_size: {pop_size}')

    return n_all, fraction_functional, n_functional, pop_size 







if __name__ == "__main__":
 
    animals = [ 
               '82913','88592', '87244', '60480',
               '97046','89841',
               '87187','88106','87245','90222',
               '94557','89622',
              ]

    print(f'Creating OV vs. Grid shuffling for {len(animals)} animal(s)')
    print(animals)
    print('\n')

    # Cutoff number of cells 
    cutoff_n_cells = 15 
    ov_cutoff_id = 'A'
    pure_only = True # Compare only "pure" cells to each other (no mixed selectivity)
    n_shuffles = 10000


    if pure_only:
        print('Comparing only "pure" cells (no mixed selectivity)')
    else:
        print('Comparing cells above threshold (can have mixed selectivity')


    n_all, fraction_functional, n_functional, pop_size = _extract_info_grid_ov(animals, 
                                                                               pure_only=pure_only,
                                                                               param_hash_id_cell='standard',
                                                                               param_hash_session='cf83e1357eefb8bd',
                                                                               ov_cutoff_id = ov_cutoff_id,
                                                                               region='MEC',
                                                                               cutoff_n_cells=cutoff_n_cells)
    
    # Build the results
    # OV ANIMALS 
    fraction_shuff_ov_ovanimal = []
    fraction_shuff_grid_ovanimal = []
    disc_shuff_ovanimal = []

    # GRID ANIMALS 
    fraction_shuff_ov_gridanimal = []
    fraction_shuff_grid_gridanimal = []
    disc_shuff_gridanimal = []

    fraction_shuff_grid = [] 
    for i in tqdm(range(n_shuffles), desc='Creating Grid/OV animals'):
        fraction_ov, fraction_grid, disc_grid_ov = _shuffle_fractions(n_all, pop_size, n_functional)
        
        if (fraction_ov>fraction_grid):
            # OV animal 
            fraction_shuff_ov_ovanimal.append(fraction_ov)
            fraction_shuff_grid_ovanimal.append(fraction_grid)
            disc_shuff_ovanimal.append(disc_grid_ov)
        elif (fraction_ov<fraction_grid):
            # Grid animal 
            fraction_shuff_ov_gridanimal.append(fraction_ov)
            fraction_shuff_grid_gridanimal.append(fraction_grid)
            disc_shuff_gridanimal.append(disc_grid_ov)

    # Export shuffled data 
    dict_export = 'data_exports/shuffled_fig_3_F'

    # OV ANIMALS 
    shuffled_ov_animal_res = {
        'frac_ov_cells'      : fraction_shuff_ov_ovanimal,
        'frac_grid_cells'    : fraction_shuff_grid_ovanimal,
        'disc_grid_ov'       : disc_shuff_ovanimal,
    }
    # GRID ANIMALS 
    shuffled_grid_animal_res = {
        'frac_ov_cells'      : fraction_shuff_ov_gridanimal,
        'frac_grid_cells'    : fraction_shuff_grid_gridanimal,
        'disc_grid_ov'       : disc_shuff_gridanimal,
    }

    with open(dict_export + f'/shuffled_ov_animals_frac_pure={pure_only}ov_cutoff_id={ov_cutoff_id}.pkl', 'wb') as f:
        pickle.dump(shuffled_ov_animal_res, f, pickle.HIGHEST_PROTOCOL)
    with open(dict_export + f'/shuffled_grid_animals_frac_pure={pure_only}ov_cutoff_id={ov_cutoff_id}.pkl', 'wb') as f:
        pickle.dump(shuffled_grid_animal_res, f, pickle.HIGHEST_PROTOCOL)


    print('Success.')