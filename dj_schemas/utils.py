### Small helpers for datajoint related functions
import sys
import copy
from datetime import datetime
from hashlib import sha512
from collections import OrderedDict

from tqdm.auto import tqdm 
import numpy as np

import datajoint as dj 

#### LOAD DATABASE #########################################
from .dj_conn import *
#### LOAD CONSTANTS ########################################
from .constants import *


def make_multi_session_object_dict(key):
    ''' 
    Loop over session / object configurations for one metasession and 
    collect session KEYs in dictionary for 
    - Base session 
    - Object session 1
    - Object session 2 (or object session 1 second object)


    This is used in the code for object vector cell identification (table OVC()) 
    and in the plotting class to identify sessions. 
    
    '''
    
    # Take care of object sessions
    # Careful, ArenaObjectPos must be fully populated, all objects must be registered
    session_dict = (Session.SessionType * \
                    ArenaObjectPos & \
                    'metasession_name = "{}"'.format(key['metasession_name'])).fetch('KEY', order_by='session_order ASC')
    
    if len(session_dict) > 2:
        raise NotImplementedError(f'{len(session_dict)} object sessions / objects found, but > 2 not supported')
    
    # Loop over sessions / objects and extract session dictionaries
    sessions = {}
    for sess in session_dict:
        if 'object1' not in sessions.keys():
            sessions['object1'] = sess
        else:
            sessions['object2'] = sess
        
    # Add base session
    sessions['base'] = (Session.SessionType & \
                        'metasession_name = "{}"'.format(key['metasession_name']) & \
                        'sessiontype = "Open Field"').fetch1('KEY')
        
    return sessions


def session_title_string(session_name):
    # Create a string for display in the title of napari window or elsewhere
    session_hash, animal_name, timestamp = (Session & f'session_name = "{session_name}"').fetch1('session_name','animal_name','timestamp')
    timestamp = datetime.strftime(timestamp, '%d.%m.%Y')
    title_string = 'Session {} | Animal {} | {}'.format(session_hash, animal_name, timestamp)
    return title_string 


##### GET FILTERED SESSIONS / CELLS ########################################################################################################################

def get_filtered_sessions(animal_name, session_type='Open Field', filter_occupancy=True, verbose=True):
    if verbose:
        print(f'{animal_name} \nSession type: {session_type}')
        if filter_occupancy:
            print(f'Cutoff occupancy [seconds]:         {OCC_SECONDS_CUTOFF}')
            print(f'Explor std cutoff:                < {EXPLOR_STD_CUTOFF}')
            print(f'Eplor ratio cutoff:               > {EXPLOR_RATIO_CUTOFF}\n')

    
    if isinstance(animal_name, str):
        animal_name = [animal_name]

    if session_type is not None:
        sessions = Session * Session.SessionType.proj() & 'sessiontype = "{}"'.format(session_type) & [f'animal_name = {animal}' for animal in animal_name]
    else:
        sessions = Session & [f'animal_name = {animal}' for animal in animal_name]

    if verbose: 
        print('{:>2} unfiltered sessions of type "{}"'.format(len(sessions), session_type))

    if filter_occupancy:
        filtered_sessions = (Session * Session.SessionType * Occupancy.proj('explor_ratio','explor_std','occupancy_time') \
                            & 'occupancy_time > {}'.format(OCC_SECONDS_CUTOFF)\
                            & 'explor_std < {}'.format(EXPLOR_STD_CUTOFF) & 'explor_ratio > {}'.format(EXPLOR_RATIO_CUTOFF)\
                            & sessions.proj()) & cell_parameter_dict
    else:
        # Supplement Tracking table since it is needed to enter filtered session entry 
        filtered_sessions = (Session * Session.SessionType * Tracking \
                            & sessions.proj()) & cell_parameter_dict

    if verbose:
        print('{:>2} filtered sessions'.format(len(filtered_sessions)))
    if (session_type is None) and verbose:
        print('Session types: {}'.format(set(filtered_sessions.fetch('sessiontype'))))

    return filtered_sessions 
        
def get_filtered_cells(filtered_sessions, corr=None, verbose=True):
    '''
    Based on a session join object, return filtered cells. 
    Cell filtering as in CONSTANTS 

    Parameters
    ----------
    filtered_sessions : datajoint join object
    corr : bool, defaults to None
           - Use "corrected" tables (i.e. RoisCorr)?
    verbose: show output 

    Returns
    -------
    filtered_cells : datajoint join object 

    '''
    animals = filtered_sessions.fetch('animal_name')

    # Find out whether to use corr or non-corr tables (FOV)
    if not len(RoisCorr.proj(signal_dataset='dataset_name') & filtered_sessions):
        if (corr is not None) and (corr): 
            raise IndexError('RoisCorr table chosen, but no cells were found')
        else:
            corr = False
    else:
        if corr is None: 
            corr = True
        else:
            corr = corr 
            
    if corr: # Search in RoisCorr table, considering the compression ratio 
        if verbose: print('Using RoisCorr and compression_ratio cutoff')
        filtered_cells = (SNR.proj(signal_dataset='dataset_name', snr='snr_df_f') * \
                            RoisCorr.proj(compression='compress_ratio',signal_dataset='dataset_name') \
                            & 'snr > {}'.format(SNR_CUTOFF)\
                            & 'compression > {}'.format(COMPRESSION_CUTOFF)\
                            & filtered_sessions \
                            & cell_parameter_dict)
    else: 
        if verbose: print('Not using RoisCorr (skipping compression_ratio cutoff)')
        filtered_cells = (SNR.proj(signal_dataset='dataset_name', snr='snr_df_f') \
                            & 'snr > {}'.format(SNR_CUTOFF)\
                            & filtered_sessions \
                            & cell_parameter_dict)
    if verbose:
        print(f'{len(filtered_cells)} cells found for {set(animals)} over {len(filtered_sessions)} sessions')

    return filtered_cells


################################################################################################################
#### FOR anatomical_distribution.py and others: Parameter hashing:

def hash_from_string(string,size=16):
    ''' Return sha512 hash from string '''
    keybytes = string.encode('utf-8')
    return sha512(keybytes).hexdigest()[:size]

def get_restriction_dict(filtered):
    ''' 
    Based on the .restriction attribute in DJ join object ('filtered'), return restriction dictionary
    This is to build a parameter hash that uniquely identifies a set of sessions / cells 
    saved in datajoint. 
    
    CAVE: The polarity is lost (ok because usually obvious)
    '''
    restriction_dict = {}
    for param in filtered.restriction:
        if isinstance(param, str):
            if '<' in param:
                param = param.split('<')
                restriction_dict[param[0].strip()] = param[-1].strip()
            elif '>' in param:
                param = param.split('>')
                restriction_dict[param[0].strip()] = param[-1].strip()
            elif '=' in param:
                param = param.split('=')
                restriction_dict[param[0].strip()] = param[-1].strip()
    return OrderedDict(sorted(restriction_dict.items()))

def get_param_hash(restriction_dict):
    ''' 
    Create hash from dictionary 
    General function. There are exception, those are handled in 
    separate functions (e.g. get_score_params_hash() below)
    '''
    restriction_dict = OrderedDict(sorted(restriction_dict.items()))
    parameters = ''.join(['' + str(v) for v in restriction_dict.values()])
    parameter_hash = hash_from_string(parameters)
    return parameter_hash

def get_score_params_hash(scoreparams_dict):
    ''' Create hash from ScoreParams dictionary '''
    scoreparams_dict = OrderedDict(sorted(scoreparams_dict.items()))
    # Skip 'scoretables' column/key since it is assumed to be too variable ... 
    parameters = ''.join(['' + str(v) for k,v in scoreparams_dict.items() if k != 'scoretables'])
    parameter_hash = hash_from_string(parameters)
    return parameter_hash

#### MAKE DJ ENTRIES FOR FILTERED SESSIONS AND CELLS 

def insert_filtered_sessions_params(filtered_sessions):
    ''' 
    Insert 1 into FilteredSessionsParams()
    Create sessions restrictions dictionary and hash and insert into datajoint 

    Parameters
    ----------
    filtered_sessions : datajoint join object
                        Filtered sessions. Make sure the whole join and not some proj()
                        is fed in 

    Returns
    -------
    session_param_dict : dict
                        Dictionary containing
                        - Parameter hash for filtered sessions, created from 
                        - Restriction dictionary
    ''' 

    restriction_dict_sessions = get_restriction_dict(filtered_sessions)
    param_hash_sessions       = get_param_hash(restriction_dict_sessions)

    session_param_dict = {
        'param_hash_session'       : param_hash_sessions,
        'restriction_dict_session' : restriction_dict_sessions
                }

    if len(FilteredSessionsParams & f'param_hash_session = "{param_hash_sessions}"'):
        print(f'FilteredSessionParams entry with hash "{param_hash_sessions}" already exists. Skipping.')
    else:
        FilteredSessionsParams.insert1(session_param_dict, ignore_extra_fields=True, skip_duplicates=True)
    return session_param_dict

def insert_filtered_cells_params(filtered_cells):
    ''' 
    Insert 1 into FilteredCellsParams()
    Create cells restrictions dictionary and hash and insert into datajoint 
    
    Parameters
    ----------
    filtered_cells    : datajoint join object
                        Filtered cells. Make sure the whole join and not some proj()
                        is fed in 

    Returns
    -------
    cell_param_dict   : dict
                        Dictionary containing
                        - Parameter hash for filtered cells, created from 
                        - Restriction dictionary and 
                        - Cell parameter dictionary

    ''' 

    restriction_dict_cells = get_restriction_dict(filtered_cells)
    ### Below "cell_parameter_dict" is loaded from .constants.py
    param_hash_cells       = get_param_hash({**restriction_dict_cells,**cell_parameter_dict})

    cell_param_dict = {
        'param_hash_id_cell'    : 'standard',
        'param_hash_cell'       : param_hash_cells,
        'restriction_dict_cell' : restriction_dict_cells,
        'parameter_dict_cell'   : cell_parameter_dict
                }
    if len(FilteredCellsParams & f'param_hash_cell = "{param_hash_cells}"'):
        print(f'FilteredCellsParams entry with hash "{param_hash_cells}" already exists. Skipping.')
    else:
        FilteredCellsParams.insert1(cell_param_dict, ignore_extra_fields=True, skip_duplicates=True)
    return cell_param_dict


def insert_filtered_sessions_cells(filtered_sessions, filtered_cells):
    ''' 
    Loop over filtered sessions and cells and create params entries and actual table entries 
    in 
    FilteredSessions()
    FilteredCells()

    This function partners with the two functions returning filtered sessions + cells:
    * get_filtered_sessions()
    * get_filtered_cells()
    => These provide the join objects that are fed into this function. 
    
    Parameters
    ----------
    filtered_sessions : datajoint join object
                        Filtered sessions. Make sure the whole join and not some proj()
                        is fed in 
    filtered_cells    : datajoint join object
                        Filtered cells. Make sure the whole join and not some proj()
                        is fed in

    ''' 
    session_param_dict = insert_filtered_sessions_params(filtered_sessions)
    cell_param_dict    = insert_filtered_cells_params(filtered_cells)

    # Insert sessions
    print(f'Processing {len(filtered_sessions)} sessions:')
    # Loop over sessions ... 
    for no , session in enumerate(filtered_sessions.proj(tracking_dataset='dataset_name')):
        sys.stdout.write(f'{no+1:<3} | Inserting cells for session_name {session["session_name"]}')
        session_dict = {
            ** session,
            'param_hash_session'       : session_param_dict['param_hash_session']
                    }
        FilteredSessions.insert1(session_dict, ignore_extra_fields=True, skip_duplicates=True)

        # Filter cells
        cells_in_session = (filtered_cells.proj(dataset_name='signal_dataset'              
                                                ) & session)
        
        # Check previous entries 
        insert_cells = True
        previous_entries = FilteredCells & session_dict & 'param_hash_id_cell = "{}"'.format(cell_param_dict['param_hash_id_cell'])
        if len(previous_entries):
            if (set(previous_entries.fetch('cell_id')) != set(cells_in_session.fetch('cell_id'))):
                sys.stdout.write(' | Entries do not match. Deleting.\n')
                previous_entries.delete()
            else: 
                sys.stdout.write(' | Entries match. Skipping.\n')
                insert_cells = False

        if insert_cells: # If none found or cells were deleted, continue to enter cells...
            for cell in cells_in_session:
                cell_dict = {
                    ** session_dict,
                    ** cell,
                    'param_hash_id_cell'  : cell_param_dict['param_hash_id_cell']
                            }
                FilteredCells.insert1(cell_dict,ignore_extra_fields=True, skip_duplicates=True)
            sys.stdout.write(f' | Inserted {len(cells_in_session)} cells.\n')

    print('Success.')
    return session_param_dict['param_hash_session'], cell_param_dict['param_hash_id_cell']

#### ADDITIONAL HELPERS 

def get_signal_indices(signal, reference):
    '''
    Given two arrays: 
    - signal
    - reference
    
    where signal is a subset of reference,
    return the indices for every entry in signal in reference.
    
    This is useful if one wants to look up for example the 
    exact times when a signal (filtered in SignalTracking)
    occurred. 
    
    e.g. 
    get_signal_indices(y_pos_signal, y_pos), where 'y_pos_signal'
    is the filtered y_pos retrieved from SignalTracking and 'y_pos' 
    is the original y_pos that was used in SignalTracking, 
    returns the indices of y_pos_signal in y_pos.    
    
    See: 
    https://stackoverflow.com/questions/8251541/
    numpy-for-every-element-in-one-array-find-the-index-in-another-array
    
    Parameter
    ---------
    signal     : np.array
                 Array of data points to look up 
                 in "reference"
    reference  : np.array
                 Reference array that "signal"
                 should be compared with
                 
    Returns
    -------
    indices    : np.array
                 Indices of "signal" in "reference"
    
    '''
    x = reference
    y = signal

    index = np.argsort(x)
    sorted_x = x[index]
    sorted_index = np.searchsorted(sorted_x, y)

    yindex = np.take(index, sorted_index, mode="clip")
    mask = x[yindex] != y

    result = np.ma.array(yindex, mask=mask)
    indices = result.data

    assert len(np.unique(indices)) == len(indices) == len(signal)

    return indices