### Scores, etc. over anatomical space
import warnings # Disable np.nanmean Runtime warning
import datajoint as dj

#### LOAD DATABASE ##########################################
from .dj_conn import *
imhotte = dj.schema(horst_imaging_db)


@imhotte 
class FilteredSessionsParams(dj.Lookup):
    definition = """
    # Restrictions saved as restriction_dict_session and uniquely identified by parameter_hash_session
    param_hash_session        : char(16)          # 16 character hash of used restrictions
    ---    
    restriction_dict_session  : blob@imgstore     # Dictionary of filters used in query
    """

@imhotte 
class FilteredCellsParams(dj.Lookup):
    definition = """
    # Parameters saved as parameter_dict_cells and restriction_dict_cell and are uniquely identified by parameter_hash_cell
    param_hash_id_cell = "standard"       : enum('standard','ov_cutoff_D','ov_cutoff_H')  #
    ---    
    param_hash_cell                       : char(16)          # 16 character hash of used restrictions AND parameters dict
    unique index (param_hash_cell)
    restriction_dict_cell                 : blob@imgstore     # Dictionary of filters used in query (Restrictions)
    parameter_dict_cell                   : blob@imgstore     # Dictionary of parameters used in query (Parameters)
    """

#### SESSION / CELL TABLES

@imhotte 
class FilteredSessions(dj.Manual):
    definition = """
    # Filtered sessions (manually selected by user)
    -> Tracking.proj(tracking_dataset = 'dataset_name')
    -> FilteredSessionsParams
    ---    
    entry_time_filt_sessions = CURRENT_TIMESTAMP : timestamp
    """

@imhotte 
class FilteredCells(dj.Manual):
    definition = """
    # Filtered cells (manually selected by user)
    -> FilteredSessions
    -> Cell
    -> FilteredCellsParams
    ---    
    entry_time_filt_cells = CURRENT_TIMESTAMP : timestamp
    """
