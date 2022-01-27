from collections import OrderedDict

#### GENERAL FILE FOR KEEPING TRACK OF CONSTANTS USED THROUGHOUT TOPOGRAPHY PROJECT

# Session level
OCC_SECONDS_CUTOFF     = 900  # Seconds, filter for "occupancy_time" in Occupancy()
EXPLOR_STD_CUTOFF      = 6    # Cutoff for "explor_std" in Occupancy() (smaller than this ...)
EXPLOR_RATIO_CUTOFF    = .7   # Cutoff of "explor_ratio" in Occupancy() (bigger than this ...)

# Single cell level
SNR_CUTOFF             = 3.5  # Signal-to-noise-ratio cutoff: "snr_df_f" in SNR()
COMPRESSION_CUTOFF     = .1   # Compression ratio cutoff. "compression" in RoisCorr()

# Cell specific 
TRACKINGPARAMS_ID      = 'A'
S_T_PARAMS_ID          = 'A'
MAP_PARAMS_ID          = 'A'
FIELD_PARAMS_ID        = 'A'
GRID_PARAMS_ID         = 'A'
BORDER_PARAMS_ID       = 'A'
ANG_PARAMS_ID          = 'A'
CHANNEL                = 'primary'
SIGNAL_TYPE            = 'spikes'
NOISE_CALC_ID          = 'A'
SPIKE_FILTER_ID        = 'A'
SHUFFLE_PARAMS_ID      = 'A'
OV_PARAMS_ID           = 'A'
OV_CUTOFF_ID           = 'A'
BV_SCORE_PARAM         = 'A'
BV_FIELD_METHOD        = 'bvs'
BV_FIELD_PARAM         = 'A'

# For binned statistics 
GLOBAL_CUTOFF_NO_CELLS         = 15 # Only calculate statistics in regions where more than n cells have been found
GLOBAL_CUTOFF_NO_STARTER_CELLS = 5  # Only calculate statistics in regions where more than n starter cells have been found

# For NN statistics 
# Reason:  GLOBAL_CUTOFF_NO_STARTER_CELLS can be too high still - there might be only a few
# functional cells remaining after filtering and I need that statistic 
NN_CUTOFF_NO_CELLS         = 15
NN_CUTOFF_NO_STARTER_CELLS = 5

MIN_NN_DIST_THRESH = 10. # [microns] Minimum distance threshold for any nearest neighbour analysis 

# Construct params dict from parameters above ( DO NOT TOUCH )
cell_parameter_dict = {
        'trackingparams_id' :    TRACKINGPARAMS_ID,
        's_t_params_id' :        S_T_PARAMS_ID,
        'map_params_id' :        MAP_PARAMS_ID,
        'field_params_id':        FIELD_PARAMS_ID,
        'grid_params_id' :       GRID_PARAMS_ID,
        'border_params_id':      BORDER_PARAMS_ID,
        'ang_params_id':         ANG_PARAMS_ID,
        'channel' :              CHANNEL,
        'signal_type':           SIGNAL_TYPE,
        'noise_calc_id':         NOISE_CALC_ID,
        'spike_filter_id':        SPIKE_FILTER_ID,
        'shuffle_params_id':      SHUFFLE_PARAMS_ID,
        'ov_params_id':          OV_PARAMS_ID,
        'ov_cutoff_id':          OV_CUTOFF_ID,
        'bvscore_params_id':     BV_SCORE_PARAM,
        'bvfield_params_id':      BV_FIELD_PARAM,   
        'bv_field_dect_method':   BV_FIELD_METHOD,
                     }    
cell_parameter_dict = OrderedDict(sorted(cell_parameter_dict.items())) # Return ordered by item