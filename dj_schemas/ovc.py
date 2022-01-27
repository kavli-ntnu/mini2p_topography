# Started 01.11.2020
# New methods for object vector cell / map calculations

import sys, os
import datajoint as dj

import numpy as np
import warnings # Disable np.nanmean Runtime warning

from scipy.ndimage import gaussian_filter
from helpers_topography.utils import corr2
from physt import special_histograms
from astropy.convolution import convolve, Gaussian2DKernel

from pointpats import PointPattern
from pointpats.centrography import weighted_mean_center

from scipy.stats import circmean, circstd

import math
from tqdm.auto import tqdm
import copy
import warnings # Disable np.nanmean Runtime warning

#### LOAD DATABASE #########################################
from .dj_conn import *
imhotte = dj.schema(horst_imaging_db)

from .utils import make_multi_session_object_dict, get_filtered_cells


@imhotte
class OVParams(dj.Lookup):
    definition = """
    # Object centered map parameters
    ov_params_id                   : char(1)        # Parameter set ID, starting with A
    ---
    bin_size_dist_ov               : float          # Bin size for distance binning in mm
    bins_angular_ov                : int            # Number of bins in 360 degrees
    sigma_time_ov                  : float          # 2D gaussian smoothing of occupancy
    sigma_signal_ov                : float          # 2D guassian smoothing of binned signal
    """
    contents = [['A', 25., 72, 2., 2.]]


@imhotte
class OVOccupancy(dj.Computed):
    definition = """
    # Object centered occupancy
    -> Tracking.OpenField
    -> SignalTrackingParams
    -> OVParams
    -> ArenaObjectPos
    ---
    occupancy_ov                   : blob@imgstore        # Smoothed 2D occupancy map [seconds], x: angles, y: distance
    mask_occ_ov                    : blob@imgstore        # Mask (where time = 0), x: angles, y: distance
    occupancy_raw_ov               : blob@imgstore        # Raw, non-smoothed 2D occupancy map,  x: angles, y: distance
    explor_ratio_ov                : double               # Exploration ratio (visited bins over all bins)
    explor_std_ov                  : double               # Exploration standard deviation (of visited bins)
    radial_edges_ov                : blob@imgstore        # Histogram edges in y (distance)
    angular_edges_ov               : blob@imgstore        # Histogram edges in x (angles)
    occupancy_time_ov              : double               # Time in seconds in occupancy
    fraction_occupancy_ov          : double               # Fraction of time in occupancy map

    """

    @property
    def key_source(self):
        return Tracking.OpenField * SignalTrackingParams * OVParams * ArenaObjectPos & Session.Apparatus

    def make(self, key):
        '''
        Get object centered occupancy maps as in
        Høydal, Ø. A. et al., doi:10.1038/s41586-019-1077-7

        These are used further in OVMap() to obtain object-centered, time-normalized firing rate maps.
        The basic layout follows the make routine of (2D) Occupancy()

        '''

        apparatus = (Apparatus.Geometry * Session.Apparatus & key).fetch1()
        params_ov = (OVParams & key).fetch1()
        params_st =  (SignalTrackingParams & key).fetch1()

        if apparatus['arena_geometry'] == 'square':
            tracking_data =  (Tracking * Tracking.OpenField & key).fetch1()
            object_data   = (ArenaObjectPos & key).fetch1()

            tracking_frametime = 1/tracking_data['sample_rate']

            speed_filter = [(tracking_data['speed'] > params_st['speed_cutoff_low'])
                & (tracking_data['speed'] < params_st['speed_cutoff_high'])][0]

            # Max possible distance in arena (one corner to the diagonally opposite one)
            max_pos_dist = np.sqrt(apparatus['arena_x_dim'] ** 2 \
                                   + apparatus['arena_y_dim'] ** 2)

            # Initialize distance and angular bins
            radial_bins = np.arange(0, max_pos_dist, params_ov['bin_size_dist_ov']) # Bin edges [mm]
            phi_bins    = params_ov['bins_angular_ov'] # Number (!) of bins in 360 degrees


            # Calculate object centered angle and distance for every tracking point
            angles, dists = get_dist_angle_obj(
                        tracking_data['x_pos'][speed_filter],
                        tracking_data['y_pos'][speed_filter],
                        object_data['obj_x_coord_calib'],
                        object_data['obj_y_coord_calib']
                                   )
            # Get polar (radial) histogram
            # https://physt.readthedocs.io/en/latest/special_histograms.html#Polar-histogram
            radial_occupancy = special_histograms.polar_histogram(
                                                            dists,
                                                            angles,
                                                            phi_bins    = phi_bins,
                                                            radial_bins = radial_bins,
                                                            transformed = True
                                                                )
            # Unpack in "numpy" style (returns structure similar to numpy.histogram)
            histogram, [yedges, xedges] = radial_occupancy.numpy_like
            # Convert to seconds
            histogram = histogram.astype(float) * tracking_frametime

            time_occupancy = np.sum(histogram)
            fraction_time_occupancy = time_occupancy / (len(tracking_data['x_pos']) * tracking_frametime)

            histogram_original = histogram.copy()
            # Smooth occupancy
            # ‘wrap’ (a b c d | a b c d | a b c d)
            # The input is extended by wrapping around to the opposite edge.
            histogram = gaussian_filter(histogram, sigma=params_ov['sigma_time_ov'], mode='wrap')
            masked_histogram = np.ma.masked_where(histogram_original < 0.001, histogram)  # arbitrary threshold of 1 ms

            # Quick statistics on coverage:
            explor_ratio = len(histogram_original[histogram_original > 0]) / len(histogram_original.ravel())
            explor_std   = np.std(histogram_original[histogram_original > 0].ravel())

            # Build entry dict
            entry_dict = {
                'occupancy_ov'     : histogram,
                'mask_occ_ov'      : masked_histogram.mask,
                'occupancy_raw_ov' : masked_histogram.data,
                'explor_ratio_ov'  : explor_ratio,
                'explor_std_ov'    : explor_std,
                'radial_edges_ov'  : yedges,
                'angular_edges_ov' : xedges,
                'occupancy_time_ov': time_occupancy,
                'fraction_occupancy_ov': fraction_time_occupancy
                        }

            self.insert1({**key, **entry_dict})
        else:
            raise NotImplementedError('Geometry "{}" not implemented yet'.format(apparatus['arena_geometry']))


@imhotte
class OVMap(dj.Computed):
    definition = """
    # Object centered ratemap (vector map)
    -> SignalTracking.proj()
    -> OVOccupancy.proj(tracking_dataset='dataset_name')
    ---
    ovmap                          : blob@imgstore        # Object centered ratemap ("vector map")
    ovmap_raw                      : blob@imgstore        # Unsmoothed (raw) 2D ratemap
    mask_ovmap                     : blob@imgstore        # Mask (where time = 0)
    binned_raw_ov                  : blob@imgstore        # Raw, binned signal
    bin_max_ov                     : blob@imgstore        # Bin with maximum signal (ovmap(bin_max) = max(ovmap))
    max_ov                         : double               # Maximum
    """

    @property
    def key_source(self):
        return super().key_source & 's_t_params_id = "A"'

    def make(self,key):
        occupancy_entry       = (OVOccupancy & key).fetch1()
        signaltracking_entry  = (SignalTracking & key).fetch1()
        params_ov             = (OVParams & key).fetch1()

        object_data   = (ArenaObjectPos & key).fetch1()

        occupancy = np.ma.array(occupancy_entry['occupancy_ov'], mask=occupancy_entry['mask_occ_ov'])
        angular_edges = occupancy_entry['angular_edges_ov'].copy()
        radial_edges  = occupancy_entry['radial_edges_ov'].copy()

        # Logic:
        # angular_edges = values in x = angles
        # radial_edges  = values in y = dists

        # Calculate object centered angle and distance for every signal tracking point
        angles, dists = get_dist_angle_obj(
                    signaltracking_entry['x_pos_signal'],
                    signaltracking_entry['y_pos_signal'],
                    object_data['obj_x_coord_calib'],
                    object_data['obj_y_coord_calib']
                )

        # Supplement 'angles' and 'dists' to retrieved signaltracking_entry
        signaltracking_entry['x_pos_signal'] = angles
        signaltracking_entry['y_pos_signal'] = dists

        # Need to rename a parameter for calc_ratemap() to work
        params_ov['sigma_signal'] = params_ov['sigma_signal_ov']

        # Get object vector (=rate) map
        ovmap_dict = calc_ratemap(occupancy, angular_edges, radial_edges, signaltracking_entry, params_ov, pad_mode='wrap')

        key['ovmap']          = ovmap_dict['ratemap'].data
        key['ovmap_raw']      = ovmap_dict['ratemap_raw'].data
        key['mask_ovmap']     = occupancy.mask  # Unchanged!
        key['binned_raw_ov']  = ovmap_dict['binned_raw']
        key['bin_max_ov']     = ovmap_dict['bin_max']
        key['max_ov']         = ovmap_dict['max']

        self.insert1(key)



@imhotte
class OVCFields(dj.Computed):
    '''
    Object vector cells

    PART A: Field based comparisons



    CAVEAT:
    This currently works only for a limited set of session types / configurations:
    - 1 baseline, 2 object sessions with one object each
    - 1 baseline, 1 object session with two objects

    '''

    definition = """
    # Object vector cell (OVC) field based calculations
    -> Ratemap.proj(base_session='session_name')
    -> ShuffleParams
    ---
    object1_session                 : varchar(16)    # Object session 1
    object2_session                 : varchar(16)    # Object session 2

    """

    ##### FIELD BASED ANALYSIS ###############################################################################################

    class Fields(dj.Part):
        definition = """
        # OVC calculated field statistics (all fields)
        -> master
        object1_field_id                : int            # Field ID of field in object session 1
        object2_field_id                : int            # Field ID of closest field to object1_field_id in object session 2
        ---
        dist_fields                     : double         # Euclidian distance between fields with object1_field_id and object2_field_id - object centered
        dist_to_object = NULL           : double         # Distance of field from object [average of field in object session 1 and 2]
        angle_to_object = NULL          : double         # Angle of field to object [average of field in object session 1 and 2]
        object1_field_base = NULL       : double         # (Object session 1 field mean rate in object session 1) / (Field mean rate in base session)
        object2_field_base = NULL       : double         # (Object session 2 field mean rate in object session 2) / (Field mean rate in base session)
        """

    #######################################################################################################################


    @property
    def key_source(self):
        '''
        Filter for metasessions for which exactly 2 objects were presented either concurrently or in two subsequent sessions.
        Return base sessions for those.
        '''
        object_sessions = Session.SessionType.proj() * ArenaObjectPos & 'sessiontype = "Open Field Object"'
        meta_object_sessions = MetaSession.aggr(object_sessions, n="count(*)") & 'n=2' # Filter out everything that is not 2 objects
        base_sessions = Session.SessionType & meta_object_sessions.proj() & 'sessiontype = "Open Field"'
        return Ratemap.proj(base_session='session_name') \
                * ShuffleParams \
                & base_sessions.proj(base_session='session_name') \
                & Shuffled.proj(base_session='session_name') \
                & 'signal_type = "spikes"'\
                & 's_t_params_id = "A"'

    def make(self, key):
        # Clean up key
        # make() of OVC()
        key_ = key.copy()

        # Get rid of some keys
        for key2pop in ['base_session', 'session_order', 'signal_dataset', 'tracking_dataset']:
            _ = key.pop(key2pop)

        session_dict = make_multi_session_object_dict(key)

        bin_dict = {}
        # Complement session dictionary
        for session, session_entry in session_dict.items():
            # Ratemap and fields
            rm, mask = (Ratemap.proj('ratemap', 'mask_rm') & session_entry & key).fetch1('ratemap', 'mask_rm')
            session_dict[session]['ratemap'] = np.ma.array(rm, mask = mask)
            # Object
            if session != 'base':
                # Take care of firing fields
                session_dict[session]['fields']  = (Ratemap.Fields & session_entry & key).fetch(order_by='field_no ASC', as_dict=True)

                # Take care of objects / positions
                obj_x, obj_y     = (ArenaObjectPos & session_entry & key).fetch1('obj_x_coord_calib','obj_y_coord_calib')
                x_edges, y_edges = (Occupancy & session_entry & key).fetch1('x_edges','y_edges')

                # ... Where is the object in ratemap "coordinates" (bins)
                bin_size_rm_x = np.mean(np.diff(x_edges))
                bin_size_rm_y = np.mean(np.diff(y_edges))

                # Save bin size for later
                bin_dict[session] = np.mean([bin_size_rm_x, bin_size_rm_y])

                obj_x_rm = ((obj_x - x_edges[0]) / bin_size_rm_x) - .5
                obj_y_rm = ((obj_y - y_edges[0]) / bin_size_rm_y) - .5

                session_dict[session]['object_x']    = obj_x
                session_dict[session]['object_y']    = obj_y
                session_dict[session]['object_x_rm'] = obj_x_rm
                session_dict[session]['object_y_rm'] = obj_y_rm


        master_dict = {
            'base_session'             :  session_dict['base']['session_name'],
            'object1_session'          :  session_dict['object1']['session_name'],
            'object2_session'          :  session_dict['object2']['session_name'],
                     }
        self.insert1({**key_,**master_dict})

        # CAVE!
        # As of November 2020 the field coordinate x/y values are "inverted"
        # both for field_peak and field_centroid

        dist_dict_list = []
        for field1 in session_dict['object1']['fields']:
            field1_no = field1['field_no']
            # Get angle and distance:
            y_diff = field1['field_centroid_x']-session_dict['object1']['object_y_rm']
            x_diff = field1['field_centroid_y']-session_dict['object1']['object_x_rm']

            object1_dist  = np.sqrt(np.square(x_diff) + np.square(y_diff))
            object1_dist *= bin_dict['object1']
            object1_angle = np.arctan2(y_diff, x_diff)
            object1_angle = (object1_angle + 2 * np.pi) % (2 * np.pi) # Make sure it's [0,2*pi]

            field1_y, field1_x = y_diff, x_diff

            # Now for the second object session (or the second object)
            field2_no = []
            dists_1_2 = [] # Distances in between fields in object session 2 and current field in object session 1

            object2_dists =  []  # Distance of field to object in object session2
            object2_angles = []  # Angles to object
            for field2 in session_dict['object2']['fields']:
                field2_no.append(field2['field_no'])
                # Get angle and distance:
                y_diff = field2['field_centroid_x']-session_dict['object2']['object_y_rm']
                x_diff = field2['field_centroid_y']-session_dict['object2']['object_x_rm']

                object2_dist  = np.sqrt(np.square(x_diff) + np.square(y_diff))
                object2_dist *= bin_dict['object2']
                object2_angle = np.arctan2(y_diff, x_diff)
                object2_angle = (object2_angle + 2 * np.pi) % (2 * np.pi) # Make sure it's [0,2*pi]
                # Save
                object2_dists.append(object2_dist)
                object2_angles.append(object2_angle)

                field2_y, field2_x = y_diff, x_diff

                # Get distance between fields (field1 and current field2)
                dist_1_2 = np.sqrt(np.square(field2_x - field1_x) + np.square(field2_y - field1_y))
                dist_1_2 *= np.mean([bin_dict['object1'], bin_dict['object2']]) # Convert to [mm]
                dists_1_2.append(dist_1_2)


            # Field number of matching field in object session 2
            if len(dists_1_2):
                field_idx_min = np.argmin(dists_1_2)
                dist_dict = {
                    'field_no_object1'    : field1_no,
                    'field_no_object2'    : field2_no[field_idx_min],
                    'dist_fields'         : np.min(dists_1_2),
                    'dist_to_object'      : np.mean([object1_dist, object2_dists[field_idx_min]]),
                    'angle_to_object'     : circmean([object1_angle, object2_angles[field_idx_min]])
                }

                dist_dict_list.append(dist_dict)
            else:
                # Just skip this field lookup
                continue

        if not len(dist_dict_list):
            # Just stop. A master entry is still being written so results are not re-calculated when missing
            return

        ################# FIELD RATE COMPARISONS #########################################################################################################################

        # Loop over fields and calculate rate change to baseline session
        # While doing that, insert into part table .Fields()

        # Transform ratemaps from masked to arrays with nans
        rm_base     = session_dict['base']['ratemap'].filled(fill_value=np.nan)
        rm_object_1 = session_dict['object1']['ratemap'].filled(fill_value=np.nan)
        rm_object_2 = session_dict['object2']['ratemap'].filled(fill_value=np.nan)

        for dist_dict in dist_dict_list:
            points_field_object1 = (Ratemap.Fields & session_dict['object1'] & key
                                    & 'field_no = {}'.format(dist_dict['field_no_object1'])).fetch1('field_coords')
            points_field_object2 = (Ratemap.Fields & session_dict['object2'] & key
                                    & 'field_no = {}'.format(dist_dict['field_no_object2'])).fetch1('field_coords')

            with warnings.catch_warnings():
                warnings.filterwarnings(action='ignore', message='Mean of empty slice')
                rate_field_object1_base = np.nanmean([rm_base[point[0], point[1]] for point in points_field_object1])
                rate_field_object2_base = np.nanmean([rm_base[point[0], point[1]] for point in points_field_object2])
                rate_field_object1      = np.nanmean([rm_object_1[point[0], point[1]] for point in points_field_object1])
                rate_field_object2      = np.nanmean([rm_object_2[point[0], point[1]] for point in points_field_object2])

            rel_base_object1 = rate_field_object1 / (rate_field_object1_base + np.finfo(float).eps) # Prevent divisions by zero
            rel_base_object2 = rate_field_object2 / (rate_field_object2_base + np.finfo(float).eps) # Prevent divisions by zero

            field_dict = {
                'object1_field_id'   :    dist_dict['field_no_object1'],
                'object2_field_id'   :    dist_dict['field_no_object2'],
                'dist_fields'        :    dist_dict['dist_fields'],
                'dist_to_object'     :    dist_dict['dist_to_object'],
                'angle_to_object'    :    dist_dict['angle_to_object'],
                'object1_field_base' :    rel_base_object1,
                'object2_field_base' :    rel_base_object2
                }

            self.Fields.insert1({**key_, **field_dict}, skip_duplicates=True)



@imhotte
class OVCScores(dj.Computed):
    '''
    Object vector cells

    PART B: Object vector map based comparisons



    CAVEAT:
    This currently works only for a limited set of session types / configurations:
    - 1 baseline, 2 object sessions with one object each
    - 1 baseline, 1 object session with two objects

    '''

    definition = """
    # Object vector cell (OVC) vector map score based calculations
    -> SignalTracking.proj(base_session='session_name')
    -> OVParams
    -> ShuffleParams
    ---
    object1_session                 : varchar(16)    # Object session 1
    object2_session                 : varchar(16)    # Object session 2
    ovscore                         : double         # Object vector score (2D correlation between OV maps)
    """

    ##### SHUFFLED OBJECT VECTOR SCORE ###############################################################################################

    class ShuffledOVScore(dj.Part):
        definition = """
        # Shuffled Object vector (OV) score and shuffling
        -> master
        ---
        shuffled_ovscores_95perc        : double          # Object vector score shuffling for cell: 95th percentile
        shuffled_ovscores_99perc        : double          # Object vector score shuffling for cell: 99th percentile
        shuffled_ovscores               : blob@imgstore   # Object vector score shuffling for cell
        """

    #######################################################################################################################

    @property
    def key_source(self):
        '''
        Filter for metasessions for which exactly 2 objects were presented either concurrently or in two subsequent sessions.
        Return base sessions for those.
        '''
        object_sessions = Session.SessionType.proj() * ArenaObjectPos & 'sessiontype = "Open Field Object"'
        meta_object_sessions = MetaSession.aggr(object_sessions, n="count(*)") & 'n=2' # Filter out everything that is not 2 objects
        base_sessions = Session.SessionType & meta_object_sessions.proj() & 'sessiontype = "Open Field"'
        return SignalTracking.proj(base_session='session_name') \
                    * OVParams \
                    * ShuffleParams \
                    & base_sessions.proj(base_session='session_name') \
                    & Shuffled.proj(base_session='session_name') \
                    & 'signal_type = "spikes"'\
                    & 's_t_params_id = "A"'


    def make(self, key):
        # Clean up key
        # make() of OVC()
        key_ = key.copy()

        # Get rid of some keys
        for key2pop in ['base_session', 'session_order', 'signal_dataset', 'tracking_dataset']:
            _ = key.pop(key2pop)

        session_dict = make_multi_session_object_dict(key)

        # Complement session dictionary
        for session, session_entry in session_dict.items():
            # Object
            if session != 'base':
                # OVMap
                ovmap, mask = (OVMap.proj('ovmap', 'mask_ovmap') & session_entry & key).fetch1('ovmap', 'mask_ovmap')
                session_dict[session]['ovmap'] = np.ma.array(ovmap, mask = mask)


        # Get score
        ovscore = corr2(session_dict['object1']['ovmap'],session_dict['object2']['ovmap'])

        master_dict = {
            'base_session'             :  session_dict['base']['session_name'],
            'object1_session'          :  session_dict['object1']['session_name'],
            'object2_session'          :  session_dict['object2']['session_name'],
            'ovscore'                  :  ovscore
                     }

        self.insert1({**key_,**master_dict})

        ##### SHUFFLING
        object1_session_shuffled_ovmaps = get_shuffled_ovmaps(session_dict['object1'], copy.deepcopy(key))
        if session_dict['object1']['session_name'] == session_dict['object2']['session_name']:
            object1_session_shuffled_ovmaps, object2_session_shuffled_ovmaps = split_list(object1_session_shuffled_ovmaps)
        else:
            object2_session_shuffled_ovmaps = get_shuffled_ovmaps(session_dict['object2'], copy.deepcopy(key))

        # Get shuffled scores
        shuffled_ovscores = []
        for ovmap1 in object1_session_shuffled_ovmaps:
            for ovmap2 in object2_session_shuffled_ovmaps:
                shuffled_ovscores.append(corr2(ovmap1, ovmap2))

        shuffled_scores_dict = {
            'shuffled_ovscores_95perc' :  np.percentile(shuffled_ovscores, 95),
            'shuffled_ovscores_99perc' :  np.percentile(shuffled_ovscores, 99),
            'shuffled_ovscores'        :  shuffled_ovscores,
                               }

        self.ShuffledOVScore.insert1({**key_,**shuffled_scores_dict})


def get_shuffled_ovmaps(session_key, shuffling_key):
    '''
    This re-implements many of the functions
    already established for the main Shuffled() table.

    Returns
    -------
    shuffled_ovmaps : list : all shuffled OV maps

    '''
    #### OVMAP related
    occupancy_entry   = (OVOccupancy & session_key).fetch1()
    params_ov         = (OVParams & session_key).fetch1()

    object_data   = (ArenaObjectPos & session_key).fetch1()
    occupancy = np.ma.array(occupancy_entry['occupancy_ov'], mask=occupancy_entry['mask_occ_ov'])

    angular_edges = occupancy_entry['angular_edges_ov'].copy()
    radial_edges  = occupancy_entry['radial_edges_ov'].copy()

    #### Spikes / Tracking / Sync
    st_params = {}
    st_params['speed_cutoff_low'], st_params['speed_cutoff_high'], time_offset = (SignalTrackingParams & session_key & shuffling_key).fetch1(
        'speed_cutoff_low', 'speed_cutoff_high', 'time_offset')

    spikes   = (FilteredSpikes.proj(signal_dataset='dataset_name', spikes='filtered_spikes')
                & session_key & shuffling_key).fetch1('spikes')
    tracking = (Tracking.OpenField * Tracking.proj(tracking_dataset='dataset_name')
                & session_key & shuffling_key).fetch1()

    center_y, center_plane, proj_mean_img = ((Projection.proj('mean_image') * Cell.Rois).proj(
        ..., signal_dataset='dataset_name') & session_key & shuffling_key).fetch1('center_y', 'center_plane', 'mean_image')

    num_planes, frame_rate_si, seconds_per_line, width_SI, height_SI = (Tif.SI & (Session & session_key & shuffling_key)).fetch1(
        'num_scanning_depths', 'framerate', 'seconds_per_line', 'width_scanimage', 'height_scanimage')


    # Special case where the SI image shape is not the same with the projection image shape
    # - suggesting the image has been cropped prior to suite2p analysis
    # - thus, the "center_y" is no longer accurate -> using the middle line for "center_y"
    if proj_mean_img.shape != (width_SI, height_SI):
        center_y = int(height_SI/2)

    sync_data_frames, sample_rate_sync, shuffling_key['sync_dataset_frames_imaging'], shuffling_key['sync_name_frames_imaging'] = \
                        (MetaSession.Setup * \
                         Setup.Sync * Sync \
                         & 'generic_name = "frames_imaging"' \
                         & session_key & shuffling_key).fetch1(
                            'sync_data', 'sample_rate', 'dataset_name', 'sync_name')

    tracking_type = (Dataset * Tracking & session_key).fetch1('datasettype')
    if tracking_type == 'DLC_tracking':
        tracking_generic = 'TrackingDLC'
    elif 'Tracking2D_2LED' in tracking_type:
        tracking_generic = 'Tracking2LED'
    else:
        raise NotImplementedError(f'Tracking dataset type {tracking_type} not implemented')

    sync_data_track = (MetaSession.Setup * Setup.Sync
                    * Sync & f'generic_name = "{tracking_generic}"' &  session_key & shuffling_key).fetch1('sync_data')

    # Sanity checks
    # 1. Compare length of tracking data and tracking sync data
    # 2. Compare length of spike data and frame sync data
    # 3. Compare last timestamp frame sync data and tracking sync data

    if len(tracking['x_pos']) != len(sync_data_track):
        raise IndexError('Mismatch between length of sync data and tracking data')
    if len(spikes) != len(sync_data_frames):
        raise IndexError('Mismatch between length of sync data and spiking data')
    if np.abs(sync_data_track[-1] - sync_data_frames[-1]) > np.mean(np.diff(sync_data_frames)):
        raise IndexError('There is more than one frame difference between the end of sync streams')

    # -> Compensate for the mismatch between timestamp at the beginning of each frame and the time it takes the laser to reach the cell body
    seconds_per_plane =  1 / (frame_rate_si * num_planes) # from Tif.SI()
    # Why is this correct? Because frame_rate_si returns the "volume" rate.

    # Careful! "samples" are floating point (real valued) timestamps for pre-synced setups since
    # sample rate = 1. for those sync data
    seconds_to_cell  = center_plane * seconds_per_plane + center_y * seconds_per_line
    samples_to_cell  = seconds_to_cell * sample_rate_sync
    samples_offset   = samples_to_cell + (time_offset * sample_rate_sync) # from 'MapParams'

    # Only spikes are shuffled at the moment
    signal_data   =  spikes
    signal_idxs   =  np.argwhere(spikes > 0).squeeze()

    # Retrieve shuffling offsets from Shuffled() table, limit by 1 (since Shuffled() depends on many more parameter sets)
    shuffling_offsets = (Shuffled & session_key & shuffling_key).fetch('shuffling_offsets', limit=1)[0]

    shuffled_ovmaps = []
    for shift in tqdm(shuffling_offsets):
        rolled_sync_data_frames = np.roll(sync_data_frames, shift)

        # Look up signal tracking
        signal_dict = calc_signal_dict(rolled_sync_data_frames, sync_data_track,
                                                tracking, signal_data, signal_idxs, samples_offset, st_params)

        # Calculate object centered angle and distance for every signal tracking point
        angles, dists = get_dist_angle_obj(
                    signal_dict['x_pos_signal'],
                    signal_dict['y_pos_signal'],
                    object_data['obj_x_coord_calib'],
                    object_data['obj_y_coord_calib']
                )

        # Supplement 'angles' and 'dists' to retrieved signaltracking_entry
        signal_dict['x_pos_signal'] = angles
        signal_dict['y_pos_signal'] = dists

        # Need to rename a parameter for calc_ratemap() to work
        params_ov['sigma_signal'] = params_ov['sigma_signal_ov']
        ovmap_dict = calc_ratemap(occupancy, angular_edges, radial_edges, signal_dict, params_ov, pad_mode='wrap')

        shuffled_ovmaps.append(ovmap_dict['ratemap'])

    return shuffled_ovmaps


##### OVC SUMMARY TABLES #########################################################################################################

@imhotte
class OVCutoffs(dj.Lookup):
    definition = """
    # Object vector cell cutoffs
    ov_cutoff_id                   : char(1)        # Parameter set ID, starting with A
    ---
    info_content_cutoff            : varchar(100)   # Information content cutoff (>). Session level: session_95, session_99
    ovscore_cutoff                 : varchar(100)   # Object vector score cutoff (>). Session level: session_95, session_99
    dist_fields_cutoff             : float          # Distance [mm] (object centered field distances) (<)
    dist_to_object_cutoff          : float          # Distance [mm] to object (>)
    object1_field_base_cutoff      : float          # Relative rate of field in object session 1 compared to base session (>)
    object2_field_base_cutoff      : float          # Relative rate of field in object session 1 compared to base session (>)
    """
    contents = [
        {
         'ov_cutoff_id'         : 'A',
         'info_content_cutoff'  : 'information_content_95',
         'ovscore_cutoff'       : 'shuffled_ovscores_95perc',
         'dist_fields_cutoff'   : 250.,
         'dist_to_object_cutoff':  40.,
         'object1_field_base_cutoff': 1.5,
         'object2_field_base_cutoff': 1.5,
        }
                ]

@imhotte
class OVC(dj.Computed):
    '''
    Object vector cells

    Final table: Summary

    '''

    definition = """
    # Object vector cell (OVC) summary table
    -> OVCScores
    -> OVCutoffs
    -> OVCFields
    ---
    object1_session                : varchar(16)    # Object session 1
    object2_session                : varchar(16)    # Object session 2
    ovscore                        : double         # Object vector score (2D correlation between OV maps)
    is_ovc                         : tinyint        # 0 - not an OVC according to cutoffs, 1 - putative OVC
    no_fields                      : int            # Number of filtered fields (matching cutoff criteria)
    mean_dist_to_object = NULL     : double         # Average distance of (filtered) fields to object [mm]
    mean_dist_fields    = NULL     : double         # Average distance between fields [mm]
    mean_angle_to_object = NULL    : double         # Circular mean of field angles to object [0, 2*pi]
    std_angle_to_object = NULL     : double         # Circular standard deviation for field angles to object [radians]
    field_ids = NULL               : blob@imgstore  # Field IDs list of dictionaries ('object1_field_id', 'object2_field_id')
    angles_to_object = NULL        : blob@imgstore  # Field angles [0, 2*pi]
    dists_to_object = NULL         : blob@imgstore  # Distances of (filtered) fields to object [mm]
    dists_fields = NULL            : blob@imgstore  # Distances of (filtered) fields to object [mm]
    """

    @property
    def key_source(self):
        # Need to constrain by CutoffsOVScore since the computation for this is implemented separately
        return OVCScores.proj() * OVCutoffs.proj() * OVCFields.proj() & CutoffsOVScore

    def make(self, key):

        ovc_cutoffs = (OVCutoffs & key).fetch1()
        ovscore     = (OVCScores & key).fetch1('ovscore')

        cell_key = {}
        for k,v in key.items():
            if k not in ['session_order','base_session','tracking_dataset','signal_dataset']:
                cell_key[k] = v
        object1_session, object2_session = (OVCScores & key).fetch1('object1_session','object2_session')

        ####### INFO CONTENT FILTER ############################################################################################
        info_labels = []
        for obj_session in [object1_session, object2_session]:
            # Check whether we are dealing with
            # - Session level cutoff
            # - Single number cutoff

            if ovc_cutoffs["info_content_cutoff"] == 'session_95':
                info_content_cutoff = _get_info_content_session_95(cell_key, obj_session)
            elif ovc_cutoffs["info_content_cutoff"] == 'session_99':
                info_content_cutoff = _get_info_content_session_99(cell_key, obj_session)
            else:
                info_content_cutoff = ovc_cutoffs["info_content_cutoff"]

            label_ = len(Ratemap.Stats * Shuffled.RatemapStats\
                            & f'information_content > {info_content_cutoff}'\
                            & f'session_name = "{obj_session}"' \
                            & cell_key)

            info_labels.append(label_)
        info_content_filter = int((np.array(info_labels) > 0).all())


        ####### OVSCORE FILTER #################################################################################################
        # Score filter
        # Check whether we are dealing with
        # - Session level cutoff
        # - Single number cutoff

        if ovc_cutoffs["ovscore_cutoff"] == 'session_95':
            ovscore_cutoff = (CutoffsOVScore & key).fetch('ovscore_95',limit=1)[0] # Because there can be multiple session params
        elif ovc_cutoffs["ovscore_cutoff"] == 'session_99':
            ovscore_cutoff = (CutoffsOVScore & key).fetch('ovscore_99',limit=1)[0]
        else:
            ovscore_cutoff = ovc_cutoffs["ovscore_cutoff"]

        score_filter = len(OVCScores * OVCScores.ShuffledOVScore.proj('shuffled_ovscores_95perc','shuffled_ovscores_99perc') \
                            & f'ovscore > {ovscore_cutoff}' \
                            & key)

        # Field filter
        field_filter = OVCFields * OVCFields.Fields  \
                            &  f'dist_fields         < {ovc_cutoffs["dist_fields_cutoff"]}' \
                            &  f'dist_to_object     > {ovc_cutoffs["dist_to_object_cutoff"]}' \
                            &  f'object1_field_base  > {ovc_cutoffs["object1_field_base_cutoff"]}' \
                            &  f'object2_field_base  > {ovc_cutoffs["object2_field_base_cutoff"]}'\
                            & key
        no_fields = len(field_filter)

        # Collect field properties
        if no_fields:
            field_ids         = []
            dists_to_object   = []
            angles_to_object  = []
            dists_fields      = []
            for field in field_filter:
                # Transcribe field ID dict
                field_id_dict = {}
                field_id_dict['object1_field_id'] = field['object1_field_id']
                field_id_dict['object2_field_id'] = field['object2_field_id']
                field_ids.append(field_id_dict)

                dists_to_object.append(field['dist_to_object'])
                angles_to_object.append(field['angle_to_object'])
                dists_fields.append(field['dist_fields'])

            mean_dist_to_object  = np.nanmean(dists_to_object)
            mean_dist_fields      = np.nanmean(dists_fields)
            mean_angle_to_object = circmean(angles_to_object, nan_policy='raise')
            std_angle_to_object  = circstd(angles_to_object, nan_policy='raise')

        else:
            field_ids              = None
            dists_to_object       = None
            angles_to_object      = None
            dists_fields           = None
            mean_dist_to_object   = None
            mean_dist_fields       = None
            mean_angle_to_object  = None
            std_angle_to_object   = None

        # Make a decision - OVC yes (1) or no (0) ?
        is_ovc = (np.array([info_content_filter, score_filter, int(no_fields>0)]) == 1).all()

        # Create dictionary
        entry_dict = {
                'object1_session'       : object1_session,
                'object2_session'       : object2_session,
                'ovscore'               : ovscore,
                'is_ovc'                : int(is_ovc),
                'no_fields'              : no_fields,
                'mean_dist_to_object'   : mean_dist_to_object,
                'mean_dist_fields'       : mean_dist_fields,
                'mean_angle_to_object'  : mean_angle_to_object,
                'std_angle_to_object'   : std_angle_to_object,
                'field_ids'              : field_ids,
                'angles_to_object'      : angles_to_object,
                'dists_to_object'       : dists_to_object,
                'dists_fields'           : dists_fields
                }
        self.insert1({**key,**entry_dict})


##############################################################################################################################
############# OVC SUMMARY TABLES HELPERS

def __get_session_key(cell_key, session_name=None):
    ''' Strip cell ID and add session name '''
    session_key = {}
    for k,v in cell_key.items():
        if k not in ['cell_id']:
            session_key[k] = v
    if session_name is not None:
        session_key['session_name'] = session_name
    return session_key

#### SPATIAL INFORMATION CONTENT #############################################################################################

def _get_info_content_session_95(cell_key, session_name):
    '''
    Retrieving from Shuffled.RatemapStats
    with STANDARD cell parameters as written in "CONSTANTS"

    '''
    info_content_key = __get_session_key(cell_key, session_name)
    filtered_cells = get_filtered_cells(Session & info_content_key, verbose=False)
    information_content_shuffles = (Shuffled.RatemapStats & info_content_key & filtered_cells).fetch('information_content_shuffles')
    percentile95 = np.nanpercentile(np.concatenate(information_content_shuffles), 95)
    return percentile95

def _get_info_content_session_99(cell_key, session_name):
    '''
    Retrieving from Shuffled.RatemapStats
    with STANDARD cell parameters as written in "CONSTANTS"

    '''
    info_content_key = __get_session_key(cell_key, session_name)
    filtered_cells = get_filtered_cells(Session & info_content_key, verbose=False)
    information_content_shuffles = (Shuffled.RatemapStats & info_content_key & filtered_cells).fetch('information_content_shuffles')
    percentile99 = np.nanpercentile(np.concatenate(information_content_shuffles), 99)
    return percentile99

#### OVSCORE  ################################################################################################################

# Deprecated since OVSCore session filtering has been implemented separately (since cutoff retrieval is expensive)
# def _get_ovscore_session_95(cell_key):
#     ovscore_key = __get_session_key(cell_key)
#     filtered_cells = get_filtered_cells(Session & ovscore_key, verbose=False)
#     # "base_name" vs. "session_name" in filtered_cells does not matter (meta session level!)
#     ovscore_shuffles = (OVCScores.ShuffledOVScore & ovscore_key & filtered_cells).fetch('shuffled_ovscores')
#     percentile95 = np.nanpercentile(np.concatenate(ovscore_shuffles), 95)
#     return percentile95

# def _get_ovscore_session_99(cell_key):
#     ovscore_key = __get_session_key(cell_key)
#     filtered_cells = get_filtered_cells(Session & ovscore_key, verbose=False)
#     # "base_name" vs. "session_name" in filtered_cells does not matter (meta session level!)
#     ovscore_shuffles = (OVCScores.ShuffledOVScore & ovscore_key & filtered_cells).fetch('shuffled_ovscores')
#     percentile99 = np.nanpercentile(np.concatenate(ovscore_shuffles), 99)
#     return percentile99



##############################################################################################################################
##############################################################################################################################

############## HELPERS

def get_dist_angle_obj(tracking_x, tracking_y, obj_x, obj_y):
    '''
    Calculate object centered distance and angle for
    each tracking point.

    '''
    angles = []
    dists = []

    for x,y in zip(tracking_x, tracking_y):
        angle = np.arctan2(y-obj_y, x-obj_x)
        angle = (angle + 2 * np.pi) % (2 * np.pi) # Make sure range is [0, 2*pi]
        angles.append(angle)

        dist = np.sqrt(np.square(x-obj_x) + np.square(y-obj_y))
        dists.append(dist)

    angles = np.array(angles)
    dists  = np.array(dists)

    return angles, dists

def find_nearest(array, value):
    ''' Find nearest element to "value" in "array" and return index of that element '''
    idx = np.searchsorted(array, value, side='left')
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1
    else:
        return idx


def calc_signal_dict(signal_sync, tracking_sync, tracking_data, signal_data, signal_idxs, samples_offset, params):
    '''
    # DUPLICATE FROM MAIN SCHEMA UNDER
    # SPATIAL_SCORES.PY!

    Look up tracking signal for every valid calcium event/spike based on sync data.

    Parameters
    ----------
    signal_sync : np.array
        Sync pulse number for every event in 'signal_data'
    tracking_sync : np.array
        Sync pulse number for every event in 'tracking_data'
    tracking_data : dict
        Fetched Tracking.OpenField() entry with:
            - x_pos
            - y_pos
            - head_angle
            - speed
    signal_data : np.array
        Calcium or spikes signal of length 'signal_sync'
    signal_idxs : np.array
        Indices of signal to keep (filter)
    samples_offset : int
        Number of samples to add to 'signal_sync' samples before
        lookup of tracking data. This shifts signal in time compared
        to tracking data and is used (1) to compensate for laser fly time
        in FOV and (2) for shuffling purposes
    params: dict
        - speed_cutoff_low: lower speed cutoff in unit 'speed' has in 'tracking_data'
        - speed_cutoff_low: upper speed cutoff in unit 'speed' has in 'tracking_data'

    Returns
    -------
    signal_dict : dict
        x_pos_signal
        y_pos_signal
        head_angle_signal
        speed_signal
        signal (amplitudes)
    '''

    tracking_keys = {'x_pos', 'y_pos', 'head_angle', 'speed'}
    signal_dict = {}

    tracking_indices_filtered = []  # Speed filtered
    signal_indices_filtered   = []  # Spikes / Fluorescence

    signal_sync          = signal_sync.astype(float)
    signal_sync          += samples_offset  # Shift in time
    signal_sync          = signal_sync[signal_idxs]

    for idx, sync_pulse in zip(signal_idxs, signal_sync): # signal_idxs either filtered or whole series (spikes vs. delta f/f)
        tracking_idx = find_nearest(tracking_sync, sync_pulse)
        if (tracking_data['speed'][tracking_idx] > params['speed_cutoff_low']) and (tracking_data['speed'][tracking_idx] < params['speed_cutoff_high']):
            tracking_indices_filtered.append(tracking_idx)
            signal_indices_filtered.append(idx)

    # Build signal tracking dictionary
    for pos_key in tracking_keys:
        signal_dict[pos_key + '_signal'] = tracking_data[pos_key][tracking_indices_filtered]

    signal_dict['signal'] = signal_data[signal_indices_filtered].squeeze()
    if not signal_dict["signal"].ndim == 1:
        # Somewhat convoluted logic: I'm not sure what circumstances call for the `.squeeze()`,
        # but where there is only a single value, that compresses it to a 0d array
        # Datajoint later converts that to a non-array floating point number during the
        # conversion to/from a blob. Therefore, enforce that the array still has at least 1 dimension
        signal_dict["signal"] = np.expand_dims(signal_dict["signal"], 0)
    return signal_dict




# DUPLICATE FROM MAIN SCHEMA UNDER
# SPATIAL_SCORES.PY!
def calc_ratemap(occupancy, x_edges, y_edges, signaltracking, params, pad_mode='symmetric'):
    '''

    # DUPLICATE FROM MAIN SCHEMA UNDER
    # SPATIAL_SCORES.PY!


    Calculate ratemap
    Parameters
    ----------
    occupancy : masked np.array
        Smoothed occupancy. Masked where occupancy low
    x_edges : np.array
        Bin edges in x
    y_edges : np.array
        Bin edges in y
    signaltracking : dict
        SignalTracking table entry
    params : dict
        MapParams table entry
    pad_mode : str
        Padding mode:
        - 'symmetric' (normal 2D ratemaps) or
        - 'wrap' (object vector maps)

    Returns
    -------
    ratemap_dict : dict
        - binned_raw : np.array: Binned raw (unsmoothed) signal
        - ratemap_raw: np masked array: Unsmoothed ratemap (mask where occupancy low)
        - ratemap    : np masked array: Smoothed ratemap (mask where occupancy low)
        - bin_max    : tuple   : (x,y) coordinate of bin with maximum signal
        - max        : float : Max of signal

    '''
    ratemap_dict = {}
    sigma_signal = params['sigma_signal']

    binned_signal = np.zeros_like(occupancy.data)
    # Add one at end to not miss signal at borders
    x_edges[-1] += 1
    y_edges[-1] += 1

    # Look up signal per bin
    for no_x in range(len(x_edges)-1):
        for no_y in range(len(y_edges)-1):
            boolean_x = (signaltracking['x_pos_signal'] >= x_edges[no_x]) & (signaltracking['x_pos_signal'] < x_edges[no_x+1])
            boolean_y = (signaltracking['y_pos_signal'] >= y_edges[no_y]) & (signaltracking['y_pos_signal'] < y_edges[no_y+1])
            extracted_signal = signaltracking['signal'][boolean_x & boolean_y]
            binned_signal[no_y, no_x] = np.nansum(extracted_signal)

    ratemap_dict['binned_raw'] = binned_signal
    binned_signal = np.ma.masked_where(occupancy.mask, binned_signal)  # Masking. This step is probably unnecessary
                                                                       # since occupancy is already masked
    ratemap_dict['ratemap_raw'] = binned_signal / occupancy


    # Use astropy.convolve to smooth padded version of the spikemap

    binned_signal = np.ma.filled(binned_signal, np.nan) # First convert masked values to nan
    kernel = Gaussian2DKernel(x_stddev=sigma_signal) # Create astropy gaussian kernel

    if pad_mode == 'wrap':
        # For radial maps like object vector maps
        # Logic: Wrap around / pad the "angular" edges and leave the "radial" edges unchanged
        pad_width = ((0,
                      0),
                     (int(5*sigma_signal),
                      int(5*sigma_signal))) # Only pad one axis (axis=0: radius, axis=1: angles)

        binned_signal_padded = np.pad(binned_signal, pad_width=pad_width, mode='wrap')  # 'wrap' for radial maps

        # There are some caveats here. For example for "extended" regions of nans that are
        # hard / impossible to interpolate. This will lead to edge effects but isn't usually
        # of concern.
        # However it does become problematic for OV (object-vector) map calculations.
        # Therefore, if pad_mode =='wrap', i.e. when radial 2D maps are fed in, convert nans to zeros and proceed.
        with warnings.catch_warnings():
                # Astropy throws a warning when nans are detected after smoothing (usually edge artefacts)
                warnings.filterwarnings(action='ignore', message=r'.*?') # Catch "any" warning
                binned_signal_smoothed = convolve(
                                            np.nan_to_num(binned_signal_padded),
                                            kernel,
                                            boundary='extend'
                                            )


    elif pad_mode == 'symmetric':
        # i.e. normal 2D ratemap. "Symmetric" is taken over from BNT function.
        pad_width = ((int(5*sigma_signal),
                      int(5*sigma_signal)),
                     (int(5*sigma_signal),
                      int(5*sigma_signal))) # Symmetric in all directions

        binned_signal_padded = np.pad(binned_signal, pad_width=pad_width, mode='symmetric')

        # Everything else but wrap (i.e. usually 'symmetric') for normal 2D ratemaps
        with warnings.catch_warnings():
                # Astropy throws a warning when nans are detected after smoothing (usually edge artefacts)
                warnings.filterwarnings(action='ignore', message=r'.*?') # Catch "any" warning
                binned_signal_smoothed = convolve(
                                            binned_signal_padded,
                                            kernel,
                                            boundary='extend'
                                            )
    else:
        raise NotImplementedError(f'Pad mode "{pad_mode}" not implemented')

    # Crop out non-padded part
    binned_signal_smoothed = binned_signal_smoothed[pad_width[0][0]:[-pad_width[0][1] if pad_width[0][1] >0 else None][0],
                                                    pad_width[1][0]:[-pad_width[1][1] if pad_width[1][1] >0 else None][0]]
    binned_signal_smoothed = np.ma.masked_where(occupancy.mask, binned_signal_smoothed)  # Masking. This step is probably unnecessary
                                                                                         # since occupancy is already masked
    masked_ratemap = binned_signal_smoothed / occupancy

    ratemap_dict['ratemap']       = masked_ratemap
    ratemap_dict['bin_max']       = np.unravel_index(masked_ratemap.argmax(), masked_ratemap.shape)
    ratemap_dict['max']           = np.max(masked_ratemap)

    return ratemap_dict

# Part of helpers:

### THIS COMES FROM HELPERS!
def split_list(a_list):
    half = len(a_list)//2
    return a_list[:half], a_list[half:]
