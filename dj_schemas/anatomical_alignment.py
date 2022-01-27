### Scores, etc. over anatomical space
import numpy as np
import datajoint as dj

import pickle
from itertools import combinations
from collections import OrderedDict
from tqdm.auto import tqdm 

from skimage.transform import warp
from skimage.feature import canny

from esda import Moran
from libpysal.weights import lat2W

#### LOAD DATABASE ###########################################################
from .dj_conn import *
imhotte = dj.schema(horst_imaging_db)
 
### CONSTANTS
from .constants import GLOBAL_CUTOFF_NO_STARTER_CELLS

### HELPERS
from helpers_topography.utils import corr2
from helpers_topography.fov_stitching import (align_session_imgs, 
                                             evaluate_unwarping, 
                                             BinningParamsFOV, 
                                             make_composite_fov, 
                                             get_fov_score_map, 
                                             shift_image,
                                             get_x_y_lim,
                                             stack_project_fov_score_maps)


@imhotte 
class AlignmentPoints(dj.Manual):
    definition = """
    # User specified alignment points between two FOVs
    -> MetaSession.proj(metasession_ref='metasession_name')
    -> MetaSession.proj(metasession_align='metasession_name')
    ---    
    projection_table        : enum('Projection','ProjectionCorr')   # Projection table the alignmet images were retrieved from  
    num_points              : int                                   # Number of user defined alignment points
    points_ref              : blob@imgstore                         # List of points specified for the reference session 
    points_align            : blob@imgstore                         # List of points specified for the aligned session 
    transform_method        : varchar(255)                          # Which transform method to use, e.g. "Euclidean" or "Affine"
    entry_time_align_points = CURRENT_TIMESTAMP : timestamp
    """

@imhotte
class FOVAlignmentParams(dj.Lookup):
    definition = """
    # Parameters for FOV alignment
    fov_align_param_id : char(1)     # Parameter set ID, starting with A
    ---
    vector_y           :  float      # How much to nudge in y 
    vector_x           :  float      # How much to nudge in x
    padding            :  float      # Overall padding

    """
    contents = [
          ['A', -200, 100, 800], # These settings center the miniscope pictures quite well 
               ]
               
@imhotte
class FOVProjParam(dj.Lookup):
    definition = """
    # Name of Projection/ProjectionCorr key used for calculating image alignment metrics
    projection_short  : char(8)      # Projection key short name
    ---
    projection_key    : varchar(50)  # Projection key (minus "_corr") used for calculating FOV alignment metrics

    """
    contents = [
          ['mean_img', 'mean_image'],
          ['mean_im2', 'mean_image_second']
               ]


@imhotte
class AlignmentFOV(dj.Computed):
    definition = """
    # Aligned FOV data
    -> AlignmentPoints  
    -> FOVAlignmentParams
    -> FOVProjParam
    ---
    padded_ref              :  blob@imgstore   # Padded / Shifted reference FOVs over planes (dictionary over planes)
    padded_align            :  blob@imgstore   # Padded / Shifted aligned FOV (warped) (dictionary over planes)
    padded_align_raw        :  blob@imgstore   # Padded / Shifted but non-warped FOV (dictionary over planes)
    ssim_original = NULL    :  double          # SSIM: Structural Similarity Index - original (average over planes)
    ssim_warped = NULL      :  double          # SSIM: Structural Similarity Index - after warping / alignment (average over planes)
    mse_original = NULL     :  double          # MSE: Mean squared error - original (average over planes)
    mse_warped = NULL       :  double          # MSE: Mean squared error - after warping / alignment (average over planes)
    transform               :  blob@imgstore   # Skimage transform object (same for all planes)
    """
    
    def make(self, key):
        '''
        Create a warping transform from user specified points. 
        Since FOV boundaries might be too tight to accomodate warping, 
        first pad and nudge everything such that there is enough space
        around the images for alignment. 
        
            
        '''        
        alignment_params = (FOVProjParam * FOVAlignmentParams & key).fetch1()
        alignment_points = (AlignmentPoints & key).fetch1()
        transform_method = alignment_points['transform_method']


        # Clarify metasession name for ref and align session 
        metasession_ref    = alignment_points['metasession_ref']
        metasession_align  = alignment_points['metasession_align']

        # Clarify user defined points
        points_align = alignment_points['points_align']
        points_ref   = alignment_points['points_ref']

        # Package in to ordered dicts and sort just to be sure 
        points_align = OrderedDict(sorted(points_align.items()))
        points_ref   = OrderedDict(sorted(points_ref.items()))

        points_align = np.array(list(points_align.values()))
        points_ref   = np.array(list(points_ref.values()))

        projection_table = alignment_points['projection_table']
        if projection_table == 'ProjectionCorr':
            projection_key = alignment_params['projection_key'] + '_corr'
            x_range_key = 'x_range_microns_full'
            y_range_key = 'y_range_microns_full'
        elif projection_table == 'Projection':
            projection_key = alignment_params['projection_key'] 
            x_range_key = 'x_range'
            y_range_key = 'y_range'
        else:
            raise NotImplementedError(f'{projection_table} not understood')        

        vector = (alignment_params['vector_y'], alignment_params['vector_x']) # Nudge everything over a bit ...
        padding = alignment_params['padding'] # ... and add overall padding
        vector_corr = tuple(np.array(vector) - padding)


        # Create transform options dict
        trans_ops = {
            'vector'         : vector,
            'padding'        : int(padding),
            'vector_corr'    : vector_corr,
            'projection_key' : projection_key,
             }

        #### CAREFUL
        #### Results could stem from different number of planes - accounted for below
        #### REFERENCE = LEFT 
        results_left = (Session * eval(projection_table) & \
                            f'metasession_name = "{metasession_ref}"').fetch(
                                                                            projection_key, 
                                                                            x_range_key, 
                                                                            y_range_key,
                                                                            'center_plane',
                                                                            as_dict=True
                                                                            )
        ##### ALIGN = RIGHT                                                                       
        results_right = (Session * eval(projection_table) & \
                            f'metasession_name = "{metasession_align}"').fetch(
                                                                            projection_key, 
                                                                            x_range_key, 
                                                                            y_range_key,
                                                                            'center_plane',
                                                                            as_dict=True
                                                                            )

        ################ SHIFT / WARP ######################################################################
        # LEFT: NON WARPED 
        padded_ref = {}
        for result in results_left:
            assert np.isnan(result[projection_key]).all() == False, f'No data found for projection key "{projection_key}"'
            projection_shifted, _, _ = align_session_imgs(
                                                        result[projection_key], 
                                                        points_ref, 
                                                        points_align, 
                                                        trans_ops, 
                                                        transform=False,
                                                        transform_method=transform_method,
                                                                )
            
            padded_ref[result['center_plane']] = projection_shifted

        # RIGHT: NON WARPED
        padded_align_raw = {}
        for result in results_right:
            assert np.isnan(result[projection_key]).all() == False, f'No data found for projection key "{projection_key}"'
            projection_shifted, _, _ = align_session_imgs(
                                                        result[projection_key], 
                                                        points_ref, 
                                                        points_align, 
                                                        trans_ops, 
                                                        transform=False,
                                                        transform_method=transform_method,
                                                                )
            
            padded_align_raw[result['center_plane']] = projection_shifted

        # RIGHT: WARPED
        padded_align = {}
        for result in results_right:
            _, projection_unwarped, tform = align_session_imgs(
                                                        result[projection_key], 
                                                        points_ref, 
                                                        points_align, 
                                                        trans_ops, 
                                                        transform=True,
                                                        transform_method=transform_method,
                                                                )
            
            padded_align[result['center_plane']] = projection_unwarped


        ################ EVALUATE UNWARPING ##############################################################
        x_range_left, y_range_left   = results_left[0][x_range_key], results_left[0][y_range_key]
        x_range_right, y_range_right = results_right[0][x_range_key], results_right[0][y_range_key]

        # Loop over results and at the end, save average / max 
        # padded align = picture unwarped 
        # padded ref   = picture_reference
        # projection left 
        # projection right 

        # Loop over results from left and right (reference and unwarped images)
        ssim_originals  = []
        ssim_warpeds    = []
        mse_originals   = []
        mse_warpeds     = []

        for result_left in results_left:
            for result_right in results_right:
                
                picture_unwarped  = padded_align[result_right['center_plane']]
                picture_reference = padded_ref[result_left['center_plane']]
                projection_left   = result_left[projection_key]
                projection_right  = result_right[projection_key]
                # ssim_original, ssim_warped, mse_original, mse_warped
                results_comparison = evaluate_unwarping(picture_unwarped, picture_reference, projection_left, projection_right,  \
                                                    y_range_left, x_range_left, y_range_right, x_range_right, verbose=True)
                
                # Append to lists ... 
                ssim_originals.append(results_comparison[0])
                ssim_warpeds.append(results_comparison[1])
                mse_originals.append(results_comparison[2])
                mse_warpeds.append(results_comparison[3])
                
        entry_dict = {
            'padded_ref'       : padded_ref, # = picture_reference
            'padded_align'     : padded_align,
            'padded_align_raw' : padded_align_raw,
            'ssim_original'    : np.nanmin(ssim_originals),
            'ssim_warped'      : np.nanmax(ssim_warpeds),
            'mse_original'     : np.nanmax(mse_originals),
            'mse_warped'       : np.nanmin(mse_warpeds),
            'transform'        : pickle.dumps(tform),
                     }
        self.insert1({**key,**entry_dict})




#######################################################################################################################################
######################### SCORE MAPS FOV


@imhotte
class ScoremapFOVParams(dj.Lookup):
    definition = """
    # Parameters for score maps
    binning_param_short : char(6)      # Parameter short name
    ---
    binning_param_set   : varchar(100) # Score map binning parameters set (identifies parameter set under BinningParamsFOV)
    bin_size_microns    : float        # Bin size in microns

    """
    # Pre-fill with some basics
    contents = [
          ['grid__', 'grid', 5.],
          ['border', 'border_bvs', 5.],
          ['m_v_l_', 'mvl', 5.],
          ['o_v_c_', 'ovc', 5.],
          ['info_c', 'info_content', 5.],
               ]

@imhotte
class ScoremapFOV(dj.Computed):
    definition = """
    # Aligned score maps
    -> MetaSession.proj(metasession_ref='metasession_name')
    -> ScoremapFOVParams
    -> FilteredCellsParams
    -> FilteredSessionsParams
    ---
    bins_x = NULL                     :   blob@imgstore   # Bins in x 
    bins_y = NULL                     :   blob@imgstore   # Bins in y
    aligned_metas                     :   blob@imgstore   # Dictionary containing all aligned sessions
    aligned_metas_shuff               :   blob@imgstore   # Dictionary containing shuffled + aligned sessions
    no_sessions                       :   int             # Number of aligned sessions
    no_shuffles                       :    blob@imgstore   # Number of shuffles
    min_no_shuffles                   :    int             # Minimum number of shuffles across sessions
    """ 

    @property
    def key_source(self):
        keys = super().key_source & AlignmentFOV
        return keys & FilteredSessions.proj(metasession_ref = 'metasession_name') \
                    & FilteredCells.proj(metasession_ref = 'metasession_name') 
        
        
    def make(self, key):
        '''
        Create collection of aligned score maps
        This inherits from MetaSession and not Session and the only reason why this works is because
        FilteredSession entries are unique per metasession - and that is only because I filtered for 
        "Open Field" sessiontype and there is usually not more than one session. 


        '''  

        # Fix entries from AlignmentFOV to single parameter set before retrieving metasessions 
        alignment_keys = {
            'fov_align_param_id'  : 'A',
            'projection_short'    : 'mean_img'
        }
        # Retrieve alignment parameters
        align_params = (FOVAlignmentParams & alignment_keys).fetch1()
        vector = (align_params['vector_y'], align_params['vector_x']) # Nudge everything over a bit ...
        padding = align_params['padding'] # ... and add overall padding
        vector_corr = tuple(np.array(vector) - padding)

        # Get parameter dictionary to constrain cell results
        parameter_dict_cell = (FilteredCellsParams & key).fetch1('parameter_dict_cell')

        # Get list of valid metasessions (difference SSIM > .025, hard cutoff)
        metasessions = AlignmentFOV.proj(diff_ssim = 'ssim_warped-ssim_original') & alignment_keys & 'diff_ssim > .00001' & key 

        # Create composite from all metasessions to extract boundaries of stitched (=composite) FOV
        _, x_min, x_max, y_max, y_min = make_composite_fov(metasessions, padding_composite=0)

        # Create binning vector
        params  = (ScoremapFOVParams & key).fetch1()
        bin_size_microns = params['bin_size_microns'] # micrometer, but NOT(!) FOV corrected (FOV_RATIO)
        bins_x = np.arange(x_min, x_max, bin_size_microns)
        bins_y = np.arange(y_min, y_max, bin_size_microns)

        # Define set of parameters 
        binning_params_class = BinningParamsFOV(bins_x=bins_x, bins_y=bins_y) # Initialize class
        binning_params = getattr(binning_params_class, params['binning_param_set'])


        ##### LOOP OVER METASESSIONS ##################################################################################################

        shuffle_iter = 1000 # Max number of iterations
        scoremap_dict = {}
        scoremap_dict_shuff = {}
        no_shuffles = [] # Collect number of actual shuffles 

        for no, metasess in enumerate(metasessions.proj()): 

            # Create new metasession key from current key             
            metasession_ref_key   = MetaSession.proj() & 'metasession_name = "{}"'.format(metasess["metasession_ref"])
            metasession_align_key = MetaSession.proj() & 'metasession_name = "{}"'.format(metasess["metasession_align"])
                
            # Transform
            tform = pickle.loads((AlignmentFOV & metasess).fetch1('transform')) 

            ######## ROI DATA #########################################################################################################
            
            # Find out which table to probe (Corr or normal)
            if len(RoisCorr & metasession_ref_key):
                corr = True
            else:
                corr = False

            if corr:
                cells_unfiltered = (RoisCorr.proj('center_x_corr', 
                                                  'center_y_corr', 
                                                  'ypix_corr', 
                                                  'xpix_corr') \
                                                  * binning_params['scoretables']) 
            else:
                # Cannot use FOV corrected data, so revert to Cell.Rois and just rename columns ...
                cells_unfiltered = (Cell.Rois.proj(center_x_corr = 'center_x',
                                                   center_y_corr = 'center_y',
                                                   ypix_corr = 'ypix',
                                                   xpix_corr = 'xpix') \
                                                   * binning_params['scoretables']) 

            ######## PROCESS POINTS AND GENERATE MAPS #################################################################################
            if no == 0:
                # Process reference once ... 
                #### REFERENCE
                cells_ref_filtered = cells_unfiltered & FilteredCells & parameter_dict_cell & metasession_ref_key & key     
                if len(cells_ref_filtered) > GLOBAL_CUTOFF_NO_STARTER_CELLS: 
                    score_map_ref, score_map_ref_shuf, _, _  = get_fov_score_map(cells_ref_filtered, binning_params, vector_corr, shuffle_iter=shuffle_iter)
                    # Save
                    scoremap_dict[metasess["metasession_ref"]]       = score_map_ref
                    scoremap_dict_shuff[metasess["metasession_ref"]] = score_map_ref_shuf
                    no_shuffles.append(len(score_map_ref_shuf))

            #### ALL OTHER
            cells_align_filtered = cells_unfiltered & FilteredCells & parameter_dict_cell & metasession_align_key & key
            if len(cells_align_filtered) > GLOBAL_CUTOFF_NO_STARTER_CELLS:
                score_map_align, score_map_align_shuff, _, _ = get_fov_score_map(cells_align_filtered, binning_params, vector_corr, tform=tform, shuffle_iter=shuffle_iter)
                # Save
                scoremap_dict[metasess["metasession_align"]]         = score_map_align
                scoremap_dict_shuff[metasess["metasession_align"]]   = score_map_align_shuff
                no_shuffles.append(len(score_map_align_shuff))


        ### BUILD ENTRY DICT

        entry_dict = {
            'bins_x'               : bins_x,
            'bins_y'               : bins_y,
            'aligned_metas'        : scoremap_dict,
            'aligned_metas_shuff'  : scoremap_dict_shuff,
            'no_sessions'          : len(scoremap_dict),
            'no_shuffles'          : no_shuffles,
            'min_no_shuffles'      : [np.min(no_shuffles) if len(scoremap_dict_shuff) else 0][0]
        }

        self.insert1({**key,**entry_dict})



@imhotte
class ScoremapFOVProjParams(dj.Lookup):
    definition = """
    # FOV Score map projection parameters
    proj_param_id       : char(1)      # Parameter ID
    ---
    proj_statistic      : varchar(100) # Score map projection statistic (as in BinningsStats())
    proj_sigma          : float        # Smoothing sigma

    """
    contents = [
          ['A', 'nanmean', 2.],
          ['B', 'nanmean', 0.],
               ]



@imhotte
class ScoremapFOVMoran(dj.Computed):
    definition = """
    # Calculate Moran's I (spat. autocorr measure) for FOV score maps 
    -> ScoremapFOV
    -> ScoremapFOVProjParams
    ---
    moran_i = NULL             :  double         # Moran's I of current FOV score map 
    moran_i_shuffles = NULL    :  blob@imgstore  # Moran's I shuffling distribution 
    moran_i_95 = NULL          :  double         # 95th percentile of Moran's I shuffling distribution 
    moran_i_99 = NULL          :  double         # 99th percentile of Moran's I shuffling distribution 
    moran_i_5 = NULL           :  double         # 5th percentile of Moran's I shuffling distribution 
    moran_i_1 = NULL           :  double         # 1st percentile of Moran's I shuffling distribution 


    """ 
    @property
    def key_source(self):
        ''' 
        For Moran's I calculation do NOT take any parameter sets (ScoremapFOVProjParams) that smooth the 
        score maps. This will lead to wrong assumptions. 
        '''
        sigma_0 = (ScoremapFOVProjParams & 'proj_sigma < 0.00001').proj() # so ... zero
        len_session = (ScoremapFOV & 'no_sessions > 0').proj()
        keys = super().key_source & sigma_0 & len_session
        return keys

    def make(self,key):
        ''' 
        Perform Moran’s I global autocorrelation statistic calculation 
        on FOV score maps. 
        More info: 
        https://pysal.org/esda/generated/esda.Moran.html
        ''' 

        params = (ScoremapFOVProjParams & key).fetch1()
        score_map_dict, score_map_dict_shuff = (ScoremapFOV & key).fetch1('aligned_metas','aligned_metas_shuff')  

        projection, _, _, _, _ = stack_project_fov_score_maps(
                                                                score_map_dict, 
                                                                statistic=params['proj_statistic'], 
                                                                sigma=0
                                                                )


        # Create W object (spatial weights class W) for a regular lattice
        # https://pysal.org/libpysal/generated/libpysal.weights.lat2W.html
        w = lat2W(projection.shape[0],projection.shape[1]) 

        moran = Moran(np.nan_to_num(projection).flatten(), w)
        real_moran_I = moran.I

        # Create shuffled distribution for Moran's I  
        mean_no_shuffles = np.mean([len(v) for v in score_map_dict_shuff.values()])
        moran_i_shuff = []

        # Fix number of shuffles 
        n_shuffles = 5000
        for _ in tqdm(range(n_shuffles)):
            
            # For each shuffling iteration, create random projection
            # by sampling from shuffling dictionary 
            # The logic for sampling is: 
            # Take a random integer (index idx) for every list of shuffled maps in score_map_dict_shuff
            # Some of them can have very low numbers of cells, and therefore total number of shuffles.
            # In order to prevent over-sampling those examples, create an index of at least mean(no_shuffles) over all 
            # shuffles. If an index is out-of-bounds just skip this map for this shuffling iteration  

            scoremap_dict_ = {}
            for k,v in score_map_dict_shuff.items():
                idx = np.random.randint(0,np.max([len(v),mean_no_shuffles])-1)
                if idx > (len(v)-1):
                    continue
                else:
                    scoremap_dict_[k] = v[idx]
                
            projection_shuff, _, _, _, _= stack_project_fov_score_maps(
                                                    scoremap_dict_, 
                                                    statistic=params['proj_statistic'], 
                                                    sigma=0
                                                    )
            moran_shuff = Moran(np.nan_to_num(projection_shuff).flatten(), w)
            moran_i_shuff.append(moran_shuff.I)


        moran_i_shuff = np.array(moran_i_shuff)

        entry_dict = {
            'moran_i'          : real_moran_I,
            'moran_i_shuffles' : moran_i_shuff,
            'moran_i_95'       : np.nanpercentile(moran_i_shuff,95),
            'moran_i_99'       : np.nanpercentile(moran_i_shuff,99),
            'moran_i_5'        : np.nanpercentile(moran_i_shuff,5),
            'moran_i_1'        : np.nanpercentile(moran_i_shuff,1),
        }

        self.insert1({**key,**entry_dict})


@imhotte
class ScoremapCorr(dj.Computed):
    definition = """
    # Correlation of FOV score maps 
    -> ScoremapFOV.proj(binning_param_A = "binning_param_short")
    -> ScoremapFOV.proj(binning_param_B = "binning_param_short")
    -> ScoremapFOVProjParams
    ---
    xcorr  = NULL                     :  double           # Cross correlation 
    xcorr_shuffles = NULL             :  blob@imgstore    # Cross correlation shuffling distribution
    xcorr_95 = NULL                   :  double           # 95th percentile shuffling distribution xcorr
    xcorr_99 = NULL                   :  double           # 99th percentile shuffling distribution xcorr
    xcorr_5 = NULL                    :  double           # 5th percentile shuffling distribution xcorr
    xcorr_1 = NULL                    :  double           # 1st percentile shuffling distribution xcorr

    """ 
    @property
    def key_source(self):
        keys = super().key_source.proj() - (super().key_source.proj()  & 'binning_param_A = binning_param_B')

        # Exclude session dictionary lengths < 1
        A_toofew = ScoremapFOV.proj(binning_param_A='binning_param_short', no_sessions='no_sessions') & keys & 'no_sessions=0'
        B_toofew = ScoremapFOV.proj(binning_param_B='binning_param_short', no_sessions='no_sessions') & keys & 'no_sessions=0'

        keys = keys - A_toofew.proj() - B_toofew.proj()

        # Exclude all "reverse" pairs. Unelegant but it works. 
        available_binnings = ScoremapFOVParams.fetch('binning_param_short')
        available_combinations = list(combinations(available_binnings, 2))
        key2exclude = []
        for comb in available_combinations:
            exc = (keys & f'binning_param_A = "{comb[1]}"' & f'binning_param_B = "{comb[0]}"').fetch('KEY')
            key2exclude.append(exc)
        key2exclude = [item for sublist in key2exclude for item in sublist]
        keys = keys - key2exclude

        return keys


    def make(self, key):
        ''' 
        Perform cross correlation analysis between projected score maps 
        and compare to shuffled maps 

        '''

        params = (ScoremapFOVProjParams & key).fetch1()
        score_map_dict_A, score_map_dict_shuff_A = (ScoremapFOV.proj(
                                                            binning_param_A='binning_param_short', 
                                                            aligned_metas="aligned_metas", 
                                                            aligned_metas_shuff="aligned_metas_shuff"
                                                            ) & key).fetch1('aligned_metas','aligned_metas_shuff')
        score_map_dict_B, score_map_dict_shuff_B = (ScoremapFOV.proj(
                                                            binning_param_B='binning_param_short', 
                                                            aligned_metas="aligned_metas", 
                                                            aligned_metas_shuff="aligned_metas_shuff"
                                                            ) & key).fetch1('aligned_metas','aligned_metas_shuff')

        # Make sure that the same (meta) sessions are compared to each other
        score_map_dict_A, score_map_dict_B             = prune_keys(score_map_dict_A, score_map_dict_B)
        score_map_dict_shuff_A, score_map_dict_shuff_B = prune_keys(score_map_dict_shuff_A, score_map_dict_shuff_B)

        projection_A, _, _, _, _ = stack_project_fov_score_maps(
                                                            score_map_dict_A, 
                                                            statistic=params['proj_statistic'], 
                                                            sigma=params['proj_sigma'],
                                                            )
        projection_B, _, _, _, _ = stack_project_fov_score_maps(
                                                            score_map_dict_B, 
                                                            statistic=params['proj_statistic'], 
                                                            sigma=params['proj_sigma'],
                                                            )


        assert projection_A.shape == projection_B.shape, 'Map A and B have different dimensions.'
        real_x_corr = corr2(projection_A, projection_B)

        # Create shuffled distribution for x corr
        mean_no_shuffles_A = np.mean([len(v) for v in score_map_dict_shuff_A.values()])
        mean_no_shuffles_B = np.mean([len(v) for v in score_map_dict_shuff_B.values()])

        x_corr_shuff = []

        # Fix number of shuffles 
        n_shuffles = 10000
        for _ in tqdm(range(n_shuffles)):

            # For each shuffling iteration, create random projection
            # by sampling from shuffling dictionary 
            # The logic for sampling is: 
            # Take a random integer (index idx) for every list of shuffled maps in score_map_dict_shuff
            # Some of them can have very low numbers of cells, and therefore total number of shuffles.
            # In order to prevent over-sampling those examples, create an index of at least mean(no_shuffles) over all 
            # shuffles. If an index is out-of-bounds just skip this map for this shuffling iteration  

            
            #### FOR A ################################################################################################
            scoremap_dict_ = {}
            for k,v in score_map_dict_shuff_A.items():
                idx = np.random.randint(0,np.max([len(v),mean_no_shuffles_A])-1)
                if idx > (len(v)-1):
                    continue
                else:
                    scoremap_dict_[k] = v[idx]

            projection_shuff_A, _, _, _, _= stack_project_fov_score_maps(
                                                    scoremap_dict_, 
                                                    statistic=params['proj_statistic'], 
                                                    sigma=params['proj_sigma'],
                                                    )
            
            
            
            #### FOR B ################################################################################################
            scoremap_dict_ = {}
            for k,v in score_map_dict_shuff_B.items():
                idx = np.random.randint(0,np.max([len(v),mean_no_shuffles_B])-1)
                if idx > (len(v)-1):
                    continue
                else:
                    scoremap_dict_[k] = v[idx]

            projection_shuff_B, _, _, _, _= stack_project_fov_score_maps(
                                                    scoremap_dict_, 
                                                    statistic=params['proj_statistic'], 
                                                    sigma=params['proj_sigma'],
                                                    )
            
            assert projection_shuff_A.shape == projection_shuff_B.shape

            x_corr_shuff.append(corr2(projection_shuff_A,projection_shuff_B))
        
        x_corr_shuff = np.array(x_corr_shuff)


        entry_dict = {
            'xcorr'            : real_x_corr,
            'xcorr_shuffles'   : x_corr_shuff,
            'xcorr_95'         : np.nanpercentile(x_corr_shuff, 95),
            'xcorr_99'         : np.nanpercentile(x_corr_shuff, 99),
            'xcorr_5'          : np.nanpercentile(x_corr_shuff, 5),
            'xcorr_1'          : np.nanpercentile(x_corr_shuff, 1),
        }
        self.insert1({**key,**entry_dict})



####################################################################################################################################
# HELPERS

def prune_keys(A, B):
    '''
    Function used in the processing of FOV scoremap xcorr 
    e.g. dj_schemas 'ScoremapCorr()' 
    and helpers_topography.plotting_helpers 'plot_summary_xcorr_scoremaps'

    Make sure that for xcorrs the same sessions are being compared.

    '''
    B_notA = set(B.keys()) - set(A.keys())
    A_notB = set(A.keys()) - set(B.keys())
    for key in B_notA:
        _ = B.pop(key)
    for key in A_notB:
        _ = A.pop(key)
        
    assert set(A) == set(B), 'Dictionaries do not have the same keys'
    return A,B


def _draw_mask_collage(collected_masks, cmap='cubehelix_r'):
    ''' 
    Helper for plotting results from 
    get_aligned_anatomical_masks()


    '''
    import seaborn as sns 
    from matplotlib import pyplot as plt 

    limits = get_x_y_lim(collected_masks, thresh=-1, padding=-120)

    sns.set(style='white',font_scale=2)

    plt.figure(figsize=(15,15))
    plt.imshow(collected_masks, cmap=cmap) 
    plt.xlim(limits[0]+50,limits[1]-50)
    plt.ylim(limits[1],limits[0])
    plt.gca().get_xaxis().set_ticks([]);plt.gca().get_yaxis().set_ticks([])
    sns.despine(left=True,bottom=True)
    return


def get_aligned_anatomical_masks(animal_name, 
                                 region, 
                                 projection_short, 
                                 center_plane=0,
                                 edge=False,
                                 plot=False):
    
    '''
    Retrieve anatomical mask entries and align them 
    with the help of AlignmentFOV() entries 

    Parameter
    ---------
    animal_name : str
                  Animal name (mLIMS)
    region      : str

    '''

    def __edge(image, thresh=1):
        ''' 
        Extract (surrounding) edge only

        '''
        image = image.copy()

        image_mask = image.copy()
        image_mask[image_mask==0] = np.min(image_mask)
        image_mask[image_mask>thresh] = image.max()
        image_mask[image_mask<thresh] = 0

        edges = canny(
            image=image_mask,
            sigma=4,
            low_threshold=0.,
            high_threshold=.99999,
            mask=None,
            use_quantiles=True,
        )
        return edges.astype(int)


    # Anatomical masks
    masks_animal   = (AnatomicalMask * AnatomicalMaskParams) & (Session & f'animal_name = "{animal_name}"') & f'center_plane={center_plane}'
    # Alignments 
    alignment_fovs = AlignmentFOV & masks_animal.proj(metasession_ref='metasession_name') & f'projection_short = "{projection_short}"'

    
    collected_masks=None
    for mask_entry in masks_animal:


        # Prepare mask
        mask = mask_entry['anat_mask'].copy()
        region_label = mask_entry[region]
        assert region_label is not None, f'"{region}" is not a valid'
        mask[mask!=region_label] = 0

        # Figure out if this is a reference or an aligned session 
        # Retrieve padding parameters and tform if available
        ref_   = alignment_fovs.proj(...,metasession_name = "metasession_ref")   & mask_entry
        align_ = alignment_fovs.proj(...,metasession_name = "metasession_align") & mask_entry

        if ref_:
            #print('This is the reference')
            tform = None
            alignment_params = (FOVAlignmentParams & ref_).fetch1()
        elif align_:
            #print('This is an aligned session')
            tform_ = align_.fetch1('transform')
            tform  = pickle.loads(tform_)
            alignment_params = (FOVAlignmentParams & align_).fetch1()

        else:
            #raise IndexError(f'No AlignmentFOV entry found for this session\n{mask_entry}')
            print((f'No AlignmentFOV entry found for this session\n{mask_entry}'))
            continue

        vector = (alignment_params['vector_y'], alignment_params['vector_x'])
        padding = alignment_params['padding']


        mask_padded = (np.pad(mask, int(padding)))
        mask_shifted = shift_image(mask_padded, vector)

        if tform is not None:
            mask_warped = warp(mask_shifted, tform, order=3, preserve_range=True)
        else:
            mask_warped = mask_shifted

        # Extract (only) edge? 
        if edge: 
            mask_warped = __edge(mask_warped)

        if collected_masks is None:
            collected_masks = mask_warped
        else: 
            collected_masks += mask_warped

    # Plot
    if plot: _draw_mask_collage(collected_masks)
    return collected_masks





def get_scoremaps_a_b(key):
    '''
    Helper for plot_summary_xcorr_scoremaps()
    under helpers_topography.plotting_helpers

    Retrieve entry of (ScoremapFOV * ScoremapFOVParams) for two scoremaps -
    `binning_param_A` and `binning_param_B` in `key`

    Parameters
    ----------
    key     :  dict
               Retrieval key for ScoremapFOV * ScoremapFOVParams
               containing both `binning_param_A` and `binning_param_B` entries

    '''
    score_maps_entry_A = ((ScoremapFOV * ScoremapFOVParams).proj(binning_param_A     = 'binning_param_short', 
                                                                 binning_param_set_A = 'binning_param_set', 
                                                                 aligned_metas       = 'aligned_metas', 
                                                                 ) & key).fetch1()

    score_maps_entry_B = ((ScoremapFOV * ScoremapFOVParams).proj(binning_param_B     = 'binning_param_short', 
                                                                 binning_param_set_B = 'binning_param_set',
                                                                 aligned_metas       = 'aligned_metas', 
                                                                 ) & key).fetch1()


    return score_maps_entry_A, score_maps_entry_B



def get_scoremaps_moran(animal, 
                        binning_param_set, 
                        proj_param_id,
                        metasession_ref=None,
                        param_hash_id_cell='standard',
                        param_hash_session='cf83e1357eefb8bd'
                        ):
    '''
    Helper for plot_summary_scoremaps() 
    under helpers_topography.plotting_helpers
    
    Retrieve ScoremapFOV entries and ScoremapFOVMoran results for a given animal. 
    
    Parameters
    ----------
    animal            : str
                        animal_name 
    binning_param_set : str
                        secondary key(binning_param_set) of ScoremapFOVParams()
    proj_param_id     : str
                        Primary key of ScoremapFOVProjParams()
    metasession_ref   : str
                        Metasession name in case there is more than one 
                        valid metasession 
    param_hash_id_cell: str
                        Primary key of FilteredCellsParams()
                        Default: 'standard'
    param_hash_session: str
                        Primary key of FilteredSessionsParams()   
                        Default: 'cf83e1357eefb8bd'

    Returns
    -------
    score_maps        : dict
                        Dictionary of retrieved ScoremapFOV entries
                        (aligned_metas)
    no_sessions       : int
                        No of sessions (len(score_maps))
    moran_is          : dict
                        Dictionary of retrieved ScoremapFOVMoran entries
                        ('moran_i_shuffles', 'moran_i', 'moran_i_95', 'moran_i_99')
    score_map_keys    : dict
                        Primary keys used to fetch from 
                        ScoremapFOV * ScoremapFOVMoran * ScoremapFOVParams
    '''
    

    sessions_animal = (Session.proj(animal='animal_name', 
                                    metasession_ref='metasession_name') 
                                    & f'animal = {animal}')
    if metasession_ref is not None: 
        sessions_animal = sessions_animal & f'metasession_ref = "{metasession_ref}"'
    assert len(sessions_animal), f'No sessions for animal {animal} (metasession_ref filter: {metasession_ref}) were found'
    
    scoremap_moran = (ScoremapFOV * ScoremapFOVMoran * ScoremapFOVParams 
                            & f'binning_param_set = "{binning_param_set}"' 
                            & f'proj_param_id = "{proj_param_id}"'
                            & f'param_hash_id_cell = "{param_hash_id_cell}"'
                            & f'param_hash_session = "{param_hash_session}"'
                            & sessions_animal.proj())
    assert len(scoremap_moran), f'No results found for animal {animal} | binning param {binning_param_set}'
    
    
    score_maps, no_sessions = scoremap_moran.fetch1('aligned_metas','no_sessions')
    moran_is                = scoremap_moran.fetch('moran_i_shuffles', 'moran_i', 'moran_i_95', 'moran_i_99', as_dict=True)[0]
    score_map_keys          = scoremap_moran.fetch('KEY')

    return score_maps, no_sessions, moran_is, score_map_keys


