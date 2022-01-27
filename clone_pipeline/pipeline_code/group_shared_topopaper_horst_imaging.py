import datajoint as dj


# -------------- group_shared_topopaper_horst_imaging --------------


schema = dj.Schema('group_shared_topopaper_horst_imaging')

vmod0 = dj.VirtualModule('vmod0', 'group_shared_topopaper_main_imaging')


@schema
class AlignmentPoints(dj.Manual):
    definition = """
    # User specified alignment points between two FOVs
    -> vmod0.MetaSession.proj(metasession_ref="metasession_name")
    -> vmod0.MetaSession.proj(metasession_align="metasession_name")
    ---
    projection_table     : enum('Projection','ProjectionCorr') # Projection table the alignmet images were retrieved from
    num_points           : int                                 # Number of user defined alignment points
    points_ref           : longblob                            # List of points specified for the reference session
    points_align         : longblob                            # List of points specified for the aligned session
    transform_method     : varchar(255)                        # Which transform method to use, e.g. "Euclidean" or "Affine"
    entry_time_align_points=CURRENT_TIMESTAMP : timestamp                    
    """


@schema
class FOVAlignmentParams(dj.Lookup):
    definition = """
    # Parameters for FOV alignment
    fov_align_param_id   : char(1)                             # Parameter set ID, starting with A
    ---
    vector_y             : float                                # How much to nudge in y
    vector_x             : float                                # How much to nudge in x
    padding              : float                                # Overall padding
    """


@schema
class FOVProjParam(dj.Lookup):
    definition = """
    # Name of Projection/ProjectionCorr key used for calculating image alignment metrics
    projection_short     : char(8)                            # Projection key short name
    ---
    projection_key       : varchar(50)                        # Projection key (minus "_corr") used for calculating FOV alignment metrics
    """


@schema
class FilteredCellsParams(dj.Lookup):
    definition = """
    # Parameters saved as parameter_dict_cells and restriction_dict_cell and are uniquely identified by parameter_hash_cell
    param_hash_id_cell="standard" : enum('standard','ov_cutoff_D','ov_cutoff_H','ov_cutoff_I') 
    ---
    param_hash_cell      : char(16)                     # 16 character hash of used restrictions AND parameters dict
    restriction_dict_cell : longblob                    # Dictionary of filters used in query (Restrictions)
    parameter_dict_cell  : longblob                     # Dictionary of parameters used in query (Parameters)
    """


@schema
class ScoremapFOVParams(dj.Lookup):
    definition = """
    # Parameters for score maps
    binning_param_short  : char(6)                      # Parameter short name
    ---
    binning_param_set    : varchar(100)                 # Score map binning parameters set (identifies parameter set under BinningParamsFOV)
    bin_size_microns     : float                         # Bin size in microns
    """


@schema
class ScoremapFOVProjParams(dj.Lookup):
    definition = """
    # FOV Score map projection parameters
    proj_param_id        : char(1)                      # Parameter ID
    ---
    proj_statistic       : varchar(100)                 # Score map projection statistic (as in BinningsStats())
    proj_sigma           : float                         # Smoothing sigma
    """


@schema
class FilteredSessionsParams(dj.Lookup):
    definition = """
    # Restrictions saved as restriction_dict_session and uniquely identified by parameter_hash_session
    param_hash_session   : char(16)                     # 16 character hash of used restrictions
    ---
    restriction_dict_session : longblob                 # Dictionary of filters used in query
    """


@schema
class PairwDistParams(dj.Lookup):
    definition = """
    # Parameters for pairwise distance analysis
    pairwise_dist_param  : char(1)                      # Param ID
    ---
    score                : varchar(50)                  # Score (column) name
    score_cutoff         : varchar(50)                  # Score cutoff value (defines starter cell population)
    scoretables          : varchar(1000)                # Datajoint tables that 'score' can be found in (e.g. 'GridScore')
    """


@schema
class AnatomicalMaskParams(dj.Lookup):
    definition = """
    # LUT for anatomical masks drawn in Napari to identify subregions in FOV
    timestamp_mask_lookup=CURRENT_TIMESTAMP : timestamp                    # Auto created timestamp
    ---
    mec_label            : tinyint                      # Medial entorhinal cortex (MEC) label name (integer)
    pas_label=null       : tinyint                      # Parasubiculum label name (integer)
    rsa_label=null       : tinyint                      # Retrosplenial agranular cortex (integer)
    prh_label=null       : tinyint                      # Perirhinal cortex
    """


@schema
class AlignmentFOV(dj.Computed):
    definition = """
    # Aligned FOV data
    -> AlignmentPoints
    -> FOVAlignmentParams
    -> FOVProjParam
    ---
    padded_ref           : longblob               # Padded / Shifted reference FOVs over planes (dictionary over planes)
    padded_align         : longblob               # Padded / Shifted aligned FOV (warped) (dictionary over planes)
    padded_align_raw     : longblob               # Padded / Shifted but non-warped FOV (dictionary over planes)
    ssim_original=null   : double                 # SSIM: Structural Similarity Index - original (average over planes)
    ssim_warped=null     : double                 # SSIM: Structural Similarity Index - after warping / alignment (average over planes)
    mse_original=null    : double                 # MSE: Mean squared error - original (average over planes)
    mse_warped=null      : double                 # MSE: Mean squared error - after warping / alignment (average over planes)
    transform            : longblob               # Skimage transform object (same for all planes)
    """


@schema
class ScoremapFOV(dj.Computed):
    definition = """
    # Aligned score maps
    -> vmod0.MetaSession.proj(metasession_ref="metasession_name")
    -> ScoremapFOVParams
    -> FilteredCellsParams
    -> FilteredSessionsParams
    ---
    bins_x=null          : longblob               # Bins in x
    bins_y=null          : longblob               # Bins in y
    aligned_metas        : longblob               # Dictionary containing all aligned sessions
    aligned_metas_shuff  : longblob               # Dictionary containing shuffled + aligned sessions
    no_sessions          : int                    # Number of aligned sessions
    no_shuffles           : longblob               # Number of shuffles
    min_no_shuffles       : int                    # Minimum number of shuffles across sessions
    """


@schema
class ScoremapFOVMoran(dj.Computed):
    definition = """
    # Calculate Moran's I (spat. autocorr measure) for FOV score maps
    -> ScoremapFOV
    -> ScoremapFOVProjParams
    ---
    moran_i=null         : double                  # Moran's I of current FOV score map
    moran_i_shuffles=null : longblob                # Moran's I shuffling distribution
    moran_i_95=null      : double                  # 95th percentile of Moran's I shuffling distribution
    moran_i_99=null      : double                  # 99th percentile of Moran's I shuffling distribution
    moran_i_5=null       : double                  # 5th percentile of Moran's I shuffling distribution
    moran_i_1=null       : double                  # 1st percentile of Moran's I shuffling distribution
    """


@schema
class ScoremapCorr(dj.Computed):
    definition = """
    # Correlation of FOV score maps
    -> ScoremapFOV.proj(binning_param_A="binning_param_short")
    -> ScoremapFOV.proj(binning_param_B="binning_param_short")
    -> ScoremapFOVProjParams
    ---
    xcorr=null           : double                       # Cross correlation
    xcorr_shuffles=null  : longblob                      # Cross correlation shuffling distribution
    xcorr_95=null        : double                       # 95th percentile shuffling distribution xcorr
    xcorr_99=null        : double                       # 99th percentile shuffling distribution xcorr
    xcorr_5=null         : double                       # 5th percentile shuffling distribution xcorr
    xcorr_1=null         : double                       # 1st percentile shuffling distribution xcorr
    """


@schema
class FilteredSessions(dj.Manual):
    definition = """
    # Filtered sessions (manually selected by user)
    -> vmod0.Tracking.proj(tracking_dataset="dataset_name")
    -> FilteredSessionsParams
    ---
    entry_time_filt_sessions=CURRENT_TIMESTAMP : timestamp                    
    """


@schema
class FilteredCells(dj.Manual):
    definition = """
    # Filtered cells (manually selected by user)
    -> FilteredSessions
    -> vmod0.Cell
    -> FilteredCellsParams
    ---
    entry_time_filt_cells=CURRENT_TIMESTAMP : timestamp                    
    """


@schema
class PairwDist(dj.Computed):
    definition = """
    # Pairwise distance statistic
    -> FilteredSessions
    -> FilteredCellsParams
    -> PairwDistParams
    -> vmod0.ProjectionCorr
    """

    class Cells(dj.Part):
        definition = """
        # Cell numbers
        -> PairwDist
        region               : char(3)                      # Brain region (3 letter abbreviation)
        ---
        n_all                : smallint                     # How many cells in total were considered?
        n_startr             : smallint                     # How many starter cells was the average calculated over (population A)?
        """

    class NN(dj.Part):
        definition = """
        # Nearest neighbour (NN) distance results per region
        -> PairwDist
        region               : char(3)                      # Brain region (3 letter abbreviation)
        ---
        mean_nn              : longblob                     # Mean NN over 1-10 neighbours
        mean_nn_shuff_all    : longblob                     # Shuffled mean NN over 1-10 neighbours taking all cells as start population
        mean_nn_shuff_ref=null : longblob                   # Shuffled mean NN over 1-10 neighbours taking only reference cells as start population
        mean_nn_csr          : longblob                     # Shuffled mean NN over 1-10 neighbours taking CSR as start population




        """

    class PairwD(dj.Part):
        definition = """
        # Pairw. distance results per region
        -> PairwDist
        region               : char(3)                      # Brain region (3 letter abbreviation)
        ---
        med_pairw_dist       : double                       # Median pairwise distance
        mean_pairw_dist=null : double                       # Mean pairwise distance
        med_pairw_dist_shuffall : double                    # Shuffled median pairwise distance taking all cells as start population
        mean_pairw_dist_shuffall : double                   # Shuffled mean pairwise distance taking all cells as start population
        med_pairw_dist_shuffref=null : double               # Shuffled median pairwise distance taking only reference cells as start population
        mean_pairw_dist_shuffref=null : double              # Shuffled mean pairwise distance taking only reference cells as start population
        med_pairw_dist_csr   : double                       # Shuffled median pairwise distance taking CSR as start population
        mean_pairw_dist_csr  : double                       # Shuffled mean pairwise distance taking CSR as start population
        """


@schema
class AnatomicalMask(dj.Manual):
    definition = """
    # Anatomical mask identifying anatomical regions in FOV
    -> vmod0.ProjectionCorr
    ---
    -> AnatomicalMaskParams
    anat_mask            : longblob               # Anatomical mask for regions in FOV
    """


@schema
class RoisCorrBrainLoc(dj.Computed):
    definition = """
    # Cell IDs and anatomical location
    -> AnatomicalMask
    -> vmod0.RoisCorr
    """

    class MEC(dj.Part):
        definition = """
        # Cells in MEC
        -> RoisCorrBrainLoc
        """

    class PAS(dj.Part):
        definition = """
        # Cells in Parasubiculum
        -> RoisCorrBrainLoc
        """

    class PRH(dj.Part):
        definition = """
        # Cells in perirhinal cortex
        -> RoisCorrBrainLoc
        """

    class RSA(dj.Part):
        definition = """
        # Cells in Retrosplenial / Agranular cortex
        -> RoisCorrBrainLoc
        """

    class Undefined(dj.Part):
        definition = """
        # Cells elsewhere (not defined by any anatomical label)
        -> RoisCorrBrainLoc
        """


@schema
class NumberCellTypes(dj.Computed):
    definition = """
    # Number of (pure) cell types per session and brain region (chose MEC and PAS)
    -> FilteredSessions
    -> FilteredCellsParams
    ---
    n_cells              : int                          # Total number of (filtered) cells in this session
    """

    class MEC(dj.Part):
        definition = """
        # Cell numbers MEC
        -> NumberCellTypes
        ---
        mec_n_cells          : int                      # MEC Total number of (filtered) cells
        mec_ovc=null         : int                      # MEC Number of OV cells
        mec_n_grid_95        : int                      # MEC Number of grid cells > 95th
        mec_n_info_95        : int                      # MEC Number of cells w info content > 95th
        mec_n_hd_95          : int                      # MEC Number of HD cells > 95th
        mec_n_border_95      : int                      # MEC Number of Border cells > 95th
        mec_n_bv_95          : int                      # MEC Number of BV cells > 95th
        mec_n_grid_99        : int                      # MEC Number of grid cells > 99th
        mec_n_info_99        : int                      # MEC Number of cells w info content > 99th
        mec_n_hd_99          : int                      # MEC Number of HD cells > 99th
        mec_n_border_99      : int                      # MEC Number of Border cells > 99th
        mec_n_bv_99          : int                      # MEC Number of BV cells > 99th
        """

    class PAS(dj.Part):
        definition = """
        # Cell numbers PAS
        -> NumberCellTypes
        ---
        pas_n_cells          : int                          # PAS Total number of (filtered) cells
        pas_ovc=null         : int                          # PAS Number of OV cells
        pas_n_grid_95        : int                          # PAS Number of grid cells > 95th
        pas_n_info_95        : int                          # PAS Number of cells w info content > 95th
        pas_n_hd_95          : int                          # PAS Number of HD cells > 95th
        pas_n_border_95      : int                          # PAS Number of Border cells > 95th
        pas_n_bv_95          : int                          # PAS Number of BV cells > 95th
        pas_n_grid_99        : int                          # PAS Number of grid cells > 99th
        pas_n_info_99        : int                          # PAS Number of cells w info content > 99th
        pas_n_hd_99          : int                          # PAS Number of HD cells > 99th
        pas_n_border_99      : int                          # PAS Number of Border cells > 99th
        pas_n_bv_99          : int                          # PAS Number of BV cells > 99th
        """


@schema
class OVParams(dj.Lookup):
    definition = """
    # Object centered map parameters
    ov_params_id         : char(1)                     # Parameter set ID, starting with A
    ---
    bin_size_dist_ov     : float                        # Bin size for distance binning in mm
    bins_angular_ov      : int                         # Number of bins in 360 degrees
    sigma_time_ov        : float                        # 2D gaussian smoothing of occupancy
    sigma_signal_ov      : float                        # 2D guassian smoothing of binned signal
    """


@schema
class OVOccupancy(dj.Computed):
    definition = """
    # Object centered occupancy
    -> vmod0.Tracking.OpenField
    -> vmod0.SignalTrackingParams
    -> OVParams
    -> vmod0.ArenaObjectPos
    ---
    occupancy_ov         : longblob               # Smoothed 2D occupancy map [seconds], x: angles, y: distance
    mask_occ_ov          : longblob               # Mask (where time = 0), x: angles, y: distance
    occupancy_raw_ov     : longblob               # Raw, non-smoothed 2D occupancy map,  x: angles, y: distance
    explor_ratio_ov      : double                 # Exploration ratio (visited bins over all bins)
    explor_std_ov        : double                 # Exploration standard deviation (of visited bins)
    radial_edges_ov      : longblob               # Histogram edges in y (distance)
    angular_edges_ov     : longblob               # Histogram edges in x (angles)
    occupancy_time_ov    : double                 # Time in seconds in occupancy
    fraction_occupancy_ov : double                # Fraction of time in occupancy map
    """


@schema
class OVMap(dj.Computed):
    definition = """
    # Object centered ratemap (vector map)
    -> vmod0.SignalTracking
    -> OVOccupancy.proj(tracking_dataset="dataset_name")
    ---
    ovmap                : longblob               # Object centered ratemap ("vector map")
    ovmap_raw            : longblob               # Unsmoothed (raw) 2D ratemap
    mask_ovmap           : longblob               # Mask (where time = 0)
    binned_raw_ov        : longblob               # Raw, binned signal
    bin_max_ov           : longblob               # Bin with maximum signal (ovmap(bin_max) = max(ovmap))
    max_ov               : double                 # Maximum
    """


@schema
class OVCFields(dj.Computed):
    definition = """
    # Object vector cell (OVC) field based calculations
    -> vmod0.Ratemap.proj(base_session="session_name")
    -> vmod0.ShuffleParams
    ---
    object1_session      : varchar(16)                  # Object session 1
    object2_session      : varchar(16)                  # Object session 2
    """

    class Fields(dj.Part):
        definition = """
        # OVC calculated field statistics (all fields)
        -> OVCFields
        object1_field_id     : int                          # Field ID of field in object session 1
        object2_field_id     : int                          # Field ID of closest field to object1_field_id in object session 2
        ---
        dist_fields          : double                       # Euclidian distance between fields with object1_field_id and object2_field_id - object centered
        dist_to_object=null  : double                      # Distance of field from object [average of field in object session 1 and 2]
        angle_to_object=null : double                      # Angle of field to object [average of field in object session 1 and 2]
        object1_field_base=null : double                    # (Object session 1 field mean rate in object session 1) / (Field mean rate in base session)
        object2_field_base=null : double                    # (Object session 2 field mean rate in object session 2) / (Field mean rate in base session)
        """


@schema
class OVCScores(dj.Computed):
    definition = """
    # Object vector cell (OVC) vector map score based calculations
    -> vmod0.SignalTracking.proj(base_session="session_name")
    -> OVParams
    -> vmod0.ShuffleParams
    ---
    object1_session      : varchar(16)                  # Object session 1
    object2_session      : varchar(16)                  # Object session 2
    ovscore              : double                       # Object vector score (2D correlation between OV maps)
    """

    class ShuffledOVScore(dj.Part):
        definition = """
        # Shuffled Object vector (OV) score and shuffling
        -> OVCScores
        ---
        shuffled_ovscores_95perc : double             # Object vector score shuffling for cell: 95th percentile
        shuffled_ovscores_99perc : double             # Object vector score shuffling for cell: 99th percentile
        shuffled_ovscores    : longblob               # Object vector score shuffling for cell

        """


@schema
class OVCutoffs(dj.Lookup):
    definition = """
    # Object vector cell cutoffs
    ov_cutoff_id         : char(1)                      # Parameter set ID, starting with A
    ---
    info_content_cutoff  : varchar(100)                 # Information content cutoff (>)
    ovscore_cutoff       : varchar(100)                 # Object vector score cutoff (>)
    dist_fields_cutoff   : float                          # Distance [mm] (object centered field distances) (<)
    dist_to_object_cutoff : float                        # Distance [mm] to object (>)
    object1_field_base_cutoff : float                     # Relative rate of field in object session 1 compared to base session (>)
    object2_field_base_cutoff : float                     # Relative rate of field in object session 1 compared to base session (>)
    """


@schema
class OVC(dj.Computed):
    definition = """
    # Object vector cell (OVC) summary table
    -> OVCScores
    -> OVCutoffs
    -> OVCFields
    ---
    object1_session      : varchar(16)                  # Object session 1
    object2_session      : varchar(16)                  # Object session 2
    ovscore              : double                       # Object vector score (2D correlation between OV maps)
    is_ovc               : tinyint                      # 0 - not an OVC according to cutoffs, 1 - putative OVC
    no_fields            : int                           # Number of filtered fields (matching cutoff criteria)
    mean_dist_to_object=null : double                   # Average distance of (filtered) fields to object [mm]
    mean_dist_fields=null : double                       # Average distance between fields [mm]
    mean_angle_to_object=null : double                  # Circular mean of field angles to object [0, 2*pi]
    std_angle_to_object=null : double                   # Circular standard deviation for field angles to object [radians]
    field_ids=null       : longblob                      # Field IDs list of dictionaries ('object1_field_id', 'object2_field_id')
    angles_to_object=null : longblob                    # Field angles [0, 2*pi]
    dists_to_object=null : longblob                     # Distances of (filtered) fields to object [mm]
    dists_fields=null    : longblob                      # Distances of (filtered) fields to object [mm]
    """


@schema
class CutoffsInfoContent(dj.Computed):
    definition = """
    # Session level info content cutoffs
    -> FilteredSessions
    -> FilteredCellsParams
    ---
    n_shuffles           : int                          # Number of shuffles (total)
    info_content_90=null : double                      # Spatial information content 90th cutoff
    info_content_95=null : double                      # Spatial information content 95th cutoff
    info_content_99=null : double                      # Spatial information content 99th cutoff
    """


@schema
class CutoffsGridscore(dj.Computed):
    definition = """
    # Session level gridscore cutoffs
    -> FilteredSessions
    -> FilteredCellsParams
    ---
    n_shuffles           : int                          # Number of shuffles (total)
    gridscore_90=null    : double                      # GridScore 90th cutoff
    gridscore_95=null    : double                      # GridScore 95th cutoff
    gridscore_99=null    : double                      # GridScore 99th cutoff
    """


@schema
class CutoffsOVScore(dj.Computed):
    definition = """
    # Session level OV score cutoffs
    -> FilteredSessions.proj(base_session="session_name")
    -> FilteredCellsParams
    ---
    n_shuffles           : int                          # Number of shuffles (total)
    ovscore_90=null      : double                      # OV Score 90th cutoff
    ovscore_95=null      : double                      # OV Score 95th cutoff
    ovscore_99=null      : double                      # OV Score 99th cutoff
    """


@schema
class CutoffsMVL(dj.Computed):
    definition = """
    # Session level MVL (head direction tuning) cutoffs
    -> FilteredSessions
    -> FilteredCellsParams
    ---
    n_shuffles           : int                          # Number of shuffles (total)
    mvl_90=null          : double                      # MVL 90th cutoff
    mvl_95=null          : double                      # MVL 95th cutoff
    mvl_99=null          : double                      # MVL 99th cutoff
    """


@schema
class CutoffsBorderScore(dj.Computed):
    definition = """
    # Session level Borderscore cutoffs
    -> FilteredSessions
    -> FilteredCellsParams
    ---
    n_shuffles           : int                          # Number of shuffles (total)
    borderscore_90=null  : double                      # Borderscore 90th cutoff
    borderscore_95=null  : double                      # Borderscore 95th cutoff
    borderscore_99=null  : double                      # Borderscore 99th cutoff
    """


@schema
class CutoffsBVScore(dj.Computed):
    definition = """
    # Session level boundary vector score cutoffs
    -> FilteredSessions
    -> FilteredCellsParams
    ---
    n_shuffles           : int                          # Number of shuffles (total)
    bvs_90=null          : double                      # BVS 90th cutoff
    bvs_95=null          : double                      # BVS 95th cutoff
    bvs_99=null          : double                      # BVS 99th cutoff
    """


@schema
class NNeighbourInterIntra(dj.Computed):
    definition = """
    # NN scores inter vs. intra score
    -> FilteredSessions
    -> FilteredCellsParams
    -> PairwDistParams.proj(pairwise_dist_param_A="pairwise_dist_param")
    -> PairwDistParams.proj(pairwise_dist_param_B="pairwise_dist_param")
    -> vmod0.ProjectionCorr
    """
    
    class Cells(dj.Part):
        definition = """
        # Cell numbers
        -> NNeighbourInterIntra
        region               : char(3)                      # Brain region (3 letter abbreviation)
        ---
        n_all                : smallint                     # How many cells in total were considered?
        n_startr_a           : smallint                     # How many starter cells was the average calculated over (population A)?
        n_startr_b           : smallint                     # How many starter cells was the average calculated over (population B)?
        """
        
    class DistAll(dj.Part):
        definition = """
        # NN distance results per region all cells - raw results
        -> NNeighbourInterIntra
        region               : char(3)                      
        ---
        nns_ab               : longblob              
        nns_ba               : longblob              
        nns_ab_shuffab       : longblob              
        nns_ba_shuffab       : longblob              
        nns_ab_shuffall      : longblob              
        nns_ba_shuffall      : longblob              
        nns_ab_csr           : longblob              
        nns_ba_csr           : longblob              
        """

    class DistSub(dj.Part):
        definition = """
        # NN distance results per region subsampled populations (size matched) - raw results
        -> NNeighbourInterIntra
        region               : char(3)                      
        ---
        nns_ab               : longblob              
        nns_ba               : longblob              
        nns_ab_shuffab       : longblob              
        nns_ba_shuffab       : longblob              
        nns_ab_shuffall      : longblob              
        nns_ba_shuffall      : longblob              
        nns_ab_csr           : longblob              
        nns_ba_csr           : longblob              
        """
        
    class NNAll(dj.Part):
        definition = """
        # Mean NN distance results per region all cells
        -> NNeighbourInterIntra
        region               : char(3)                      
        ---
        nns_ab               : double                       
        nns_ba               : double                       
        nns_ab_shuffab       : double                       
        nns_ba_shuffab       : double                       
        nns_ab_shuffall      : double                       
        nns_ba_shuffall      : double                       
        nns_ab_csr           : double                       
        nns_ba_csr           : double                       
        """
        
    class NNSub(dj.Part):
        definition = """
        # Mean NN distance results per region subsampled populations (size matched)
        -> NNeighbourInterIntra
        region               : char(3)                      
        ---
        nns_ab               : double                       
        nns_ba               : double                       
        nns_ab_shuffab       : double                       
        nns_ba_shuffab       : double                       
        nns_ab_shuffall      : double                       
        nns_ba_shuffall      : double                       
        nns_ab_csr           : double                       
        nns_ba_csr           : double                       
        """
    class RatioAll(dj.Part):
        definition = """
        # Inter to intra NN distances per region all cells
        -> NNeighbourInterIntra
        region               : char(3)                      
        ---
        ratio_ab             : double                       
        ratio_ba             : double                       
        ratio_ab_shuffab     : double                       
        ratio_ba_shuffab     : double                       
        ratio_ab_shuffall    : double                       
        ratio_ba_shuffall    : double                       
        ratio_ab_csr         : double                       
        ratio_ba_csr         : double                       
        """
        
    class RatioSub(dj.Part):
        definition = """
        # Inter to intra NN distances per region subsampled populations (size matched)
        -> NNeighbourInterIntra
        region               : char(3)                      
        ---
        ratio_ab             : double                       
        ratio_ba             : double                       
        ratio_ab_shuffab     : double                       
        ratio_ba_shuffab     : double                       
        ratio_ab_shuffall    : double                       
        ratio_ba_shuffall    : double                       
        ratio_ab_csr         : double                       
        ratio_ba_csr         : double                       
        """
        
    class Unprocessed(dj.Part):
        definition = """
        # Inter to intra NN distances per region (and cell numbers)
        -> NNeighbourInterIntra
        """
