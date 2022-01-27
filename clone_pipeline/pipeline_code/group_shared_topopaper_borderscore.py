import datajoint as dj


# -------------- group_shared_topopaper_borderscore --------------


schema = dj.Schema('group_shared_topopaper_borderscore')

vmod0 = dj.VirtualModule('vmod0', 'group_shared_topopaper_main_imaging')


@schema
class BVFieldParams(dj.Lookup):
    definition = """
    # Boundary vector score (BVS) field detection params
    bvfield_params_id    : char(1)                      # Parameter set ID, starting with A
    ---
    std_detect           : float                        # Number of standard deviations over median for field detection
    std_include          : float                        # Number of standard deviations over median for field inclusion
    min_bin              : smallint                    # Minimum no of bins for field inclusion
    """


@schema
class BVField(dj.Computed):
    definition = """
    # Boundary vector score (BVS) fields
    -> vmod0.Ratemap
    -> BVFieldParams
    ---
    no_fields=null       : smallint                     # Total number of detected fields
    fields_map=null      : longblob                     # Map of detected field
    """

    class Fields(dj.Part):
        definition = """
        # Field detection results
        -> BVField
        field_no             : int                          # Field number
        ---
        field_coords         : longblob                     # Coordinates of all bins in the firing field
        field_centroid_x     : double                       # Field centroid x coordinate
        field_centroid_y     : double                       # Field centroid y coordinate
        field_area           : int                          # Area in number of bins
        field_bbox           : longblob                     # Field bounding box
        """


@schema
class BVScoreParams(dj.Lookup):
    definition = """
    # Boundary vector score (BVS) analysis params
    bvscore_params_id    : char(1)                        # Parameter set ID, starting with A
    ---
    r_factor             : float                           # r-factor (weighing contribution of extra fields in diminishing score)
    barwidth_max         : tinyint                        # Maximum bar width
    """


@schema
class BVScoreFieldMethod(dj.Lookup):
    definition = """
    # Field detection method for calculating boundary vector score (bvs)
    bv_field_dect_method : enum('opexebo','bvs')           # Specifies how fields were extracted
    """


@schema
class BVScore(dj.Computed):
    definition = """
    # Boundary vector score (BVS)
    -> BVField
    -> BVScoreParams
    -> BVScoreFieldMethod
    ---
    bvs=null             : double                       # Boundary vector score (BVS)
    orientation=null     : enum('vertical','horizontal') 
    """

    class ScoreX(dj.Part):
        definition = """
        # Score x (bars spanning X)
        -> BVScore
        ---
        score_x              : double                       # Maximum score for bars spanning x (horizontal bars)
        bar_width            : tinyint                      # Barwidth at maximum
        ypos                 : smallint                     # Bar position at maximum (center of bar)
        ypos_rel             : float                         # relative bar position at maximum (center of bar)
        bar_map              : longblob                     # barMap (streak of ones) at maximum
        """

    class ScoreY(dj.Part):
        definition = """
        # Score y (bars spanning Y)
        -> BVScore
        ---
        score_y              : double                       # Maximum score for bars spanning y (vertical bars)
        bar_width            : tinyint                      # Barwidth at maximum
        xpos                 : smallint                     # Bar position at maximum (center of bar)
        xpos_rel             : float                         # relative bar position at maximum (center of bar)
        bar_map              : longblob                     # barMap (streak of ones) at maximum
        """


@schema
class ShuffledBVS(dj.Computed):
    definition = """
    # Shuffling table for boundary vector score (BVS)
    -> vmod0.FilteredSpikes.proj(signal_dataset="dataset_name")
    -> vmod0.Tracking.OpenField.proj(tracking_dataset="dataset_name")
    -> vmod0.ShuffleParams
    -> vmod0.SignalTrackingParams
    -> vmod0.MapParams
    -> BVFieldParams
    -> BVScoreParams
    -> BVScoreFieldMethod
    ---
    -> [nullable] vmod0.Sync.proj(sync_dataset_frames_imaging="dataset_name",sync_name_frames_imaging="sync_name")
    number_shuffles      : int                          # Total number of shuffles (can vary from expected number)
    shuffling_offsets    : longblob                     # Shuffling offset
    """

    class BVS(dj.Part):
        definition = """
        # Shuffled bvs
        -> ShuffledBVS
        ---
        bvs_99               : double                  # BVS 99th percentile
        bvs_95               : double                  # BVS 95th percentile
        bvs_shuffles          : longblob                # Individual shuffles BVS
        """


