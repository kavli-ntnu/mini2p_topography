import datajoint as dj


# -------------- group_shared_topopaper_main_imaging --------------


schema = dj.Schema('group_shared_topopaper_main_imaging')

vmod0 = dj.VirtualModule('vmod0', 'group_shared_topopaper_mlims')


@schema
class Suite2pyOps(dj.Manual):
    definition = """
    # Suite2P Python options (ops)
    ops_name             : varchar(100)                 # Name for option set
    ---
    ops                  : longblob                     # Suite2P python ops
    entry_time=CURRENT_TIMESTAMP : timestamp            # Auto created timestamp
    """


@schema
class Scope(dj.Manual):
    definition = """
    scope_name           : varchar(50)                  # Scope name
    ---
    entry_time_scope=CURRENT_TIMESTAMP : timestamp                    # Auto created timestamp
    scope_notes=null     : varchar(5000)                
    """

    class Picture(dj.Part):
        definition = """
        # Photos of the scope body
        scope_photo_id       : int unsigned auto_increment  # Photo ID
        -> Scope
        ---
        scope_photo          : longblob               # Photo of the scope body
        entry_time_scope_photo=CURRENT_TIMESTAMP : timestamp                    # Auto created timestamp
        """


@schema
class CueCard(dj.Lookup):
    definition = """
    cue_card             : varchar(32)                  
    ---
    card_content         : varchar(100)                 
    """


@schema
class SignalTrackingParams(dj.Lookup):
    definition = """
    # Parameters for SignalTracking and Occupancy tables
    s_t_params_id        : char(1)                      # Parameter set ID, starting with A
    ---
    speed_cutoff_low     : float                        # Speed filter low cutoff in mm/s
    speed_cutoff_high    : float                        # Speed filter high cutoff in mm/s
    time_offset          : float                        # Temporal offset of signal (spikes/calcium) vs. tracking in seconds
    s_t_params_hash=""   : varchar(32)                  
    """


@schema
class User(dj.Manual):
    definition = """
    username             : varchar(255)                 # NTNU username
    ---
    email                : varchar(255)                 # email address
    first_name           : varchar(255)                 # first name
    last_name            : varchar(255)                 # last name
    entry_time=CURRENT_TIMESTAMP : timestamp                    # current timestamp
    is_active            : tinyint                      # active
    is_staff             : tinyint                      # staff status
    is_superuser         : tinyint                      # superuser status
    """


@schema
class ArenaObject(dj.Lookup):
    definition = """
    obj_name             : varchar(32)                  
    ---
    obj_geometry         : varchar(16)                  # e.g. cube, cylinder
    obj_width            : float                        # Object width in cm
    obj_length           : float                        # Object length in cm
    obj_height           : float                        # Object height in cm
    obj_desc=""          : varchar(500)                 # Object description
    """

    class Picture(dj.Part):
        definition = """
        # Photos of the object
        obj_photo_id         : int unsigned auto_increment  # Photo ID
        -> ArenaObject
        ---
        obj_photo            : longblob               # Photo of the object
        entry_time_obj_photo=CURRENT_TIMESTAMP : timestamp                    # Auto created timestamp
        """


@schema
class ShuffleParams(dj.Lookup):
    definition = """
    # Shuffling table parameters
    shuffle_params_id    : char(1)                      # Parameter set ID, starting with A
    ---
    number_shuffles      : int                          # Expected number of shuffles
    margin_seconds       : float                        # Margin in seconds at start, end of sync samples
    break_seconds        : float                        # Break in seconds between sync samples (-> non-continuous sampling of sync stream)
    shuffle_params_hash="" : varchar(32)                  
    """


@schema
class NoiseCalcParams(dj.Lookup):
    definition = """
    # Noise Calculation parameters
    noise_calc_id        : char(1)                      # Parameter set ID, starting with A
    ---
    seconds_back         : float                        # How many seconds to go back from each extracted spike for noisy episode extraction
    seconds_forw         : float                        # How many seconds to go forward from each extracted spike for noisy episode extraction
    noise_calc_params_hash="" : varchar(32)                  
    """


@schema
class FileType(dj.Lookup):
    definition = """
    filetype             : varchar(50)                  # Name for filetype
    ---
    description          : varchar(500)                 # Short description of filetype
    """


@schema
class Repository(dj.Lookup):
    definition = """
    # A physical location for storing files
    repository_name      : varchar(100)                 # Short name for repository
    ---
    repository_type      : enum('Local','NTNUNetwork','Amazon') # Local or network
    path                 : varchar(1000)                # Path without trailing slash
    """


@schema
class DLCTrackingType(dj.Lookup):
    definition = """
    dlc_tracking_type    : varchar(32)                  
    """


@schema
class DLCModel(dj.Manual):
    definition = """
    dlc_model            : varchar(32)                  # lab-friendly model name
    ---
    dlc_task             : varchar(32)                  
    dlc_date             : varchar(16)                  
    dlc_iteration        : int                          # iteration/version of this model
    dlc_snapshotindex    : int                          # which snapshot index used for prediction (if -1 then use the latest snapshot)
    dlc_shuffle          : int                          # which shuffle of the training dataset used for training the network (typically 1)
    dlc_trainingsetindex : int                          # which training set fraction used to generate the model (typically 0)
    dlc_scorer           : varchar(128)                 # scorer/network name for a particular shuffle, training fraction etc.
    dlc_cfg_template     : longblob                     # dictionary of the config yaml needed to run the deeplabcut.analyze_videos()
    dlc_model_description="" : varchar(1000)               
    """

    class ModelPath(dj.Part):
        definition = """
        -> DLCModel
        ---
        dlc_project_path     : varchar(255)                 # relative path of the DLC project path, used for the project_path var in the config.yaml
        -> Repository
        """


@schema
class DLCTrackingProcessingMethod(dj.Lookup):
    definition = """
    dlc_tracking_processing_method : varchar(16)                  # e.g 2points_leftrightear
    ---
    method_description="" : varchar(1000)                
    dlc_tracking_processing_params : longblob                     
    -> DLCTrackingType.proj(applicable_tracking_type="dlc_tracking_type")
    function_to_invoke   : varchar(32)                  # DLC processing method in the "helpers.dlc_tracking" module
    param_hash           : uuid                        
    """


@schema
class DLCProcessingMethod(dj.Manual):
    definition = """
    -> DLCModel
    -> DLCTrackingProcessingMethod
    ---
    method_desc=""       : varchar(1000)                # description for this model-method combination
    """


@schema
class GridParams(dj.Lookup):
    definition = """
    # Grid analysis parameters
    grid_params_id       : char(1)                      # Parameter set ID, starting with A
    ---
    min_orientation      : int                          # Minimum orientation between adjacent fields in autocorr [degrees]
    grid_params_hash=""  : varchar(32)                  
    """


@schema
class FieldParams(dj.Lookup):
    definition = """
    # Field detection parameters
    field_params_id      : char(1)                      # Parameter set ID, starting with A
    ---
    init_thresh          : float                        # Initial threshold for field detection
    min_bins             : int                          # Minimum number of bins in fields
    fraction_min_peak    : float                        # min_peak (as fraction of ratemap maximum)
    fraction_min_mean    : float                        # min_mean (as fraction of ratemap maximum)
    search_method        : varchar(10)                  # Peak initialisation method ('default','sep', ...)
    field_params_hash="" : varchar(32)                  
    """


@schema
class ImplantSpec(dj.Lookup):
    definition = """
    # Implant specifications / properties
    implant_name         : varchar(100)                 # Short name of implant
    ---
    fov_scaler_x         : double                       # Additional scaler in x that ImagingFOV() calibration does not account for
    fov_scaler_y         : double                       # Additional scaler in y that ImagingFOV() calibration does not account for
    fov_scaler           : double                       # Average of scaler in x and y
    implant_desc         : varchar(1000)                # Description of implant
    """


@schema
class Implant(dj.Manual):
    definition = """
    # Implants
    -> vmod0.Animal
    -> ImplantSpec
    brain_region         : char(3)                      # Desired implant location - brain structure (e.g. 'ENT')
    hemisphere           : enum('left','right')         # Hemisphere (left/right)
    ---
    ml_reference=null    : varchar(100)                 # Anatomical reference structure for medio-lateral positioning (e.g. 'sagittal sinus')
    implant_ml=null      : float                        # Desired implant location (midpoint) medio-lateral [mm]
    dv_reference=null    : varchar(100)                 # Anatomical reference structure for dorso-ventral positioning (e.g. 'brain surface')
    implant_dv=null      : float                        # Desired implant location (midpoint) dorso-ventral [mm]
    ap_reference=null    : varchar(100)                 # Anatomical reference structure for anterior-posterior positionig (e.g. 'transverse sinus')
    implant_ap=null      : float                        # Desired implant location (midpoint) anterior-posterior [mm]
    angle_ap=null        : float                        # Angle of implant in anterior-posterior direction
    angle_ap_dir=null    : enum('forwards','backwards') # In the anterior-posterior direction: Where did the implant point towards - Forwards or backwards (mouse centered coordinates)?
    angle_ml=null        : float                        # Angle of implant in medio-lateral direction
    angle_ml_dir=null    : enum('left','right')         # In the medio-lateral direction: Where did the implant point towards - left or right (mouse centered coordinates)?
    implant_comment=null : varchar(1000)                # Space for comments about the implant / surgery
    """


@schema
class RateParams(dj.Lookup):
    definition = """
    # Smoothed spike rates parameters
    rate_id              : char(1)                      # Parameter set ID, starting with A
    ---
    spike_rate_sigma     : float                        # Gaussian smoothing sigma for filtered spikes
    comment=""           : varchar(255)                 # Comment
    rate_params_hash=""  : varchar(32)                  
    """


@schema
class ApparatusCategory(dj.Lookup):
    definition = """
    category             : varchar(40)                  # Arena / Apparatus category (openfield, linear track, etc.)
    """


@schema
class Apparatus(dj.Lookup):
    definition = """
    apparatus            : varchar(32)                  
    -> ApparatusCategory
    ---
    apparatus_note=""    : varchar(1000)                
    """

    class Geometry(dj.Part):
        definition = """
        -> Apparatus
        ---
        arena_geometry       : enum('square','rectangle','circle','linear','treadmill') 
        arena_x_dim=null     : float                        # Arena x extend in mm
        arena_y_dim=null     : float                        # Arena y extend in mm
        arena_radius=null    : float                        # Arena radius (for circular environments) in mm
        arena_wall_height=null : float                        # Arena wall height in mm
        """

    class Treadmill(dj.Part):
        definition = """
        -> Apparatus
        ---
        treadmill_radius     : float                        # Treadmill radius in mm
        treadmill_width      : float                        # Treadmill belt/wheel width in mm
        """


@schema
class DatasetType(dj.Lookup):
    definition = """
    datasettype          : varchar(50)                  # Name for Dataset Type
    ---
    description          : varchar(500)                 # Short description of Dataset Type
    """


@schema
class Dataset(dj.Manual):
    definition = """
    # Dataset
    dataset_name         : char(16)                     # 16 character hash
    ---
    -> DatasetType
    """


@schema
class SessionType(dj.Lookup):
    definition = """
    sessiontype          : varchar(50)                  # Name for Session Type
    ---
    description          : varchar(500)                 # Short description of session type
    """


@schema
class TrackedBodyPart(dj.Lookup):
    definition = """
    body_part            : varchar(32)                  
    """


@schema
class PhysicalFile(dj.Manual):
    definition = """
    # Single file
    -> Dataset
    -> FileType
    order                : int                          # When among the other files in the dataset?
    ---
    -> Repository
    file_path            : varchar(1000)                # File path on repository
    """


@schema
class SpikeFilterParams(dj.Lookup):
    definition = """
    # Spike filtering parameters
    spike_filter_id      : char(1)                      # Parameter set ID, starting with A
    ---
    spike_std            : float                        # Number of standard deviations over mean for spike threshold calculation.
    comment=""           : varchar(255)                 # Comment
    spike_filter_params_hash="" : varchar(32)                  
    """


@schema
class MapParams(dj.Lookup):
    definition = """
    # Ratemap / Occupancy map parameters
    map_params_id        : char(1)                      # Parameter set ID, starting with A
    ---
    bin_size             : float                        # Bin size in mm
    sigma_time           : float                        # 2D gaussian smoothing of occupancy
    sigma_signal         : float                        # 2D guassian smoothing of binned signal
    map_params_hash=""   : varchar(32)                  
    """


@schema
class AngularParams(dj.Lookup):
    definition = """
    # Angular rate / occupancy parameters
    ang_params_id        : char(1)                      # Parameter set ID, starting with A
    ---
    bins_angular         : int                          # Number of bins in 360 deg
    sigma_time_ang       : float                        # 1D gaussian smoothing of occupancy
    sigma_signal_ang     : float                        # 1D gaussian smoothing of signal
    ang_params_hash=""   : varchar(32)                  
    """


@schema
class ExperimentType(dj.Lookup):
    definition = """
    # Experiment category that defines what to do with the data
    experiment_type      : varchar(100)                 
    """

@schema
class Setup(dj.Manual):
    definition = """
    # Imaging setup parameters
    setup_name           : varchar(50)                  # Setup name
    setup_timestamp      : datetime                     # Timestamp of setup configuration
    ---
    -> ExperimentType
    """

    class Labview(dj.Part):
        definition = """
        # Labview hardware parameters
        -> Setup
        ---
        labview_config       : longblob               # LabVIEW config file
        """

    class Scanimage(dj.Part):
        definition = """
        # Scanimage hardware parameters
        -> Setup
        ---
        machine_file         : longblob               # ScanImage config file
        """

    class Sync(dj.Part):
        definition = """
        # Synchronisation hardware parameters
        -> Setup
        sync_name            : varchar(50)                  # Name of sync stream
        ---
        generic_name         : varchar(50)                  # What is it called in e.g. Sync? E.g. '2P frames'
        polarity=null        : enum('positive','negative')  # Rising / Falling edge?
        master=null          : varchar(50)                  # Name of master clock (generic name)
        tolerance=null       : float                         # How many seconds difference are allowed between sync stream ends
        """

    class Tracking1D(dj.Part):
        definition = """
        # Linear tracking hardware parameters
        -> Setup
        ---
        encoder_steps=null   : int                          # Steps for full turn rotary encoder
        """

    class Tracking2D(dj.Part):
        definition = """
        # 2D tracking hardware parameters
        -> Setup
        ---
        header               : smallint                     # How many header lines in csv?
        timebase=null        : double                       # If time column is in samples what is the timebase?
        time                 : varchar(20)                  # Name of the time column
        led1_x               : varchar(20)                  # Name of the x coordinate of LED1
        led1_y               : varchar(20)                  # Name of the y coordinate of LED1
        led2_x               : varchar(20)                  # Name of the x coordinate of LED2
        led2_y               : varchar(20)                  # Name of the y coordinate of LED2
        """

    class Treadmill(dj.Part):
        definition = """
        -> Setup
        ---
        calibration          : float                         # Conversion factor steps to mm (1: No calibration supplied)
        """


@schema
class MetaSession(dj.Manual):
    definition = """
    metasession_name     : char(16)                         # Meta session name (hash)
    ---
    entry_time=CURRENT_TIMESTAMP : timestamp                # Auto created timestamp
    """

    class Note(dj.Part):
        definition = """
        -> MetaSession
        note_hash            : char(16)                     # Hash of note content
        ---
        note                 : varchar(5000)                # Note
        entry_time_note=CURRENT_TIMESTAMP : timestamp       # Auto created timestamp
        """

    class Scope(dj.Part):
        definition = """
        -> MetaSession
        ---
        -> Scope
        """

    class Suite2pyOps(dj.Part):
        definition = """
        -> MetaSession
        ---
        -> Suite2pyOps
        """

    class Setup(dj.Part):
        definition = """
        -> MetaSession
        ---
        -> Setup
        """


@schema
class ImagingFOVRaw(dj.Manual):
    definition = """
    # Raw imaging FOV parameters (uncalibrated and uncorrected)
    -> Setup
    -> Scope
    zoom                 : decimal(5,3)                   # Zoom level
    ---
    grid_spacing_microns : float                           # Spacing in microns of the imaged grid
    grid_picture=null    : longblob                       # Picture of the scanned grid
    original_points      : longblob                       # Coordinates of points drawn by user on image (intersections in grid)
    transformed_points   : longblob                       # Coordinates of transformed points (what it should look like after correction)
    entry_time_imaging_fov=CURRENT_TIMESTAMP : timestamp  # Auto created timestamp
    """

    class Beads(dj.Part):
        definition = """
        # Photos of beads
        beads_picture_id     : int unsigned auto_increment  # Beads picture ID
        -> ImagingFOVRaw
        ---
        beads_size           : float                        # Size of beads in microns
        beads_picture=null   : longblob               # Picture of scanned beads
        """


@schema
class ImagingFOV(dj.Computed):
    definition = """
    # Corrected and calibrated imaging FOV parameters
    -> ImagingFOVRaw
    ---
    calibration_fov      : double                       # Calibration px to microns
    grid_spacing_px      : double                       # Grid spacing in pixels (derived from transformed_points)
    grid_picture_unwarped=null : longblob               # Picture of the scanned grid
    affine_transform=null : longblob                     # Affine transform pickle dump (skimage. transform. PiecewiseAffineTransform)
    """


@schema
class Intervention(dj.Lookup):
    definition = """
    intervention      : varchar(50)               # Name for Intervention
    ---
    description      : varchar(500)               # Short description of intervention type
    """


@schema
class Session(dj.Manual):
    definition = """
    -> MetaSession
    session_order        : int                          # Order of session within meta sessions (zero index!)
    session_name         : char(16)                     # Session name: Hash of animal_id, datasource_id, timestamp and combined 'yes'/'no' label
    ---
    -> vmod0.Animal
    animal_name          : varchar(255)                 # Animal name in mlims
    timestamp            : timestamp                    
    combined             : enum('yes','no','-')         
    timeseries_name=""   : varchar(255)                 # Timeseries name [e.g. MUnit_0]
    -> ExperimentType
    -> User
    """

    class Apparatus(dj.Part):
        definition = """
        -> Session
        -> Apparatus
        """

    class Data(dj.Part):
        definition = """
        # Raw Dataset
        -> Session
        -> Dataset
        """

    class Intervention(dj.Part):
        definition = """
        -> Session
        -> Intervention
        """

    class SessionType(dj.Part):
        definition = """
        -> Session
        -> SessionType
        ---
        note=""              : varchar(5000)                # Note for this session
        """

    class ArenaObject(dj.Part):
        definition = """
        -> Session.Apparatus
        -> ArenaObject
        obj_x_coord          : decimal(6,2)                 # x coordinate of object in pixels
        obj_y_coord          : decimal(6,2)                 # y coordinate of object in pixels
        ---
        obj_note=""          : varchar(1000)                
        """

    class CueCard(dj.Part):
        definition = """
        -> Session.Apparatus
        -> CueCard
        card_position        : enum('north','east','south','west') # Cue card position (in tracked picture)
        ---
        card_note=""         : varchar(1000)                
        """


@schema
class SessionDLC(dj.Manual):
    definition = """
    -> Session
    ---
    -> DLCProcessingMethod
    """


@schema
class TrackingRaw(dj.Imported):
    definition = """
    # Tracking data
    -> Session.Data
    ---
    start_time           : datetime                       # Registered system start time
    stop_time            : datetime                       # Registered system stop time
    timestamps=null      : longblob                       # Timestamps in seconds - synced timestamps
    timestamps_sys=null  : longblob                       # Timestamps in seconds - explicitly derived from system clock
    sample_rate          : double                         # Sampling rate
    sync                 : smallint                       # Was a sync signal considered? (0=no / 1=yes)
    is_calibrated=0      : tinyint                        # Was calibration (conversion px -> mm) applied?
    entry_time_tracking_raw=CURRENT_TIMESTAMP : timestamp # Auto created timestamp
    """

    class DLCPart(dj.Part):
        definition = """
        # Deeplabcut tracking data
        -> TrackingRaw
        -> TrackedBodyPart
        ---
        -> DLCModel
        dlc_training_iteration : int                          # DLC training iteration of the model used for prediction (snapshot number)
        bodypart_x_pos       : longblob               
        bodypart_y_pos       : longblob               
        bodypart_likelihood  : longblob               
        """

    class Linear(dj.Part):
        definition = """
        # Linear tracking data - constrained movement
        -> TrackingRaw
        ---
        linear_type          : enum('optical','rotary')     # Type of linear tracking setup
        file_version=null    : varchar(5)                   # File version for 'rotary' type (not applicable for 'optical')
        pos                  : longblob               # Absolute position zeroed to first sample (steps)
        rel_pos              : longblob               # Relative position change from one frame to the next (steps)
        speed                : longblob               # Speed (steps/s)
        lap_sensor=null      : smallint                     # Was a lap sensor in place? (0=no / 1=yes)
        lap_pos=null         : longblob               # Absolute position reset each lap (steps)
        lap_indices=null     : longblob               # Beam break indices (Lap start indices)
        lap_indices_var=null : float                        # Variance of lap start positions
        motor_indices=null   : longblob               # Episodes (indices) during which motor was engaged
        motor_onset_indices=null : longblob               # Indices when motor was started
        time_offset=null     : double                       # Time offset in seconds (first sync pulse)
        polarity=null        : smallint                     # Is position increasing (1=positive) or decreasing (-1=negative)
        """

    class OpenField(dj.Part):
        definition = """
        # LED tracking data - unconstrained movement
        -> TrackingRaw
        ---
        x_pos                : longblob               # Tracking x position (px)
        y_pos                : longblob               # Tracking y position (px)
        boundary_left        : double                       # Extracted boundary towards left (px)
        boundary_right       : double                       # Extracted boundary towards right (px)
        boundary_bottom      : double                       # Extracted boundary towards bottom (px)
        boundary_top         : double                       # Extracted boundary towards top (px)
        speed                : longblob               # Speed (px/s)
        head_angle           : longblob               # Head angle [0,2*pi]
        """


@schema
class Tif(dj.Imported):
    definition = """
    # Tif stack info
    -> Session.Data
    ---
    motion_cor           : int                          # Motion corrected? 1-yes / 0-no
    num_frames           : int                          # Number of frames (sum across channels)
    width                : int                          # Image width in pixels
    height               : int                          # Image height in pixels
    datatype             : varchar(50)                  # Datatype
    config               : varchar(50)                  # Configuration
    bigtiff              : int                          # Is it a bigtiff? 1-yes,0-no
    mean_image           : longblob               # Mean image of at most first 50 images
    """

    class SI(dj.Part):
        definition = """
        # Scanimage tif stack info
        -> Tif
        ---
        num_frames_scanimage : int                          # Number of frames (channels accounted for!)
        num_channels         : int                          # Number of channels in SI tif
        num_fields           : int                          # Number of fields in SI tif
        num_scanning_depths=1 : int                          # Number of scanning depths (planes) in this SI tif
        is_bidirectional     : int                          # Bidirection scan?
        scanner_frequency    : double                       # Scanner frequency
        seconds_per_line     : double                       # Seconds per line
        framerate            : double                       # Frame rate
        volume_rate          : double                       # Volume scan rate
        spatial_fill_fraction : double                       # Spatial fill fraction
        temporal_fill_fraction : double                       # Temporal fill fraction
        scanner_type         : varchar(50)                  # Scanner type
        width_scanimage      : int                          # Image width in pixels
        height_scanimage     : int                          # Image height in pixels
        width_scanimage_micro : double                       # Image width in microns
        height_scanimage_micro : double                       # Image height in microns
        zoom                 : double                       # Scanning zoom level
        version              : varchar(10)                  # Scanimage version
        header=null          : longblob               # Scanimage header information
        epoch=CURRENT_TIMESTAMP : timestamp                    # Start time
        power                : longblob               # Beam percentage
        lines_per_frame      : int                          # Lines per frame
        pixels_per_line      : int                          # Pixels per line
        scan_pixel_time_mean : double                       # Average scan pixel time
        filenames=null       : longblob               # List of filenames padded with '+++'
        timestamps_sys=null  : longblob               # List of system timestamps for every frame
        aux_trigger_0=null   : longblob               # List of AUX Trigger signal 0
        aux_trigger_1=null   : longblob               # List of AUX Trigger signal 1
        aux_trigger_2=null   : longblob               # List of AUX Trigger signal 2
        aux_trigger_3=null   : longblob               # List of AUX Trigger signal 3
        frame_numbers=null   : longblob               # (Sequential) frame numbers
        """


@schema
class Sync(dj.Imported):
    definition = """
    # Synchronisation data
    -> Session.Data
    sync_name            : varchar(50)                  # Name of sync stream
    ---
    num_indices=null     : int unsigned                 # Total number of indices
    sync_data=null       : longblob               # Indices over all subsessions
    sample_rate=null     : double                       # Sample rate
    lengths=null         : longblob               # Sample lengths of all subsessions
    sample_boundaries=null : longblob               # Array of sample boundaries for each subsession
    timestamp=CURRENT_TIMESTAMP : timestamp                    # Timestamp of acquisition
    entry_time=CURRENT_TIMESTAMP : timestamp                    # Auto created timestamp
    """

    class Mesc(dj.Part):
        definition = """
        # MESC digital I/O data
        -> Sync
        ---
        curve_name           : varchar(50)                  # Name of the curve in MESC
        curve_format_version : varchar(50)                  # Curve format version
        """

    class SI(dj.Part):
        definition = """
        # Scanimage Tif synchronisation data
        -> Sync
        ---
        si_version           : varchar(50)                  # Scanimage version
        imaging_sys_name     : varchar(255)                 # Imaging system name
        """

    class Unprocessed(dj.Part):
        definition = """
        # Unprocessed Sync Files
        -> Sync
        """

    class Ws(dj.Part):
        definition = """
        # Wavesurfer digital I/O data
        -> Sync
        ---
        di_ch_num            : int                          # Number of digital input channel in original file
        di_ch_status         : longblob               # Status (1-active,0-inactive) of all input channels in original file
        device               : varchar(50)                  # Device name (NI hardware) that was used for recording
        version              : varchar(50)                  # Wavesurfer version
        yoked                : tinyint unsigned             # Yoked to scanimage?
        """


@schema
class ImagingAnalysis(dj.Imported):
    definition = """
    # Meta table for all imaging analysis output. Supports up to suite2p v.0.10.0
    -> Session.Data
    ---
    plane=null           : int                                # Plane number (in multiplane imaging)
    order                : int                                # When among the other files in the dataset?
    filetype             : varchar(50)                         # Name for filetype
    entry_time_imaginganalysis=CURRENT_TIMESTAMP : timestamp  # Auto created timestamp
    is_tif_preproc=0     : smallint                     
    """

    class Suite2py(dj.Part):
        definition = """
        # Suite2P python parameters
        -> ImagingAnalysis
        ---
        ops_diameter=null    : int                          # Cell diameter in pixels
        ops_fs=null          : double                       # Image rate (per plane)
        ops_tau=null         : double                       # Timescale of sensor (in seconds)
        ops_nplanes=null     : int                          # Number of imaging planes in tif
        ops_nchannels=null   : int                          # Number of channels in tif
        ops_functional_chan=null : int                      # Index of functional channel (1-based)
        ops_sparse_mode=null : tinyint                      # Run sparse mode?
        ops_spatial_scale=null : tinyint                    # Spatial scale (sparse mode) 0: multi-scale
        ops_spatscale_pix=null : float                       # Spatial scale in pixels
        ops_preclassify=null : double                       # Threshold applied to filter out ROIs before signal extraction based on classifier
        ops_save_mat=null    : tinyint                      # Whether matlab export was saved
        ops_frames_per_folder=null : longblob               # Frames per folder
        ops_delete_bin=null  : int                          # Kept .bin file?
        ops_do_registration=null : int                      # Registration run?
        ops_do_phasecorr=null : int                         # Phasecorrelation for motion correction?
        ops_bidiphase=null   : float                         # User defined bidirectional (line to line) offset
        ops_do_bidiphase=null : tinyint                     # Whether to perform bidirectional offset estimation
        ops_pre_smooth=null  : float                         # Smoothing window before high-pass filtering before registration
        ops_spatial_hp=null  : float                         # Window for spatial high-pass filtering before registration
        ops_nbinned=null     : int                          # Maximum number of binned frames for cell detection
        ops_roidetect=null   : tinyint                      # Run ROI detection?
        ops_subpixel=null    : double                       
        ops_nimg_init=null   : int                          # Number of images included in motion correction reference image
        ops_batch_size=null  : int                          # Number of images processed in batch for motion correction
        ops_maxregshift=null : double                       # Maximum allowed reg shift as fraction of frame max
        ops_pad_fft=null     : tinyint                      # Pad FFT?
        ops_align_by_chan=null : int                        # Align by which cahnnel? (1-based)
        ops_reg_tif=null     : int                          # Export reg tif channel 1?
        ops_reg_tif_chan2=null : int                        # Export reg tif channel 2?
        ops_nonrigid=null    : int                          # Non-rigid motion correction applied?
        ops_block_size=null  : longblob                     # Block size for non rigid motion correction in Y and X
        ops_snr_thresh=null  : double                       # If SNR of any block (non rigid) is below this threshold it gets smmoothed until it is above
        ops_smooth_sigma=null : double                       
        ops_baseline=null    : varchar(10)                  # Baselining mode (constant, maximin, percentile, etc.)
        ops_win_baseline=null : double                      # Window for maximin
        ops_sig_baseline=null : double                      # Smoothing constant for gaussian during deconvolution
        ops_prctile_baseline=null : double                       
        ops_nsvd_for_roi=null : int                         # Maximum number of SVD components to keep
        ops_navg_frames_svd=null : int                      # Maximum number of binned frames for SVD
        ops_inner_neuropil_radius=null : float               # inner_neuropil_radius - number of pixels in donut between ROI and Neuropil
        ops_outer_neuropil_radius=null : float               # outer_neuropil_radius
        ops_tile_factor=null : float                         # Finer (>1) or coarser (<1) tiles for neuropil estimation during cell detection
        ops_allow_overlap=null : int                        # Allow fluorescence extraction from shared pixels of adjacent ROIs or exclude those pixels
        ops_min_neuropil_pixels=null : int                  # Minimum number of pixels in Neuropil
        ops_ratio_neuropil=null : double                    # Ratio of neuropil basis size to cell radius
        ops_ratio_neuropil_to_cell=null : double                       
        ops_high_pass=null   : double                       # ROI detection - running mean subtraction window length
        ops_connected=null   : int                          # Whether or not to require of ROIs to be fully connected (0 for dendrites)
        ops_max_iterations=null : int                       # Maximum number of iterations for ROI detection
        ops_smooth_masks=null : int                         # Whether to smooth masks in final step of cell detection
        ops_threshold_scaling=null : double                 # ROI detection - adjust auto threshold by this scaler
        ops_max_overlap=null : double                       # Maximum allowed overlap ROIs
        ops_neucoeff=null    : double                       # Neuropil coefficient
        ops_filelist=null    : longblob                      # Tif file list (natsorted)
        ops_mesoscan=null    : int                          
        ops_h5py_key=null    : varchar(10)                  # If hdf5 files were analyzed
        ops_chan2_thres=null : float                         # Minimum for detection of brightness on channel 2
        ops_chan2_max=null   : float                         # Maximum for NEGATIVE detection of brightness on channel 2
        ops_1preg=null       : tinyint                      # Whether to perform high-pass filtering and tapering
        ops_spatial_taper=null : float                       # Slope of taper mask at the edges
        ops_th_badframes=null : float                        # Threshold for bad frame rejection
        ops_suite2p_version=null : varchar(50)              # Suite2p python version (first version: '0.9.0')
        ops_force_sktiff=null : tinyint                     # Whether or not to use scikit-image for tiff reading
        ops_frames_include=null : int                       # Number of frames that were included in processing (-1 = all frames)
        ops_combined=null    : tinyint                      # Whether a combined output was generated 
        ops_aspect=null      : float                         # (um/pixels in X) / (um/pixels in Y) (for correct aspect ratio in GUI)
        ops_bidi_corrected=null : tinyint                   # Was bidirectional phase offset corrected? 
        ops_two_step_registration=null : tinyint            # Was two step registration process used?
        ops_smooth_sigma_time=null : float                   # Gaussian smoothing in time (Standard deviation in time frames)
        ops_spikedetect=null : tinyint                      # Run spike detection? 
        ops_tiff_list=null   : longblob                     # Tif file list (natsorted) - like 'ops_filelist' but without folder hierarchy
        ops_input_format=null : varchar(15)                 # Input file format
        ops_frames_per_file=null : longblob                  # Array of frames per file
        ops_date_proc=null   : datetime                     # Timestamp of processing time
        ops_multiplane_parallel=null : double               # Comment from s2p - whether or not to run on server
        ops_ignore_flyback=null : longblob                   # list of planes to ignore
        ops_anatomical_only=null : float                     # use cellpose masks from mean image (no functional segmentation)
        ops_cellprob_threshold=null : float                  # cellprob_threshold for cellpose (if anatomical_only > 1)
        ops_flow_threshold=null : float                       # flow_threshold for cellpose (if anatomical_only > 1)
        ops_denoise=null     : tinyint                      # denoise binned movie for cell detection in sparse_mode
        ops_soma_crop=null   : tinyint                      # crop dendrites for cell classification stats like compactness
        ops_lam_percentile=null : float                      # percentile of lambda within area to ignore when excluding cell pixels for neuropil extraction
        """


@schema
class Projection(dj.Imported):
    definition = """
    # Motion corrected stack projections
    -> ImagingAnalysis
    center_plane         : int                          # imaging plane of this projection
    ---
    mean_image           : longblob                     # Average projection of motion corrected stack
    mean_image_enhanced  : longblob                     # Average projection of motion corrected stack - contrast enhanced
    mean_image_ref       : longblob                     # Reference image for motion correction
    max_img              : longblob                     # Maximum projection of motion corrected stack
    vcorr_image          : longblob                     # Correlation map
    vmap_image           : longblob                     # Vmap
    mean_image_second    : longblob                     # Average projection of second channel
    mean_image_second_corr : longblob                   # Corrected average projection of second channel
    x_range              : longblob                     # X axis positioning of cropped images (after motion correction) in relation to original frame size
    y_range              : longblob                     # Y axis positioning of cropped images (after motion correction) in relation to original frame size
    entry_time_projection=CURRENT_TIMESTAMP : timestamp # Auto created timestam
    """

    class Motion(dj.Part):
        definition = """
        # Results of motion correction
        -> Projection
        ---
        x_mocor=null         : longblob               # X axis rigid motion offset for every frame
        y_mocor=null         : longblob               # Y axis rigid motion offset for every frame
        pxshift_pc_rigid=null : longblob              # Pixel shift over PCs - rigid
        pxshift_pc_nonrigid=null : longblob           # Pixel shift over PCs - nonrigid
        pxshift_pc_nonrigidmax=null : longblob        # Pixel shift over PCs - nonrigidmax
        max_pxshift_pc_rigid=null : double            # Maximum pixel shift over PCs - rigid
        max_pxshift_pc_nonrigid=null : double         # Maximum pixel shift over PCs - nonrigid
        max_pxshift_pc_nonrigidmax=null : double      # Maximum pixel shift over PCs - nonrigidmax
        image_pc_first=null  : longblob                # Image PCs of first images in stack [PCs x height x width ]
        image_pc_last=null   : longblob               # Image PCs of last images in stack [PCs x height x width ]
        pc_time=null         : longblob               # Shifts PC over time (magnitude over time)
        block_coords_x=null  : longblob               # X axis ranges of motion correction blocks
        block_coords_y=null  : longblob               # Y axis ranges of motion correction blocks
        block_mocor_x=null   : longblob               # Nonrigid motion correction x offsets for each block over frames
        block_mocor_y=null   : longblob               # Nonrigid motion correction y offsets for each block over frames
        zdrift               : longblob               # ZDrift over frames
        bad_frames           : longblob               # Boolean for every frame (bad=False / ok=True)
        two_step_registration=null : smallint         # Whether or not to run registration twice (for low SNR data)
        smooth_sigma_time=null : double               # How many frames to smooth in time before getting phase correlations (motion correction)
        """


@schema
class ProjectionCorr(dj.Computed):
    definition = """
    # Unwarped projections (FOVs)
    -> ImagingFOV
    -> Projection
    ---
    width_microns_full   : double                 # Full, uncropped width of corrected FOV in microns
    height_microns_full  : double                 # Full, uncropped height of corrected FOV in microns
    width_microns_eff    : double                 # Effective width of corrected FOV in microns that remains after cropping / motion correction
    height_microns_eff   : double                 # Effective height of corrected FOV in microns that remains after cropping / motion correction
    mean_image_corr      : longblob               # Average projection of motion corrected stack (unwarped)
    mean_image_enhanced_corr : longblob           # Average projection of motion corrected stack - contrast enhanced (unwarped)
    mean_image_ref_corr  : longblob               # Reference image for motion correction (unwarped)
    max_img_corr=null    : longblob               # Maximum projection of motion corrected stack (unwarped)
    vcorr_image_corr=null : longblob              # Correlation map (unwarped)
    vmap_image_corr=null : longblob               # Vmap (unwarped)
    mean_image_second_corr=null : longblob        # Average projection of second channel (unwarped)
    mean_image_second_corr_corr=null : longblob   # Corrected average projection of second channel (unwarped)
    x_range_microns_full : longblob               # Full, uncropped x coordinates of FOV in microns after unwarping
    y_range_microns_full : longblob               # Full, uncropped y coordinates of FOV in microns after unwarping
    x_range_microns_eff  : longblob               # Effective x coordinates of FOV in microns after unwarping that remains after cropping / motion correction
    y_range_microns_eff  : longblob               # Effective y coordinates of FOV in microns after unwarping that remains after cropping / motion correctio
    """


@schema
class Cell(dj.Imported):
    definition = """
    # Extracted cells
    -> ImagingAnalysis
    cell_id              : int                          # Cell ID
    ---
    entry_time=CURRENT_TIMESTAMP : timestamp            # Auto created timestamp
    """

    class Rois(dj.Part):
        definition = """
        # First (functional) channel ROI properties
        -> Cell
        ---
        npix=null            : double                 # Number of pixels in ROIs
        center_x             : double                 # center x coordinate
        center_y             : double                 # center y coordinate
        xpix                 : longblob               # x coordinates in pixels
        ypix                 : longblob               # y coordinates in pixels
        lambda               : longblob               # Lambda (weights) for every pixel in ROI
        mrs=null             : double                 # Mean root square error (distance) to center point
        mrs0=null            : double                 # Mean root square error (distance) to center point
        compact=null         : double                 # ROI compactness
        footprint=null       : double                 # Footprint of ROIs
        isoverlap            : longblob               # Overlap of every pixel with other ROIs?
        aspect_ratio=null    : double                 # Aspect ratio
        ipix_neuropil        : longblob               # Index of neuropil pixels
        skew=null            : double                 # Skewness of fluorescence signal
        std=null             : double                 # Standard deviation of fluorescence signal
        ellipse              : longblob               # Ellipse around ROI (x/y)
        radius=null          : double                 # ROI radius in pixels
        noiselevel=null      : double                 # Suite2P noiselevel - suite2p matlab only
        cell_prob=null       : double                 # Cell probability (Based on classifier)
        redcell=null         : tinyint                # Is red cell? (0-no / 1-yes)
        redcell_prob=null    : double                 # Probability of being a red cell
        center_plane=0       : smallint                     
        neuropil_mask=null   : longblob               # Neuropil mask for each cell 
        cell_solidity=null   : double                 # Solidity of every ROI.
        npix_norm=null       : double                 # Number of pixels per ROI, normalised to max
        npix_norm_no_crop=null : double               # Number of pixels per ROI, normalised to max, without cropping
        soma_crop=null       : longblob               # Boolean mask for ROI pixels, =1 where soma
        npix_soma=null       : double                 # Number of pixels per ROI, soma crop
        """

    class Traces(dj.Part):
        definition = """
        # Traces
        -> Cell
        channel              : enum('primary','secondary')  # the channel that this trace comes from (ROI masks are always primary)
        ---
        samples=null         : double                       # Number of samples in time
        neuropil_coef=null   : double                       # Extracted neuropil coefficient
        fluo                  : longblob                     # Raw fluorescence trace
        neuropil             : longblob                     # Neuropil fluorescence trace
        fcorr                : longblob                     # Neuropil corrected fluorescence trace (fluo-neuropil_coef*neuropil)
        df_f                 : longblob                     # Delta F/F trace
        """

    class Spikes(dj.Part):
        definition = """
        # Spikes
        -> Cell.Traces
        ---
        spikes               : longblob                     # Unfiltered, deconvolved spikes
        """


@schema
class FilteredSpikes(dj.Computed):
    definition = """
    # Spike filtering
    -> Cell.Spikes
    -> SpikeFilterParams
    ---
    filtered_spikes       : longblob                        # Filtered, deconvolved spikes
    spike_threshold      : double                          # Calculated spike threshold
    entry_time=CURRENT_TIMESTAMP : timestamp               # Auto created timestam
    """


@schema
class SNR(dj.Computed):
    definition = """
    # Signal to noise ratio of fluorescence timeseries
    -> Cell.Traces
    -> FilteredSpikes
    -> NoiseCalcParams
    ---
    snr_df_f=null        : double                       # SNR over dF/F
    snr_f_corr=null      : double                       # SNR over neuropil corrected traces (FCorr)
    noise_indices=null   : longblob                     # Indices of episodes classified as noise
    noise_df_f=null      : longblob                     # Noise over dF/F
    noise_f_corr=null    : longblob                     # Noise over neuropil corrected traces (FCorr)
    noise_std=null       : double                       # Standard Dev of dF/F noise
    noise_std_half1=null : double                       # Standard Dev of dF/F noise - 1st half
    noise_std_half2=null : double                       # Stanrard Dev of dF/F noise - 2nd half
    noise_med_half1=null : double                       # Median of dF/F noise - 1st half
    noise_med_half2=null : double                       # Median of dF/F noise - 2nd half
    noise_ratio=null     : double                       # Ratio of noise median in first vs. second half (dF/F)
    noise_slope=null     : double                       # Linear regression slope of samples vs. noise (dF/F)
    t_value=null         : double                       # T value: Result of ttest_ind (scipy) for first half noise vs. second half noise (dF/F)
    p_value=null         : double                       # p value: Result of ttest_ind (scipy) for first half noise vs. second half noise (dF/F)
    entry_time=CURRENT_TIMESTAMP : timestamp            # Auto created timestam
    """


@schema
class RoisCorr(dj.Computed):
    definition = """
    # FOV corrected ROIs
    -> ImagingFOV
    -> Cell.Rois
    ---
    npix_corr            : double                 # Number of pixels in ROIs
    center_x_corr        : double                 # center x coordinate in microns
    center_y_corr        : double                 # center y coordinate in microns
    xpix_corr            : longblob               # x coordinates in microns
    ypix_corr            : longblob               # y coordinates in microns
    lambda_corr          : longblob               # Lambda (weights) for every pixel in ROI
    compress_ratio       : double                 # Compression ratio (unwarped/original)
    center_plane=null    : int                    # Imaging plane of this ROI - identical to center_plane in Cell.Roi
    """


@schema
class TrackingParams(dj.Lookup):
    definition = """
    # Open field tracking preprocessing parameters
    trackingparams_id    : char(1)                     # Parameter set ID, starting with A
    ---
    sigma_tracking       : float                        # Gaussian smoothing for tracking positions
    sigma_speed          : float                        # Gaussian smoothing for speed
    sigma_angle          : float                        # Gaussian smoothing for head direction
    rotate               : smallint                    # Auto-rotate arena coordinates? (0-no / 1-yes)
    tracking_params_hash="" : varchar(32)                  
    """


@schema
class Tracking(dj.Computed):
    definition = """
    # Tracking data
    -> TrackingRaw
    -> TrackingParams
    ---
    sample_rate          : double                       # Sampling rate
    entry_time_tracking=CURRENT_TIMESTAMP : timestamp   # Auto created timestamp
    """

    class Linear(dj.Part):
        definition = """
        # Linear tracking data - constrained movement
        -> Tracking
        ---
        linear_type          : enum('optical','rotary')     # Type of linear tracking setup
        calibration          : double                       # Conversion factor steps to mm (1: No calibration supplied)
        pos                  : longblob                     # Absolute position zeroed to first sample (steps or mm)
        rel_pos              : longblob                     # Relative position change from one frame to the next (steps or mm)
        lap_pos=null         : longblob                     # Absolute position reset each lap (steps)
        speed                : longblob                     # Speed (steps/s or mm/s)
        timestamps=null      : longblob                     # Timestamps in seconds - synced timestamps (identical to TrackingRaw timestamps)
        """

    class OpenField(dj.Part):
        definition = """
        # LED tracking data - unconstrained movement
        -> Tracking
        ---
        calibration          : double                 # Conversion factor px to mm (1: No calibration supplied)
        x_pos                : longblob               # Tracking x position (px or mm)
        y_pos                : longblob               # Tracking y position (px or mm)
        boundary_left        : double                 # Extracted boundary towards left (px or mm)
        boundary_right       : double                 # Extracted boundary towards right (px or mm)
        boundary_bottom      : double                 # Extracted boundary towards bottom (px or mm)
        boundary_top         : double                 # Extracted boundary towards top (px or mm)
        speed                : longblob               # Speed (mm or px/s)
        head_angle           : longblob               # Head angle [0,2*pi]
        timestamps           : longblob               # Timestamps in seconds - synced timestamps (identical to TrackingRaw timestamps)
        """

    class DLCPart(dj.Part):
        definition = """
        # Deeplabcut tracking data
        -> Tracking
        -> TrackingRaw.DLCPart
        ---
        bodypart_x_pos       : longblob               
        bodypart_y_pos       : longblob               
        bodypart_likelihood  : longblob               
        """


@schema
class Occupancy(dj.Computed):
    definition = """
    # Occupancy maps open field
    -> Tracking.OpenField
    -> SignalTrackingParams
    -> MapParams
    ---
    occupancy            : longblob               # Smoothed 2D occupancy map [seconds]
    mask_occ             : longblob               # Mask (where time = 0)
    occupancy_raw        : longblob               # Raw, non-smoothed 2D occupancy map
    explor_ratio         : double                 # Exploration ratio (visited bins over all bins)
    explor_std           : double                 # Exploration standard deviation (of visited bins)
    x_edges              : longblob               # Histogram edges in x
    y_edges              : longblob               # Histogram edges in y
    occupancy_time       : double                 # Time in seconds in occupancy
    fraction_occupancy   : double                 # Fraction of time in occupancy ma
    """


@schema
class ArenaObjectPos(dj.Computed):
    definition = """
    # Calibrated position of objects
    object_hash          : varchar(16)                  # Object hash
    ---
    -> Session.ArenaObject
    -> Tracking.OpenField
    obj_x_coord_calib=null : float                       # Object x coord [mm]
    obj_y_coord_calib=null : float                       # Object y coord [mm]
    """


@schema
class AngularOccupancy(dj.Computed):
    definition = """
    # Angular occupancy from tracking angle
    -> Tracking.OpenField
    -> AngularParams
    -> SignalTrackingParams
    ---
    angular_occupancy    : longblob               # Smoothed 1D occupancy [seconds]
    angular_occupancy_raw : longblob              # Raw, non-smoothed 1D occupancy
    angle_edges          : longblob               # Histogram bin edges in radians
    angle_centers        : longblob               # Histogram bin centers in radian
    """


@schema
class SignalTracking(dj.Computed):
    definition = """
    # Tracking entry for every spike / calcium event
    -> FilteredSpikes.proj(signal_dataset="dataset_name")
    -> Tracking.OpenField.proj(tracking_dataset="dataset_name")
    -> SignalTrackingParams
    signal_type          : enum('df_f','spikes')        # Signal type (Fluorescence vs. Spikes)
    ---
    -> [nullable] Sync.proj(sync_dataset_frames_imaging="dataset_name",sync_name_frames_imaging="sync_name")
    signal               : longblob               # Signal (spikes/calcium) amplitudes
    x_pos_signal         : longblob               # Tracking x position for signal (center of 2 LEDs)
    y_pos_signal         : longblob               # Tracking y position for signal (center of 2 LEDs)
    speed_signal         : longblob               # Speed in pixels / sec for signal
    head_angle_signal    : longblob               # Head angle [0,2*pi] for signa
    """


@schema
class AngularRate(dj.Computed):
    definition = """
    # 1D ratemaps head direction
    -> SignalTracking
    -> AngularOccupancy.proj(tracking_dataset="dataset_name")
    ---
    angular_rate         : longblob               # Smoothed 1D tuning curve
    angular_rate_raw     : longblob               # Unsmoothed (raw) 1D tuning curve
    binned_raw_angular   : longblob               # Raw, binned signal
    bin_max_angular=null : double                 # Bin (center, in radians [0, 2*pi]) with maximum signal
    max_angular=null     : double                 # Maximum
    """

    class Stats(dj.Part):
        definition = """
        # Tuning curve statistics
        -> AngularRate
        ---
        mvl=null             : double                       # Mean resultant vector length
        angular_mean=null    : double                       # Angular mean in radians ([0, 2*pi])
        angular_var=null     : double                       # Angular variance
        angular_std=null     : double                       # Angular standard deviation
        rayleigh_p=null      : double                       # Rayleigh test two-tailed p-value
        rayleigh_z=null      : double                       # Rayleigh test value of the z-statistic
        """


@schema
class Ratemap(dj.Computed):
    definition = """
    # 2D ratemaps open field
    -> SignalTracking
    -> Occupancy.proj(tracking_dataset="dataset_name")
    -> FieldParams
    ---
    ratemap              : longblob               # Smoothed 2D ratemap
    ratemap_raw          : longblob               # Unsmoothed (raw) 2D ratemap
    mask_rm              : longblob               # Mask (where time = 0)
    binned_raw           : longblob               # Raw, binned signal
    bin_max              : longblob               # Bin with maximum signal (ratemap(bin_max) = max(ratemap))
    max                  : double                 # Maximum
    no_fields=null       : smallint                # Total number of detected fields
    fields_map=null      : longblob                # Map of detected field
    """

    class Fields(dj.Part):
        definition = """
        # Field detection results
        -> Ratemap
        field_no             : int                    # Field number
        ---
        field_coords         : longblob               # Coordinates of all bins in the firing field
        field_peak_x         : double                 # Field peak firing x coordinate
        field_peak_y         : double                 # Field peak firing y coordinate
        field_centroid_x     : double                 # Field centroid x coordinate
        field_centroid_y     : double                 # Field centroid y coordinate
        field_area           : int                    # Area in number of bins
        field_bbox           : longblob               # Field bounding box
        field_mean_rate      : double                 # Field mean value (rate)
        field_peak_rate      : double                 # Field peak value (rate)
        field_map            : longblob               # Field map
        """

    class Stats(dj.Part):
        definition = """
        # Ratemap statistics
        -> Ratemap
        ---
        information_rate=null : double                      # Information rate
        information_content=null : double                   # Information content
        sparsity=null        : double                       # Sparsity
        selectivity=null     : double                       # Selectivity
        coherence=null       : double                       # Coherence
        """


@schema
class GridScore(dj.Computed):
    definition = """
    # Grid score and stats
    -> Ratemap
    -> FieldParams
    -> GridParams
    ---
    gridscore=null       : double                 # Grid score
    acorr=null           : longblob               # Spatial autocorrelatio
    """

    class Stats(dj.Part):
        definition = """
        # Grid statistics
        -> GridScore
        ---
        gs_spacings          : longblob               # Spacing of three adjacent fields closest to center in autocorr (in [bins])
        gs_spacing=null      : double                 # Nanmean of 'spacings' in [bins]
        gs_positions         : longblob               # [y,x] coordinates of six fields closest to center
        gs_orientations      : longblob               # Orientation of three adjacent fields closest to center in autocorr (in [radians])
        gs_orientations_std=null : double             # Standard deviation of orientations modulo 60 (in [radians])
        gs_orientation=null  : double                 # Orientation of grid in [radians] (mean of fields of 3 main axes)
        gs_ellipse_theta=null : double                # Ellipse theta (corrected according to previous BNT standard) in [radians]
        gs_ellipse_aspect_ratio=null : double         # Ellipse aspect ratio (major radius / minor radius)
        gs_ellipse           : longblob               # Ellipse fit returning [x coordinate, y coordinate, major radius, minor radius, theta]
        """


@schema
class BorderParams(dj.Lookup):
    definition = """
    # Borderscore analysis parameters
    border_params_id     : char(1)                      # Parameter set ID, starting with A
    ---
    search_width         : int                          # Pixel number from border that should be checked
    walls                : varchar(4)                   # Definition of walls that the borderscore is calculated for (TLBR)
    border_params_hash="" : varchar(32)                  
    """


@schema
class BorderScore(dj.Computed):
    definition = """
    # Border score and stats
    -> Ratemap
    -> BorderParams
    ---
    borderscore=null     : double                       # Border Score
    """


@schema
class Shuffled(dj.Computed):
    definition = """
    # Shuffling table
    -> FilteredSpikes.proj(signal_dataset="dataset_name")
    -> Tracking.OpenField.proj(tracking_dataset="dataset_name")
    -> ShuffleParams
    -> SignalTrackingParams
    -> Occupancy.proj(tracking_dataset="dataset_name")
    -> MapParams
    -> FieldParams
    -> AngularOccupancy.proj(tracking_dataset="dataset_name")
    -> AngularParams
    -> GridParams
    -> BorderParams
    ---
    -> [nullable] Sync.proj(sync_dataset_frames_imaging="dataset_name",sync_name_frames_imaging="sync_name")
    number_shuffles      : int                    # Total number of shuffles (can vary from expected number)
    shuffling_offsets    : longblob               # Shuffling offset
    """

    class AngularRateStats(dj.Part):
        definition = """
        # Shuffled tuning curve stats
        -> Shuffled
        ---
        mvl_99               : double                # Mean vector length 99th percentile
        mvl_95               : double                # Mean vector length 95th percentile
        mvl_shuffles         : longblob               # Individual shuffles mean vector length
        """

    class BorderScore(dj.Part):
        definition = """
        # Shuffled  borderscore stats
        -> Shuffled
        ---
        borderscore_99       : double                # Borderscore 99th percentile
        borderscore_95       : double                # Borderscore 95th percentile
        borderscore_shuffles : longblob               # Individual shuffles borderscore
        """

    class Fields(dj.Part):
        definition = """
        # Shuffled fields stats
        -> Shuffled
        ---
        fields_peak_x        : longblob               # Field peak firing x coordinate
        fields_peak_y        : longblob               # Field peak firing y coordinate
        fields_centroid_x    : longblob               # Field centroid x coordinate
        fields_centroid_y    : longblob               # Field centroid y coordinate
        fields_area          : longblob               # Area in number of bins
        """

    class GridScore(dj.Part):
        definition = """
        # Shuffled gridscore stats
        -> Shuffled
        ---
        gridscore_99         : double                # Gridscore 99th percentile
        gridscore_95         : double                # Gridscore 95th percentile
        gridscore_shuffles   : longblob               # Individual shuffles gridscore
        """

    class RatemapStats(dj.Part):
        definition = """
        # Shuffled ratemap stats
        -> Shuffled
        ---
        information_rate_99  : double                # Information rate 99th percentile
        information_rate_95  : double                # Information rate 95th percentile
        information_rate_shuffles : longblob          # Individual shuffles information rate
        information_content_99 : double              # Information content 99th percentile
        information_content_95 : double              # Information content 95th percentile
        information_content_shuffles : longblob       # Individual shuffles information content
        sparsity_99          : double                # Sparsity 99th percentile
        sparsity_95          : double                # Sparsity 95th percentile
        sparsity_shuffles    : longblob               # Individual shuffles sparsity
        selectivity_99       : double                # Selectivity 99th percentile
        selectivity_95       : double                # Selectivity 95th percentile
        selectivity_shuffles : longblob               # Individual shuffles selectivity
        coherence_99         : double                # Coherence 99th percentile
        coherence_95         : double                # Coherence 95th percentile
        coherence_shuffles   : longblob               # Individual shuffles coherence
        """
