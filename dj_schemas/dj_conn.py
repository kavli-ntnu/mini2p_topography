import datajoint as dj 


# Update as necessary, or make use of a local or global datajoint configuration
# The values given below are the defaults if running a docker container as noted in the readme
# dj.config["database.user"] = "root"
# dj.config["database.password"] = "simple"
# dj.config["database.host"] = "localhost"
# dj.config["database.port"] = 3306

#### DATABASE CONNECTION #########################################

mlims_db         = 'group_shared_topopaper_mlims'
main_imaging_db  = 'group_shared_topopaper_main_imaging' # group_imaging_1b # group_shared_topopaper_main_imaging
horst_imaging_db = 'group_shared_topopaper_horst_imaging' # user_horsto_imaging # group_shared_topopaper_horst_imaging
border_db        = 'group_shared_topopaper_borderscore' # `user_horsto_borderscore` # group_shared_topopaper_borderscore

# Load mlims
mlims = dj.schema(mlims_db)
mlims.spawn_missing_classes()
# Load main imaging
main_imaging = dj.schema(main_imaging_db)
main_imaging.spawn_missing_classes()
# Load horst imaging
horst_imaging = dj.schema(horst_imaging_db)
horst_imaging.spawn_missing_classes()
# Load borderscore
borderscore = dj.schema(border_db)
borderscore.spawn_missing_classes()

