# DB_HOST=datajoint.it.ntnu.no DB_USER=user DB_PASS=pass ./dump.sh
mysqldump -h $DB_HOST -u$DB_USER -p$DB_PASS --databases group_shared_topopaper_mlims group_shared_topopaper_main_imaging group_shared_topopaper_horst_imaging group_shared_topopaper_borderscore --max_allowed_packet=1G> dump.sql
