# DB_USER=root DB_PASS=simple ./load.sh
docker exec -it sqldump_db_1 bash -c "mysql -u$DB_USER -p$DB_PASS < /tmp/dump.sql"
