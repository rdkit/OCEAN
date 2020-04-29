curwd=$(pwd)
db_filename="${curwd}/ocean_db_chembl20_data.csv"
#psql -c "copy ocean_db from 'ocean_db_chembl20_data.csv' with delimiter ',' CSV HEADER" chembl17
psql -c "copy ocean_db from '${db_filename}' with delimiter ',' CSV HEADER" chembl_20
