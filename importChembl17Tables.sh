psql -c "copy target_dictionary from '/home/chembl/OCEAN/chembl17_target_dictionary.csv' with delimiter ',' CSV HEADER" chembl_17
psql -c "copy target_components from '/home/chembl/OCEAN/chembl17_target_components.csv' with delimiter ',' CSV HEADER" chembl_17
psql -c "copy component_class from '/home/chembl/OCEAN/chembl17_component_class.csv' with delimiter ',' CSV HEADER" chembl_17
psql -c "copy protein_family_classification from '/home/chembl/OCEAN/chembl17_protein_family_classification.csv' with delimiter ',' CSV HEADER" chembl_17
