DROP TABLE ocean_db;

create table ocean_db as (
    select result.target_chembl_id,
           result.target_pref_name,
           result.organism,
           result.molregno,
           result.molecule_chembl_id,
           result.avg as standard_value,
           result.standard_relation,
           result.standard_units,
           result.canonical_smiles
    from (
             select subQ.target_chembl_id,
                    subQ.target_pref_name,
                    subQ.organism,
                    subQ.molregno,
                    subQ.molecule_chembl_id,
                    subQ.canonical_smiles,
                    subQ.standard_relation,
                    subQ.standard_units,
                    count(*)                                                   as cnt,
                    stddev(subQ.standard_value)                                as std,
                    avg(subQ.standard_value)                                   as avg,
                    min(subQ.standard_value)                                   as min,
                    max(subQ.standard_value)                                   as max,
                    avg(subQ.standard_value) - 3 * stddev(subQ.standard_value) as allowed_min,
                    avg(subQ.standard_value) + 3 * stddev(subQ.standard_value) as allowed_max
             from (
                      select td.chembl_id        as target_chembl_id,
                             td.pref_name        as target_pref_name,
                             td.organism,
                             md.molregno,
                             md.chembl_id        as molecule_chembl_id,
                             act.standard_value,
                             act.standard_relation,
                             act.standard_units,
                             act.standard_type,
                             cs.canonical_smiles as canonical_smiles
                      from target_dictionary td
                               join assays a on td.tid = a.tid
                               join activities act on a.assay_id = act.assay_id
                               join molecule_dictionary md on act.molregno = md.molregno
                               join compound_properties cp on md.molregno = cp.molregno
                               join compound_structures cs on md.molregno = cs.molregno
                      where standard_relation = '='
                        and standard_units = 'nM'
                        and standard_type in ('IC50', 'Ki')
                        and target_type in ('SINGLE PROTEIN', 'PROTEIN COMPLEX')
                        and organism = 'Homo sapiens'
                  ) as subQ
             group by subQ.molregno,
                      subQ.molecule_chembl_id,
                      subQ.target_chembl_id,
                      subQ.target_pref_name,
                      subQ.organism,
                      subQ.standard_relation,
                      subQ.standard_units,
                      subQ.canonical_smiles
         ) as result
    where
	result.avg < 10000 and
	(
		result.std is null or
		(
			result.min >= result.allowed_min and
			result.max <= result.allowed_max
		)
	)
);

CREATE UNIQUE INDEX OCEAN_DB_IDX1 ON ocean_db (TARGET_CHEMBL_ID, MOLREGNO)   ;
CREATE INDEX OCEAN_DB_IDX2 ON ocean_db (TARGET_CHEMBL_ID)   ;
CREATE INDEX OCEAN_DB_IDX3 ON ocean_db (MOLREGNO)   ;
CREATE INDEX OCEAN_DB_IDX4 ON ocean_db (TARGET_CHEMBL_ID, STANDARD_VALUE)   ;
CREATE INDEX OCEAN_DB_IDX5 ON ocean_db (CANONICAL_SMILES) ;

delete from ocean_db where target_chembl_id in (select chembl_id from target_dictionary where target_type not in ('SINGLE PROTEIN','PROTEIN COMPLEX'));
delete from ocean_db where target_chembl_id in (select insel.target_chembl_id from (select target_chembl_id,count(molregno) as cs from ocean_db group by target_chembl_id) as insel where insel.cs<10)
