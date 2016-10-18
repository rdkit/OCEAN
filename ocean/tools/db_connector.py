import os,sys
import ocean.settings


class DB_connector:

    def __init__(self,db):
        self.db = db
        if not db in ocean.settings.DATABASES:
            raise Exception('Database {} not in settings.DATABASES {}'.format(db,ocean.settings.DATABASES))
        self.current_db = ocean.settings.DATABASES[db]
        if 'oracle' in self.current_db['ENGINE']:
            import cx_Oracle
            self.db_con = cx_Oracle.connect(self.current_db['USER'],
                                            self.current_db['PASSWORD'],
                                            self.current_db['HOST'])
            self.db_type = 'oracle'
        elif 'postgre' in self.current_db['ENGINE']:
            import psycopg2
            self.db_con = psycopg2.connect(database=self.current_db['NAME'],
                                           user=self.current_db['USER'],
                                           password=self.current_db['PASSWORD'],
                                           host=self.current_db['HOST'])
            self.db_type = 'postgre'
        else:
            raise Exception('unknown ENGINE {}'.format(self.current_db['ENGINE']))

        self.cursor = self.db_con.cursor()
        self.filter = ""
        self.cursor.arraysize = 64
        self.verbose = False

    def close(self):
        self.db_con.close()

    def getTargetName(self,chembl_id):
        return self.getTargetNameForChemblID(chembl_id)

    def getAllMolFiles(self):
        self.cursor.execute("select molregno,molfile from compound_structures")
        iterator = self.cursor.__iter__()
        names = [x[0] for x in self.cursor.description]
        return names,iterator

    def getCompoundCount(self):
        result = self.cursor.execute("select count(*) from (select min(molregno) from %s group by molregno)" % ocean.settings.OCEAN_DB_TABLE)
        return result.fetchone()[0]

    def getTargetCount(self):
        result = self.cursor.execute("select count(*) from (select min(target_chembl_id) from %s group by target_chembl_id)" % ocean.settings.OCEAN_DB_TABLE)
        return result.fetchone()[0]

    def getRndCompounds(self,count):
        if count>1000: raise Exception("can fetch max 1000 entries")
        import random
        count_entries = self.getCompoundCount()
        print count_entries

        rstr = "(" + str(random.randint(1,count_entries))

        for rnd in range(1,count):
            rstr += ","+ str(random.randint(1,count_entries))
        rstr += ")"

        query = "select * from (select MOLREGNO,rownum rn from compound_structures) where rn in %s" % rstr
        self.cursor.execute(query)
        iterator = self.cursor.__iter__()
        names = [x[0] for x in self.cursor.description]
        return names,iterator

    def getCmpdsForTarget2(self,target_id):
        query =  "select molregno,target_pref_name,organism from %s where target_chembl_id='%s' %s"
        final_query = query % (ocean.settings.OCEAN_DB_TABLE,target_id,self.filter)
        self.cursor.execute(final_query)
        iterator = self.cursor.__iter__()
        target_pref_name = None
        organism = None
        compound_molregnos = []
        for entry in iterator:
            compound_molregnos.append(entry[0])
            if target_pref_name is None:
                target_pref_name = entry[1]
                organism = entry[2]
        return compound_molregnos,target_pref_name,organism

    def getModelDataset(self,tableName):
        query =  "select molecule_chembl_id,standard_value,canonical_smiles from %s"
        final_query = query % (tableName)
        self.cursor.execute(final_query)
        iterator = self.cursor.__iter__()
        result = []
        for entry in iterator:
            result.append((entry[0],(entry[1],entry[2])))
        return result

    def getCmpdsForTarget_with_IC50_dict(self,target_id):
        query =  "select molregno,standard_value from %s where target_chembl_id='%s' %s"
        final_query = query % (ocean.settings.OCEAN_DB_TABLE,target_id,self.filter)
        self.cursor.execute(final_query)
        iterator = self.cursor.__iter__()
        result = {cmpd[0]:{"value":cmpd[1]} for cmpd in iterator}
        return result

    def getMolregnos(self):
        query =  "select min(molregno) from %s group by molregno" % ocean.settings.OCEAN_DB_TABLE
        final_query = query
        self.cursor.execute(final_query)
        iterator = self.cursor.__iter__()
        result = [cmpd[0] for cmpd in iterator]
        return result

    def getMolregnosAndSmiles(self):
        query =  "select /*+ PARALLEL({0},4) */ molregno,canonical_smiles from {0} group by molregno,canonical_smiles".format(ocean.settings.OCEAN_DB_TABLE)
        final_query = query
        self.cursor.execute(final_query)
        iterator = self.cursor.__iter__()
        result = [(cmpd[0],cmpd[1]) for cmpd in iterator]
        return result

    def getInfoForTargetAndCompound(self,molregno,target_id):
        if molregno == -1:
            query = "select molregno,molecule_chembl_id,target_pref_name,organism,canonical_smiles,standard_value from %s where target_chembl_id='%s' %s"
            final_query = query % (ocean.settings.OCEAN_DB_TABLE,target_id,self.filter)
        else:
            query =  "select molregno,molecule_chembl_id,target_pref_name,organism,canonical_smiles,standard_value from %s where target_chembl_id='%s' and molregno=%d %s"
            final_query = query % (ocean.settings.OCEAN_DB_TABLE,target_id,molregno,self.filter)
        self.cursor.execute(final_query)
        iterator = self.cursor.__iter__()
        if molregno == -1:
            result = [(hits[0],hits[1],hits[2],hits[3],hits[4],hits[5]) for hits in iterator]
        else:
            result = [(hits[0],hits[1],hits[2],hits[3],hits[4],hits[5]) for hits in iterator][0]
        return result

    def getMolregnoForSmiles(self,smiles):
        query = "select molregno from compound_structures where canonical_smiles='%s'"
        final_query = query % (smiles)
        self.cursor.execute(final_query)
        iterator = self.cursor.__iter__()
        result = [hits[0] for hits in iterator]
        if len(result)>0: return result[0]
        else: return None


    def getSmilesForMolregno(self,molregno):
        query = "select canonical_smiles from compound_structures where molregno='%d'"
        final_query = query % (molregno)
        self.cursor.execute(final_query)
        iterator = self.cursor.__iter__()
        result = [hits[0] for hits in iterator][0]
        return result

    def getValueForMolregno(self,molregno,target_id):
        query =  "select standard_value from %s where target_chembl_id='%s' and molregno='%d' %s"
        final_query = query % (ocean.settings.OCEAN_DB_TABLE,target_id,molregno,self.filter)
        self.cursor.execute(final_query)
        iterator = self.cursor.__iter__()
        result = [cmpd[0] for cmpd in iterator][0]
        return result

    def getChemblIDforMolregno(self,molregno):
        query = "select chembl_id from molecule_dictionary where molregno='%d'"
        final_query = query % (molregno)
        self.cursor.execute(final_query)
        iterator = self.cursor.__iter__()
        result = [hits[0] for hits in iterator][0]
        return result

    def getTargetNameForChemblID(self,chemblID):
        query = "select pref_name from target_dictionary where chembl_id = '%s'"
        final_query = query % (chemblID)
        self.cursor.execute(final_query)
        iterator = self.cursor.__iter__()
        result = [hits[0] for hits in iterator][0]
        return result

    def getOrganismForChemblID(self,chemblID):
        query = "select organism from target_dictionary where chembl_id = '%s'"
        final_query = query % (chemblID)
        self.cursor.execute(final_query)
        iterator = self.cursor.__iter__()
        result = [hits[0] for hits in iterator][0]
        return result

    def getTargetNameAndOrganismForChemblID(self,chemblID):
        query = "select pref_name,organism from target_dictionary where chembl_id = '%s'"
        final_query = query % (chemblID)
        self.cursor.execute(final_query)
        iterator = self.cursor.__iter__()
        result = [(hits[0],hits[1]) for hits in iterator][0]
        return result

    def getProteinClassifications(self):
        query = "select td.chembl_id,pfc.l1,pfc.l2 from target_dictionary td,target_components tc, component_class cc, protein_family_classification pfc where td.tid = tc.tid and tc.component_id = cc.component_id and cc.protein_class_id = pfc.protein_class_id"
        final_query = query
        self.cursor.execute(final_query)
        iterator = self.cursor.__iter__()
        result = [(entry[0],entry[1],entry[2]) for entry in iterator]
        return result

    def setFilter(self,sql_filter):
        self.filter=sql_filter
        return

    def execute(self,query):
        self.cursor.execute(query)
        data = self.cursor.fetchall()
        names = [x[0] for x in self.cursor.description]
        return names,data

    def getTargets(self):
        if self.db_type == 'postgre':
            query = "select target_chembl_id,pa.cs from (select target_chembl_id,count(molregno) cs from %s group by target_chembl_id) as pa where cs>%d" % (
            ocean.settings.OCEAN_DB_TABLE, ocean.settings.CMPD_COUNT_CUTOFF - 1)
        else:
            if self.filter:
                if self.filter == "parallel":
                    query = "select target_chembl_id,cs from (select /*+ PARALLEL({0},4) */ target_chembl_id,count(molregno) cs from {0} group by target_chembl_id) where cs>{1}".format(
                        ocean.settings.OCEAN_DB_TABLE, ocean.settings.CMPD_COUNT_CUTOFF - 1)
                else:
                    query = "select distinct target_chembl_id from %s where %s" % (
                    ocean.settings.OCEAN_DB_TABLE, self.filter)

            else:
                query = "select min(target_chembl_id) from %s group by target_chembl_id" % ocean.settings.OCEAN_DB_TABLE

        final_query = query
        print "query is", final_query
        self.cursor.execute(final_query)

        result = [t[0] for t in self.cursor.__iter__()]
        return result

if os.path.exists('ocean/tools/custom_db_connector.py'):
    import custom_db_connector
    cdc = {k:v for k,v in custom_db_connector.__dict__.items() if not k.startswith('__')}
    current_locals = locals()
    for entry,value in cdc.items():
        if not entry in current_locals.keys():
            current_locals.update({entry:value})
            print >> sys.stderr, "monkey patch custom db_connector",entry,value