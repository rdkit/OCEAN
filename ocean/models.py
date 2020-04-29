from django.conf       import settings
from django.db         import models

import numpy as np
# import settings
from ocean import settings

from ocean.tools.db_connector import DB_connector

import os,sys
import pickle

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))

class FP_Manager:

    def __init__(self):
        self.fp_vault = None
        self.verbose = False
        self.targets = {}
        self.clear()
        pass

    def clear(self):
        self.fp_vault = {x: {} for x in settings.DATASOURCES}

    def __str__(self):
        template = """\
Class: {0}
DataSources:"""
        for source in sorted(self.fp_vault.keys(),key=lambda x: x.name):
            template += "\n\t{0} ({1})".format(source.name,", ".join(["FP {0}: {1} Items".format(x,len(self.fp_vault[source][x])) for x in sorted(self.fp_vault[source].keys())]))
        return template.format(self.DataSource.__name__)

    def fpAvailable(self,fp_id):
        for datasource,data in self.fp_vault.items():
            if not fp_id in data or len(data[fp_id]) == 0:
                return False
        return True

    def loadFPsFromDataSource(self,datasource,fp_id):
        # import pickle
        if datasource in settings.DATASOURCES:
            if not fp_id in self.fp_vault[datasource]:
                self.fp_vault[datasource][fp_id] = {}
            else:
                self.fp_vault[datasource][fp_id].clear()

            ds = DataSources.objects.filter(name=datasource.name)
            if ds.count()==0:
                print(f"Datasource {datasource.name} does not exist", file=sys.stderr)
                print(f"create Datasource {datasource.name}", file=sys.stderr)
                new_ds = DataSources(name=datasource.name)
                new_ds.save()
                ds = DataSources.objects.filter(name=datasource.name)

            all_fps = All_FPS.objects.filter(fp_id=fp_id,datasource=ds[0])
            i = 0
            for entry in all_fps:
                i+=1
                # if type(entry.fp_dict) is buffer:
                #     tmp_dict = pickle.loads(str(entry.fp_dict))
                # else:
                tmp_dict = pickle.loads(entry.fp_dict)
                self.fp_vault[datasource][fp_id].update(tmp_dict)
            if i == 0:
                print(f"Datasource {datasource.name} is empty!", file=sys.stderr)
            if self.verbose:
                print(f"loaded {len(self.fp_vault[datasource][fp_id])} Entries of {datasource}, FP {fp_id} into fp_vault")
        else:
            raise Exception("unknown DataSource")

    def loadEverything(self):
        for source in DataSources.objects.all():
            x = All_FPS.objects.filter(datasource=source).values('fp_id').distinct()
            for fp in x:
                fp_id = fp['fp_id']
                print("load source.name", source.name)
                print("this is datasource", settings.DATASOURCES[source.name])
                print("of datasources", settings.DATASOURCES)
                self.loadFPsFromDataSource(settings.DATASOURCES[source.name], fp_id)
        pass

    def __getitem__(self, item):
        return self.fp_vault[item]


class Target_Compounds:
    vault = {}
    vault_set = set()
    targets = {}
    counts = {}
    pass

    @staticmethod
    def getTargets(datasource):
        if not datasource in Target_Compounds.targets:
            raise Exception("invalid datasource {0}".format(datasource))

        if not datasource in Target_Compounds.targets:
            return None
        else:
            return Target_Compounds.targets[datasource]

    @staticmethod
    def setTargets(datasource,targets):
        print(f"set {len(targets)} targets for {datasource}")
        Target_Compounds.targets[datasource] = set(targets)

    @staticmethod
    def _fill_targets_():
        for datasource in settings.DATASOURCES:
            if datasource is settings.DATASOURCES.CHEMBL:
                chembl_connection = DB_connector(settings.CHEMBL_VERSION)
                chembl_connection.setFilter("parallel")
                t = set(chembl_connection.getTargets())
                Target_Compounds.setTargets(datasource,t)
                chembl_connection.close()
            else:
                pass
            Target_Compounds.counts[datasource.name] = len(t)
        pass

    @staticmethod
    def _fill_vault_(drop_cmpds = set()):
        PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
        for datasource in settings.DATASOURCES:
            if not datasource in Target_Compounds.vault:
                Target_Compounds.vault[datasource] = {}
            print(f"Collect Target-Compounds of Dataset {datasource.name}")
            if datasource is settings.DATASOURCES.CHEMBL:
                chembl_connection = DB_connector(settings.CHEMBL_VERSION)
                chembl_connection.setFilter("and standard_value < %d" % settings.CMPD_NM_CUTOFF)
                for target in Target_Compounds.targets[datasource]:
                    target_cmpds,target_name,target_organism = chembl_connection.getCmpdsForTarget2(target)
                    chembl_entry = (datasource,target)
                    if not chembl_entry in Target_Compounds.vault_set:
                        Target_Compounds.vault_set.add(chembl_entry)
                    cmpd_set_old = set(target_cmpds)
                    cmpd_set = cmpd_set_old.difference(drop_cmpds)
                    Target_Compounds.vault[datasource][target] = {'cmpds':cmpd_set,
                                                                  'desc':target_name}
                chembl_connection.close()
            else:
                pass
        pass

    @staticmethod
    def clear():
        Target_Compounds.targets = {}
        Target_Compounds.vault = {}

    @staticmethod
    def fill():
        Target_Compounds.clear()
        Target_Compounds._fill_targets_()
        if settings.drop_compounds is None:
            Target_Compounds._fill_vault_()
        else:
            Target_Compounds._fill_vault_(drop_cmpds = settings.drop_compounds)

class DataSources(models.Model):
    id = models.AutoField('id',primary_key=True,help_text="<< Primary Key >>")
    name = models.CharField('name',max_length=50)

    @staticmethod
    def getDict():
        ds = DataSources.objects.all()
        datasources = {d.id:str(d.name) for d in ds}
        return datasources

    @staticmethod
    def get(id):
        d = DataSources.objects.get(id=id)
        return d.name

    @staticmethod
    def getDataSourceActivities(datasource):
        if datasource is settings.DATASOURCES.CHEMBL:
            data = DataSources.getDataSourceChEMBL()
        else:
            raise Exception("unknown DataSource")
        return data

    @staticmethod
    def getDataSourceChEMBL():
        chembl_connection = DB_connector(settings.CHEMBL_VERSION)
        idAndSmiles = chembl_connection.getMolregnosAndSmiles()
        chembl_connection.close()
        return idAndSmiles

class ProcessManager:
    pm = []
    datasources = set()

    def __init__(self, process, pipeProcess, pipeOutside, datasource):
        self.process = process
        self.pipeProcess = pipeProcess
        self.pipeOutside = pipeOutside
        self.datasource = datasource
        ProcessManager.pm.append(self)

    def send(self, msg):
        self.pipeOutside.send(msg)

    def recv(self):
        return self.pipeOutside.recv()

    def start(self):
        self.process.start()

    def join(self):
        self.process.join()

    def is_alive(self):
        return self.process.is_alive()

    @staticmethod
    def kill_all():
        print("kill all processes")
        for p in ProcessManager.pm:
            p.send('die')

        ProcessManager.join_all()
        print("done")
        ProcessManager.pm = []

    @staticmethod
    def start_all():
        for p in ProcessManager.pm:
            if not p.is_alive():
                p.start()
            else:
                pass

    @staticmethod
    def join_all():
        for p in ProcessManager.pm:
            p.join()

    @staticmethod
    def send_all(msg, datasources):
        for p in ProcessManager.pm:
            if p.datasource in datasources:
                p.send(msg)

    @staticmethod
    def recv_all(datasources):
        result = []
        for p in ProcessManager.pm:
            if p.datasource in datasources:
                result += p.recv()
        return result


class Ocean_hit():

    def __init__(self,target,p_value=None,e_value=None):
        self.target = target
        self.compounds = None
        self.compounds_valid = None
        self.compounds_count = None
        self.tclist = None
        self.tc_avr = None
        self.p_value = p_value
        self.e_value = e_value
        self.comparedTo = None
        self.targetName = None
        self.organism = None
        self.classification = None
        self.datasourceName = None
        self.targetLink = None

    def __str__(self):
        if self.p_value is None or self.e_value is None:
            return "%s: %d Compounds, tc_average: %f" % (self.target,self.compounds_valid,self.tc_avr)
        else:
            return "%s: %d Compounds, tc_average: %f, p_value: %f, e_value: %f" % (self.target,self.compounds_valid,self.tc_avr,self.p_value,self.e_value)

    def __unicode__(self):
        return "%s" % self.__str__()

    def __repr__(self):
        return "<%s>" % self.__str__()

    def __eq__(self, other):
        if type(other)== str:
            return self.target == other
        if type(other)==Ocean_hit:
            return type(self) == type(other) and \
               self.target == other.target and \
               self.compounds == other.compounds and \
               self.tclist == other.tclist and \
               self.tc_avr == other.tc_avr

    def __lt__(self, other):
        if self.e_value is None:
            return self.tc_avr < other.tc_avr
        else:
            return self.e_value < other.e_value

    def setTClist(self,tclist):
        self.tclist = tclist
        self.tc_avr = np.mean(tclist)
        self.compounds_valid = len(tclist)

    def setCompounds(self,compounds):
        result = []
        for i in range(len(compounds)):
            l = list(compounds[i])
            l.append(self.tclist[i])
            result.append(l)
        self.compounds = result
        self.compounds_count = len(result)

    def setComparedTo(self,smiles):
        self.comparedTo = smiles

    def setTargetName(self,targetName):
        self.targetName = targetName

    def setOrganism(self,organism):
        self.organism = organism

    def setClassification(self,classification):
        self.classification = classification

    def setThreshold(self,threshold):
        self.threshold = threshold

    def setActive(self,active):
        self.active = active

    def setTargetLink(self):
        id = self.target
        self.targetLink = settings.DATASOURCE_LINK_TARGET[self.datasourceName].format(id)

    def setDataSource(self,datasourceName):
        self.datasourceName = datasourceName


class Rnd_set_comparison(models.Model):
    id = models.AutoField('id',primary_key=True,help_text="<< Primary Key >>")
    setsize = models.IntegerField('setsize')
    fp = models.IntegerField('fp')
    threshold = models.FloatField('value', help_text="Property Value", default=0.0)
    rawscore = models.FloatField('rawscore',default=0.0)
    datasource = models.ForeignKey(DataSources, models.CASCADE)

    def __unicode__(self):
        return "id:%d, size:%d, fp:%d, threshold:%f, rawscore:%f" % (self.id,self.setsize,self.fp,self.threshold,self.rawscore)


class Protein_Classifications:
    classifications = {}
    classification_targets_set = {}

    @staticmethod
    def getClassifications():
        chembl_connection = DB_connector(settings.CHEMBL_VERSION)
        prot_classifications = chembl_connection.getProteinClassifications()
        chembl_connection.close()

        for chembl_id,l1,l2 in prot_classifications:
            if l2 is None:
                Protein_Classifications.classifications[chembl_id] = l1
            else:
                Protein_Classifications.classifications[chembl_id] = l2

        Protein_Classifications.classification_targets_set = set(Protein_Classifications.classifications.keys())

    @staticmethod
    def get(target):
        if target in Protein_Classifications.classification_targets_set:
            return Protein_Classifications.classifications[target]
        else:
            return None

    @staticmethod
    def __str__():
        template = "<{0} - {1} Entries>"
        return template.format("Protein_Classifications",len(Protein_Classifications.classification_targets_set))


class FP_Parameter(models.Model):
    id = models.AutoField('id',primary_key=True,help_text="<< Primary Key >>")
    fp_id = models.IntegerField('fp_id')
    datasource = models.ForeignKey(DataSources, models.CASCADE)
    threshold = models.FloatField('value', default=0.0)
    formula_raw_mean = models.CharField('formula_raw_mean',max_length=100)
    formula_raw_stddev = models.CharField('formula_raw_stddev',max_length=100)
    chisquare_mean = models.FloatField('chisquare_mean', default=0.0)
    chisquare_evd = models.FloatField('chisquare_evd', default=0.0)
    date = models.DateTimeField(auto_now_add=True, blank=True)

    vault = {}

    @staticmethod
    def get_from_vault(ds,fp,thresh):
        key = (ds,fp,thresh)
        if not key in FP_Parameter.vault:
            FP_Parameter.vault[key] = FP_Parameter._loadSpecific_(ds,fp,thresh)
        return FP_Parameter.vault[key]

    @staticmethod
    def clear():
        FP_Parameter.vault.clear()

    @staticmethod
    def _loadSpecific_(ds,fp,thresh):
        entries = FP_Parameter.objects.filter(datasource=DataSources.objects.get(name=ds),
                                              fp_id=fp,threshold=thresh)
        if entries.count()==0:
            raise Exception("Didn't found datasource{0},fp{1},thresh{2} in FP_Parameter".format(ds,fp,thresh))
        y = None
        for entry in entries:
            y = entry
            break

        return {'formula_raw_mean':lambda x: eval(y.formula_raw_mean),
                'formula_raw_stddev':lambda x: eval(y.formula_raw_stddev)}

    @DeprecationWarning
    @staticmethod
    def load():
        FP_Parameter.vault = {}
        all_entries = FP_Parameter.objects.all()

        for entry in all_entries:
            ds = settings.DATASOURCES[entry.datasource.name]
            if not ds in FP_Parameter.vault:
                FP_Parameter.vault[ds] = {}
            if not entry.fp_id in FP_Parameter.vault[ds]:
                FP_Parameter.vault[ds][entry.fp_id] = {}
            if not entry.threshold in FP_Parameter.vault[ds][entry.fp_id]:
                FP_Parameter.vault[ds][entry.fp_id][entry.threshold] = {}

            mean_func_str = str(entry.formula_raw_mean)
            mean_func = lambda x: eval(mean_func_str)
            mean_func.func_name = mean_func_str
            std_func_str = str(entry.formula_raw_stddev)
            stddev_func = lambda x: eval(std_func_str)
            stddev_func.func_name = std_func_str

            FP_Parameter.vault[ds][entry.fp_id][entry.threshold] = {'formula_raw_mean':entry.formula_raw_mean,
                                                                    'formula_raw_stddev':entry.formula_raw_stddev,
                                                                    'chisquare_mean':entry.chisquare_mean,
                                                                    'chisquare_evd':entry.chisquare_evd,
                                                                    'date':entry.date,
                                                                    'mean_func':mean_func,
                                                                    'stddev_func':stddev_func}


class All_FPS(models.Model):
    id = models.AutoField('id',primary_key=True,help_text="<< Primary Key >>")
    datasource = models.ForeignKey(DataSources, models.CASCADE)
    fp_id = models.IntegerField('fp_id')
    fp_dict = models.BinaryField('fp_dict')

    def __unicode__(self):
        return self.__str__()
    def __str__(self):
        return str(self.id)


class Target(models.Model):
    chembl_id = models.CharField('chembl_id',primary_key=True,help_text="<< Primary Key >>",max_length=100)
    avr_fp = models.BinaryField('avr_fp')
    num_cmpds = models.IntegerField('num_cmpds')
    name = models.CharField('name',help_text="Name of Target",max_length=100)
    compound_dict = models.BinaryField('avr_fp')
    def __unicode__(self):
        return self.chembl_id


if os.path.exists('ocean/custom_models.py'):
    import custom_models
    mcm = {k:v for k,v in custom_models.__dict__.items() if not k.startswith('__')}
    current_locals = locals()
    for entry,value in mcm.items():
        if not entry in current_locals.keys():
            current_locals.update({entry:value})
            print("monkey patch custom model", entry, value, file=sys.stderr)
