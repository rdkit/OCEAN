from django.shortcuts import render_to_response
from django.http import HttpResponse

from rdkit import Chem
from tools import ocean_kit
import base64
try:
    from rdkit.Chem.Draw import SimilarityMaps as SM
except:
    try:
        import SimilarityMaps as SM
    except:
        import ocean.SimilarityMaps as SM
import cPickle
import time
from tools.score_calculator import Calculator
from ocean.models import *
from multiprocessing import Pool, Process, Queue, Lock, Manager, Pipe
from xml.etree.ElementTree import Element,SubElement,tostring
from rdkit.Chem import Draw
from collections import deque
import random

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
manager = Manager()

FP_VALIDATION_PROCESS = []
FP_MANAGER = FP_Manager() # for fps from database
TARGET_CMPDS = manager.dict()
PROTEIN_CLASSIFICATIONS = {}
TARGET_THRESHOLDS = {}
COMPOUND_DEPICTIONS = deque(maxlen=1000)
COMPOUND_DEPICTIONS_lock = Lock()

FP_VAULT = {} # for created fingerprints

searchBuffer = deque(maxlen=30)
searchBuffer_lock = Lock()

def loadFPS(fp=None,datasource=None):

    if datasource is None:
        for datasource in settings.DATASOURCES:
            if fp is None:
                fp = settings.SCORING_PARAMS[datasource.name]['FP']
            print "recursive call loadFPS for fp {} and datasource {}".format(fp, datasource)
            loadFPS(fp,datasource)
    else:
        dn = datasource if type(datasource) is str else datasource.name
        if fp is None:
            fp = settings.SCORING_PARAMS[dn]['FP']
        print "load FP {} for datasource {}".format(fp, datasource)
        FP_MANAGER.loadFPsFromDataSource(datasource,fp)

def init():
    print "databases",settings.DATABASES
    print >>sys.stderr,"1: Load Fingerprints"
    loadFPS()

    print >>sys.stderr,"3: Load Protein-Classifications"
    Protein_Classifications.getClassifications()

    print >>sys.stderr,"4: Load Target-Compound Relationships"
    Target_Compounds.fill()

def createFPforEntry(data):
    id,smiles,fp=data[0],data[1],data[2]
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return id,None
    fingerprint = settings.FP_METHODS[fp](mol)
    return id,fingerprint

def createAllFPSforAllMolregnos(request):
    print "createAllFPSforMolregnos requested"

    result_dict = {}
    fp = None
    datasource = None

    if u'fp' in request.GET and request.GET.get(u'fp')!='':
        fp = request.GET.get(u'fp')
    if u'datasource' in request.GET and request.GET.get(u'datasource')!='':
        datasource = request.GET.get(u'datasource')
    
    if u'fp' in request.POST and request.POST.get(u'fp')!='':
        fp = request.POST.get(u'fp')
    if u'datasource' in request.POST and request.POST.get(u'datasource')!='':
        datasource = request.POST.get(u'datasource')

    if fp is None or datasource is None:
        return HttpResponse('no fp or no datasource defined')

    fp = int(fp)
    datasource = str(datasource)
    datasource = settings.DATASOURCES[datasource]

    All_FPS.objects.filter(fp_id = fp).filter(datasource=DataSources.objects.get(name=datasource.name)).delete()
    print "fp is",fp

    datasources = [x.name for x in settings.DATASOURCES]

    for ds in datasources:
        current_sources = DataSources.objects.filter(name=ds)
        if current_sources.count()==0:
            print "no DataSource entry in DB for {} .. create entry".format(ds)
            new_entry = DataSources(name=ds)
            new_entry.save()

    dict_pack_size = 20000
    pool = Pool(processes=settings.PARALLEL_PROCESSES)
    print "add Fingerprints {0} to DataSource {1}".format(fp,datasource.name)

    idAndSmiles = DataSources.getDataSourceActivities(datasource)

    new_idAndSmiles = [(x[0],getLargestFragment(x[1]),fp) for x in idAndSmiles]
    i2 = 0
    for datapoint in pool.imap_unordered(createFPforEntry, new_idAndSmiles, 50):
        if i2 % 10000 == 0:
            print i2

        id,fingerprint = datapoint
        if fingerprint is not None:
            result_dict[id] = fingerprint

        if i2 % dict_pack_size == 0 and i2>0:
            print "pack dict into pickle and save into db"
            pickle = cPickle.dumps(result_dict,2)
            entry = All_FPS(fp_dict=pickle,fp_id = fp,datasource=DataSources.objects.get(name=datasource.name))
            entry.save()
            result_dict.clear()

        i2 += 1

    # we have to save the last (smaller than dict_pack_size) chunk too!
    if len(result_dict)>0:
        pickle = cPickle.dumps(result_dict,2)
        entry = All_FPS(fp_dict=pickle,fp_id = fp,datasource=DataSources.objects.get(name=datasource.name))
        entry.save()
        result_dict.clear()

    print "save last fp_dict"

    print "load FPS of Datasource {0} into Memory again".format(datasource.name)
    loadFPS(fp,datasource)

    report = "Saved {0} Items of FP {1} for DataSource {2}".format(i2,fp,datasource.name)

    return None if request is None else HttpResponse(report)


def main_marvin(request):
    datasources = map(lambda x: x.name, settings.DATASOURCES)
    return render_to_response("byCompounds_marvin.html", {"datasources": datasources})

def main_smiles(request):
    datasources = map(lambda x: x.name, settings.DATASOURCES)
    return render_to_response("byCompounds_smiles.html", {"datasources": datasources})

def findNeighbourhoodForCompounds(request):
    req = request.POST.copy()
    req.update(request.GET)

    start = time.time()

    smiles = []
    if 'smiles' in req and req.get('smiles') != '':
        if 'hit_profile' in req:
            smiles.append(base64toSmiles(str(req.get('smiles'))))
        else:
            try:
                m = Chem.MolFromSmiles(str(req.get('smiles')))
                if m is None:
                    raise Exception()
                smi = str(req.get('smiles'))
            except:
                smi_decoded = base64toSmiles(str(req.get('smiles')))
                m = Chem.MolFromSmiles(smi_decoded)
                if m is None:
                    raise Exception("Couldn't Decode Smiles: %s" % str(req.get('smiles')))
                smi = smi_decoded
            smiles.append(smi)
    elif 'molBlock' in req:
        mb = req.get(u'molBlock')
        m = Chem.MolFromMolBlock(str(mb))
        sm = Chem.MolToSmiles(m,isomericSmiles=True,canonical=True)
        smiles.append(sm)
    else:
        return HttpResponse("Please give a Smiles or Draw Molecule with MarvinJS.")

    datasource = req[u'datasource']
    datasource = settings.DATASOURCES[datasource]

    print "smiles: %s" % str(smiles)

    print "datasource: {0}".format(datasource.name)

    print "############"

    searchBuffer_lock.acquire()
    searchBuffer_dict = dict(list(searchBuffer))
    bufferEntry = (smiles[0],datasource.name)

    fp = settings.SCORING_PARAMS[datasource.name]['FP']
    threshold = settings.SCORING_PARAMS[datasource.name]['THRESHOLD']

    if bufferEntry in searchBuffer_dict:
        ranked_targets = searchBuffer_dict[bufferEntry]
    else:
        ranked_targets = getRankedTargetList(smiles[0], verbose=False, datasources=datasource, fp=fp, fp_threshold=threshold)
        searchBuffer.append((bufferEntry,ranked_targets))
    searchBuffer_lock.release()

    class_counter = {}
    for entry in ranked_targets[:30]:
        if entry.classification in class_counter:
            class_counter[entry.classification] += 1
        else:
            class_counter[entry.classification] = 1

    class_counter_sorted = sorted(class_counter.items(), key=lambda x: x[1], reverse=True)

    if "hit_profile" in req:
        # TODO replace XML with simple json
        totals = sum(map(lambda x: x[1],class_counter_sorted))
        root = Element('root')
        for i in range(len(class_counter_sorted[:4])):
            child_prop = SubElement(root,"e_%d" % i)
            percent = int(round((float(class_counter_sorted[i][1]) / float(totals))*100))
            child_prop.text = "%s: %s %%" % (str(class_counter_sorted[i][0]),str(percent))
        root = '<?xml version="1.0" encoding="UTF-8"?>' + tostring(root)
        return HttpResponse(root,content_type="text/xml")

    print "ranked by compare {0} Targets".format(len(ranked_targets))

    b64smiles = "".join(smilesTobase64(smiles[0]).splitlines())
    print "b64_smiles [%s]" % b64smiles
    print "request finished:",time.time() - start
    return render_to_response("searchResults.html",{"results" : ranked_targets,
                                                    "smiles"  : smiles[0],
                                                    "b64smiles": b64smiles,
                                                    "class_counter_sorted": class_counter_sorted,
                                                    "datasource":datasource.name})

def getTC(smiles,molecule_id,fp,datasource):
    smiles = getLargestFragment(smiles)
    if not datasource in FP_VAULT:
        FP_VAULT[datasource] = {}
    if not fp in FP_VAULT[datasource]:
        FP_VAULT[datasource][fp] = {}
    if smiles in FP_VAULT[datasource][fp]:
        fpThis = FP_VAULT[datasource][fp][smiles]
    else:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print >>sys.stderr, "couldn't create molecule of smiles",smiles
        fpThis = settings.FP_METHODS[fp](mol)
        FP_VAULT[datasource][fp][smiles] = fpThis

    fpThat = FP_MANAGER[datasource][fp][molecule_id]
    tc = ocean_kit.get_tc(fpThis,fpThat)
    return tc

def base64toSmiles(base64string):
    smiles = base64.decodestring(base64string)
    return smiles

def smilesTobase64(smiles):
    base64string = base64.encodestring(smiles)
    return base64string

def getCmpdsForTargetOfDatasource(target,datasource):
    if datasource is settings.DATASOURCES.CHEMBL:
        chembl_connection = DB_connector(settings.CHEMBL_VERSION)
        data = chembl_connection.getInfoForTargetAndCompound(-1, target)
        chembl_connection.close()
        return data
    else:
        raise Exception("unknown Datasource {}, try to MonkeyPatch this function and add own Datasource".format(datasource))

def getCmpdsForTarget(request):
    if 'target' in request.POST:
        target = str(request.POST['target'])
        smiles = str(request.POST['smiles'])
        smiles = base64toSmiles(smiles)
    elif u'target' in request.POST:
        target = str(request.POST[u'target'])
        smiles = str(request.POST[u'smiles'])
        smiles = base64toSmiles(smiles)

    datasource = request.GET[u'datasource'] if u'datasource' in request.GET else request.POST[u'datasource']
    datasource = settings.DATASOURCES[datasource]

    fp_tmp = settings.SCORING_PARAMS[datasource.name]['FP']
    data = getCmpdsForTargetOfDatasource(target,datasource)

    # TODO replace XML with simple json
    root = Element('root')
    attributes = ["molregno","molecule_chembl_id","target_pref_name","organism","canonical_smiles","standard_value","tc","b64smiles"]
    for entry in data:
        if entry[0] in settings.drop_compounds:
            continue
        child = SubElement(root,'entry')
        for i in range(len(attributes)):
            child_prop = SubElement(child,str(attributes[i]))
            if attributes[i] == "molecule_chembl_id":
                link = settings.DATASOURCE_LINK_COMPOUND[datasource.name].format(entry[i])
                # print "link is",link
                child_prop.text = "<a href='{0}'>{1}</a>".format(link,entry[i])
            elif attributes[i] == "molregno":
                continue
            elif attributes[i] == "tc":
                child_prop.text = str(getTC(smiles,entry[0],fp_tmp,datasource))
            elif attributes[i] == "b64smiles":
                b64smiles = "".join(smilesTobase64(entry[4]).splitlines())
                child_prop.text = str(b64smiles)
            else:
                child_prop.text = str(entry[i])

    root = '<?xml version="1.0" encoding="UTF-8"?>' + tostring(root)
    return HttpResponse(root,content_type="text/xml")

def calcOceanHits_parallel(pipe, target_jobs, fp, fp_thresh, datasource):
    target_job_info = []
    for target,target_cmpds in target_jobs.items():
        fps = [FP_MANAGER[datasource][fp].get(molecule_id) for molecule_id in target_cmpds]
        valid_cmpds = sum(map(lambda x: x is not None,fps))
        target_classification = Protein_Classifications.get(target)
        entry = (target,valid_cmpds,target_classification,fps)
        target_job_info.append(entry)
    while True:
        try:
            data = pipe.recv()
        except KeyboardInterrupt:
            time.sleep(0.8*random.random())
            pipe.close()
            print "Close Thread"
            return

        if data == 'die':
            break
        smiles,query_fp = data
        result = []
        for job in target_job_info:
            target,valid_cmpds,target_classification,fps = job
            ocean_hit = Ocean_hit(target)
            ocean_hit.setTClist(ocean_kit.get_tc_list(query_fp,fps,drop_identical=settings.VALIDATING_PROCESS))
            ocean_hit.setComparedTo(smiles)
            ocean_hit.setTargetName(target)
            ocean_hit.setClassification(target_classification)
            ocean_hit.setDataSource(datasource.name)
            ocean_hit.setTargetLink()
            result.append(ocean_hit)
        pipe.send(result)
    pipe.close()
    return

def getRankedTargetList(smiles,
                        targetList=None,
                        threads=settings.PARALLEL_PROCESSES,
                        cutoff=settings.CMPD_COUNT_CUTOFF,
                        cmpd_nm_cutoff=settings.CMPD_NM_CUTOFF,
                        verbose=False,
                        fp=None,
                        fp_threshold=None,
                        datasources=None):

    if datasources is None:
        datasources = list(settings.DATASOURCES)
    else:
        if not type(datasources) is list:
            datasources = [datasources]

    if datasources is None:
        datasources = list(settings.DATASOURCES)
    else:
        if not type(datasources) is list:
            datasources = [datasources]

    if len(ProcessManager.pm) == 0 or True:
        for datasource in settings.DATASOURCES:
            if datasource not in datasources or datasource in ProcessManager.datasources:
                continue
            if not fp in FP_MANAGER[datasource]:
                print "Fingerprint %d is not in FP_MANAGER " % fp, FP_MANAGER[datasource].keys()
                loadFPS(datasource=datasource, fp=fp)
                print "Fingerprint-Set %d loaded into %s" % (fp,datasource.name)

            if targetList is None:
                targetList = Target_Compounds.getTargets(datasource)
            ti = []
            total_cmpds = 0
            for target in targetList:
                target_info = Target_Compounds.vault[datasource].get(target)
                target_cmpds = target_info['cmpds']
                ti.append((target,len(target_cmpds)))
                total_cmpds += len(target_cmpds)

            unallocated_cmpds = total_cmpds
            cmpds_per_thread = float(total_cmpds) / settings.PARALLEL_PROCESSES
            ti.sort(key=lambda x: x[1],reverse=True)
            thread_targets = []
            thread_chunk = []
            thread_chunk_size = 0
            for ti_entry in ti:
                if ti_entry[1] < cutoff:
                    print "target {0} has to few cmpds {1}".format(ti_entry[0], ti_entry[1])
                thread_chunk.append(ti_entry[0]) #add target_id to chunk
                thread_chunk_size += ti_entry[1] #increads current chunk_size
                unallocated_cmpds -= ti_entry[1]
                if thread_chunk_size >= cmpds_per_thread:
                    thread_targets.append(thread_chunk)
                    thread_chunk = []
                    thread_chunk_size = 0
            if thread_chunk_size > 0:
                thread_targets.append(thread_chunk)
                thread_chunk = []
                thread_chunk_size = 0

            for thread_chunk in thread_targets:
                process_job = {}
                for target in thread_chunk:
                    target_cmpds = Target_Compounds.vault[datasource].get(target)['cmpds']
                    process_job[target] = target_cmpds
                (p1, p2) = Pipe(True)
                process = Process(target=calcOceanHits_parallel, args=(p1, process_job, fp, fp_threshold, datasource))
                p = ProcessManager(process, p1, p2, datasource)
            ProcessManager.datasources.add(datasource)

            print "start {0} Processes".format(len(ProcessManager.pm))
            ProcessManager.start_all()

    smiles = getLargestFragment(smiles)
    query_fp = [settings.FP_METHODS[fp](Chem.MolFromSmiles(smiles))]

    tii = time.time()
    ProcessManager.send_all((smiles,query_fp),datasources)
    result = ProcessManager.recv_all(datasources)
    print "duration",time.time()-tii

    tmp_cache = {}
    fp_tmp = None
    thresh_tmp = None

    for entry in result:
        entry.datasource = settings.DATASOURCES[entry.datasourceName]
        if fp_tmp is None and fp is None:
            fp_tmp = settings.SCORING_PARAMS[entry.datasource.name]['FP']
        else:
            fp_tmp = fp
        if thresh_tmp is None and fp_threshold is None:
            thresh_tmp = settings.SCORING_PARAMS[entry.datasource.name]['THRESHOLD']
        else:
            thresh_tmp = fp_threshold
        calc_functions = FP_Parameter.get_from_vault(entry.datasourceName,fp_tmp,thresh_tmp)


        c = Calculator(entry.tclist,
                       Target_Compounds.counts[entry.datasourceName],
                       calc_functions['formula_raw_mean'],
                       calc_functions['formula_raw_stddev'],
                       thresh_tmp)

        c.calculate()
        scores = c.result
        entry.p_value = c.result[-2]
        entry.e_value = c.result[-1]
        entry.e_valuestr = "{:.3g}".format(entry.e_value)
        entry.setTargetName(Target_Compounds.vault[entry.datasource].get(entry.target)['desc'])

    result.sort()
    return result


##
## WE USE THIS TO COLLECT CHEMBL_DATA FOR OCEAN_DB:
##
# create table tmp_transforms2 as
# (select target_chembl_id,target_pref_name,organism,molregno,molecule_chembl_id,avg as standard_value,standard_relation,standard_units,canonical_smiles
# from
#  (select target_chembl_id,target_pref_name,organism,molregno,molecule_chembl_id,canonical_smiles,standard_relation,standard_units,count(*) as cnt,stddev(standard_value) as std,avg(standard_value) as avg,min(standard_value) as min, max(standard_value) as max, avg(standard_value)-3*stddev(standard_value) as allowed_min, avg(standard_value)+3*stddev(standard_value) as allowed_max
#  from (select /*+ PARALLEL(target_dictionary,4) PARALLEL(assays,4) PARALLEL(molecule_dictionary,4) PARALLEL(compound_properties,4) PARALLEL(compound_structures,4) */ td.chembl_id as target_chembl_id,td.pref_name as target_pref_name,td.organism,md.molregno,md.chembl_id as molecule_chembl_id,act.standard_value,act.standard_relation,act.standard_units,act.standard_type,cs.canonical_smiles as canonical_smiles from target_dictionary td join assays a on td.tid=a.tid join activities act on a.assay_id=act.assay_id join molecule_dictionary md on act.molregno=md.molregno join compound_properties cp on md.molregno=cp.molregno join compound_structures cs on md.molregno=cs.molregno where standard_relation='=' and standard_units='nM' and standard_type in ('IC50','Ki') and target_type in ('SINGLE PROTEIN','PROTEIN COMPLEX') and organism='Homo sapiens')
#  group by molregno,molecule_chembl_id,target_chembl_id,target_pref_name,organism,standard_relation,standard_units,canonical_smiles)
# where std<=avg and min>=allowed_min and max<=allowed_max and avg<10000);
# create unique index tmp_t_idx1 on tmp_transforms (target_chembl_id,molregno);
# create index tmp_t_idx2 on tmp_transforms (target_chembl_id);
# create index tmp_t_idx3 on tmp_transforms (molregno);
# create index tmp_t_idx4 on tmp_transforms (target_chembl_id,standard_value);
# create index tmp_t_idx5 on tmp_transforms (canonical_smiles);
#
# we want only targets ('SINGLE PROTEIN','PROTEIN COMPLEX')
# delete from tmp_transforms where target_chembl_id in (select chembl_id from target_dictionary where target_type not in ('SINGLE PROTEIN','PROTEIN COMPLEX'));
#
# delete from tmp_transforms2 where target_chembl_id in (select target_chembl_id from (select target_chembl_id,count(molregno) as cs from tmp_transforms group by target_chembl_id) where cs<10)

def getLargestFragment(smiles):
    largest = sorted(smiles.split('.'),key=lambda x: len(x),reverse=True)[0]
    return largest

def calcOceanStatistics(request):
    fp = None
    datasource = None

    if u'fp' in request.GET and request.GET.get(u'fp')!='':
        fp = request.GET.get(u'fp')
    if u'datasource' in request.GET and request.GET.get(u'datasource')!='':
        datasource = request.GET.get(u'datasource')

    if u'fp' in request.POST and request.POST.get(u'fp')!='':
        fp = request.POST.get(u'fp')
    if u'datasource' in request.POST and request.POST.get(u'datasource')!='':
        datasource = request.POST.get(u'datasource')

    if fp is None or datasource is None:
        return HttpResponse('no fp or no datasource defined')

    fp = int(fp)
    datasource = str(datasource)
    datasource = settings.DATASOURCES[datasource]

    recalc = True if request.GET.get(u'recalc') == 'True' else False

    print "statistical recalculation requested for FP {0} of DataSource {1}".format(fp,datasource.name)
    cpd_fps = {}
    i = 0

    loadFPS(fp,datasource)
    ocean_kit.calc_ocean_parameter(FP_MANAGER,fp,datasource,recalc=recalc)
    FP_Parameter.clear()
    return HttpResponse("statistical recalculation for FP {0} of DataSource {1} is done.".format(fp,datasource.name))

def png_for_smiles(request):
    width = 200
    height = 200
    smiles = base64toSmiles(str(request.GET['smiles']))
    if 'width' in request.GET:
        width = int(str(request.GET['width']))
    if 'height' in request.GET:
        height = int(str(request.GET['height']))

    smiles = getLargestFragment(smiles)

    key = "%s[w=%d][h=%d]" % (smiles,width,height)
    COMPOUND_DEPICTIONS_dict = dict(list(COMPOUND_DEPICTIONS))
    COMPOUND_DEPICTIONS_lock.acquire()
    if key in COMPOUND_DEPICTIONS_dict:
        response = COMPOUND_DEPICTIONS_dict[key]
    else:
        mol = Chem.MolFromSmiles(smiles)
        img = Draw.MolToImage(mol,size=(width,height))
        response = HttpResponse(content_type="image/png")
        img.save(response,"PNG")
        COMPOUND_DEPICTIONS.append((key,response))

    COMPOUND_DEPICTIONS_lock.release()
    return response

if os.path.exists('ocean/custom_views.py'):
    import custom_views

    mcv = {k:v for k,v in custom_views.__dict__.items() if not k.startswith('__')}
    current_locals = locals()
    for entry,value in mcv.items():
        if not entry in current_locals.keys() or \
                        entry in ['init',
                                  'getCmpdsForTargetOfDatasource']:
            current_locals.update({entry:value})
            print >> sys.stderr, "monkey patch custom view",entry,value

init()