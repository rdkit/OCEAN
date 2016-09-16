from rdkit import Chem

try:
    from ocean.tools.db_connector import DB_connector
except:
    from tools.db_connector import DB_connector

from rdkit.DataStructs import TanimotoSimilarity as tc
import sys

try:
    from rdkit.Chem.Draw import SimilarityMaps as SM
except:
    import ocean.SimilarityMaps as SM
import numpy as np
from ocean.models import Rnd_set_comparison,FP_Parameter,DataSources

import random
from score_calculator import Calculator
try:
    import ocean.settings as settings
except:
    import settings

import os
from scipy.stats import genextreme as ge
from scipy.stats import norm
from django.db import transaction
from multiprocessing import Pool
from math import isinf
from math import isnan

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))

def get_tc(fp0,fp1):
    return tc(fp0,fp1)

def getRandList(data,count):
    result = []
    datacount = len(data)
    while len(result) < count:
        i = random.randint(0,datacount-1)
        result.append(data[i])
    return result

def createRandLists(min,max,iter,molregnos):
    randomLists = []
    for i in range(min,max+1,iter):
        randomLists.append(getRandList(molregnos,i))
    return randomLists

def get_tc_list(fplist0,fplist1,drop_identical=False):
    i = 0
    tcList = [float()]*(len(fplist0)*len(fplist1))
    for fp1 in fplist0:
        if fp1 is None:
            tcList.pop()
            continue
        for fp2 in fplist1:
            if fp2 is None:
                tcList.pop()
                continue
            val = tc(fp1,fp2)
            if drop_identical and val==1.0:
                tcList.pop()
            else:
                tcList[i] = val
                i += 1
    return tcList

def get_tc_list_para(fplist):
    fplist0,fplist1 = fplist
    i = 0
    tcList = [float()]*(len(fplist0)*len(fplist1))
    for fp1 in fplist0:
        if fp1 is None:
            tcList.pop()
            continue
        for fp2 in fplist1:
            if fp2 is None:
                tcList.pop()
                continue
            tcList[i] = tc(fp1,fp2)
            i += 1
    return len(fplist0),np.asarray(tcList)

def get_tc_list_selfsimilarity(fplist):
    tcList = []
    l_fplist = len(fplist)
    for ii in range(l_fplist):
        for jj in range(0,ii):
            tcList.append(tc(fplist[ii],fplist[jj]))
    return tcList

def calc_ocean_parameter(FP_MANAGER, fp, datasource, recalc=False):
    """
    http://www.jamesphoughton.com/2013/08/making-gif-animations-with-matplotlib.html
    """
    print "calcOceanStatistics function start"
    db_ocean = DB_connector("default")  # chembl
    cursor = db_ocean.cursor

    ds = DataSources.objects.get(name=datasource.name)
    if recalc:
        print "delete rnd set items for fp",fp
        Rnd_set_comparison.objects.all().filter(fp=fp).filter(datasource=ds).delete()
        print "done"
    print "delete parameter entries for fp",fp
    FP_Parameter.objects.all().filter(fp_id=fp).filter(datasource=ds).delete()
    print "done"

    if not recalc and Rnd_set_comparison.objects.all().filter(fp=fp).filter(datasource=ds).count()==0:
        return "no entries for fp %d, try ?recalc=True" % fp

    repeats = settings.CALC_OCEAN_PARAMETER_REPEATS

    start = settings.CALC_OCEAN_PARAMETER_START
    end = settings.CALC_OCEAN_PARAMETER_END
    steps = settings.CALC_OCEAN_PARAMETER_STEPS

    thresh_start = settings.CALC_OCEAN_PARAMETER_THRESH_START
    thresh_end   = settings.CALC_OCEAN_PARAMETER_THRESH_END
    thresh_steps = settings.CALC_OCEAN_PARAMETER_THRESH_STEPS

    animatedGif = True

    try:
        from PIL import Image
        from images2gif import writeGif
    except:
        print >> sys.stderr, "Couldn't import Image from PIL or writeGif from images2gif, so plotting is deactivated now"
        animatedGif = False

    plotting = True
    try:
        import matplotlib.pyplot as plt
    except:
        plotting = True
        animatedGif = False

    processes = settings.PARALLEL_PROCESSES
    if recalc: walker = Pool(processes=processes)

    thresh_list = np.arange(thresh_start,thresh_end,thresh_steps)
    molecule_ids = np.asarray(FP_MANAGER[datasource][fp].keys())

    ds = DataSources.objects.get(name=datasource.name)
    for runde in range(repeats):
        if not recalc: continue

        print "runde %d" % runde
        result = {}
        rand_lists1 = createRandLists(start,end,steps,molecule_ids)
        rand_lists2 = createRandLists(start,end,steps,molecule_ids)

        tasks = [([FP_MANAGER[datasource][fp].get(x1) for x1 in rand_lists1[i]],[FP_MANAGER[datasource][fp].get(x2) for x2 in rand_lists2[i]]) for i in range(len(rand_lists2))]

        if processes>1:
            np.random.shuffle(tasks)
            result2 = {}
            for data_entry in walker.imap_unordered(get_tc_list_para,tasks,20):
                result2[data_entry[0]] = data_entry[1]
                print "addet %d of %d" % (len(result2),len(tasks))
        else:
            result2 = {}
            while (len(tasks)>0):
                task = tasks.pop()
                score = get_tc_list_para(task)
                result2[score[0]] = score[1]
                print "addet %d of %d" % (len(result2),len(tasks))

        print "create %d Result-Objects for DB-Table rnd_set_comparison" % (len(thresh_list) * len(result2))
        with transaction.atomic():
            buffer = []
            for threshold in thresh_list:
                for key,value in result2.iteritems():
                    raw_score = np.sum(value[value>=threshold])
                    item = (key**2,fp,threshold,raw_score)
                    buffer.append(item)
            print "created %d buffered items" % len(buffer)

            for w,x,y,z in buffer:
                obj = Rnd_set_comparison(setsize=w,fp=x,threshold=y,rawscore=z,datasource=ds)
                obj.save()

    figures = []

    data_cache = {}

    min_mean = None
    max_mean = None
    min_stddev = None
    max_stddev = None

    for threshold in thresh_list:
        if db_ocean.db_type=='postgre':
            query = "select setsize,threshold, round(stddev_pop(rawscore)::numeric,2) as stddev_pop,round(avg(rawscore)::numeric,2) as mean from ocean_rnd_set_comparison where fp=%d and threshold=%f and datasource_id=%d group by setsize,threshold order by setsize" % (fp,threshold,ds.id)
        else:
            query = "select setsize,threshold,round(stddev(rawscore),2) as stddev,round(avg(rawscore),2) as mean from ocean_rnd_set_comparison where fp=%d and threshold=%f and datasource_id=%d group by setsize,threshold order by setsize" % (
            fp, threshold, ds.id)
        cursor.execute(query)

        x_data = []
        stddev_data = []
        mean_data = []
        for result in cursor.fetchall():
            x_data.append(float(result[0]))
            mean_data.append(float(result[3]))
            stddev_data.append(float(result[2]))

        if min_mean is None:
            if len(mean_data) > 0:
                min_mean,max_mean = min(mean_data),max(mean_data)
            if len(stddev_data) > 0:
                min_stddev,max_stddev = min(stddev_data),max(stddev_data)
        else:
            if len(mean_data) > 0:
                min_mean, max_mean = min([min_mean,min(mean_data)]), max([max_mean,max(mean_data)])
            if len(stddev_data) > 0:
                min_stddev, max_stddev = min([min_stddev, min(stddev_data)]), max([max_stddev, max(stddev_data)])

        data_cache[threshold] = (x_data,mean_data,stddev_data)

    skip_3_to_6 = True

    for threshold in thresh_list:
        x_data,mean_data,stddev_data = data_cache[threshold]

        if len(x_data) == 0 or len(mean_data)==0 or len(stddev_data)==0:
            continue
        if plotting:
            plt.clf()

        if plotting:
            if skip_3_to_6:
                fig,(r0,r1,r2,r6) = plt.subplots(nrows=4,figsize=(12,14))
            else:
                fig,(r0,r1,r2,r3,r4,r5,r6) = plt.subplots(nrows=7,figsize=(6,14))

        raw_mean_func = Calculator.getRawScoreExpFunction(x_data,mean_data)
        print "\nmean function for threshold: %f is [%s]" % (threshold,raw_mean_func.func_name)

        exp_mean_data = [raw_mean_func(en) for en in x_data]
        if plotting:
            r0.plot(np.array(x_data),np.array(mean_data),linewidth=1.0)
            r0.plot(x_data,exp_mean_data,alpha=0.5,linewidth=2.5)
            r0.set_title("Mean, Threshold: %.2f" % threshold)

            r0.set_ylim((min_mean,max_mean))
            r1.set_ylim((min_stddev,max_stddev))
            r2.set_xlim((-1,1.5))
            r2.set_ylim((0,2.5))
        new_std_function = Calculator.getRawScoreStdDevExpFunction(x_data,stddev_data)

        print "stddev function for threshold: %f is [%s]" % (threshold,new_std_function.func_name)

        newdata2 = new_std_function(x_data)

        if plotting:
            r1.plot(x_data,stddev_data)
            r1.plot(x_data, newdata2, alpha=0.8, linewidth=2.0)
            r1.set_title("StdDev")

        z_Scores = Calculator.getZScores(x_data,mean_data,raw_mean_func,new_std_function)

        histo_bins = 50
        counts,bin_edges = np.histogram(z_Scores,histo_bins,normed=True)
        bin_centres = (bin_edges[:-1] + bin_edges[1:])/2.


        if plotting:
            n,bins,patches = r2.hist(z_Scores,bins=histo_bins,normed=True,alpha=0.5)
            r2.set_title("z-Scores")

        e_val_function = Calculator.getZScoreDistExpFunction(z_Scores)
        e_val_data_x = np.linspace(min(z_Scores),max(z_Scores),num=500)
        e_val_data = [e_val_function(entry) for entry in e_val_data_x]
        if plotting:
            if not skip_3_to_6: r3.plot(e_val_data_x,e_val_data,alpha=0.5)

        c=-0.1
        for c in [-0.05]:
            x_ls = np.linspace(ge.ppf(0.01,c),ge.ppf(0.99,c),100)
            if plotting:
                if not skip_3_to_6: r4.plot(x_ls,ge.pdf(x_ls,c),linewidth=1.6-c*4)

        (shape_evd,loc_evd,scale_evd) = ge.fit(z_Scores)

        loc_norm,scale_norm = norm.fit(z_Scores)
        x = ge.pdf(bin_centres,shape_evd,loc=loc_evd,scale=scale_evd)

        if plotting:
            evd_plot, = r2.plot(bin_centres,x,'b',color='black',label='Extreme Value Distribution')

        ndist = norm.pdf(bin_centres,loc=loc_norm,scale=scale_norm)
        if plotting:
            norm_plot, = r2.plot(bin_centres,ndist,'b',color="red",label='Normal Distribution')
            r2.legend([evd_plot,norm_plot],['Extreme Value Distribution','Normal Distribution'],loc=1)

        def getDecNpArray(value):
            return np.asarray(value).astype(float)

        expected_evd = getDecNpArray(x)
        expected_norm = getDecNpArray(ndist)
        observed = getDecNpArray(counts)

        def normalizedChisquare(observed,expected):
            if len(observed) != len(expected): raise Exception("len of observed and expected has to be the same")

            zipped = zip(observed,expected)
            fun = lambda input: ((input[0]-input[1])**2 / (input[0]+input[1]))
            result = sum(map(fun,zipped))

            return result

        chisq_mean = normalizedChisquare(observed,expected_norm)
        chisq_evd = normalizedChisquare(observed,expected_evd)

        print "chisquare_norm",chisq_mean
        print "chisquare_evd",chisq_evd

        #django doesn't like inf or -inf in float-fields of oracle database, so we change it..
        if isinf(chisq_mean) or isnan(chisq_mean):
            print "chisquare_norm seems to be inf or nan (%s), change to -1.0" % str(chisq_mean)
            chisq_mean = -1.0
        if isinf(chisq_evd) or isnan(chisq_evd):
            print "chisquare_evd seems to be inf or nan (%s), change to -1.0" % str(chisq_evd)
            chisq_evd = -1.0

        if plotting:
            if not skip_3_to_6: n,bins,patches = r5.hist(z_Scores,bins=histo_bins,normed=True,alpha=0.75)#,bins=20)

            if not skip_3_to_6:
                import matplotlib.mlab as mlab
                y = mlab.normpdf(bins,loc_evd,scale_evd)

        fp_parameter = FP_Parameter(fp_id=fp,
                                    threshold=threshold,
                                    formula_raw_mean=raw_mean_func.func_name,
                                    formula_raw_stddev=new_std_function.func_name,
                                    chisquare_mean=chisq_mean,
                                    chisquare_evd=chisq_evd,
                                    datasource=ds)
        fp_parameter.save()
        if plotting:
            if not skip_3_to_6: r5.plot(bins,y)

        if threshold==thresh_list[-1]:      #this is last round
            print "last round"

            query = "select threshold,chisquare_mean,chisquare_evd from ocean_fp_parameter where fp_id=%d and datasource_id=%d order by threshold" % (fp,ds.id)
            cursor.execute(query)
            data_chi2_mean = []
            data_chi2_evd = []
            x_chidata = []
            for val in cursor.fetchall():
                x_chidata.append(float(val[0]))
                data_chi2_mean.append(float(val[1]))
                data_chi2_evd.append(float(val[2]))

            print x_chidata,data_chi2_mean,data_chi2_evd

            if plotting:
                if not skip_3_to_6: r6.plot(x_chidata,data_chi2_mean,'o')
                if not skip_3_to_6: r6.plot(x_chidata,data_chi2_evd,'.')
                chi2_mean, = r6.plot(x_chidata,data_chi2_mean,'o')
                chi2_evd, = r6.plot(x_chidata,data_chi2_evd,'.')
                r6.legend([chi2_mean,chi2_evd],['ChiSquare Normal Distribution','ChiSquare Extreme Value Distribution'],loc=1)

        def fitfunc(p,x):
            if p[0]==0:
                return np.exp(-np.exp(-x))*np.exp(-x)
            else:
                print p[0],type(x)
                return np.exp(-(1-p[0]*x)**(1/p[0]))*(1-p[0]*x)**(1/p[0]-1)
        errfunc = lambda p,x,y: (y-fitfunc(p,x))

        init = [0.2]

        bins = bins[:-1]
        bins = np.array(bins)
        n = np.array(n)

        if plotting:
            plt.tight_layout()
            filename = "%f.png" % threshold
            plt.savefig(filename)
            figures.append(filename)

    if animatedGif:
        file_names = figures
        print "d",file_names
        images = [Image.open(fn) for fn in file_names]
        writeGif("animation_mean_stddev.gif",images,duration=0.5)
        for image in images:
            image.close()