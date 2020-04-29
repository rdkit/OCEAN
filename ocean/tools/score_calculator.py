from scipy import stats
import numpy as np
from math import exp,pi,sqrt,e
from decimal import Decimal
import ocean.settings


class Calculator:
    threshold = None
    expected_rawscore_stddev_exponent = None
    expected_rawscore_mean_function = None
    expected_rawscore_stddev_function = None

    def __init__(self,data,target_count,raw_mean_func,raw_mean_std_func,threshold):
        self.byTClist = None
        self.tcList = None
        self.target_count = target_count

        self.threshold = threshold

        self.expected_rawscore_mean_function = raw_mean_func
        self.expected_rawscore_stddev_function = raw_mean_std_func

        self.byTClist = True
        self.tcList = data

        return

    @staticmethod
    def getRawScoreExpFunction(x_data,y_data):
        slope,intercept,r_val,p_val,std_err = stats.linregress(x_data,y_data)
        fun = lambda x: slope*x+intercept
        fun.func_name = "(%f * x) + %f" % (slope,intercept)
        return fun

    @staticmethod
    def getRawScoreStdDevExpFunction(x_data,y_data):
        from scipy.optimize import curve_fit
        def fitfunc2(x,a,b):
            return a*(x**b)

        fp2,fc2 = curve_fit(fitfunc2,x_data,y_data)

        fun2 = lambda x: fp2[0] * (x**fp2[1])
        fun2.func_name = "%f * (x ** %f)" % (fp2[0],fp2[1])

        return fun2

    @staticmethod
    def getZScoreDistExpFunction(y_data):
        from scipy.stats import genextreme as ge
        fit_params = ge.fit(y_data)
        mean = fit_params[1]
        sigma = fit_params[2]
        shape = fit_params[0]
        return lambda x: ge.pdf(x,shape,loc=mean,scale=sigma)

    @staticmethod
    def getZScoreDistExpFunction2(x_data,y_data):
        from scipy.optimize import curve_fit
        from math import exp
        def fitfunc(x,c):
            if c==0:
                return exp(-exp(-x))*exp(-x)
            else:
                return exp(-(1-c*x)**(1/c))*(1-c*x)**(1/c-1)
        fp,fc = curve_fit(fitfunc,x_data,y_data)
        print(fp, fc)
        return lambda x: fp[1]*np.log(x*fp[0]+fp[2])+fp[3]

    @staticmethod
    def getZScores(x_data,raw_data,exp_raw_func,exp_std_func):
        exp_raw = np.asfarray([exp_raw_func(entry) for entry in x_data])
        exp_std = np.asfarray([exp_std_func(entry) for entry in x_data])
        result = (np.asfarray(raw_data) - exp_raw) / exp_std
        result = sorted(result,reverse=True)
        return result

    def __str__(self):
        return str(self.names)+str(self.result)

    def calculate(self):
        pass
        self.names = ["n","exp_raw","raw","exp_z","z","p_value"]
        if self.byTClist:
            n = len(self.tcList)

        exp_raw = Decimal(self.expected_rawscore_mean_function(n))

        exp_raw_stddev = Decimal(self.expected_rawscore_stddev_function(n))
        raw = Decimal(self._raw_score_()[0])

        exp_z = exp_raw_stddev
        try:
            z = (raw - exp_raw) / exp_raw_stddev
            p_value = -Decimal(e)**((-z*Decimal(pi))/(Decimal(sqrt(6))-Decimal(0.577215665)))
            p_value = -(p_value+(p_value**2) / 2 + (p_value**3)/6)

            e_value = p_value * self.target_count
        except Exception as err:
            z = None
            p_value = Decimal(0.0)
            e_value = Decimal(0.0)
        self.result = [n,exp_raw,raw,exp_z,z,p_value,e_value]

    def _listTCs_(self,pair):
        from rdkit.DataStructs import TanimotoSimilarity as tc
        result = []
        for x in pair.set1:
            for y in pair.set2:
                result.append(tc(x,y))
        return result
    
    def _raw_score_(self):
        result = []
        if self.byTClist:
            pair_score = sum([value for value in self.tcList if value > self.threshold]) #possible that this is only important if i want to determine the expectation-values?!
            return [pair_score]
        return result

    def _mean_(self):
        return