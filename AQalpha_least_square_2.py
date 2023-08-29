# A,Q,alpha least square optimization
# third optimization step for REBO2 parameters fitting

import bond_order_calculation
import Bbeta_least_square
import REBO
import numpy as np
from math import exp, sin, cos, pi, sqrt
from lmfit import Model
import matplotlib.pyplot as plt


class AQalpha_least_square_minimization(object):

    def __init__(self, objREBO):
        self.training_value2l = []
        self.obj = objREBO
        return

    '''def energy_func(self,x,self.obj,rindex,A,Q,alpha):
        self.obj.p['A'],self.obj.p['Q'],self.obj.p['alpha']=A,Q,alpha
        etot=self.obj.energy_calculation(x,rindex)
        return etot'''

    def training_value2(self, E0, rij, rindex):
        bij = self.obj.p['bij'][rindex]
        return (E0 + bij * self.obj.VA(rij))

    def set_training_values2(self):
        # r is the list defined in initial
        self.training_value2l = []
        for i in range(0, 6):
            # for i in range(0, self.obj.r.shape[0]):
            self.training_value2l.append(
                self.training_value2(self.obj.E[i], self.obj.r[i], i))
        return

    @staticmethod
    def exponential1(x, A, alpha):
        results = []
        for x1 in x:
            '''if x1 < 2.5:
                f_2deriv = 0
            if x1 > 2.5 and x1 < 3.05:
                Rij_1, Rij_2 = 2.5, 3.05
                f_2deriv = -1. / 2 * (pi / (Rij_2 - Rij_1))**2 * cos(pi * (x1 - Rij_1) /
                                                                    (Rij_2 - Rij_1))'''
            Q = 0
            f_2deriv = 0
            results.append(A * exp(-alpha * x1))
        return np.asarray(results)

    def AQalpha_ls_min(self):
        xdata = self.obj.r[0:6]
        # ydata_E=np.asarray(self.obj.E)
        ydata_k = np.asarray(self.training_value2l)
        mod = Model(self.exponential1, prefix='f1_')
        mod.set_param_hint('f1_alpha', value=self.obj.p[
                           'alpha'], min=0)
        mod.set_param_hint('f1_A', value=self.obj.p['A'], min=0)
        params = mod.make_params()
        results = mod.fit(ydata_k, params, x=xdata)
        #f = open(self.obj.outfile, 'a')
        #f.write("AQalpha_analysis \n")
        # f.write(results.fit_report())
        # f.close()
        '''plt.plot(xdata, ydata_k, 'bo', label='data')
        plt.plot(xdata, results.init_fit, 'k--', label='init_fit')
        plt.plot(xdata, results.best_fit, 'r-', label='best_fit')
        plt.legend()
        plt.show()'''
        # print(results.fit_report())
        self.obj.p['A'], self.obj.p[
            'alpha'] = float(results.best_values['f1_A']),\
            float(results.best_values['f1_alpha'])
        return
