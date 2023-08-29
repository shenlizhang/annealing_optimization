# least_square_minimization for B1,2,3 and beta1,2,3
# first optimization step for REBO2 parameters fitting

import numpy as np
import REBO
from math import sin, cos, pi, exp
from lmfit import Model
import matplotlib.pyplot as plt


class Bbeta_least_square_minimization(object):

    def __init__(self, objREBO):
        self.training_value = []
        self.obj = objREBO
        return

    def get_numerator(self, E0, rij):
        return (self.obj.VR(rij) - E0)

    def get_denominator(self, rij):
        termd = self.obj.fR(rij)
        return termd

    def training_value1(self, E0, rij):
        energy = self.get_numerator(E0, rij)
        force = self.get_denominator(rij)
        # term = energy * self.obj.fij(rij) / \
        #   (self.obj.f_derivative(rij) * energy - force * self.obj.fij(rij))
        term = -energy / force
        return float(term)

    def set_training_values1(self):
        # r is the list defined in initial
        self.training_value = []
        for i in range(0, 6):
            # for i in range(0, self.obj.r.shape[0]):
            self.training_value.append(
                self.training_value1(self.obj.E[i], self.obj.r[i]))
        return

    @staticmethod
    def exponential(x, B, beta):
        return B * np.exp(-beta * x)

    @staticmethod
    def deriv_exponential(x, B, beta):
        return B * beta * np.exp(-beta * x)

    def Bbeta_ls_min(self):
        x = self.obj.r[0:6]
        y = np.asarray(self.training_value)
        mod = Model(self.exponential, prefix='f1_') / \
            Model(self.deriv_exponential, prefix='f4_')
        mod.set_param_hint('f1_beta', value=self.obj.p['beta1'], min=0)
        params = mod.make_params()
        params.add('f1_B', value=self.obj.p['B1'], min=0)
        params.add('f4_B', expr='f1_B')
        params.add('f4_beta', expr='f1_beta')
        results = mod.fit(y, params, x=x)  # , method='dogbox')
        #f = open(self.obj.outfile, 'a')
        #f.write("Bbeta_analysis \n")
        # f.write(results.fit_report())
        # f.close()
        # print(results.fit_report())
        '''plt.plot(x, y, 'bo', label='data')
        plt.plot(x, results.init_fit, 'k--', label='init_fit')
        plt.plot(x, results.best_fit, 'r-', label='best_fit')
        plt.legend()
        plt.show()'''
        self.obj.p['B1'] = float(results.best_values['f1_B'])
        self.obj.p['beta1'] = float(results.best_values['f1_beta'])
        return
