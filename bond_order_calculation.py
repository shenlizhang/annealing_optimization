# calculate bond order parameter
# second optimization step for REBO2 parameters fitting

import REBO
import Bbeta_least_square
from math import exp


class Bond_order_calculation(object):

    def __init__(self, objREBO, objBbeta):
        self.obj = objREBO
        self.objBbeta = objBbeta
        return

    def normalization(self):
        B1, B2, B3 = self.obj.p['B1'], self.obj.p['B2'], self.obj.p['B3']
        beta1, beta2, beta3 = self.obj.p[
            'beta1'], self.obj.p['beta2'], self.obj.p['beta3']
        rn=2.22
        term1 = self.objBbeta.get_denominator(rn)
        b=B1 * beta1 * exp(-beta1 * rn) + B2 * beta2 * \
                exp(-beta2 * rn) + B3 * beta3 * exp(-beta3 * rn)
        k=-term1/b
        self.obj.p['B1']=B1*k
        self.obj.p['B2']=B2*k
        self.obj.p['B3']=B3*k
        return

    def bond_order_calculation(self):
        self.obj.p['bij']=[]
        B1, B2, B3 = self.obj.p['B1'], self.obj.p['B2'], self.obj.p['B3']
        beta1, beta2, beta3 = self.obj.p[
            'beta1'], self.obj.p['beta2'], self.obj.p['beta3']
        for i in range(0, len(self.obj.r)):
            rn = self.obj.r[i]
            term1 = self.objBbeta.get_denominator(rn)
            term2 = self.obj.f_derivative(rn)
            term3 = self.obj.fij(rn)
            a = B1 * exp(-beta1 * rn) + B2 * \
                exp(-beta2 * rn) + B3 * exp(-beta3 * rn)
            b = B1 * beta1 * exp(-beta1 * rn) + B2 * beta2 * \
                exp(-beta2 * rn) + B3 * beta3 * exp(-beta3 * rn)
            #print rn, term1,term2, term3
            bij = term1 / (term2 * a - term3 * b)
            self.obj.p['bij'].append(bij)
            print bij
        #print self.obj.p
        return

    '''def bond_order_calculation(self):
        self.obj.p['bij']=[]
        B1, B2, B3 = self.obj.p['B1'], self.obj.p['B2'], self.obj.p['B3']
        beta1, beta2, beta3 = self.obj.p[
            'beta1'], self.obj.p['beta2'], self.obj.p['beta3']
        for i in range(0, len(self.obj.r)):
            rn = self.obj.r[i]
            E=self.obj.E[i]
            VR=self.obj.VR(rn)
            VA=self.obj.VA(rn)
            bij = (VR-E)/VA
            self.obj.p['bij'].append(bij)
            print bij
        #print self.obj.p
        return'''
