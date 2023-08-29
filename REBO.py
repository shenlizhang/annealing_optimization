# python script for REBO 2 optimization

import numpy as np
from math import sin, cos, pi, sqrt, exp


class REBO_fit(object):

    def __init__(self, r, E, k, type1, type2, type3=''):
        P_GeGe_t1 = dict(A=1769, Q=0, alpha=2.4451,
                         B1=419.23, B2=0, B3=0, beta1=1.7047, beta2=0, beta3=0,
                         R_1=2.8, R_2=3.1, bij=[])
        '''P_GeGe_t1 = dict(A=0.20698542724759195, Q=3722.292135406406, alpha=1.0364045307513337,
                         B1=193.89666613262503, B2=1.288514053021563e-07, B3=481.0879924731308,
                         beta1=1.3082291094987513, beta2=9.655525184611948, beta3=1.3082290923491513,
                         R_1=2.8, R_2=3.1, bij=[1.0, 0.99267487139592347, 0.99100576076619074,
                                                0.94582102535346879, 0.88477163328350417, 0.87262208944312769])'''
        '''P_GeCl_t1 = dict(A=3580.803, Q=0, alpha=3.2269,
                         bij=[],
                         beta1=1.8545, beta2=0, beta3=0,
                         B1=336.259, B2=0, B3=0, R_1=2.53, R_2=2.83)'''
        P_GeCl_t1 = dict(A=521.8337881789156, Q=0, alpha=2.3248392567132847,
                         bij=[0.99999999999999989, 0.98669870896449574, 0.93523941112754183,
                              0.89841393527772528, 0.88646387005424243, 0.87467275612619122],
                         beta1=0.9857847009534442, beta2=0, beta3=0,
                         B1=62.96734218343653, B2=0, B3=0, R_1=2.53, R_2=2.83)
        '''P_GeCl_t1 = dict(A=521.8337881789156, Q=0, alpha=2.3248392567132847,
                         bij=[1.02333519301, 1.02077058912, 0.615357230431,
                              1.01870596652, 0.883710022841, 1.00882899704],
                         beta1=0.9857847009534442, beta2=0, beta3=0,
                         B1=62.96734218343653, B2=0, B3=0, R_1=2.53, R_2=2.83)'''
        '''P_SiSi_t1 = dict(A=90.1964, Q=15.6614, alpha=2.13083,
                         B1=92.74551, B2=255.329, B3=-3.4026, beta1=1.72687, beta2=1.64617, beta3=132.454,
                         R_1=2.5, R_2=3.05, bij=[])
        P_CC_t1 = dict(A=10953.54416, Q=0.313, alpha=4.7465,
                       B1=12388.79198, B2=17.5674, B3=30.7149, beta1=4.720, beta2=1.433, beta3=1.383,
                       R_1=1.7, R_2=2.0, bij=[])
        P_CC_t1 = dict(A=10, Q=30, alpha=1,
                       B1=300, B2=2, B3=100, beta1=2, beta2=10, beta3=10,
                       R_1=1.7, R_2=2.0, bij=[])'''
        #self.outfile = file
        '''P_GeGeGe_t34 = dict(c=0, d=0.1702, h=-0.43884, 'alpha'=0, 'beta'=0)
		P_GeGeCl_t34 = dict('c'=0.0027508, 'd'=0.38211, 'alpha'=6.8708, 'beta'=3)
		P_ClClCl_t34 = dict('c'=4.0, 'd'=0.0, 'h'=10, 'alpha'=6, 'beta'=1)
		P_ClGeCl_t34 = dict('c'=0.0027508, 'd'=0.38211, 'alpha'=6.8708, 'beta'=3)
		P_ClGeGe_t34 = dict('c'=0, 'd'=0.1702, 'h'=-0.43884, 'alpha'=0, 'beta'=0)
		P_GeClGe_t34 = dict('c'=1.28, 'd'=1.000, 'h'=-1.000, 'alpha'=3, 'beta'=1)
		P_GeClCl_t34 = dict('c'=4.0, 'd'=0.0, 'h'=10, 'alpha'=6, 'beta'=1)'''

        self.r, self.E, self.k = r, E, k
        self.typei, self.typej, self.typek = type1, type2, type3

        if self.typei == 1 and self.typej == 2:
            self.p = P_GeCl_t1
        if self.typei == 1 and self.typej == 1:
            self.p = P_GeGe_t1
        '''if self.typei != self.typej:
			self.p = P_GeCl_t1
		if self.typei == 2 and self.typej == 2:
			self.p = P_ClCl_t1
		if typek == 1:
			if self.typei == 1 and self.typej == 1:
				self.p = P_GeGeGe_t34
			if self.typei == 1 and self.typej == 2:
				self.p = P_ClGeGe_t34
			if self.typei == 2 and self.typej == 1:
				self.p = P_GeClGe_t34
			if self.typei == 2 and self.typej == 2:
				self.p = P_ClClCl_t34
		if typek == 2:
			if self.typei == 1 and self.typej == 1:
				self.p = P_GGeCl_t34
			if self.typei == 2 and self.typej == 2:
				self.p = P_ClClCl_t34
			if self.typei == 1 and self.typej == 2:
				self.p = P_ClGeCl_t34
			if self.typei == 2 and self.typej == 1:
				self.p = P_GeClCl_t34'''
        return

    def fij(self, rij):
        Rij_1, Rij_2 = self.p['R_1'], self.p['R_2']
        if rij < Rij_1:
            return 1
        if Rij_1 <= rij <= Rij_2:
            results = 1. / 2 * \
                (1 + cos(pi * (rij - Rij_1) / (Rij_2 - Rij_1)))
            return results
        if rij > Rij_2:
            return 0

    def f_derivative(self, rij):
        Rij_1, Rij_2 = self.p['R_1'], self.p['R_2']
        if Rij_1 <= rij <= Rij_2:
            results = -1. / 2 * pi / (Rij_2 - Rij_1) * sin(pi * (rij - Rij_1) /
                                                           (Rij_2 - Rij_1))
            return results
        else:
            return 0

    def f_2nd_derivative(self, rij):
        Rij_1, Rij_2 = self.p['R_1'], self.p['R_2']
        if Rij_1 <= rij <= Rij_2:
            results = -1. / 2 * (pi / (Rij_2 - Rij_1))**2 * cos(pi * (rij - Rij_1) /
                                                                (Rij_2 - Rij_1))
            return results
        else:
            return 0

    def VR(self, rij):
        A, Q, alpha = self.p['A'], self.p['Q'], self.p['alpha']
        return self.fij(rij) * (1 + Q / rij) * A * exp(-alpha * rij)

    def VA(self, rij):
        B1, B2, B3, beta1, beta2, beta3 = self.p['B1'], self.p[
            'B2'], self.p['B3'], self.p['beta1'], self.p['beta2'], self.p['beta3']
        f = self.fij(rij)
        return f * (B1 * exp(-beta1 * rij) + B2 * exp(-beta2 * rij) + B3 * exp(-beta3 * rij))

    def energy_calculation(self, rij, rindex):
        E = self.VR(rij) - self.p['bij'][rindex] * self.VA(rij)
        return E

    def fA(self, rij):
        B1, B2, B3, beta1, beta2, beta3 = self.p['B1'], self.p[
            'B2'], self.p['B3'], self.p['beta1'], self.p['beta2'], self.p['beta3']
        f = self.fij(rij)
        df = self.f_derivative(rij)
        '''fa = df * (B1 * exp(-beta1 * rij) + B2 * exp(-beta2 * rij) + B3 * exp(-beta3 * rij)) -\
            f * (B1 * beta1 * exp(-beta1 * rij) + B2 * beta2 *
                 exp(-beta2 * rij) + B3 * beta3 * exp(-beta3 * rij))'''
        fa = -f * (B1 * beta1 * exp(-beta1 * rij) + B2 * beta2 *
                   exp(-beta2 * rij) + B3 * beta3 * exp(-beta3 * rij))
        return fa

    def fR(self, rij):
        A, Q, alpha = self.p['A'], self.p['Q'], self.p['alpha']
        f = self.fij(rij)
        df = self.f_derivative(rij)
        '''fr = df * (1 + Q / rij) * A * exp(-alpha * rij) - f * (Q / rij**2) * A * \
            exp(-alpha * rij) - f * alpha * \
            (1 + Q / rij) * A * exp(-alpha * rij)'''
        fr = - f * (Q / rij**2) * A * exp(-alpha * rij) - f * \
            alpha * (1 + Q / rij) * A * exp(-alpha * rij)
        return fr

    def force_constant_R(self, rij):
        A, Q, alpha = self.p['A'], self.p['Q'], self.p['alpha']
        f_2deriv = self.f_2nd_derivative(rij)
        f_deriv = self.f_derivative(rij)
        f = self.fij(rij)
        '''VR_2deriv = f_2deriv * (1 + Q / rij) * A * exp(-alpha * rij) - 2 * f_deriv * Q / rij**2 * A * exp(-alpha * rij) +\
            2 * f * Q / rij**3 * A * exp(-alpha * rij) + 2 * alpha * f * Q / rij**2 * A * exp(-alpha * rij) -\
            2 * alpha * f_deriv * (1 + Q / rij) * A * exp(-alpha * rij) + \
            alpha**2 * f * (1 + Q / rij) * A * exp(-alpha * rij)'''
        VR_2deriv = 2 * f * Q / rij**3 * A * exp(-alpha * rij) + 2 * alpha * f * Q / rij**2 * \
            A * exp(-alpha * rij) + alpha**2 * f * \
            (1 + Q / rij) * A * exp(-alpha * rij)
        return VR_2deriv

    def force_constant_A(self, rij, rindex):
        bij = self.p['bij'][rindex]
        B1, B2, B3, beta1, beta2, beta3 = self.p['B1'], self.p[
            'B2'], self.p['B3'], self.p['beta1'], self.p['beta2'], self.p['beta3']
        f_2deriv = self.f_2nd_derivative(rij)
        f_deriv = self.f_derivative(rij)
        f = self.fij(rij)
        a = B1 * exp(-beta1 * rij) + B2 * \
            exp(-beta2 * rij) + B3 * exp(-beta3 * rij)
        b = B1 * beta1 * exp(-beta1 * rij) + B2 * beta2 * \
            exp(-beta2 * rij) + B3 * beta3 * exp(-beta3 * rij)
        c = B1 * beta1**2 * exp(-beta1 * rij) + B2 * beta2**2 * \
            exp(-beta2 * rij) + B3 * beta3**2 * exp(-beta3 * rij)
        #VA_2deriv = -bij * (f_2deriv * a - 2 * f_deriv * b + f * c)
        VA_2deriv = -bij * f * c
        return VA_2deriv

    def force_constant_calculation(self, rij, rindex):
        VR_2deriv = self.force_constant_R(rij)
        VA_2deriv = self.force_constant_A(rij, rindex)
        return VR_2deriv + VA_2deriv

    def force_calculation(self, rij, rindex):
        bij = self.p['bij'][rindex]
        force = self.fR(rij) - bij * self.fA(rij)
        return force
