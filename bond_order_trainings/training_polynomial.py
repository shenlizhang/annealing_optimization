# python script to train three body terms
# Shenli, 02/23/2018
from math import pi, cos
from lmfit import Model
import numpy as np
import matplotlib.pyplot as plt


def fik(rik, R1, R2):
    fvalue = []
    for item in rik:
        if item < R1:
            fvalue.append(1)
        if R1 <= item <= R2:
            results = 1. / 2 * \
                (1 + cos(pi * (item - R1) / (R2 - R1)))
            fvalue.append(results)
        if item > R2:
            fvalue.append(0)
    fvalue = np.asarray(fvalue)
    return fvalue


def gfunc(theta, a1, a2, a3):
    x = cos(theta)
    f = a1 * x**2 + a2 * x + a3
    return f


def ksi(theta, rik, rij, lamda, a1, a2, a3, R1, R2):
    # this is just for one molecule; rij is one number; rik could be an array
    s = []
    for i in range(0, rik.shape[0]):
        #print fik(rik[i], R1, R2)
        #print gfunc(theta[i] / 180 * pi, a1, a2, a3)
        s1 = np.dot(fik(rik[i], R1, R2), np.exp(lamda) * np.ones(len(rik[i])) *
                    gfunc(theta[i] / 180 * pi, a1, a2, a3))
        s.append(s1)
    return np.asarray(s)


def bij(theta, rik, rij, lamda, a1, a2, a3, R1, R2):
    bij = (1 + ksi(theta, rik, rij, lamda, a1, a2, a3, R1, R2))**(-1. / 2)
    print "bij", bij
    #print "lamda, a1, a2, a3", lamda, a1, a2, a3
    return bij

# here the training sets are for i=Ge, j=Cl, k=Cl
# GeCl2 and GeCl4; Cl-Ge-Cl bonds
val_bij = np.array([0.8984, 0.9352, 1])
theta = np.array([101.09, 108.69, 109.47])
rij = np.array([2.3, 2.27, 2.22])
rik = np.array([[2.3], [2.27, 2.27], [2.22, 2.22, 2.22]])

# here the training sets are for i=Ge,j=Ge,k=Cl
# beta and n should use Ge-Ge parameters
'''val_bij = np.array([0.9507, 0.9507])
theta = np.array([109.78, 109.78])
rij = np.array([2.50, 2.50])
rik = np.array([[2.23, 2.23, 2.23], [2.23, 2.23, 2.23]])'''

# here the training sets are for i=Ge,j=Cl, k=Ge
'''val_bij = np.array([0.8747, 0.8747])
theta = np.array([104.21, 104.21])
rij = np.array([2.32, 2.32])
rik = np.array([[2.49, 2.49, 2.49], [2.49, 2.49, 2.49]])'''

mod = Model(bij, independent_vars=['theta', 'rik', 'rij'])
mod.set_param_hint('R1', value=2.8, vary=False)
mod.set_param_hint('R2', value=3.1, vary=False)
# here the lamda3 can't be trained if rij=rik
mod.set_param_hint('lamda', value=0, vary=False)
mod.set_param_hint('a1',value=-1)
mod.set_param_hint('a2', value=0)
mod.set_param_hint('a3', value=0)
results = mod.fit(val_bij, a1=-1, a2=0, a3=0, lamda=0,
                  theta=theta,  rik=rik, rij=rij)
print results.best_values
#plt.plot(theta, results.init_fit, 'k--', label='init_fit')
plt.plot(theta, results.best_fit, 'b*-.', label='best_fit')
plt.plot(theta, val_bij, 'ro', label='traing_data')
plt.ylim([0.88, 1.01])
plt.xlabel('j-i-k angles')
plt.ylabel('$b_{ij}$')
plt.legend(loc='best')
plt.savefig("withK_polynominal/fig_ClGeCl_withK.png")
