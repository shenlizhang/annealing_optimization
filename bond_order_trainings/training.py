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


def gfunc(theta, gamma, c, d, h):
    f = gamma * (1 + c**2 / d**2 - c**2 / (d**2 + (cos(theta) - h)**2))
    return f


def ksi(theta, rik, rij, gamma, c, d, h, lamda3, m, R1, R2):
    # this is just for one molecule; rij is one number; rik could be an array
    s = []
    for i in range(0, rik.shape[0]):
        s1 = np.dot(fik(rik[i], R1, R2), np.exp(
            np.power(lamda3 * (rij[i] - rik[i]), m))) * gfunc(theta[i] / 180 * pi, gamma, c, d, h)
        s.append(s1)
    return np.asarray(s)


def bij(beta, n, theta, rik, rij, gamma, c, d, h, lamda3, m, R1, R2):
    bij = (1 + (beta * ksi(theta, rik, rij, gamma, c,
                           d, h, lamda3, m, R1, R2))**n)**(1. / (-2 * n))
    print "bij", bij
    print "beta, gamma, c, d, h, lamda3", beta, gamma, c, d, h, lamda3
    return bij

# here the training sets are for i=Ge, j=Cl, k=Cl
# GeCl2 and GeCl4; Cl-Ge-Cl bonds
val_bij = np.array([0.7567, 0.84, 1])
theta = np.array([101.09, 108.69, 109.47])
rij = np.array([2.3, 2.27, 2.22])
rik = np.array([[2.3],[2.27, 2.27],[2.22, 2.22, 2.22]])

# here the training sets are for i=Ge,j=Ge,k=Cl
# beta and n should use Ge-Ge parameters
'''val_bij = np.array([0.9507])
theta = np.array([109.78])
rij = np.array([2.50])
rik = np.array([[2.23, 2.23, 2.23]])'''

# here the training sets are for i=Ge,j=Cl, k=Ge
'''val_bij = np.array([0.7978])
theta = np.array([104.21])
rij = np.array([2.32])
rik = np.array([[2.49, 2.49, 2.49]])'''

mod = Model(bij, independent_vars=['theta', 'rik', 'rij'])
mod.set_param_hint('m', value=1, vary=False)
mod.set_param_hint('R1', value=2.8, vary=False)
mod.set_param_hint('R2', value=3.1, vary=False)
# here the lamda3 can't be trained if rij=rik
mod.set_param_hint('lamda3', min=0)
mod.set_param_hint('c', min=0)
mod.set_param_hint('d', min=0)
mod.set_param_hint('h', min=0)
mod.set_param_hint('beta', min=0)
mod.set_param_hint('n', min=0)
results = mod.fit(val_bij, gamma=4, c=0.1, d=1, h=1,lamda3=2, beta=1, n=1,
                  theta=theta,  rik=rik, rij=rij, method='dogbox')
print results.best_values
plt.plot(theta, results.init_fit, 'k--', label='init_fit')
plt.plot(theta, results.best_fit, 'r-', label='best_fit')
plt.plot(theta, val_bij, 'ro', label='traing_data')
plt.xlabel('j-i-k angles')
plt.ylabel('$b_{ij}$')
plt.legend(loc='best')
plt.savefig("3rd_try/fig_ClGeGe.png")
