# Shenli, 04/05/2018
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


def gfunc2(theta, a1, a2, a3, a4, a5, a6):
    x = np.cos(theta)
    f = a1 * np.power(x, 5) + a2 * np.power(x, 4) + \
        a3 * np.power(x, 3) + a4 * np.power(x, 2) + a5 * x + a6
    return f


def ksi2(theta, rik, rij, lamda, a1, a2, a3, a4, a5, a6, R1, R2):
    # this is just for one molecule; rij is one number; rik could be an array
    s = []
    for i in range(0, rik.shape[0]):
        s1 = np.dot(fik(rik[i], R1[i], R2[i]), np.exp(lamda) * np.ones(len(rik[i])) *
                    gfunc2(theta[i] * pi / 180, a1, a2, a3, a4, a5, a6))
        s.append(s1)
    return np.asarray(s)


def bij2(theta, rik, rij, lamda, a1, a2, a3, a4, a5, a6, R1, R2):
    bij = (1 + ksi2(theta, rik, rij, lamda, a1, a2,
                    a3, a4, a5, a6, R1, R2))**(-1. / 2)
    return bij

par = np.array([0.025095542734161067, -0.020490382789339085, 0.0080300114984608456,
                0.16657760597198434, 0.13141982857766421, 0.02778701715262941])
theta = np.array([101.09, 108.69, 109.47, 109.78, 104.21])
rij = np.array([2.3, 2.27, 2.22, 2.50, 2.32])
rik = np.array([[2.3], [2.27, 2.27], [2.22, 2.22, 2.22],
                [2.23, 2.23, 2.23], [2.49, 2.49, 2.49]])
R1 = [2.284, 2.284, 2.284, 2.284, 2.8]
R2 = [2.584, 2.584, 2.584, 2.584, 3.1]
val_bij = bij2(theta, rik, rij, 0, par[0], par[1], par[
              2], par[3], par[4], par[5], R1, R2)
print val_bij
