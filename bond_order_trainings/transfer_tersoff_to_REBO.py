# Shenli, 04/04/2018
from math import pi, cos
from lmfit import Model
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


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
    f = gamma * (1 + c**2 / d**2 - c**2 / (d**2 + (np.cos(theta) - h)**2))
    return f


def ksi(theta, rik, rij, gamma, c, d, h, lamda3, m, R1, R2):
    # this is just for one molecule; rij is one number; rik could be an array
    s = []
    for i in range(0, rik.shape[0]):
        s1 = np.dot(fik(rik[i], R1, R2), np.exp(
            np.power(lamda3 * (rij[i] - rik[i]), m))) * gfunc(theta[i] / 180 * pi, gamma, c, d, h)
        s.append(s1)
    print ksi
    return np.asarray(s)


def bij1(beta, n, theta, rik, rij, gamma, c, d, h, lamda3, m, R1, R2):
    bij = (1 + (beta * ksi(theta, rik, rij, gamma, c,
                           d, h, lamda3, m, R1, R2))**n)**(1. / (-2 * n))
    # print "bij", bij
    # print "beta, gamma, c, d, h, lamda3", beta, gamma, c, d, h, lamda3
    print ksi(theta, rik, rij, gamma, c,
              d, h, lamda3, m, R1, R2)**n, beta**n
    return bij


def gfunc2(theta, a1, a2, a3, a4, a5, a6):
    x = np.cos(theta)
    f = a1 * np.power(x, 5) + a2 * np.power(x, 4) + \
        a3 * np.power(x, 3) + a4 * np.power(x, 2) + a5 * x + a6
    # f = a1 * np.power(x, 2) + a2 * np.power(x, 1) + a3
    return f


def ksi2(theta, rik, rij, lamda, a1, a2, a3, a4, a5, a6, R1, R2):
    # this is just for one molecule; rij is one number; rik could be an array
    s = []
    for i in range(0, rik.shape[0]):
        # print fik(rik[i], R1, R2)
        # print gfunc(theta[i] / 180 * pi, a1, a2, a3)
        s1 = np.dot(fik(rik[i], R1, R2), np.exp(lamda) * np.ones(len(rik[i])) *
                    gfunc2(theta[i] * pi / 180, a1, a2, a3, a4, a5, a6))
        s.append(s1)
    return np.asarray(s)


def bij2(theta, rik, rij, lamda, a1, a2, a3, a4, a5, a6, R1, R2):
    bij = (1 + ksi2(theta, rik, rij, lamda, a1, a2,
                    a3, a4, a5, a6, R1, R2))**(-1. / 2)
    # print "bij", bij
    # print "lamda, a1, a2, a3", lamda, a1, a2, a3
    return bij

# calculate Ge-Ge-Ge bij values with provided Tersoff parameters
'''theta = np.linspace(0, 180, 100) * pi / 180
# theta1 = np.linspace(0, 109.47, 20) * pi / 180
# theta2 = np.linspace(109.47, 120, 20) * pi / 180
# theta3 = np.linspace(120, 180, 20) * pi / 180
# gtheta = gfunc(theta, 1, 1.0643e5, 15.652, -0.43884)
# gtheta1 = gfunc(theta1, 1, 1.0643e5, 15.652, -0.43884)
# gtheta2 = gfunc(theta2, 1, 1.0643e5, 15.652, -0.43884)
# gtheta3 = gfunc(theta3, 1, 1.0643e5, 15.652, -0.43884)
gtheta = gfunc(theta * pi / 180, 4, 0, 1, 1)
mod = Model(gfunc2, independent_vars=['theta'])
results = mod.fit(gtheta, a1=0, a2=0, a3=0, a4=0,
                  a5=1, a6=1, a7=1, theta=theta)
results1 = mod.fit(gtheta1, a1=1, a2=1, a3=1, a4=1,
                   a5=1, a6=1, a7=1, theta=theta1)
results2 = mod.fit(gtheta2, a1=1, a2=1, a3=1, a4=1,
                   a5=1, a6=1, a7=1, theta=theta2)
results3 = mod.fit(gtheta3, a1=1, a2=1, a3=1, a4=1,
                   a5=1, a6=1, a7=1, theta=theta3)
print results.best_values  # , results2.best_values, results3.best_values
plt.plot(theta * 180 / pi, gtheta, '-', linewidth='2', label="Tersoff")'''
'''plt.plot(theta1 * 180 / pi, results1.best_fit, 'gs-', label="fitted REBO")
plt.plot(theta2 * 180 / pi, results2.best_fit, 'gs-')
plt.plot(theta3 * 180 / pi, results3.best_fit, 'gs-')
plt.plot(theta * 180 / pi, results.best_fit, 'gs-', label="fitted REBO")
plt.xlim([0, 180])
plt.xlabel("j-i-k angle (degree)")
plt.ylabel(r'$g(\theta)$')
plt.title(r'Ge-Ge-Ge G($\theta$) function')
plt.legend(loc='best')
plt.savefig("gtheta_Cl_fitting.png")'''
theta = np.array([180, 179, 120, 109.5, 90, 70.53, 60])
rik = np.array([[2.665], [2.665], [2.44, 2.44], [2.447, 2.447, 2.447], [2.665, 2.665,
                                                                        2.665, 2.665], [2.857, 2.857], [2.939, 2, 939]])
rij = np.array([2.665, 2.665, 2.44, 2.447, 2.665, 2.857, 2.939])
val_bij = bij1(9.0166e-7, 0.75627, theta, rik, rij, 1,
               1.0643e5, 15.652, -0.43884, 0.0, 3.0, 2.8, 3.1)
# print val_bij
val2_bij = val_bij ** (-2) - 1
# print val2_bij
numbers = np.dot(np.array([1, 1, 2, 3, 4, 2, 2]), fik(rij, 2.8, 3.1))
gtheta = np.divide(val2_bij, numbers)
# print gtheta
fc1 = interp1d(theta * pi / 180, gtheta, kind='cubic')
new_theta = np.linspace(60, 180, num=121, endpoint=True)
y = fc1(new_theta * pi / 180)

theta1=new_theta[60:121]
theta2=new_theta[50:61]
#theta1 = theta[0:3]
#theta2 = theta[2:4]
theta3 = theta[3:7]
#gtheta1 = gtheta[0:3]
#gtheta2 = gtheta[2:4]
gtheta1=y[60:121]
gtheta2=y[50:61]
gtheta3 = gtheta[3:7]
mod = Model(gfunc2, independent_vars=['theta'])
results1 = mod.fit(gtheta1, a1=0, a2=0, a3=0, a4=3,
                   a5=1, a6=1, theta=theta1 * pi / 180, method="dogbox")
results2 = mod.fit(gtheta2, a1=0, a2=0, a3=0, a4=3,
                   a5=1, a6=1, theta=theta2 * pi / 180, method="dogbox")
results3 = mod.fit(gtheta3, a1=0, a2=0, a3=0, a4=3,
                   a5=1, a6=1, theta=theta3 * pi / 180, method="dogbox")
print results1.best_values, results2.best_values, results3.best_values
'''y1 = mod.eval(results1.best_values, theta=np.linspace(120, 180, 50) * pi / 180)
y2 = mod.eval(results2.best_values, theta=np.linspace(
    109.5, 120, 50) * pi / 180)
y3 = mod.eval(results3.best_values,
              theta=np.linspace(60, 109.5, 50) * pi / 180)'''

#plt.plot(new_theta, y, 'b-', label="spline fit")
plt.plot(theta, gtheta, linestyle='none', marker='s', linewidth='2',
         label=r"target g($\theta$) in REBO")
plt.plot(theta1, results1.best_fit, 'b-', label="fitted")
plt.plot(theta2, results2.best_fit, 'b-')
plt.plot(theta3, results3.best_fit, 'b-')
plt.xlim([50, 180])
plt.xlabel("j-i-k angle (degree)")
plt.ylabel(r'$g(\theta)$')
plt.title(r'j-Ge-k G($\theta$) function')
plt.legend(loc='best')
plt.savefig("REBO_gtheta_3.png")
