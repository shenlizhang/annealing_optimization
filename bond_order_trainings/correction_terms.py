# python script to calculate and then do a bicubic interpolation for
# correction terms Hij

import numpy as np
#from scipy.interpolate import interp1d
from scipy import interpolate
import matplotlib.pyplot as plt

bij_old = np.array([1, 1, 0.9979, 0.9989, 0.99865, 1, 1])
bij_new = np.array([0.8837, 0.8837, 1.0187, 0.6154, 1.0233, 1, 1])

h = bij_new**(-2) - bij_old
print h
h_new = np.ma.masked_where(h, h)
print h_new

x = np.array([0, 0.1, 1, 2, 3, 3.9, 4])
#fc1 = interp1d(x, h, kind='cubic')
fc1 = interpolate.splrep(x, h, k=3, s=0)
print fc1

xnew = np.linspace(0, 4, num=100, endpoint=True)
ynew = interpolate.splev(xnew, fc1, der=0)
#plt.plot(x, h, 'o', xnew, fc1(xnew), '-')
plt.plot(x, h_new, 'o', xnew, ynew, '-')
plt.legend(['data', 'cubic'], loc='best')
plt.xlabel('$N^{Cl}_{ij}$')
plt.ylabel('H value')
plt.title("$H(N^{Ge}_{ij}=0$) cubic interpolation ")
plt.savefig('withK_polynominal/h_cubic_interp_2.png')
