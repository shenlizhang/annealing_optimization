# check and plot data using optimized parameters
import annealing_test
import numpy as np
import REBO
import Bbeta_least_square
import bond_order_calculation
import AQalpha_least_square
import matplotlib.pyplot as plt
import math
# Ge-Ge parameters
'''r1 = np.array([2.41, 2.44, 2.447, 2.665, 2.857, 2.939])
E1 = np.array([3.47, 2.06, 1.925, 1.18, 0.853, 0.565]) * (-1)
k1 = np.array([174.92, 155.97, 148.8, 89.9, 70.4, 69.8]) * 6.242 * 10**(-2)'''

'''r1 = np.array([1.941, 2.17, 2.24, 2.3, 2.351, 2.4,  2.635, 2.8])
E1 = np.array([8.650, 3.874, 3.13, 2.681, 2.315,  2.015, 1.048, 0.691]) * (-1)
k1 = np.array([580460, 270446, 222321.2, 180642.1, 161241.3,
               142697, 95834.4, 80232.4]) * 6.242 * 10**(-5)'''

r1 = np.array([2.22, 2.23, 2.27, 2.30, 2.31, 2.32])
E1 = np.array([4.23, 4.21, 1.47, 4.16, 3.28, 4.08]) * (-1)
k1 = np.array([163.41, 117.72, 109, 141.26, 122.60, 131.67]) * 6.242 * 10**(-2)

REBO_obj1 = REBO.REBO_fit(r1, E1, k1, 1, 2)
Bbeta_obj1 = Bbeta_least_square.Bbeta_least_square_minimization(REBO_obj1)
bond_obj1 = bond_order_calculation.Bond_order_calculation(
    REBO_obj1, Bbeta_obj1)
AQalpha_obj1 = AQalpha_least_square.AQalpha_least_square_minimization(
    REBO_obj1)
GeGe_optimizer = annealing_test.Potential_minimization(
    0.01, REBO_obj1, Bbeta_obj1, AQalpha_obj1, bond_obj1)

# 1. check the error function value is the same as during optimization process
#bond_obj1.normalization()
#bond_obj1.bond_order_calculation()
#e = GeGe_optimizer.energy()
#print e

# 2. plot the fitted training data curve
# Bbeta part
B1, B2, B3 = REBO_obj1.p['B1'], REBO_obj1.p['B2'], REBO_obj1.p['B3']
beta1, beta2, beta3 = REBO_obj1.p[
    'beta1'], REBO_obj1.p['beta2'], REBO_obj1.p['beta3']
A, alpha = REBO_obj1.p['A'], REBO_obj1.p['alpha']
x1 = r1[0:6]
'''Bbeta_obj1.set_training_values1()
y1 = np.asarray(Bbeta_obj1.training_value)
y2 = (Bbeta_obj1.exponential(x1, B1, beta1) + Bbeta_obj1.exponential(x1, B2, beta2)
      + Bbeta_obj1.exponential(x1, B3, beta3)) / (Bbeta_obj1.deriv_exponential(x1, B1, beta1)
                                                  + Bbeta_obj1.deriv_exponential(x1, B2, beta2) +
                                                  Bbeta_obj1.deriv_exponential(x1, B3, beta3))
print y2
plt.plot(x1, y1, 'bo', label='Bbeta training data')
plt.plot(x1, y2, 'b-', label='best_fit')
plt.legend()
plt.xlabel("Bond Length ($\AA$)")
plt.ylabel("Training function ($\AA$)")
plt.savefig("C_Bbeta_training_data_best_param.png")
plt.close()

# AQalpha part
A, Q, alpha = REBO_obj1.p['A'], REBO_obj1.p['Q'], REBO_obj1.p['alpha']
AQalpha_obj1.set_training_values2()
y3 = np.asarray(AQalpha_obj1.training_value2l)
y4 = AQalpha_obj1.exponential1(x1, A, Q, alpha) + AQalpha_obj1.exponential2(
    x1, A, Q, alpha) + AQalpha_obj1.exponential3(x1, A, Q, alpha)
print y4
plt.plot(x1, y3, 'go', label='AQalpha training data')
plt.plot(x1, y4, 'g-', label='best_fit')
plt.legend()
plt.xlabel("Bond Length ($\AA$)")
plt.ylabel("Training function (ev/$\AA^2$)")
plt.savefig("C_AQalpha_training_data_best_param.png")
plt.close()'''

# 3. plot bond energy and force constant versus bond length
'''e, k = [], []
#e = []
for i in range(0, r1.shape[0]):
    e.append(REBO_obj1.energy_calculation(r1[i], i))
    k.append(REBO_obj1.force_constant_calculation(r1[i], i))

e, k = np.asarray(e), np. asarray(k)
#e = np.asarray(e)
#print e, k/ 6.242 / 10**(-2)
#print np.dot((k1-k),(k1-k))
plt.plot(r1, E1, 'bo', label='training data')
plt.plot(r1, e, 'rs-.', label='fitted curve')
plt.xlim([2.2, 2.33])
plt.xlabel("Bond Length ($\AA$)")
plt.ylabel("Bond Energy (eV)")
plt.legend(loc='best')
plt.savefig("analysis_GeCl_6trainingdata_withK/bond energy_adjusted.png")
plt.close()'''
x = np.arange(6)
bar_width = 0.35
rnew = []
for i in range(0, r1.shape[0]):
    bij = REBO_obj1.p['bij'][i]
    rn = math.log(bij * B1 * beta1 / (A * alpha)) / (beta1 - alpha)
    rnew.append(rn)
print rnew
b1 = plt.plot(x, rnew, 'ro--', label='fitted value')
b2 = plt.plot(x, r1, 'bs--', label='DFT value')
plt.xlabel('Molecules')
plt.ylabel('bond length ($\AA$)')
plt.ylim([2.0, 2.7])
plt.xticks(x, ('$GeCl_4$', '$Ge_2Cl_6$',
                           '$GeCl_3$', '$GeCl_2$', '$GeCl$', '$Ge_4H_9Cl$'))
plt.legend(loc='best')
plt.savefig("analysis_GeCl_6trainingdata_withK/bond_length_adjusted.png")
'''plt.plot(r1, k1/6.242/10**(-2), 'bo', label='training data')
plt.plot(r1, k/6.242/10**(-2), 'rs-.', label='fitted curve')
plt.xlim([2.2, 2.33])
plt.xlabel("Bond Length ($\AA$)")
plt.ylabel("Force constant (N/$m^2$)")
#plt.ylabel("Force constant (eV/$\AA^2$)")
plt.legend(loc='best')
plt.savefig("analysis_GeCl_6trainingdata_withK/Force constant_adjusted.png")
plt.close()'''
