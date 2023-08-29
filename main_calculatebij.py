# main script to run annealing optimization
import annealing_test
import numpy as np
import REBO
import Bbeta_least_square
import bond_order_calculation
import AQalpha_least_square
import matplotlib.pyplot as plt
# Ge-Cl parameters
r1 = np.array([2.22, 2.23, 2.27, 2.30, 2.31, 2.32])
E1 = np.array([4.23, 4.21, 1.47, 4.16, 3.28, 4.08]) * (-1)
'''r1=np.array([2.50])
E1=np.array([2.91])*(-1)'''
k1 = np.array([163.41, 117.72, 109, 141.26, 122.60, 131.67]) * 6.242 * 10**(-2)


# initialize GeGe REBO object
#outfile = 'analysis/parameters_Ge_T5000_300_0.001_200steps_test_3.txt'
REBO_obj1 = REBO.REBO_fit(r1, E1, k1, 1, 2)
Bbeta_obj1 = Bbeta_least_square.Bbeta_least_square_minimization(REBO_obj1)
bond_obj1 = bond_order_calculation.Bond_order_calculation(
    REBO_obj1, Bbeta_obj1)
#bond_obj1.normalization()
bond_obj1.bond_order_calculation()
