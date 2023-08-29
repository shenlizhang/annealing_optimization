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
k1 = np.array([163.41, 117.72, 109, 141.26, 122.60, 131.67]) * 6.242 * 10**(-2)


# initialize GeGe REBO object
#outfile = 'analysis/parameters_Ge_T5000_300_0.001_200steps_test_3.txt'
REBO_obj1 = REBO.REBO_fit(r1, E1, k1, 1, 2)
Bbeta_obj1 = Bbeta_least_square.Bbeta_least_square_minimization(REBO_obj1)
bond_obj1 = bond_order_calculation.Bond_order_calculation(
    REBO_obj1, Bbeta_obj1)
AQalpha_obj1 = AQalpha_least_square.AQalpha_least_square_minimization(
    REBO_obj1)
# Bbeta_obj1.set_training_values1()
# Bbeta_obj1.Bbeta_ls_min()
# bond_obj1.normalization()
# bond_obj1.bond_order_calculation()
# initialize annealing optimization
'''print "B1 before optimization:", REBO_obj1.p['B1']
Bbeta_obj1.set_training_values1()
print "B1 after set value:", REBO_obj1.p['B1']
Bbeta_obj1.Bbeta_ls_min()
print "B1 after optimization:", REBO_obj1.p['B1']'''
GeCl_optimizer = annealing_test.Potential_minimization(
    0.05, REBO_obj1, Bbeta_obj1, AQalpha_obj1, bond_obj1)

# GeGe_optimizer.local_optimization()
'''VA = []
VR = []
for i in range(0, r1.shape[0]):
    va = REBO_obj1.VA(r1[i])
    vr = REBO_obj1.VR(r1[i])
    d_va = REBO_obj1.fA(r1[i])
    d_vr = REBO_obj1.fR(r1[i])
    VA.append(va / d_va)
    VR.append((vr - E1[i]) / d_vr)
    print(vr - E1[i]) / va, d_vr / d_va
print "before_energy VA: ", VA
print "before_energy VR: ", VR'''

'''df_value=[]
f_value=[]
df2_value=[]
rtest=np.linspace(2.1,3.5,100)
for item in rtest:
	f_value.append(REBO_obj1.fij(item))
	df2_value.append(REBO_obj1.f_2nd_derivative(item))
	df_value.append(REBO_obj1.f_derivative(item))
plt.plot(rtest, np.asarray(f_value),'ko',label='f')
plt.plot(rtest, np.asarray(df_value),'rs',label='df')
plt.plot(rtest, np.asarray(df2_value),'b*',label='df2')
plt.legend()
plt.show()'''
GeCl_optimizer.steps = 600
GeCl_optimizer.Tmax = 30000.0
GeCl_optimizer.Tmin = 50
print "initial state:", GeCl_optimizer.state
#GeGe_optimizer.energy()
state, e = GeCl_optimizer.anneal(REBO_obj1)
print "final state:", GeCl_optimizer.state,e
