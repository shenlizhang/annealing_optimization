# annealing algorithm test
from numpy import random
import copy
from simanneal import Annealer
import REBO2_whole
import json
import time
import datetime
import math


def round_figures(x, n):
    """Returns x rounded to n significant figures."""
    return round(x, int(n - math.ceil(math.log10(abs(x)))))


def time_string(seconds):
    """Returns time in seconds as a string formatted HHHH:MM:SS."""
    s = int(round(seconds))  # round to nearest second
    h, s = divmod(s, 3600)   # get hours and remainder
    m, s = divmod(s, 60)     # split remainder into minutes and seconds
    return '%4i:%02i:%02i' % (h, m, s)


class Potential_minimization(Annealer):

    def __init__(self, delta, objREBO, objBbeta, objAQalpha, objbond):
        """, but only A, Q and alpha will be changed
        delta is the percentage to change the picked variable value
        """
        self.state = objREBO.p
        self.delta = delta
        self.obj1 = objREBO
        self.obj2 = objBbeta
        self.obj3 = objAQalpha
        self.obj4 = objbond
        # used parameters are A, Q, alpha
        super(Potential_minimization, self).__init__(self.state)
        return

    def move(self):
        """modify a randomly chosen parameter in A,Q,alpha by plus/minus delta percent."""
        a = random.randint(0, 1)
        parlist = ['A', 'alpha']
        b = random.choice([-1, 1], 1)  # randomly choose if add/minus delta
        self.state[parlist[a]] = float(
            self.state[parlist[a]] * (1 + b * self.delta))
        return

    def energy(self):
        """calculate the least-squares sum of energy and force constant"""
        e_diff, f_diff = 0, 0
        # for i in range(0, len(self.obj1.r)):
        for i in range(0, 6):
            e_diff += (self.obj1.energy_calculation(
                self.obj1.r[i], i) - self.obj1.E[i])**2
            f_diff += (self.obj1.force_calculation(
                self.obj1.r[i], i))**2
        # print "e diff and k diff: \n", e_diff * 81 + k_diff
        #f = open("analysis/energy_GeGe_T5000_300_0.001_200steps_3.txt", 'a')
        #f.write("\n energy \n")
        #f.write("%.5e \n" % (e_diff * 81 + k_diff))
        # f.close()
        return e_diff * 80 + f_diff

    def local_optimization(self):
        self.obj2.set_training_values1()
        self.obj2.Bbeta_ls_min()
        self.obj4.normalization()
        self.obj4.bond_order_calculation()
        # print self.obj1.p['B1'], self.obj1.p['bij']
        self.obj3.set_training_values2()
        self.obj3.AQalpha_ls_min()
        self.obj2.set_training_values1()
        self.obj2.Bbeta_ls_min()
        self.obj4.normalization()
        self.obj4.bond_order_calculation()
        e1 = self.energy()
        f = open("analysis_GeCl_6trainingdata/local_energy_T30000_3_0.05_600steps.txt", 'a')
        f.write("%.5e \n" % e1)
        f.close()
        '''f = open(self.obj1.outfile, 'a')
        f.write("\n optimized state parameters \n")
        f.write(json.dumps(self.state))
        f.close()'''

    def anneal(self, objREBO):
        """Minimizes the energy of a system by simulated annealing.

        Parameters
        state : an initial arrangement of the system

        Returns
        (state, energy): the best state and energy found.
        """
        step = 0
        self.start = time.time()

        # Precompute factor for exponential cooling from Tmax to Tmin
        if self.Tmin <= 0.0:
            raise Exception('Exponential cooling requires a minimum "\
                "temperature greater than zero.')
        Tfactor = -math.log(self.Tmax / self.Tmin)

        # Note initial state
        self.local_optimization()
        T = self.Tmax
        E = self.energy()
        prevState = self.copy_state(self.state)
        prevEnergy = E
        self.best_state = self.copy_state(self.state)
        self.best_energy = E
        trials, accepts, improves = 0, 0, 0
        if self.updates > 0:
            updateWavelength = self.steps / self.updates
            self.update(step, T, E, None, None)

        # Attempt moves to new states
        while step < self.steps and not self.user_exit:
            step += 1
            T = self.Tmax * math.exp(Tfactor * step / self.steps)
            self.move()
            E1 = self.energy()
            f = open(
                "analysis_GeCl_6trainingdata/random_energy_T30000_3_0.05_600steps_AQalpha.txt", 'a')
            f.write("%.5e \n" % E1)
            f.close()
            self.local_optimization()  # change parameter dicitonary : p

            E = self.energy()
            dE = E - prevEnergy
            trials += 1

            if dE > 0.0 and math.exp(-dE / T) < random.random():
                # Restore previous state
                self.state = self.copy_state(prevState)
                E = prevEnergy
                # print "\n first branch %.4e \n" %E
                objREBO.p = self.state
                self.obj1.p = self.state
                self.obj2.obj.p = self.state
                self.obj3.obj.p = self.state
                self.obj4.obj.p = self.state

            else:
                # Accept new state and compare to best state
                # print "\n 2nd branch %.4e \n" %E
                accepts += 1
                if dE < 0.0:
                    improves += 1
                prevState = self.copy_state(self.state)
                prevEnergy = E
                if E < self.best_energy:
                    self.best_state = self.copy_state(self.state)
                    self.best_energy = E

            if self.updates > 1:
                if (step // updateWavelength) > ((step - 1) // updateWavelength):
                    self.update(
                        step, T, E, accepts / trials, improves / trials)
                    trials, accepts, improves = 0, 0, 0
            f = open(
                "analysis_GeCl_6trainingdata/accept_energy_T30000_3_0.05_600steps_AQalpha.txt", 'a')
            f.write("%.5e \n" % E)
            f.close()

        self.state = self.copy_state(self.best_state)
        objREBO.p = self.state
        self.obj1.p = objREBO.p
        self.obj2.obj.p = objREBO.p
        self.obj3.obj.p = objREBO.p
        self.obj4.obj.p = objREBO.p
        if self.save_state_on_exit:
            self.save_state()

        # Return best state and energy
        return self.best_state, self.best_energy
