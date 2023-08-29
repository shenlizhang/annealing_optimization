# main scrpit for optimization

import REBO
import Bbeta_least_square
import bond_order_calculation
import AQalpha_least_square
import numpy as np


class REBO2_optimization(REBO.REBO_fit,
                         Bbeta_least_square.Bbeta_least_square_minimization,
                         bond_order_calculation.Bond_order_calculation,
                         AQalpha_least_square.AQalpha_least_square_minimization):

    def __init__(self, r, E, k, type1, type2, type3=''):
        #REBO.REBO_fit.__init__(self, r, E, k, type1, type2, type3)
        self.REBO_obj1 = REBO.REBO_fit(r, E, k, type1, type2, type3)
       # Bbeta_least_square.Bbeta_least_square_minimization.__init__(
        #    self, self.REBO_obj1)
        self.Bbeta_obj1 = Bbeta_least_square.Bbeta_least_square_minimization(
            self.REBO_obj1)
        bond_order_calculation.Bond_order_calculation.__init__(
            self, self.REBO_obj1, self.Bbeta_obj1)
        AQalpha_least_square.AQalpha_least_square_minimization.__init__(
            self, self.REBO_obj1)
        #AQalpha_least_square.AQalpha_least_square(REBO_obj1)
        return
