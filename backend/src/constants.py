from math import pow

DELTATIME = 0.005

SIGMA = 1
EPSILON = 2
SIGMA2 = pow(SIGMA, 2)
SIGMA6 = pow(SIGMA, 6)
SIGMA12 = pow(SIGMA6, 2)
RMIN = pow(2., (1. / 6)) * SIGMA