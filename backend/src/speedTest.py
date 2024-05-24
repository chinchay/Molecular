from cProfile import Profile
from pstats import SortKey, Stats

import numpy as np
from atom import Atom
from system import System
from dynamics import updateLennarJonesAcc, applyThermalBath, computeMechanicalEnergy, computeDistance, computeKinetic
from constants import RMIN

atom1 = Atom()
atom1.pos = np.asarray([0, 0, 0], dtype=float)
atom1.vel = np.asarray([0, 0, 0], dtype=float)
atom1.acc = np.asarray([0, 0, 0], dtype=float)
atom1.mass = 100


atom2 = Atom()
delta = 0.1
x = RMIN + delta
atom2.pos = np.asarray([x, 0, 0], dtype=float)
atom2.vel = np.asarray([0, 0, 0], dtype=float)
atom2.acc = np.asarray([0, 0, 0], dtype=float)
atom2.mass = 0.5

system = System([atom1, atom2])

initialEnergy = computeMechanicalEnergy(system)
desiredKinetic = initialEnergy / 2

n = 201

with Profile() as profile:

    isInsideBall = False
    for i in range(n):
        atom1.move()
        atom2.move()

        updateLennarJonesAcc(atom1, atom2)

        _, distance = computeDistance(atom1, atom2)

        if abs(distance - RMIN) < 0.01:
            if not isInsideBall:
                isInsideBall = True

                # kinetic = computeKinetic(system)
                # print("Applying thermal bath. Kinetic energy = ", kinetic)

                applyThermalBath(system, desiredKinetic)           
        else:
            isInsideBall = False
    #
    (
        Stats(profile).strip_dirs().sort_stats(SortKey.CUMULATIVE).print_stats()
    )