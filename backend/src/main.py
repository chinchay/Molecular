import numpy as np
import matplotlib.pyplot as plt

from atom import Atom
from system import System
from dynamics import updateLennarJonesAcc, applyThermalBath, computeMechanicalEnergy, computeDistance, computeKinetic
# from . import constants
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

n = 1201
fig, ax = plt.subplots()
ax.set_xlim(-1, 4), ax.set_xticks([])
scatter1 = ax.scatter([], [], color="b")
scatter2 = ax.scatter([], [], color="r")

nPeriod = 200

isInsideBall = False
for i in range(n):
    atom1.move()
    atom2.move()

    [x1, y1, z1] = atom1.pos
    [x2, y2, z2] = atom2.pos

    x_values = [x1] 
    y_values = [y1] 
    scatter1.set_offsets(list(zip(x_values, y_values)))


    x_values2 = [x2] 
    y_values2 = [y2] 
    scatter2.set_offsets(list(zip(x_values2, y_values2)))

    plt.pause(0.01)

    updateLennarJonesAcc(atom1, atom2)

    _, distance = computeDistance(atom1, atom2)

    energy = computeMechanicalEnergy(system)

    if abs(distance - RMIN) < 0.01:
        if not isInsideBall:
            isInsideBall = True

            kinetic = computeKinetic(system)
            print("Applying thermal bath. Kinetic energy = ", kinetic)

            applyThermalBath(system, desiredKinetic)           
    else:
        isInsideBall = False

plt.show()

