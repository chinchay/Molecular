import numpy as np
import matplotlib.pyplot as plt

from atom import Atom
from system import System
from dynamics import applyThermalBath, computeMechanicalEnergy, computeDistance, computeKinetic, updateAcc
# from . import constants
from constants import RMIN

atom1 = Atom()
atom1.pos = np.asarray([0, 0, 0], dtype=float)
atom1.vel = np.asarray([0, 0, 0], dtype=float)
atom1.acc = np.asarray([0, 0, 0], dtype=float)
atom1.mass = 100


atom2 = Atom()
delta = 0.1
d = RMIN + delta
atom2.pos = np.asarray([d, 0, 0], dtype=float)
atom2.vel = np.asarray([0, 0, 0], dtype=float)
atom2.acc = np.asarray([0, 0, 0], dtype=float)
atom2.mass = 0.5

atom3 = Atom()
delta = 0.1
d = RMIN + delta
atom3.pos = np.asarray([0, d, 0], dtype=float)
atom3.vel = np.asarray([0, 0, 0], dtype=float)
atom3.acc = np.asarray([0, 0, 0], dtype=float)
atom3.mass = 0.5

atom4 = Atom()
delta = 0.1
d = RMIN + delta
atom4.pos = np.asarray([-d, 0, 0], dtype=float)
atom4.vel = np.asarray([0, 0, 0], dtype=float)
atom4.acc = np.asarray([0, 0, 0], dtype=float)
atom4.mass = 0.5


system = System([atom1, atom2, atom3, atom4])

initialEnergy = computeMechanicalEnergy(system)
desiredKinetic = initialEnergy / 2

n = 3201
fig, ax = plt.subplots()
ax.set_xlim(-4, 4), ax.set_xticks([])
ax.set_ylim(-4, 4), ax.set_xticks([])
scatter1 = ax.scatter([], [], color="b")
scatter2 = ax.scatter([], [], color="r")
scatter3 = ax.scatter([], [], color="g")
scatter4 = ax.scatter([], [], color="k")

# nPeriod = 200

isInsideBall = False
for i in range(n):
    system.move()

    [x1, y1, z1] = atom1.pos
    [x2, y2, z2] = atom2.pos
    [x3, y3, z3] = atom3.pos
    [x4, y4, z4] = atom4.pos

    x_values = [x1] 
    y_values = [y1] 
    scatter1.set_offsets(list(zip(x_values, y_values)))

    x_values2 = [x2] 
    y_values2 = [y2] 
    scatter2.set_offsets(list(zip(x_values2, y_values2)))

    x_values3 = [x3] 
    y_values3 = [y3] 
    scatter3.set_offsets(list(zip(x_values3, y_values3)))

    x_values4 = [x4] 
    y_values4 = [y4] 
    scatter4.set_offsets(list(zip(x_values4, y_values4)))

    plt.pause(0.01)

    # updateLennarJonesAcc(atom1, atom2)
    updateAcc(system)

    _, distance = computeDistance(atom1, atom2) #<<< include other atoms!

    # energy = computeMechanicalEnergy(system)

    if abs(distance - RMIN) < 0.01:
        if not isInsideBall:
            isInsideBall = True

            kinetic = computeKinetic(system)
            print("Applying thermal bath. Kinetic energy = ", kinetic)

            applyThermalBath(system, desiredKinetic)           
    else:
        isInsideBall = False

plt.show()

