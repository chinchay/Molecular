import numpy as np
import matplotlib.pyplot as plt

from atom import Atom
from dynamics import updateLennarJonesAcc
# from . import constants
from constants import RMIN

atom1 = Atom()
atom1.pos = np.asarray([0, 0, 0])
atom1.vel = np.asarray([0, 0, 0])
atom1.acc = np.asarray([0, 0, 0])
atom1.mass = 100


atom2 = Atom()
delta = 0.1
x = RMIN + delta
atom2.pos = np.asarray([x, 0, 0])
atom2.vel = np.asarray([0, 0, 0])
atom2.acc = np.asarray([0, 0, 0])
atom2.mass = 0.5

n = 1200
fig, ax = plt.subplots()
ax.set_xlim(-1, 4), ax.set_xticks([])
scatter1 = ax.scatter([], [], color="b")
scatter2 = ax.scatter([], [], color="r")

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

plt.show()

