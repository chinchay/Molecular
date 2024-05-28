from cProfile import Profile
from pstats import SortKey, Stats

import numpy as np
from atoms import Atoms
from dynamics import applyThermalBath, computeMechanicalEnergy, computeDistance, computeKinetic, updateAcc, calculate
from constants import RMIN

x = 6
y = 6
z = 6
x2 = 12
y2 = 12
nAtoms = 30
delta = 0.1
d = RMIN + delta
atoms = Atoms(nAtoms)
atoms.positions = np.asarray(
    [
        [ 0,  0, 0],
        [ d,  0, 0],
        [ 0,  d, 0],
        [-d,  0, 0],
        [ 0, -d, 0],

        [ 0,  y, 0],
        [ d,  y, 0],
        [ 0,  y + d, 0],
        [-d,  y, 0],
        [ 0,  y - d, 0],

        [ x,  0, 0],
        [ x + d,  0, 0],
        [ x,  d, 0],
        [ x - d,  0, 0],
        [ x, -d, 0],
        ########################

        [ 0,  0, z],
        [ d,  0, z],
        [ 0,  d, z],
        [-d,  0, z],
        [ 0, -d, z],

        [ 0,  y2, z],
        [ d,  y2, z],
        [ 0,  y2 + d, z],
        [-d,  y2, z],
        [ 0,  y2 - d, z],

        [ x2,  0, z],
        [ x2 + d,  0, z],
        [ x2,  d, z],
        [ x2 - d,  0, z],
        [ x2, -d, z],
    ],
    dtype=float
)
nAtoms, _ = atoms.positions.shape
atoms.velocities = np.zeros((nAtoms, 3), dtype=float)
atoms.accelerations = np.zeros((nAtoms, 3), dtype=float)

for i in range(nAtoms):
    if i in [0, 5, 15, 20, 25]:
        atoms.masses[i] = 100
    else:
        atoms.masses[i] = 0.5
#


calculate(atoms)
kineticEnergy = computeKinetic(atoms)
potentialEnergy = np.sum(atoms.potentialEnergies)
initialEnergy = kineticEnergy + potentialEnergy
desiredKinetic = initialEnergy / 2

n = 20000 #200000

with Profile() as profile:

    isInsideBall = False
    for i in range(n):
        atoms.move()

        # updateLennarJonesAcc(atom1, atom2)
        updateAcc(atoms)

        distance = np.linalg.norm(atoms.positions[0] - atoms.positions[1])

        if abs(distance - RMIN) < 0.01:
            if not isInsideBall:
                isInsideBall = True

                # kinetic = computeKinetic(system)
                # print("Applying thermal bath. Kinetic energy = ", kinetic)

                applyThermalBath(atoms, desiredKinetic)  
        else:
            isInsideBall = False
    #
    (
        Stats(profile).strip_dirs().sort_stats(SortKey.CUMULATIVE).print_stats()
    )