import numpy as np
import matplotlib.pyplot as plt

from atoms import Atoms
from dynamics import applyThermalBath, computeMechanicalEnergy, computeDistance, computeKinetic, updateAcc, calculate
# from . import constants
from constants import RMIN

delta = 0.1
d = RMIN + delta
atoms = Atoms(5)
atoms.positions = np.asarray(
    [
        [ 0,  0, 0],
        [ d,  0, 0],
        [ 0,  d, 0],
        [-d,  0, 0],
        [ 0, -d, 0],
    ],
    dtype=float
)
nAtoms, _ = atoms.positions.shape
atoms.velocities = np.zeros((nAtoms, 3), dtype=float)
atoms.accelerations = np.zeros((nAtoms, 3), dtype=float)
atoms.masses[0, 0] = 100
atoms.masses[:, 0] = 0.5

calculate(atoms)
kineticEnergy = computeKinetic(atoms)
potentialEnergy = np.sum(atoms.potentialEnergies)
initialEnergy = kineticEnergy + potentialEnergy
desiredKinetic = initialEnergy / 2

n = 320#1
fig, ax = plt.subplots()
ax.set_xlim(-4, 4), ax.set_xticks([])
ax.set_ylim(-4, 4), ax.set_xticks([])
scatter1 = ax.scatter([], [], color="b")
scatter2 = ax.scatter([], [], color="r")
scatter3 = ax.scatter([], [], color="g")
scatter4 = ax.scatter([], [], color="k")
scatter5 = ax.scatter([], [], color="k")

# nPeriod = 200

isInsideBall = False
for i in range(n):
    atoms.move()

    [x1, y1, z1] = atoms.positions[0]
    [x2, y2, z2] = atoms.positions[1]
    [x3, y3, z3] = atoms.positions[2]
    [x4, y4, z4] = atoms.positions[3]
    [x5, y5, z5] = atoms.positions[4]

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

    x_values5 = [x5] 
    y_values5 = [y5] 
    scatter5.set_offsets(list(zip(x_values5, y_values5)))

    plt.pause(0.01)

    # updateLennarJonesAcc(atom1, atom2)
    updateAcc(atoms)

    distance = np.linalg.norm(atoms.positions[0] - atoms.positions[1])

    # energy = computeMechanicalEnergy(system)

    if abs(distance - RMIN) < 0.01:
        if not isInsideBall:
            isInsideBall = True

            kinetic = computeKinetic(atoms)
            print("Applying thermal bath. Kinetic energy = ", kinetic)

            applyThermalBath(atoms, desiredKinetic)           
    else:
        isInsideBall = False

plt.show()


# https://github.com/libAtoms/pymatnest/blob/master/ns_run.py
            # from ase.constraints import FixAtoms
            # c = FixAtoms(indices=[atom.index for atom in at if atom.index < movement_args['keep_atoms_fixed']])
            # at.set_constraint(c)
            # # ...
            # pot.lmp.command('fix 1 mobile nve')
            # pot.lmp.command('fix 1 all ns/gmc {} {} '.format(rng.int_uniform(1, 900000000), Emax))
            # pot.restart_lammps(at)



# lennearJones calculator: https://gitlab.com/ase/ase/-/blob/master/ase/calculators/lj.py?ref_type=heads


        # for ii in range(natoms):
        #     neighbors, offsets = self.nl.get_neighbors(ii)
        #     cells = np.dot(offsets, cell)

        #     # pointing *towards* neighbours
        #     distance_vectors = positions[neighbors] + cells - positions[ii]

        #     r2 = (distance_vectors ** 2).sum(1)
        #     c6 = (sigma ** 2 / r2) ** 3
        #     c6[r2 > rc ** 2] = 0.0
        #     c12 = c6 ** 2

        #     if smooth:
        #         cutoff_fn = cutoff_function(r2, rc**2, ro**2)
        #         d_cutoff_fn = d_cutoff_function(r2, rc**2, ro**2)

        #     pairwise_energies = 4 * epsilon * (c12 - c6)
        #     pairwise_forces = -24 * epsilon * (2 * c12 - c6) / r2  # du_ij





# from ase.build import molecule
# from ase.calculators.lj import LennardJones # https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html
# from ase import units
# from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
# from ase.md.verlet import VelocityVerlet

# atoms = molecule("CH")
# calc = LennardJones()
# atoms.calc = calc
# # pot = atoms.get_potential_energy()


# # Set the momenta corresponding to T=300K
# MaxwellBoltzmannDistribution(atoms, temperature_K=300)

# dyn = VelocityVerlet(atoms, 5 * units.fs)

# def printenergy(a):
#     """Function to print the potential, kinetic and total energy"""
#     epot = a.get_potential_energy() / len(a)
#     ekin = a.get_kinetic_energy() / len(a)
#     print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
#           'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))

# # Now run the dynamics
# printenergy(atoms)
# for i in range(20):
#     dyn.run(10)
#     printenergy(atoms)




# from ase import units
# from ase.lattice.cubic import FaceCenteredCubic
# from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
# from ase.md.verlet import VelocityVerlet

# from ase.build import molecule
# from ase.io.trajectory import Trajectory

# # Use Asap for a huge performance increase if it is installed
# use_asap = False

# if use_asap:
#     from asap3 import EMT
#     size = 10
# else:
#     from ase.calculators.emt import EMT
#     size = 3

# # Set up a crystal
# # atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
# #                           symbol='Cu',
# #                           size=(size, size, size),
# #                           pbc=True)

# atoms = molecule("CH")

# # Describe the interatomic interactions with the Effective Medium Theory
# atoms.calc = EMT()

# # Set the momenta corresponding to T=300K
# MaxwellBoltzmannDistribution(atoms, temperature_K=1)

# # We want to run MD with constant energy using the VelocityVerlet algorithm.
# dyn = VelocityVerlet(atoms, 5 * units.fs)  # 5 fs time step.


# def printenergy(a):
#     """Function to print the potential, kinetic and total energy"""
#     epot = a.get_potential_energy() / len(a)
#     ekin = a.get_kinetic_energy() / len(a)
#     print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
#           'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))


# # Now run the dynamics
# printenergy(atoms)
# for i in range(20):
#     dyn.run(10)
#     printenergy(atoms)



# # We also want to save the positions of all atoms after every 100th time step.
# traj = Trajectory('moldyn3.traj', 'w', atoms)
# dyn.attach(traj.write, interval=50)

# # Now run the dynamics
# dyn.run(5000)

# # https://wiki.fysik.dtu.dk/ase/tutorials/md/md.html#md-tutorial
# # ase gui moldyn3.traj
# # using 2 atoms only, but the forces are not linear! why ?