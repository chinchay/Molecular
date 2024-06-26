import numpy as np
import matplotlib.pyplot as plt

from atom import Atom
from system import System
from dynamics import applyThermalBath, computeMechanicalEnergy, computeDistance, computeKinetic, updateAcc
# from . import constants
from constants import RMIN

listAtom = []
atom = Atom()
atom.pos = np.asarray([0, 0, 0], dtype=float)
atom.vel = np.asarray([0, 0, 0], dtype=float)
atom.acc = np.asarray([0, 0, 0], dtype=float)
atom.mass = 100
listAtom.append(atom)

atom = Atom()
delta = 0.1
d = RMIN + delta
atom.pos = np.asarray([d, 0, 0], dtype=float)
atom.vel = np.asarray([0, 0, 0], dtype=float)
atom.acc = np.asarray([0, 0, 0], dtype=float)
atom.mass = 0.5
listAtom.append(atom)

atom = Atom()
delta = 0.1
d = RMIN + delta
atom.pos = np.asarray([0, d, 0], dtype=float)
atom.vel = np.asarray([0, 0, 0], dtype=float)
atom.acc = np.asarray([0, 0, 0], dtype=float)
atom.mass = 0.5
listAtom.append(atom)

atom = Atom()
delta = 0.1
d = RMIN + delta
atom.pos = np.asarray([-d, 0, 0], dtype=float)
atom.vel = np.asarray([0, 0, 0], dtype=float)
atom.acc = np.asarray([0, 0, 0], dtype=float)
atom.mass = 0.5
listAtom.append(atom)

atom = Atom()
delta = 0.1
d = RMIN + delta
atom.pos = np.asarray([0, -d, 0], dtype=float)
atom.vel = np.asarray([0, 0, 0], dtype=float)
atom.acc = np.asarray([0, 0, 0], dtype=float)
atom.mass = 0.5
listAtom.append(atom)

# ################################################################################
y = 6
atom = Atom()
atom.pos = np.asarray([0, y, 0], dtype=float)
atom.vel = np.asarray([0, 0, 0], dtype=float)
atom.acc = np.asarray([0, 0, 0], dtype=float)
atom.mass = 100
listAtom.append(atom)

atom = Atom()
delta = 0.1
d = RMIN + delta
atom.pos = np.asarray([d, y, 0], dtype=float)
atom.vel = np.asarray([0, 0, 0], dtype=float)
atom.acc = np.asarray([0, 0, 0], dtype=float)
atom.mass = 0.5
listAtom.append(atom)

atom = Atom()
delta = 0.1
d = RMIN + delta
atom.pos = np.asarray([0, y + d, 0], dtype=float)
atom.vel = np.asarray([0, 0, 0], dtype=float)
atom.acc = np.asarray([0, 0, 0], dtype=float)
atom.mass = 0.5
listAtom.append(atom)

atom = Atom()
delta = 0.1
d = RMIN + delta
atom.pos = np.asarray([-d, y, 0], dtype=float)
atom.vel = np.asarray([0, 0, 0], dtype=float)
atom.acc = np.asarray([0, 0, 0], dtype=float)
atom.mass = 0.5
listAtom.append(atom)

atom = Atom()
delta = 0.1
d = RMIN + delta
atom.pos = np.asarray([0, y - d, 0], dtype=float)
atom.vel = np.asarray([0, 0, 0], dtype=float)
atom.acc = np.asarray([0, 0, 0], dtype=float)
atom.mass = 0.5
listAtom.append(atom)

################################################################################
x = 6
atom = Atom()
atom.pos = np.asarray([x, y, 0], dtype=float)
atom.vel = np.asarray([0, 0, 0], dtype=float)
atom.acc = np.asarray([0, 0, 0], dtype=float)
atom.mass = 100
listAtom.append(atom)

atom = Atom()
delta = 0.1
d = RMIN + delta
atom.pos = np.asarray([x + d, y, 0], dtype=float)
atom.vel = np.asarray([0, 0, 0], dtype=float)
atom.acc = np.asarray([0, 0, 0], dtype=float)
atom.mass = 0.5
listAtom.append(atom)

atom = Atom()
delta = 0.1
d = RMIN + delta
atom.pos = np.asarray([x, y + d, 0], dtype=float)
atom.vel = np.asarray([0, 0, 0], dtype=float)
atom.acc = np.asarray([0, 0, 0], dtype=float)
atom.mass = 0.5
listAtom.append(atom)

atom = Atom()
delta = 0.1
d = RMIN + delta
atom.pos = np.asarray([x - d, y, 0], dtype=float)
atom.vel = np.asarray([0, 0, 0], dtype=float)
atom.acc = np.asarray([0, 0, 0], dtype=float)
atom.mass = 0.5
listAtom.append(atom)

atom = Atom()
delta = 0.1
d = RMIN + delta
atom.pos = np.asarray([x, y - d, 0], dtype=float)
atom.vel = np.asarray([0, 0, 0], dtype=float)
atom.acc = np.asarray([0, 0, 0], dtype=float)
atom.mass = 0.5
listAtom.append(atom)

system = System(listAtom)

# for atom in system.listAtom:
#     print(atom.pos)

initialEnergy = computeMechanicalEnergy(system)
desiredKinetic = initialEnergy / 2

n = 320#1
fig, ax = plt.subplots()
ax.set_xlim(-4, 8), ax.set_xticks([])
ax.set_ylim(-4, 8), ax.set_xticks([])

nAtoms = len(system.listAtom)
listScatter = [ax.scatter([], [], color="b") for _ in range(nAtoms)] 

isInsideBall = False
for i in range(n):
    system.move()

    for j in range(nAtoms):
        [x, y, z] = system.listAtom[j].pos
        x_values = [x] 
        y_values = [y]
        listScatter[j].set_offsets(list(zip(x_values, y_values)))


    plt.pause(0.01)

    # updateLennarJonesAcc(atom1, atom2)
    updateAcc(system)

    _, distance = computeDistance(listAtom[0], listAtom[1]) #<<< include other atoms!

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