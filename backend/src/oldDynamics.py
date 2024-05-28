from atom import Atom
from atoms import Atoms
from system import System
import numpy as np
from constants import SIGMA, SIGMA2, SIGMA12, SIGMA6, EPSILON
from math import pow, sqrt
from calculator import Calculator
from lennardJones import LennardJones
from typing import List

# class Dynamics:
#     def __init__(self, system: System):
#         self.system = system



def computeDistance(atom1: Atom, atom2: Atom):
    dR21 = atom2.pos - atom1.pos

    # distance = np.linalg.norm(dR21)
    # this is more efficient:
    distance = sqrt((dR21[0] ** 2) + (dR21[1] ** 2) + (dR21[2] ** 2))
    
    return dR21, distance

def computeLennardJonesForce(distance: float):
    repulsive = SIGMA12 / pow(distance, 13)
    attractive = 0.5 * SIGMA6 / pow(distance, 7)
    force = 48 * EPSILON * (repulsive - attractive)
    return force

# def updateLennarJonesAcc(atom1: Atom, atom2: Atom):
#     dR21, distance = computeDistance(atom1, atom2)
#     force = computeLennardJonesForce(distance)
#     u21 = dR21 / distance
#     force21 = force * u21
#     force12 = -force21
#     atom2.acc = force21 / atom2.mass
#     atom1.acc = force12 / atom1.mass

def computeForceVector(atom1: Atom, atom2: Atom):
    dR21, distance = computeDistance(atom1, atom2)
    force = computeLennardJonesForce(distance)
    u21 = dR21 / distance
    force21 = force * u21
    return force21

def updateResultantForces(system: System):
    system.zeroResultantForces()
    nAtoms = len(system.listAtom)
    for i in range(nAtoms):
        atomI = system.listAtom[i]
        for j in range(i + 1, nAtoms):
            atomJ = system.listAtom[j]
            forceJI = computeForceVector(atomI, atomJ)
            forceIJ = -forceJI
            atomJ.resultantForce += forceJI
            atomI.resultantForce += forceIJ


def updateAcc(system: System):
    updateResultantForces(system)
    for atom in system.listAtom:
        atom.acc = atom.resultantForce / atom.mass

# def updateAcc(atoms: Atoms):
#     calculate(atoms) # compute energies and forces
#     atoms.accelerations = atoms.forces / atoms.masses


def computeKinetic(system: System):
    kinetic = 0
    for atom in system.listAtom:
        vel = atom.vel
        mass = atom.mass
        kinetic += np.dot(vel, vel) * mass
    kinetic *= 0.5
    return kinetic

# def computeKinetic(atoms: Atoms):
#     squaredVelocityperAtom = np.sum(atoms.velocities ** 2, axis=1)
#     kinetic = 0.5 * np.sum(atoms.masses * squaredVelocityperAtom)
#     return kinetic

def computeLennardJonesPot(atom1: Atom, atom2: Atom):
    dR21, distance = computeDistance(atom1, atom2)
    x = SIGMA / distance
    pot = 4 * SIGMA * (pow(x, 12) - pow(x, 6))
    return pot

def computePotential(system: System):
    pot = 0
    nAtoms = len(system.listAtom)
    for i in range(nAtoms):
        atomI = system.listAtom[i]
        for j in range(i + 1, nAtoms):
            atomJ = system.listAtom[j]
            dPot = computeLennardJonesPot(atomI, atomJ)
            pot += dPot
    return pot

def getNeighbors(positions: np.ndarray, indices: List[int], index: int):
    # nAtoms, _ = positions.shape
    # listNeighborIndx = np.delete(np.arange(nAtoms), index)
    listNeighborIndx = np.delete(indices, index)
    return listNeighborIndx

def computeForceLJ(relatPos, c6, c12, squaredDist):
    pairwiseForcesOverDist = -24 * EPSILON * ((2 * c12) - c6) / squaredDist
    #
    # Compute the axis force components for each pair of atoms
    # `pairwiseForcesOverDist` multiplied by `relatPos` has units
    # of forces. Notice that `relatPos / dist` is a unit vector!
    forcePerAxisPerPair = pairwiseForcesOverDist[:, np.newaxis] * relatPos
    force = np.sum(forcePerAxisPerPair, axis=0)
    return force

def computePotentialLJ(c6, c12):
    pairwiseEnergies = 4 * EPSILON * (c12 - c6)
    energy = np.sum(pairwiseEnergies)
    return energy

def computeContributions(squaredDist):
    c6 = (SIGMA2 / squaredDist) ** 3
    c12 = c6 ** 2
    return c6, c12
def computeSquaredPosAndRelatPos(positions, listNeighborIndx, i):
    relatPos = positions[listNeighborIndx] - positions[i]
    squaredDist = np.sum(relatPos ** 2, axis=1)
    return relatPos, squaredDist

def calculate(atoms: Atoms):
    positions = atoms.positions
    nAtoms, _ = positions.shape
    indices = np.arange(nAtoms)
    for i in range(nAtoms):
        listNeighborIndx = getNeighbors(positions, indices, i)
        # relatPos = positions[listNeighborIndx] - positions[i]
        # squaredDist = np.sum(relatPos ** 2, axis=1)

        relatPos, squaredDist = computeSquaredPosAndRelatPos(positions, listNeighborIndx, i)


        # c6 = (SIGMA2 / squaredDist) ** 3
        # c12 = c6 ** 2

        c6, c12 = computeContributions(squaredDist)
        atoms.forces[i] = computeForceLJ(relatPos, c6, c12, squaredDist)
        atoms.potentialEnergies[i] = computePotentialLJ(c6, c12)

def computeMechanicalEnergy(system: System):
    kinetic = computeKinetic(system)
    potential = computePotential(system)
    energy = kinetic + potential
    return energy

# def computeMechanicalEnergy(atoms: Atoms):
#     kinetic = computeKinetic(atoms)
#     potential = computePotential(atoms)
#     energy = kinetic + potential
#     return energy


# This is what should be
# def applyThermalBath(system: System, desiredKinetic: float):
#     currentKinetic = computeKinetic(system)
#     lamb = sqrt(desiredKinetic / currentKinetic)
#     for atom in system:
#         atom.vel *= lamb

# # This is what worked at this level, having only two atoms
def applyThermalBath(system: System, desiredMechanicalEnergy: float):
    currentKinetic = computeKinetic(system)
    try:
        lamb = sqrt(abs(desiredMechanicalEnergy / currentKinetic))
    except:
        print("Error with kinetic as denominotar = ", currentKinetic)
        return

    for atom in system.listAtom:
        atom.vel *= lamb

# # This is what worked at this level, having only two atoms
# def applyThermalBath(atoms: Atoms, desiredMechanicalEnergy: float):
#     currentKinetic = computeKinetic(atoms)
#     try:
#         lamb = sqrt(abs(desiredMechanicalEnergy / currentKinetic))
#     except:
#         print("Error with kinetic as denominotar = ", currentKinetic)
#         return
    
#     atoms.velocities *= lamb
