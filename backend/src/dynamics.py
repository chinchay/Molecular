from atom import Atom
from system import System
import numpy as np
from  constants import SIGMA, SIGMA12, SIGMA6, EPSILON
from math import pow, sqrt

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

def updateLennarJonesAcc(atom1: Atom, atom2: Atom):
    dR21, distance = computeDistance(atom1, atom2)
    force = computeLennardJonesForce(distance)
    u21 = dR21 / distance
    force21 = force * u21
    force12 = -force21
    atom2.acc = force21 / atom2.mass
    atom1.acc = force12 / atom1.mass

def computeKinetic(system: System):
    kinetic = 0
    for atom in system.listAtom:
        vel = atom.vel
        mass = atom.mass
        kinetic += np.dot(vel, vel) * mass
    kinetic *= 0.5
    return kinetic

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

def computeMechanicalEnergy(system: System):
    kinetic = computeKinetic(system)
    potential = computePotential(system)
    energy = kinetic + potential
    return energy

# This is what should be
# def applyThermalBath(system: System, desiredKinetic: float):
#     currentKinetic = computeKinetic(system)
#     lamb = sqrt(desiredKinetic / currentKinetic)
#     for atom in system:
#         atom.vel *= lamb

# This is what worked at this level, having only two atoms
def applyThermalBath(system: System, desiredMechanicalEnergy: float):
    currentKinetic = computeKinetic(system)
    try:
        lamb = sqrt(abs(desiredMechanicalEnergy / currentKinetic))
    except:
        print("Error with kinetic as denominotar = ", currentKinetic)
        return

    for atom in system.listAtom:
        atom.vel *= lamb