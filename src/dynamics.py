from atom import Atom
import numpy as np
from  constants import SIGMA12, SIGMA6, EPSILON
from math import pow

def computeDistance(atom1: Atom, atom2: Atom):
    dR21 = atom2.pos - atom1.pos
    distance = np.linalg.norm(dR21)
    return dR21, distance

def computeLennardJonesForce(distance: float):
    repulsive = SIGMA12 / pow(distance, 13)
    attractive = 0.5 * SIGMA6 / pow(distance, 7)
    force = 48 * EPSILON * (repulsive - attractive)
    return force

def updateLennarJonesAcc(atom1: Atom, atom2: Atom):
    dR21, distance = computeDistance(atom1, atom2)
    absForce = computeLennardJonesForce(distance)
    u21 = dR21 / distance
    force21 = absForce * u21
    force12 = -force21
    atom2.acc = force21 / atom2.mass
    atom1.acc = force12 / atom1.mass