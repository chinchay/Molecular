import numpy as np
from constants import DELTATIME

class Atoms:
    def __init__(self, nAtoms: int) -> None:
        self.positions = np.zeros((nAtoms, 3), dtype=float)
        self.velocities = np.zeros((nAtoms, 3), dtype=float)
        self.accelerations = np.zeros((nAtoms, 3), dtype=float)
        self.masses = np.ones((nAtoms, 1), dtype=float)
        # self.internalResultantForces = None
        self.forces = np.zeros((nAtoms, 3), dtype=float)
        self.potentialEnergies = np.zeros((nAtoms, 3), dtype=float)
    
    def move(self):
        self.positions += self.velocities * DELTATIME
        self.velocities += self.accelerations * DELTATIME
