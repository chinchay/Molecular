import numpy as np
# from . import constants
from constants import DELTATIME

class Atom:
    def __init__(self):
        self.pos = np.zeros(3, dtype=float)
        self.vel = np.zeros(3, dtype=float)
        self.acc = np.zeros(3, dtype=float)
        self.mass = 1
        self.resultantForce = np.zeros(3, dtype=float)

    # @property
    # def pos(self):
    #     return self._pos
    
    # @pos.setter
    # def pos(self, newPos):
    #     self._pos = newPos.copy()
    
    # @property
    # def vel(self):
    #     return self._vel
    
    # @vel.setter
    # def vel(self, newVel):
    #     self._vel = newVel
    
    # @property
    # def acc(self):
    #     return self._acc

    # @acc.setter
    # def acc(self, newAcc):
    #     self._acc = newAcc
    
    # @property
    # def mass(self):
    #     return self._mass
    
    # @mass.setter
    # def mass(self, newMass):
    #     self._mass = newMass

    def move(self):
        # self.pos += self.vel * DELTATIME
        # self.vel += self.acc * DELTATIME

        dPos = self.vel * DELTATIME
        dVel = self.acc * DELTATIME
        self.pos += dPos
        self.vel += dVel

        pass

