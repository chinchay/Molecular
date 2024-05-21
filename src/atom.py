import numpy as np
# from . import constants
from constants import DELTATIME

class Atom:
    def __init__(self):
        self._pos = np.zeros(3, dtype=float)
        self._vel = np.zeros(3)
        self._acc = np.zeros(3)
        self._mass = 1

    @property
    def pos(self):
        return self._pos
    
    @pos.setter
    def pos(self, newPos):
        self._pos = newPos.copy()
    
    @property
    def vel(self):
        return self._vel
    
    @vel.setter
    def vel(self, newVel):
        self._vel = newVel
    
    @property
    def acc(self):
        return self._acc

    @acc.setter
    def acc(self, newAcc):
        self._acc = newAcc
    
    @property
    def mass(self):
        return self._mass
    
    @mass.setter
    def mass(self, newMass):
        self._mass = newMass

    def move(self):
        # self.pos += self.vel * constants.DELTATIME
        # self.vel += self.acc * constants.DELTATIME

        dPos = self.vel * DELTATIME
        dVel = self.acc * DELTATIME
        # print("dPos = ", dPos, " dAcc = ", dAcc)
        self.pos = self.pos + dPos
        self.vel = self.vel + dVel
