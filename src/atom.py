import numpy as np
from . import constants

class Atom:
    def __init__(self):
        self._pos = np.zeros(3)
        self._vel = np.zeros(3)
        self._acc = np.zeros(3)
        self._mass = 1

    @property
    def pos(self):
        return self._pos
    
    @pos.setter
    def pos(self, newPos):
        self._pos = newPos
    
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

    def move(self):
        self.pos += self.vel * constants.DELTATIME
        self.vel += self.acc * constants.DELTATIME
