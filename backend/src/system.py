from atom import Atom
# from calculator import Calculator
from lennardJones import LennardJones
from typing import List

class System:
    def __init__(self, listAtom: List[Atom]) -> None:
        self.listAtom = listAtom
        self.calc = LennardJones() #Calculator()

        # self.listAtom = self.getListAtom()

    # def getListAtom(self):
    #     listAtom = []
    #     for i in range(2):
    #         atom = Atom()
    #         listAtom.append(atom)

    def zeroResultantForces(self):
        for atom in self.listAtom:
            atom.resultantForce[:] = 0

    def move(self):
        for atom in self.listAtom:
            atom.move()

    def calculate(self):
        calc = LennardJones(self.listAtom)
        calc.calculate()
    



