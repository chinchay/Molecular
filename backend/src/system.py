from atom import Atom
from typing import List

class System:
    def __init__(self, listAtom: List[Atom]) -> None:
        self.listAtom = listAtom

        # self.listAtom = self.getListAtom()

    # def getListAtom(self):
    #     listAtom = []
    #     for i in range(2):
    #         atom = Atom()
    #         listAtom.append(atom)



