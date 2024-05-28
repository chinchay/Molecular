from calculator import Calculator
from atom import Atom
import numpy as np
from constants import SIGMA
from typing import List

class LennardJones(Calculator):
    def __init__(self) -> None:
        Calculator.__init__(self)
        pass

    def computeDistance(atom1: Atom, atom2: Atom):
        dR21 = atom2.pos - atom1.pos

    def calculate(self, listAtom: List[Atom]):
        nAtoms = len(listAtom)
        positions = np.asarray([ listAtom[i].pos for i in range(nAtoms) ])
        indices = np.arange(nAtoms)
        
        for i in range(nAtoms):
            atomI = listAtom[i]
            neighbors = np.delete(indices, i)
            distance_vectors = positions[neighbors] - positions[i]
            r2 = np.sum(distance_vectors ** 2, axis=1)
            c6 = (SIGMA ** 2 / r2) ** 3
            # c6[r2 > rc ** 2] = 0.0
            c12 = c6 ** 2

            for j in range(i + 1, nAtoms):
                atomJ = listAtom[j]
                forceJI = computeForceVector(atomI, atomJ)
                forceIJ = -forceJI
                atomJ.resultantForce += forceJI
                atomI.resultantForce += forceIJ

        # for i in range(nAtoms):
        #     # neighbors, offsets = self.nl.get_neighbors(i)
        #     # cells = np.dot(offsets, cell)

        #     # pointing *towards* neighbours
        #     distance_vectors = positions[neighbors] + cells - positions[i]

        #     r2 = (distance_vectors ** 2).sum(1)
        #     c6 = (sigma ** 2 / r2) ** 3
        #     c6[r2 > rc ** 2] = 0.0
        #     c12 = c6 ** 2

        #     if smooth:
        #         cutoff_fn = cutoff_function(r2, rc**2, ro**2)
        #         d_cutoff_fn = d_cutoff_function(r2, rc**2, ro**2)

        #     pairwise_energies = 4 * epsilon * (c12 - c6)
        #     pairwise_forces = -24 * epsilon * (2 * c12 - c6) / r2  # du_ij


        pass