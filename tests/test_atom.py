import sys
sys.path.append("..")
import pytest

import numpy as np

from src.atom import Atom



def testPos():
    atom = Atom()
    assert np.linalg.norm(atom.pos) == 0

    pos = np.asarray([3, 4, 0])
    atom.pos = pos
    assert np.linalg.norm(atom.pos) == 5
