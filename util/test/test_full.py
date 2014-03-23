import numpy as np
from util.full import matrix

def assert_(this, ref):
    print this
    print ref
    assert np.allclose(this, ref)

def test_diag():
    ref = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    this = matrix.diag([1,1,1])
    assert_(this, ref)
    
