from util.unformatted import FortranBinary
import numpy as np

def test_1():
    """Read int, float

      integer, parameter :: n = 3
      double precision x(n)
      x = (/ 1.0D0, 2.0D0, 3.0D0 /)
      open(1, file='fort.1', status='new', form='unformatted')
      write(1) n
      write(1) x
      close(1)
      end
    """
    fb = FortranBinary('fort.1')
    # first record is int 3
    fb.readrec()
    n = fb.readbuf(1, 'i')[0]
    assert n == 3
    # first record is float 1. 2. 3.
    fb.readrec()
    xref = (1., 2., 3.)
    x = fb.readbuf(n, 'd')
    assert np.allclose(x, xref)

def test_2():
    """Find and read label

      character*5 lab
      integer n
      lab = 'LABEL'
      n = 0
      open(1, file='fort.2', status='new', form='unformatted')
      write(1) n
      write(1) lab
      close(1)
      end
    """
    fb = FortranBinary('fort.2')
    rec  = fb.find('LABEL')

    assert rec == 'LABEL'
