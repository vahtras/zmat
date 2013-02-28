import os
from util.unformatted import FortranBinary
import numpy as np

def setup():
    global tdir
    tdir, tfile = os.path.split(__file__)

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
    ffile = os.path.join(tdir, 'fort.1')
    fb = FortranBinary(ffile)
    # first record is int 3
    next(fb)
    n = fb.readbuf(1, 'i')[0]
    assert n == 3
    # first record is float 1. 2. 3.
    next(fb)
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
    ffile = os.path.join(tdir, 'fort.2')
    fb = FortranBinary(ffile)
    rec  = fb.find('LABEL')

    assert rec == 'LABEL'

def test_3():
    """Iteration protocol

      integer, parameter :: n = 3
      double precision x(n), y(n)
      x = (/ 1.0D0, 2.0D0, 3.0D0 /)
      y = (/ 5.0D0, 6.0D0, 7.0D0 /)
      open(3, file='fort.3', status='new', form='unformatted')
      write(1) x
      write(1) y
      close(1)
      end
    """
    ffile = os.path.join(tdir, 'fort.3')
    fb = FortranBinary(ffile)
    # first record is int 3
    x=[]
    for rec in fb:
        x += list(fb.readbuf(3, 'd'))
    xref = (1., 2., 3.,  5., 6., 7.)
    print x, xref
    assert np.allclose(x, xref)

if __name__ == "__main__":
    setup()
    test_3()
