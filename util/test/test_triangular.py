from util.full import triangular
def test_init():
    ft = triangular.init([1.,  0.,  1.])
    assert ft.sshape == (2, 2)
    assert ft[0, 0] == 1. 
    assert ft[1, 0] == 0. 
    assert ft[1, 1] == 1. 
