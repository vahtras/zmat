import math 
from numpy.testing import assert_allclose
from util import full
from ..zmat import Atom, Mol

tmpdir = '/tmp'

def setup():
    Atom.angular.clear()

def test_first():
    obj = Atom("H")
    assert obj.label == "H"
    assert obj.R is None
    assert obj.A is None
    assert obj.D is None
    assert obj.refs == []
    assert obj.charge == 1.0
    assert str(obj) ==  """       1.0    1
H         0.0000000000         0.0000000000         0.0000000000
"""

def test_second():
    obj = Atom("He 1 1.0")
    assert obj.label == "He"
    assert obj.R == "1.0"
    assert obj.r == 1.0
    assert obj.A is None
    assert obj.a is None
    assert obj.D is None
    assert obj.d is None
    assert obj.refs == [0]
    assert obj.charge == 2.0
    assert str(obj) ==  """       2.0    1
He        0.0000000000         0.0000000000         0.0000000000
"""

def test_second_par():
    obj = Atom("He 1 R")
    assert obj.label == "He"
    assert obj.R == "R"
    assert obj.r == None
    assert obj.A is None
    assert obj.a is None
    assert obj.D is None
    assert obj.d is None
    assert obj.refs == [0]
    assert obj.charge == 2.0
    assert str(obj) ==  """       2.0    1
He        0.0000000000         0.0000000000         0.0000000000
"""

def test_third():
    obj = Atom("Li 2 1.0 1 90")
    assert obj.label == "Li"
    assert obj.R == "1.0"
    assert obj.r == 1.0
    assert obj.A == "90"
    assert_allclose(obj.a, math.pi/2)
    assert obj.D is None
    assert obj.d is None
    assert obj.refs == [1, 0]
    assert obj.charge == 3.0
    assert str(obj) ==  """       3.0    1
Li        0.0000000000         0.0000000000         0.0000000000
"""

def test_third_par_R():
    obj = Atom("Li 2 R 1 90")
    assert obj.label == "Li"
    assert obj.R == "R"
    assert obj.r == None
    assert obj.A == "90"
    assert_allclose(obj.a, math.pi/2)
    assert obj.D is None
    assert obj.d is None
    assert obj.refs == [1, 0]
    assert obj.charge == 3.0
    assert str(obj) ==  """       3.0    1
Li        0.0000000000         0.0000000000         0.0000000000
"""

def test_third_par_A():
    obj = Atom("Li 2 1.0 1 A")
    assert obj.label == "Li"
    assert obj.R == "1.0"
    assert obj.r == 1.0
    assert obj.A == "A"
    assert obj.a == None
    assert obj.D is None
    assert obj.d is None
    assert obj.refs == [1, 0]
    assert obj.charge == 3.0
    assert str(obj) ==  """       3.0    1
Li        0.0000000000         0.0000000000         0.0000000000
"""
    #assert "A" in Atom.angular
    assert Atom.angular == set(['A'])

def test_third_par_RA():
    obj = Atom("Li 2 R 1 A")
    assert obj.label == "Li"
    assert obj.R == "R"
    assert obj.r == None
    assert obj.A == "A"
    assert obj.a == None
    assert obj.D is None
    assert obj.d is None
    assert obj.refs == [1, 0]
    assert obj.charge == 3.0
    assert str(obj) ==  """       3.0    1
Li        0.0000000000         0.0000000000         0.0000000000
"""
    assert "A" in Atom.angular
    print Atom.angular
    assert Atom.angular == set(['A'])

def test_general():
    obj = Atom("B 3 1.0 2 90 1 180")
    assert obj.label == "B"
    assert obj.R == "1.0"
    assert obj.r == 1.0
    assert obj.A == "90"
    assert_allclose(obj.a, math.pi/2)
    assert obj.D == "180"
    assert_allclose(obj.d,  math.pi)
    assert obj.refs == [2, 1, 0]
    assert obj.charge == 5.0
    assert str(obj) ==  """       5.0    1
B         0.0000000000         0.0000000000         0.0000000000
"""

def test_general_par_D():
    Atom.angular.clear()
    obj = Atom("B 3 1.0 2 90 1 D")
    assert obj.label == "B"
    assert obj.R == "1.0"
    assert obj.r == 1.0
    assert obj.A == "90"
    assert_allclose(obj.a, math.pi/2)
    assert obj.D == "D"
    assert obj.d is None
    assert obj.refs == [2, 1, 0]
    assert obj.charge == 5.0
    assert str(obj) ==  """       5.0    1
B         0.0000000000         0.0000000000         0.0000000000
"""
    print Atom.angular
    assert Atom.angular == set(['D'])


def zwrite(filename, string):
    f=open(os.path.join(tmpdir,filename), 'w')
    f.write(string)
    f.close()

def test_Mol_read_1():
    m = Mol(["H"])
    assert m.zmat == ["H"]
    assert m.params == {}
    assert len(m.atomlist) == 1
    assert "H" in m.atomtypes

def test_Mol_read_2():
    m = Mol(["H", "H 1 0.7"])
    assert m.zmat == ["H", "H 1 0.7"] 
    assert m.params == {}
    assert len(m.atomlist) == 2
    assert len(m.atomtypes) == 1
    assert "H" in m.atomtypes

def test_Mol_read_3():
    m = Mol([
"H",
"H 1 R",
"Variables:",
"R = 0.7"
    ])
    print m.zmat
    assert m.zmat == ["H", "H 1 R"] 
    assert m.params == {"R" : 0.7}
    assert len(m.atomlist) == 2
    assert len(m.atomtypes) == 1
    assert len(m.atomtypes["H"]) == 2
    assert "H" in m.atomtypes

def test_Mol_read_4():
    m = Mol([
"O",
"H 1 R",
"H 1 R 2 A",
"Variables:",
"R = 0.7",
"A = 120.0"
    ])
    print m.zmat
    assert m.zmat == ["O", "H 1 R", "H 1 R 2 A"] 
    print m.params["R"]; print m.params["A"]
    assert m.params == {"R" : 0.7, "A" : 120*math.pi/180}
    assert len(m.atomlist) == 3
    assert len(m.atomtypes) == 2
    assert "O" in m.atomtypes
    assert "H" in m.atomtypes
    assert len(m.atomtypes["O"]) == 1
    assert len(m.atomtypes["H"]) == 2

def test_Mol_read_4B():
    m = Mol([
"H",
"O 1 R",
"H 2 R 1 A",
"Variables:",
"R = 0.7",
"A = 120.0"
    ])
    print m.zmat
    assert m.zmat == ["H", "O 1 R", "H 2 R 1 A"] 
    print m.params["R"]; print m.params["A"]
    assert m.params == {"R" : 0.7, "A" : 120*math.pi/180}
    assert len(m.atomlist) == 3
    assert len(m.atomtypes) == 2
    assert "O" in m.atomtypes
    assert "H" in m.atomtypes
    assert len(m.atomtypes["O"]) == 1
    assert len(m.atomtypes["H"]) == 2

def test_Mol_read_5():
    m = Mol([
"N",
"X 1 1.0",
"H 1 R 2 A",
"H 1 R 2 A 3 120",
"H 1 R 2 A 4 120",
"Variables:",
"R=1.0058",
"A=110.2596"
    ])
    print m.zmat
    assert m.zmat == [
        "N", 
        "X 1 1.0", 
        "H 1 R 2 A", 
        "H 1 R 2 A 3 120",
        "H 1 R 2 A 4 120"
        ]
    print m.params["R"]; print m.params["A"]
    assert m.params == {"R" : 1.0058, "A" : 110.2596*math.pi/180}
    assert len(m.atomlist) == 5
    assert len(m.atomtypes) == 3
    assert "N" in m.atomtypes
    assert "H" in m.atomtypes
    assert len(m.atomtypes["N"]) == 1
    assert len(m.atomtypes["H"]) == 3

def test_update_cartesian_1():
    m = Mol(["H"])
    m.update_cartesian()
    assert_allclose(m.atomlist[0].coor.norm2(), 0)

def test_update_cartesian_2():
    m = Mol(["H", "H 1 0.7"])
    m.update_cartesian()
    H1, H2 = m.atomlist
    R12 = (H1.coor - H2.coor).norm2()
    assert_allclose(R12, 0.7)

def test_update_cartesian_3():
    m = Mol([
"H",
"H 1 R",
"Variables:",
"R = 0.7"
    ])
    m.update_cartesian()
    H1, H2 = m.atomlist
    R12 = (H1.coor - H2.coor).norm2()
    assert_allclose(R12, 0.7)


def test_Mol_linear():
    """Linear H3 molecule"""
    zm = """H
H 1 R
H 2 R 1 180
 Variables: 
R=0.925"""
    lines = [line.strip() for line in zm.split('\n')]
    print lines
    m = Mol(lines)
    m.update_cartesian()
    h1, h2, h3 = m.atomlist
    H1, H2, H3 = h1.coor, h2.coor, h3.coor

    H12 = (H1 - H2).norm2()
    H23 = (H2 - H3).norm2()
    assert_allclose(H12, 0.925)
    assert_allclose(H23, 0.925)

from ..sample_molecules import molinp

def assert_zmat_geometry(molecule, desired):
    print molecule
    lines = zinit(molecule)
    m = Mol(lines)
    m.update_cartesian()

    actual = []

    atoms = len(m.atomlist)

    if atoms > 1:
        b, a = m.atomlist[:2]
        actual.append(a.coor.dist(b.coor))

    if atoms > 2:
        a = m.atomlist[2]
        b, c = [m.atomlist[i] for i in a.refs]
        actual.append(a.coor.dist(b.coor))
        actual.append(a.coor.angle3d(b.coor, c.coor))

    for a in m.atomlist[3:]:
        b, c, d = [m.atomlist[i] for i in a.refs]
        actual.append(a.coor.dist(b.coor))
        actual.append(a.coor.angle3d(b.coor, c.coor))
        actual.append(a.coor.dihedrald(b.coor, c.coor, d.coor))


    for d, a in zip(desired, actual):
        print d, a
        assert_allclose(d, a)


def test_h2h():
    assert_zmat_geometry('h2h',[0.9246, 0.9246, 180])

def test_h():
    assert_zmat_geometry('h',[])

def test_h2():
    assert_zmat_geometry('h2',[0.7361])

def test_h2o():
    assert_zmat_geometry('h2o', [0.95, 0.95, 105])

def test_nh2():
    assert_zmat_geometry('nh2', [1.0191, 1.0191, 104.27])

def test_n2h():
    assert_zmat_geometry('n2h', [1.041, 1.169, 116.566])

def test_ph2():
    assert_zmat_geometry('ph2', [1.4149, 1.4149, 92.4849])

def test_nh3():
    assert_zmat_geometry('nh3', [
        1.0, 
        1.0058, 110.2596,
        1.0058, 110.2596, 120,
        1.0058, 110.2596, 120]
        )

def test_ph3():
    assert_zmat_geometry('ph3', [
        1.0, 
        1.410, 122.0756, 
        1.410, 122.0756, 120,
        1.410, 122.0756, 120
        ])

def test_h3():
    assert_zmat_geometry('h3', [0.925, 0.925, 180])

def test_h2h():
    assert_zmat_geometry('h2h', [0.9246, 0.9246, 180])

def test_ch3():
    assert_zmat_geometry('ch3', [1.074, 1.074, 120, 1.074, 120, 180])

def test_ch4():
    assert_zmat_geometry('ch4', [
        1.085,
        1.085, 109.4712206,
        1.085, 109.4712206, 120,
        1.085, 109.4712206, -120]
        )

def test_ch3oh():
    assert_zmat_geometry('ch3oh', [
        1.408,
        1.090, 112.020,
        1.084, 106.854, 241.4-360,
        0.953, 109.807, 298.6-360,
        1.090, 112.020, 122.8
        ]
    )


def zinit(case):
    lines = [line.strip() for line in molinp[case].split('\n')]
    # remove empty lines
    while '' in lines:
        lines.remove('')
    print "zinit:lines",lines
    return lines

def test_ch2oh():
    assert_zmat_geometry('ch2oh', [
        1.074,
        1.078, 120.906,
        1.356, 113.472, 150.2,
        0.953, 110.744, 182.6-360
        ]
    )

def test_c2h6():
    assert_zmat_geometry('c2h6',
       [1.522, 
        1.087, 111.295, 
        1.087, 111.295, 120,
        1.087, 111.295, -120,
        1.087, 111.295, 60,
        1.087, 111.295, 180,
        1.087, 111.295, -60
        ]
    )



def test_c2h5():
    assert_zmat_geometry('c2h5',
        [1.485, 
        1.094, 111.744,
        1.077, 120.892, -85.534,
        1.077, 120.892, 85.534,
        1.088, 111.661, 33.871,
        1.088, 111.661, -33.871,
        ]
    )

def test_n2h2():
    assert_zmat_geometry('n2h2', [1.023, 
        1.227, 107.343,
        1.023, 107.343, 180
        ])

def test_TS01():
    assert_zmat_geometry('TS01', [1.371, 1.114, 180])

def test_TS02():
    assert_zmat_geometry('TS02', [
        0.963,
        1.264, 100.451,
        0.848, 166.210, 0.0
        ]
    )

def test_TS03():
    assert_zmat_geometry('TS03', [
        1.080,
        1.080, 114.677,
        1.080, 114.677, 135.8,
        1.381, 103.565, 247.9 - 360,
        2.742,  53.912, 112.1
        ]
    )

def test_TS04():
    assert_zmat_geometry('TS04', [
        1.082,
        1.082, 113.227,
        1.244, 104.276, 116.7,
        1.251, 174.761, 239.5 - 360,
        0.952, 100.548, 0.0,
        1.082, 113.227, 229.4 - 360
        ]
    )

def test_TS05():
    assert_zmat_geometry('TS05', [
        1.374,
        1.087, 115.089,
        1.081, 109.674, 230.6 - 360,
        0.958, 110.524, 314.3 - 360,
        1.330, 110.196, 116.8,
        3.032,  45.235, 243.3 - 360
        ]
    )

def test_TS06():
    assert_zmat_geometry('TS06', [
        0.925,
        0.925, 180
        ]
    )


def test_TS07():
    assert_zmat_geometry('TS07', [
        0.959,
        1.248, 108.803,
        1.151, 150.779, 305.5 - 360,
        1.011, 110.384, 67.3,
        1.011, 108.054, 309.5 - 360
        ]
    )


def test_TS08():
    assert_zmat_geometry('TS08', [
        1.4265,
        1.3981, 180,
        1.0796, 101.66, 0,
        1.0796, 101.66, 120,
        1.0796, 101.66, -120
        ]
    )

def test_TS09():
    assert_zmat_geometry('TS09', [
        1.508,
        1.0866, 110.65,
        1.084, 111.01, 119.87,
        1.084, 111.01, -119.87,
        1.1670, 107.01, 180,
        1.0834, 113.45, 116.51,
        1.0834, 113.45, -116.51,
        1.3601, 175.42, 180,
        0.9675, 95.90, 0
        ]
    )

def test_TS10():
    assert_zmat_geometry('TS10', [
        1.440,
        0.772, 180
        ]
    )

def test_TS11():
    assert_zmat_geometry('TS11', [
        1.179,
        1.289, 360 - 181.05,
        1.082, 103.9, 180,
        1.081, 104.5, 120,
        1.081, 104.5, -120
        ]
    )

def test_TS12():
    assert_zmat_geometry('TS12', [
        1.2587,
        1.4851, 170.61,
        1.4152, 92.92, 45,
        1.4152, 92.92, -45
        ]
    )

def test_TS13():
    assert_zmat_geometry('TS13', [
        1.4850
        ]
    )

def test_TS14():
    assert_zmat_geometry('TS14', [
        1.2145,
        0.8941, 180
        ]
    )

def test_TS15():
    assert_zmat_geometry('TS15', [
        1.06,
        1.25, 105.8,
        1.17, 106.2, 0,
        1.20, 173.5, 0
        ]
    )

def test_TS16():
    assert_zmat_geometry('TS16', [
        1.2385,
        1.4046, 175.26,
        1.3372, 95.27, 0
        ]
    )

def test_TS17():
    assert_zmat_geometry('TS17', [
        2.44,
        2.66, 131.4
        ]
    )

def test_TS18():
    assert_zmat_geometry('TS18', [
        1.4447,
        1.0769, 100.77,
        1.0773, 104.51, 119.42,
        1.0773, 104.51, -119.42,
        1.2091, 175.14, 0,
        1.0244, 103.66, 0
        ]
    )

def test_TS19():
    assert_zmat_geometry('TS19', [
        1.3960,
        1.5118, 105.85,
        1.0820, 103.28, 113.50,
        1.0820, 103.28, -113.50,
        1.2326, 177.05, 0,
        1.0181, 99.33, 0,
        1.0887, 111.76, 0,
        1.0850, 107.77, 60.19,
        1.0850, 107.77, -60.19
        ]
    )

def test_TS20():
    assert_zmat_geometry('TS20', [
        1.517,
        1.254, 111.088,
        1.313, 175.404, 0,
        1.093, 111.088, 60,
        1.093, 111.088, -60,
        1.096, 111.088, 180,
        1.092, 114.157, 60,
        1.092, 114.157, -60,
        1.024, 100.007, 60,
        1.024, 100.007, -60,
        ]
    )


def test_zmat2her():
    zref = """O
H 1 R
H 1 R 2 A
Variables:
R=0.95
A=105
"""
    href="""BASIS
<your basis set label here>
Generated from yo.zmat by molconvert
====================================
Atomtypes=2 Units=Angtrom Nosymmetry
Charge=1.0 Atoms=2
H   0.950000   0.000000   0.000000
H  -0.245878   0.917630   0.000000
Charge=8.0 Atoms=1
O   0.000000   0.000000   0.000000
"""
    from ..molconvert import main
    z=open('yo.zmat', 'w'); z.write(zref); z.close()
    main('yo.zmat', 'yo.mol')
    h=open('yo.mol', 'r').read()
    assert href == h


if __name__ == "__main__":
    test_TS07()
