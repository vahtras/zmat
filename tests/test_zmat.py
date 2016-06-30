import math 
import unittest
import numpy.testing
from util import full
from ..zmat import Atom, Mol
from .sample_molecules import molinp

tmpdir = '/tmp'

def zinit(case):
    lines = [line.strip() for line in molinp[case].split('\n')]
    # remove empty lines
    while '' in lines:
        lines.remove('')
    print "zinit:lines",lines
    return lines


class AtomTest(unittest.TestCase):

    def setUp(self):
        Atom.atomlist = []

    def tearDown(self):
        pass

    def assert_allclose(self, *args, **kwargs):
        numpy.testing.assert_allclose(*args, **kwargs)

    def assert_zmat_geometry(self, molecule, desired):
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
            self.assert_allclose(d, a)

    def test_first(self):
        obj = Atom("H")
        self.assertTupleEqual(
            (obj.label, obj.R, obj.A, obj.D, obj.refs, obj.charge),
            ("H", None, None, None, [], 1.0)
            )

    def test_first_coor(self):
        obj = Atom("H")
        self.assert_allclose(obj.coor, (0, 0, 0))

    def test_second(self):
        Atom("H")
        obj = Atom("He 1 1.0")
        self.assertTupleEqual(
            (obj.label, obj.R, obj.r, obj.A, obj.a, obj.D, obj.d, obj.refs, obj.charge),
            ("He", "1.0", 1.0, None, None, None, None, [0], 2.0)
            )

    def test_second_coor(self):
        Atom("H")
        obj = Atom("He 1 1.0")
        obj.update_cartesian()
        self.assert_allclose(obj.coor, (1, 0, 0))

    def test_second_par_undef(self):
        Atom("H")
        obj = Atom("He 1 R")
        self.assertTupleEqual(
            (obj.label, obj.R, obj.r, obj.A, obj.a, obj.D, obj.d, obj.refs, obj.charge),
            ("He", "R", None, None, None, None, None, [0], 2.0)
            )

    def test_second_par(self):
        Atom("H")
        obj = Atom("He 1 R")
        Atom.params.update({'R': 1.0})
        self.assertTupleEqual(
            (obj.label, obj.R, obj.r, obj.A, obj.a, obj.D, obj.d, obj.refs,
                obj.charge),
            ("He", "R", 1.0, None, None, None, None, [0], 2.0)
            )

    def test_second_par_cor(self):
        Atom("H")
        obj = Atom("He 1 R")
        Atom.params.update({'R': 1.0})
        obj.update_cartesian()
        self.assert_allclose(obj.coor, (1, 0, 0))

    def test_third(self):
        Atom("H")
        Atom("He 1 R")
        obj = Atom("Li 2 1.0 1 90")
        self.assertTupleEqual(
            (obj.label, obj.R, obj.r, obj.A, obj.a, obj.D, obj.d, obj.refs,
                obj.charge),
            ("Li", "1.0", 1.0, "90", math.pi/2, None, None, [1, 0], 3.0)
            )

    def test_third_coor(self):
        Atom("H")
        Atom("He 1 R")
        obj = Atom("Li 2 1.0 1 90")
        obj.params["R"] = 1.0
        obj.update_cartesian()
        self.assert_allclose(obj.coor, (1, 1, 0))

    def test_third_par_R(self):
        Atom("H")
        Atom("He 1 R")
        obj = Atom("Li 2 R 1 90")
        self.assertTupleEqual(
            (obj.label, obj.R, obj.r, obj.A, obj.a, obj.D, obj.d, obj.refs,
                obj.charge),
            ("Li", "R", None, "90", math.pi/2, None, None, [1, 0], 3.0)
            )

    def test_third_par_R(self):
        Atom("H")
        Atom("He 1 R")
        obj = Atom("Li 2 1.0 1 90")
        obj.params["R"] = 1.0
        obj.update_cartesian()
        self.assert_allclose(obj.coor, (1, 1, 0))

    def test_third_par_A(self):
        Atom("H")
        Atom("He 1 R")
        obj = Atom("Li 2 1.0 1 A")
        obj.params["R"] = 1.0
        obj.params["A"] = math.pi/2
        obj.update_cartesian()
        self.assert_allclose(obj.coor, (1, 1, 0))

    def test_third_par_RA(self):
        Atom("H")
        Atom("He 1 R")
        obj = Atom("Li 2 R 1 A")
        obj.params.update({"R": 1.0, "A": math.pi/2})
        obj.update_cartesian()
        self.assert_allclose(obj.coor, (1, 1, 0))

    def test_general(self):
        Atom("H")
        Atom("He 1 1.0")
        Atom("Li 2 1.0 1 90")
        obj = Atom("B 3 1.0 2 90 1 180")
        self.assertTupleEqual(
            (obj.label, obj.R, obj.r, obj.A, obj.a, obj.D, obj.d, obj.refs,
                obj.charge),
            ("B", "1.0", 1.0, "90", math.pi/2, "180", math.pi, [2, 1, 0], 5.0)
            )

    def test_general_coor(self):
        Atom("H")
        Atom("He 1 1.0")
        Atom("Li 2 1.0 1 90")
        obj = Atom("B 3 1.0 2 90 1 180")
        obj.update_cartesian()
        self.assert_allclose(obj.coor, (2, 1, 0), atol=1e-7)

    def test_general_par(self):
        Atom("H")
        Atom("He 1 R")
        Atom("Li 2 R 1 A")
        obj = Atom("B 3 R 2 A 1 D")
        obj.params.update({'R': 1.0, 'A': math.pi/2, 'D': math.pi})
        self.assertTupleEqual(
            (obj.label, obj.R, obj.r, obj.A, obj.a, obj.D, obj.d, obj.refs,
                obj.charge),
            ("B", "R", 1.0, "A", math.pi/2, "D", math.pi, [2, 1, 0], 5.0)
            )

    def test_Mol_read_1(self):
        m = Mol(["H"])
        self.assertListEqual(m.zmat, ["H"])
        self.assertDictEqual(m.params, {})
        self.assertEqual(len(m.atomlist), 1)
        self.assertIn("H", m.atomtypes)

    def test_Mol_read_2(self):
        m = Mol(["H", "H 1 0.7"])
        self.assertListEqual(m.zmat, ["H", "H 1 0.7"] )
        self.assertDictEqual(m.params, {})
        self.assertEqual(len(m.atomlist), 2)
        self.assertEqual(len(m.atomtypes), 1)
        self.assertIn("H", m.atomtypes)

    def test_Mol_read_3(self):
        m = Mol([
    "H",
    "H 1 R",
    "Variables:",
    "R = 0.7"
        ])
        self.assertListEqual(m.zmat, ["H", "H 1 R"])
        self.assertDictEqual(m.params, {"R" : 0.7})
        self.assertEqual(len(m.atomlist), 2)
        self.assertEqual(len(m.atomtypes), 1)
        self.assertEqual(len(m.atomtypes["H"]), 2)
        self.assertIn("H", m.atomtypes)


    def test_Mol_read_3(self):
        m = Mol([
    "H",
    "H 1 R",
    "Variables:",
    "R = 0.7"
        ])
        self.assertListEqual(m.zmat, ["H", "H 1 R"])
        self.assertDictEqual( m.params, {"R" : 0.7})
        self.assertEqual(len(m.atomlist), 2)
        self.assertEqual(len(m.atomtypes), 1)
        self.assertEqual(len(m.atomtypes["H"]), 2)
        self.assertIn("H", m.atomtypes)

    def test_Mol_read_4(self):
        m = Mol([
    "O",
    "H 1 R",
    "H 1 R 2 A",
    "Variables:",
    "R = 0.7",
    "A = 120.0"
        ])
        self.assertListEqual(m.zmat, ["O", "H 1 R", "H 1 R 2 A"])
        self.assertDictEqual(m.params, {"R" : 0.7, "A" : 120*math.pi/180})
        self.assertEqual(len(m.atomlist), 3)
        self.assertEqual(len(m.atomtypes), 2)
        self.assertIn("O", m.atomtypes)
        self.assertIn("H", m.atomtypes)
        self.assertEqual(len(m.atomtypes["O"]), 1)
        self.assertEqual(len(m.atomtypes["H"]), 2)

    def test_Mol_read_4B(self):
            m = Mol([
        "H",
        "O 1 R",
        "H 2 R 1 A",
        "Variables:",
        "R = 0.7",
        "A = 120.0"
            ])
            self.assertListEqual(m.zmat, ["H", "O 1 R", "H 2 R 1 A"])
            self.assertDictEqual(m.params, {"R" : 0.7, "A" : 120*math.pi/180})
            self.assertEqual(len(m.atomlist), 3)
            self.assertEqual(len(m.atomtypes), 2)
            self.assertIn("O", m.atomtypes)
            self.assertIn("H", m.atomtypes)
            self.assertEqual(len(m.atomtypes["O"]), 1)
            self.assertEqual(len(m.atomtypes["O"]), 1)

    def test_Mol_read_5(self):
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
        self.assertListEqual(m.zmat, [
            "N", 
            "X 1 1.0", 
            "H 1 R 2 A", 
            "H 1 R 2 A 3 120",
            "H 1 R 2 A 4 120"
            ])
        self.assertDictEqual(m.params, {"R" : 1.0058, "A" :
            110.2596*math.pi/180})
        self.assertEqual(len(m.atomlist), 5)
        self.assertEqual(len(m.atomtypes), 3)
        self.assertIn("N", m.atomtypes)
        self.assertIn("H", m.atomtypes)
        self.assertEqual(len(m.atomtypes["N"]), 1)
        self.assertEqual(len(m.atomtypes["H"]), 3)

    def test_update_cartesian_1(self):
        m = Mol(["H"])
        m.update_cartesian()
        self.assert_allclose(m.atomlist[0].coor.norm2(), 0)

    def test_update_cartesian_2(self):
        m = Mol(["H", "H 1 0.7"])
        m.update_cartesian()
        H1, H2 = m.atomlist
        R12 = (H1.coor - H2.coor).norm2()
        self.assert_allclose(R12, 0.7)

    def test_update_cartesian_3(self):
        m = Mol([
    "H",
    "H 1 R",
    "Variables:",
    "R = 0.7"
        ])
        m.update_cartesian()
        H1, H2 = m.atomlist
        R12 = (H1.coor - H2.coor).norm2()
        self.assert_allclose(R12, 0.7)


    def test_Mol_linear(self):
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
        self.assert_allclose(H12, 0.925)
        self.assert_allclose(H23, 0.925)




    def test_h2h(self):
        self.assert_zmat_geometry('h2h',[0.9246, 0.9246, 180])

    def test_h(self):
        self.assert_zmat_geometry('h',[])

    def test_h2(self):
        self.assert_zmat_geometry('h2',[0.7361])

    def test_h2o(self):
        self.assert_zmat_geometry('h2o', [0.95, 0.95, 105])

    def test_nh2(self):
        self.assert_zmat_geometry('nh2', [1.0191, 1.0191, 104.27])

    def test_n2h(self):
        self.assert_zmat_geometry('n2h', [1.041, 1.169, 116.566])

    def test_ph2(self):
        self.assert_zmat_geometry('ph2', [1.4149, 1.4149, 92.4849])

    def test_nh3(self):
        self.assert_zmat_geometry('nh3', [
            1.0, 
            1.0058, 110.2596,
            1.0058, 110.2596, 120,
            1.0058, 110.2596, 120]
            )

    def test_ph3(self):
        self.assert_zmat_geometry('ph3', [
            1.0, 
            1.410, 122.0756, 
            1.410, 122.0756, 120,
            1.410, 122.0756, 120
            ])

    def test_h3(self):
        self.assert_zmat_geometry('h3', [0.925, 0.925, 180])

    def test_h2h(self):
        self.assert_zmat_geometry('h2h', [0.9246, 0.9246, 180])

    def test_ch3(self):
        self.assert_zmat_geometry('ch3', [1.074, 1.074, 120, 1.074, 120, 180])

    def test_ch4(self):
        self.assert_zmat_geometry('ch4', [
            1.085,
            1.085, 109.4712206,
            1.085, 109.4712206, 120,
            1.085, 109.4712206, -120]
            )

    def test_ch3oh(self):
        self.assert_zmat_geometry('ch3oh', [
            1.408,
            1.090, 112.020,
            1.084, 106.854, 241.4-360,
            0.953, 109.807, 298.6-360,
            1.090, 112.020, 122.8
            ]
        )



    def test_ch2oh(self):
        self.assert_zmat_geometry('ch2oh', [
            1.074,
            1.078, 120.906,
            1.356, 113.472, 150.2,
            0.953, 110.744, 182.6-360
            ]
        )

    def test_c2h6(self):
        self.assert_zmat_geometry('c2h6',
           [1.522, 
            1.087, 111.295, 
            1.087, 111.295, 120,
            1.087, 111.295, -120,
            1.087, 111.295, 60,
            1.087, 111.295, 180,
            1.087, 111.295, -60
            ]
        )



    def test_c2h5(self):
        self.assert_zmat_geometry('c2h5',
            [1.485, 
            1.094, 111.744,
            1.077, 120.892, -85.534,
            1.077, 120.892, 85.534,
            1.088, 111.661, 33.871,
            1.088, 111.661, -33.871,
            ]
        )

    def test_n2h2(self):
        self.assert_zmat_geometry('n2h2', [1.023, 
            1.227, 107.343,
            1.023, 107.343, 180
            ])

    def test_TS01(self):
        self.assert_zmat_geometry('TS01', [1.371, 1.114, 180])

    def test_TS02(self):
        self.assert_zmat_geometry('TS02', [
            0.963,
            1.264, 100.451,
            0.848, 166.210, 0.0
            ]
        )

    def test_TS03(self):
        self.assert_zmat_geometry('TS03', [
            1.080,
            1.080, 114.677,
            1.080, 114.677, 135.8,
            1.381, 103.565, 247.9 - 360,
            2.742,  53.912, 112.1
            ]
        )

    def test_TS04(self):
        self.assert_zmat_geometry('TS04', [
            1.082,
            1.082, 113.227,
            1.244, 104.276, 116.7,
            1.251, 174.761, 239.5 - 360,
            0.952, 100.548, 0.0,
            1.082, 113.227, 229.4 - 360
            ]
        )

    def test_TS05(self):
        self.assert_zmat_geometry('TS05', [
            1.374,
            1.087, 115.089,
            1.081, 109.674, 230.6 - 360,
            0.958, 110.524, 314.3 - 360,
            1.330, 110.196, 116.8,
            3.032,  45.235, 243.3 - 360
            ]
        )

    def test_TS06(self):
        self.assert_zmat_geometry('TS06', [
            0.925,
            0.925, 180
            ]
        )


    def test_TS07(self):
        self.assert_zmat_geometry('TS07', [
            0.959,
            1.248, 108.803,
            1.151, 150.779, 305.5 - 360,
            1.011, 110.384, 67.3,
            1.011, 108.054, 309.5 - 360
            ]
        )


    def test_TS08(self):
        self.assert_zmat_geometry('TS08', [
            1.4265,
            1.3981, 180,
            1.0796, 101.66, 0,
            1.0796, 101.66, 120,
            1.0796, 101.66, -120
            ]
        )

    def test_TS09(self):
        self.assert_zmat_geometry('TS09', [
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

    def test_TS10(self):
        self.assert_zmat_geometry('TS10', [
            1.440,
            0.772, 180
            ]
        )

    def test_TS11(self):
        self.assert_zmat_geometry('TS11', [
            1.179,
            1.289, 360 - 181.05,
            1.082, 103.9, 180,
            1.081, 104.5, 120,
            1.081, 104.5, -120
            ]
        )

    def test_TS12(self):
        self.assert_zmat_geometry('TS12', [
            1.2587,
            1.4851, 170.61,
            1.4152, 92.92, 45,
            1.4152, 92.92, -45
            ]
        )

    def test_TS13(self):
        self.assert_zmat_geometry('TS13', [
            1.4850
            ]
        )

    def test_TS14(self):
        self.assert_zmat_geometry('TS14', [
            1.2145,
            0.8941, 180
            ]
        )

    def test_TS15(self):
        self.assert_zmat_geometry('TS15', [
            1.06,
            1.25, 105.8,
            1.17, 106.2, 0,
            1.20, 173.5, 0
            ]
        )

    def test_TS16(self):
        self.assert_zmat_geometry('TS16', [
            1.2385,
            1.4046, 175.26,
            1.3372, 95.27, 0
            ]
        )

    def test_TS17(self):
        self.assert_zmat_geometry('TS17', [
            2.44,
            2.66, 131.4
            ]
        )

    def test_TS18(self):
        self.assert_zmat_geometry('TS18', [
            1.4447,
            1.0769, 100.77,
            1.0773, 104.51, 119.42,
            1.0773, 104.51, -119.42,
            1.2091, 175.14, 0,
            1.0244, 103.66, 0
            ]
        )

    def test_TS19(self):
        self.assert_zmat_geometry('TS19', [
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

    def test_TS20(self):
        self.assert_zmat_geometry('TS20', [
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


    def test_zmat2her(self):
        zref = """O
H 1 R
H 1 R 2 A
Variables:
R=0.95
A=105
"""
        href="""BASIS
<your basis set label here>
Generated from /tmp/yo.zmat by molconvert
====================================
Atomtypes=2 Units=Angtrom Nosymmetry
Charge=1.0 Atoms=2
H   0.950000   0.000000   0.000000
H  -0.245878   0.917630   0.000000
Charge=8.0 Atoms=1
O   0.000000   0.000000   0.000000
"""
        from ..molconvert import main
        z=open('/tmp/yo.zmat', 'w'); z.write(zref); z.close()
        main('/tmp/yo.zmat', '/tmp/yo.mol')
        h=open('/tmp/yo.mol', 'r').read()
        assert href == h


if __name__ == "__main__":
    unittest.main()
