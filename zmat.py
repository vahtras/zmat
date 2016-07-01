"""Module for processing a molecule input file in ZMAT format"""

import math
from numpy import allclose
from util import full
DEBUG = False
DEG2RAD = math.pi/180
ELEMENTS = [
    "X",
    "H", "He",
    "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar"
    ]

class Atom(object):
    """Input line from ZMAT section defines an atom instance
       E I R J A K D
       E element name of current atom
       R distance relative to atom number I
       A angle relative to atoms I,J: J-I-E
       D dihedral relative to atoms I,J,K: K-J-I-E
         Data members:
             label: E (string)
             R (string)
             A
             D
             refs: refer to previous atoms, e.g. [I], [I, J], [I, J, K]
             coor: cartesian coordinates: are only set to zero here
    """
    angular = set()
    params = {None: None}
    atomlist = []


    def __init__(self, line):
        """Member parameters describe relation to other atoms
        """
        words = line.split()
       
        self.R = None
        self.A = None
        self.D = None

        self.refs = []

        if words: #what happens if empty, nothing
            self.label = words[0]

            try:
                self.R = words[2]
                self.update_params(self.R)

                self.A = words[4]
                self.update_params(self.A, angular=True)

                self.D = words[6]
                self.update_params(self.D, angular=True)
            except IndexError:
                pass
                
            self.refs = [int(i) - 1 for i in words[1::2]]
            self.coor = full.matrix(3)
            self.atomrefs = [Atom.atomlist[i] for i in self.refs]
            Atom.atomlist.append(self)

    def update_params(self, word, angular=False):
        try:
            Atom.params[word] = float(word)
            if angular:
                Atom.params[word] *= DEG2RAD
        except ValueError:
            Atom.params[word] = None
            if angular:
                Atom.angular.update([word])

    @property
    def charge(self):
        return float(ELEMENTS.index(self.label))

    @property
    def r(self):
        return Atom.params[self.R]

    @property
    def a(self):
        return Atom.params[self.A]

    @property
    def d(self):
        return Atom.params[self.D]
          

    def update_cartesian(self):
        """This aims to generate cartesian coordinates from zmat values"""
        update_cartesian(self.atomlist, self.params)

    def isbonded(self, other):
        return self.atomrefs[0] is other

def first_index(tags, lines):
    indices = [lines.index(tag) if tag in lines else len(lines) for tag in tags]
    return min(indices)

def read_params(tag, lines):
    param_values = {}
    if tag in lines:
        start = lines.index(tag) + 1
        for i in range(start, len(lines)):
            try:
                P, V = lines[i].split('=')
            except ValueError:
                break
            v = float(V)
            if P.strip() in Atom.angular:
                v *= math.pi/180
            param_values.update({P.strip(): v})
    return param_values


class Mol():
    """Molecule class holing all zmat data"""

    def __init__(self, lines):
        """An instance of class Mol is created for an input of a list of strings
        in Z-matrix format"""
        self.atomtypes = {}
        #
        zstart = 0
        zend = first_index(('Variables:', 'Constants:'), lines)

        
        self.zmat = lines[zstart:zend]
        self.atomlist = [Atom(line) for line in self.zmat]
        for atom in self.atomlist:
            if atom.label in self.atomtypes:
                self.atomtypes[atom.label].append(atom)
            else:
                self.atomtypes[atom.label] = [atom]
        #
        # Extract Parameter section
        #
        self.params = {}
        self.params.update(read_params("Variables:", lines))
        self.params.update(read_params("Constants:", lines))

    def update_cartesian(self):
        """This aims to generate cartesian coordinates from zmat values"""
        #
        # First special cases up to three atoms
        #
        update_cartesian(self.atomlist, self.params)
        return


def update_cartesian(atomlist, params):
    """This aims to generate cartesian coordinates from zmat values"""
    #
    # First special cases up to three atoms
    #
    if len(atomlist) == 1:
        pass
    elif len(atomlist) == 2:
        a, b = atomlist[:]
        assert b.refs[0] == 0
        b.coor[0] = params.get(b.R, b.r)
    elif len(atomlist) == 3:
        a, b, c = atomlist
        rb, rc = b.coor, c.coor
        #
        b.coor[0] = params.get(b.R, b.r)
        #
        if c.isbonded(a):
        #
        #      b
        #     /
        # c--a
        #
            Rca = params.get(c.R, c.r)
            cab = params.get(c.A, c.a)
            rc[0] = Rca*math.cos(cab)
            rc[1] = Rca*math.sin(cab)
        elif c.isbonded(b):
        #
        # c is bonded to b
        #
            Rcb = params.get(c.R, c.r)
            cba = params.get(c.A, c.a)
            rc[0] = rb[0] - Rcb*math.cos(cba)
            rc[1] = Rcb*math.sin(cba)
        else:
            raise SystemExit(1)
    else: 
        #
        # General case now.
        #
        # First three
        #
        a, b, c = atomlist[:3]
        b.coor[0] = params.get(b.R, b.r)
        if c.isbonded(a):
            CA = params.get(c.R, c.r)
            CAB = params.get(c.A, c.a)
            c.coor[0] = CA*math.cos(CAB)
            c.coor[1] = CA*math.sin(CAB)
        else:
            CB = params.get(c.R, c.r)
            CBA = params.get(c.A, c.a)
            c.coor[0] = b.coor[0] - CB*math.cos(CBA)
            c.coor[1] = CB*math.sin(CBA)
        for a in atomlist[3:]:
            b, c, d = a.atomrefs
            A, B, C, D = a.coor, b.coor, c.coor, d.coor
            AB = params.get(a.R, a.r)
            ABC = params.get(a.A, a.a)
            ABCD = params.get(a.D, a.d)
            if ABCD is None:
            # one known case where this fails, negated paramter
                if a.D[0] == "-":
                    ABCD = -params.get(a.D[1:], a.d)
            #
            # Translate A
            #
            # Initial setup (from origin)
            if allclose(A,  [0, 0, 0]):
                # Translate along CB
                A[:] = B + (AB/(B-C).norm2()) * (B-C)
                # Rotate in BCD plane
                n = ((D-C).cross(B-C))
                ABC0 = A.angle3(B, C)
#  rot      >>>>>>  A.rot(ABC-ABC0, n)  #rotate a around B
                A[:] = B + (A-B).rot(ABC-ABC0, n)
                # Dihedral rotation
                ABCD0 = A.dihedral(B, C, D)
                #A.rot(ABCD - ABCD0, B-C)
                A[:] = B + (A-B).rot(ABCD - ABCD0, B-C)
            else:
            # Update from preious origin)
                A[:] = B + AB*(A - B)/(A - B).norm2()
            #
            # Rotate A-B in the ABC plane:
            #
            # Current angle
            #
                ABC0 = A.angle3(B, C)
            #
            # Normal #If parallel, after initial x translation
                      #more cases?
            #
                eps = 1e-7
                if abs(ABC0) < eps:
                    # if AB andj
                    N = full.init([0., 1., 0.])
                else:
                    N = (A-B).cross(C-B)

                A.rot(ABC-ABC0, N, B)
            #
            # Current dihedral
            #

                ABCD0 = A.dihedral(B, C, D)
                print "ABCD0", ABCD0
                A.rot(ABCD - ABCD0, B-C)
            


