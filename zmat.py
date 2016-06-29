"""Module for processing a molecule input file in ZMAT format"""

import math
from numpy import allclose
from util import full
DEBUG = False
DEG2RAD = math.pi/180

class Atom:
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
             charge: charge of element E (float)
             coor: cartesian coordinates: are only set to zero here
    """
    angular = set()
    params = {None: None}
    atomlist = []
    charge = ["X",
         "H", "He",
         "Li", "Be", "B", "C", "N", "O", "F", "Ne",
         "Na", "Mg","Al","Si", "P", "S","Cl", "Ar"]


    def __init__(self, line):
        """Member parameters describe relation to other atoms
        """
        words = line.split()
        lw = len(words)
       
        self.R = None
        self.A = None
        self.D = None

        self.refs = []

        if words: #what happens if empty, nothing
            self.label = words[0]
            self.charge = float(Atom.charge.index(self.label))

            if lw > 2: 
                self.R = words[2]
                try:
                    Atom.params[self.R] = float(self.R)
                except ValueError:
                    Atom.params[self.R] = None
                
            if lw > 4: 
                self.A = words[4]
                try:
                    Atom.params[self.A] = float(self.A)*DEG2RAD
                except ValueError:
                    Atom.params[self.A] = None
                    Atom.angular.update([self.A])

            if lw > 6: 
                self.D = words[6]
                try:
                    Atom.params[self.D] = float(self.D)*DEG2RAD
                except ValueError:
                    Atom.params[self.D] = None
                    Atom.angular.update([self.D])

            self.refs = [int(i) - 1 for i in words[1:lw:2]]
            self.coor = full.matrix(3)
            self.atomrefs = [Atom.atomlist[i] for i in self.refs]
            Atom.atomlist.append(self)

    @property
    def r(self):
        return Atom.params[self.R]

    @property
    def a(self):
        return Atom.params[self.A]

    @property
    def d(self):
        return Atom.params[self.D]
          
    def __str__(self):
        """Return atom line in mol style"""
        retstr = "      %4s    1\n" % self.charge
        #retstr = ""
        retstr += "%-2s" % self.label
        retstr += "%20.10f %20.10f %20.10f\n" % (
            self.coor[0], self.coor[1], self.coor[2]
            )
        if DEBUG:
            print self.label, self.charge, self.R, self.A, self.D, self.refs
        return retstr

    def update_cartesian(self):
        """This aims to generate cartesian coordinates from zmat values"""
        #
        # First special cases up to three atoms
        #
        if len(self.atomlist) == 1:
            return 0
        elif len(self.atomlist) == 2:
            a, b = self.atomlist[:]
            rb = b.coor
            
            assert b.refs[0] == 0
            #
            # Retrievs BA parametrized distance if defined 
            # the constant value in the ZMAT input
            #
            rb[0] = self.params.get(b.R, b.r)
        elif len(self.atomlist) == 3:
            a, b, c = self.atomlist
            rb, rc = b.coor, c.coor
            #
            rb[0] = self.params.get(b.R, b.r)
            #
            if c.refs[0] == 0: #want: is a
            #
            # c is bonded to a
            #
                Rca = self.params.get(c.R, c.r)
                cab = self.params.get(c.A, c.a)
                rc[0] = Rca*math.cos(cab)
                rc[1] = Rca*math.sin(cab)
            elif c.refs[0] == 1:
            #
            # c is bonded to b
            #
                Rcb = self.params.get(c.R, c.r)
                cba = self.params.get(c.A, c.a)
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
            a, b, c = self.atomlist[:3]
            b.coor[0] = self.params.get(b.R, b.r)
            if c.refs[0] == 0:
                CA = self.params.get(c.R, c.r)
                CAB = self.params.get(c.A, c.a)
                c.coor[0] = CA*math.cos(CAB)
                c.coor[1] = CA*math.sin(CAB)
            else:
                CB = self.params.get(c.R, c.r)
                CBA = self.params.get(c.A, c.a)
                c.coor[0] = b.coor[0] - CB*math.cos(CBA)
                c.coor[1] = CB*math.sin(CBA)
            for a in self.atomlist[3:]:
                #b, c, d = [self.atomlist[i] for i in a.refs] BUG
                b = self.atomlist[a.refs[0]]
                c = self.atomlist[a.refs[1]]
                d = self.atomlist[a.refs[2]]
                A, B, C, D = a.coor, b.coor, c.coor, d.coor
                AB = self.params.get(a.R, a.r)
                ABC = self.params.get(a.A, a.a)
                ABCD = self.params.get(a.D, a.d)
                if ABCD is None:
                # one known case where this fails, negated paramter
                    if a.D[0] == "-":
                        ABCD = -self.params.get(a.D[1:], a.d)
                #print "AB", AB
                #print "ABC", ABC
                #print "ABCD", ABCD
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

class Mol():
    """Molecule class holing all zmat data"""

    def __init__(self, lines):
        """An instance of class Mol is created for an input of a list of strings
        in Z-matrix format"""
        self.atomtypes = {}
        #
        zstart = 0
        if 'Variables:' in lines:
            zend = lines.index('Variables:')
        elif 'Constants:' in lines:
            zend = lines.index('Constants:')
        else:
            zend = len(lines)

        
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
        pstart = zend + 1
        for pline in lines[pstart:]:
            P, V = pline.split('=')
            P = P.strip()
            if P in Atom.angular:
                self.params[P] = float(V)*math.pi/180
            else:
                self.params[P] = float(V)

    def __str__(self):
        retstr = "\n".join(self.zmat) + "\n\n"
        for (k, v) in zip(self.params.keys(), self.params.values()):
            if k in Atom.angular:
                retstr += "%s = %f\n" % (k, v*180/math.pi)
            else:
                retstr += "%s = %f\n" % (k, v)
        #for a in self.atomlist: retstr += str(a)
        return retstr


    def update_cartesian(self):
        """This aims to generate cartesian coordinates from zmat values"""
        #
        # First special cases up to three atoms
        #
        if len(self.atomlist) == 1:
            return 0
        elif len(self.atomlist) == 2:
            a, b = self.atomlist[:]
            rb = b.coor
            
            assert b.refs[0] == 0
            #
            # Retrievs BA parametrized distance if defined 
            # the constant value in the ZMAT input
            #
            rb[0] = self.params.get(b.R, b.r)
        elif len(self.atomlist) == 3:
            a, b, c = self.atomlist
            rb, rc = b.coor, c.coor
            #
            rb[0] = self.params.get(b.R, b.r)
            #
            if c.refs[0] == 0: #want: is a
            #
            # c is bonded to a
            #
                Rca = self.params.get(c.R, c.r)
                cab = self.params.get(c.A, c.a)
                rc[0] = Rca*math.cos(cab)
                rc[1] = Rca*math.sin(cab)
            elif c.refs[0] == 1:
            #
            # c is bonded to b
            #
                Rcb = self.params.get(c.R, c.r)
                cba = self.params.get(c.A, c.a)
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
            a, b, c = self.atomlist[:3]
            b.coor[0] = self.params.get(b.R, b.r)
            if c.refs[0] == 0:
                CA = self.params.get(c.R, c.r)
                CAB = self.params.get(c.A, c.a)
                c.coor[0] = CA*math.cos(CAB)
                c.coor[1] = CA*math.sin(CAB)
            else:
                CB = self.params.get(c.R, c.r)
                CBA = self.params.get(c.A, c.a)
                c.coor[0] = b.coor[0] - CB*math.cos(CBA)
                c.coor[1] = CB*math.sin(CBA)
            for a in self.atomlist[3:]:
                #b, c, d = [self.atomlist[i] for i in a.refs] BUG
                b = self.atomlist[a.refs[0]]
                c = self.atomlist[a.refs[1]]
                d = self.atomlist[a.refs[2]]
                A, B, C, D = a.coor, b.coor, c.coor, d.coor
                AB = self.params.get(a.R, a.r)
                ABC = self.params.get(a.A, a.a)
                ABCD = self.params.get(a.D, a.d)
                if ABCD is None:
                # one known case where this fails, negated paramter
                    if a.D[0] == "-":
                        ABCD = -self.params.get(a.D[1:], a.d)
                #print "AB", AB
                #print "ABC", ABC
                #print "ABCD", ABCD
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
                


