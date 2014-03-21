"""Script to convert from zmat format to other formats
 
    Supported output formats:
        xyz: standard format
        pot: as Dalton potential input file
        mol: as Dalton molecule input file
"""

import argparse
import os
import zmat

class UnsupportedMoltype(Exception): pass
# make the class return error output?

FROMS = ['.zmat', ]
TOS = ['.xyz', '.mol', '.pot']


def main(*args, **kwargs):
    """Convert molecule file foo to bar"""
    foo, bar = args
    #
    # Determine type of from file by extension
    #
    _, from_ext  = os.path.splitext(foo)
    if from_ext not in FROMS: raise(UnsupportedMoltype)

    #
    # Determine output file by extension
    #
    _, to_ext = os.path.splitext(bar)
    if to_ext not in TOS: 
        print to_ext, TOS
        raise(UnsupportedMoltype)


    #Process input file
    if from_ext == '.zmat':
        lines = [line.strip() for line in open(foo).readlines()]
        mol = zmat.Mol(lines)
        mol.update_cartesian()
        #
        # Filer out fake centers X in zmat input
        alist = [atom for atom in mol.atomlist if atom.label != "X"]

        #Open output flie
        with open(bar, 'w') as outf:
            if to_ext  == '.xyz':
                outf.write("%d\nGenerated from %s by molconvert\n" % (len(alist), foo, ))
                for atom in alist:
                    x, y, z = atom.coor[:]
                    outf.write("%s %10.6f %10.6f %10.6f\n"%(atom.label, x, y, z))
                    print("%s %10.6f %10.6f %10.6f\n"%(atom.label, x, y, z))
            elif to_ext == '.pot':
                outf.write("AA\n%d 1 0 1\n" % (len(alist),))
                for atom in alist:
                    x, y, z = atom.coor[:]
                    outf.write(
                        "1 %10.6f %10.6f %10.6f 0 0 0 0 %g\n" % 
                        (x, y, z, kwargs[atom.label])
                    )
            elif to_ext == '.mol':
                outf.write("""BASIS
<your basis set label here>
Generated from %s by molconvert
====================================
Atomtypes=%d Units=Angtrom Nosymmetry
""" % (foo, len(mol.atomtypes))
                    )
                for atype in mol.atomtypes.keys():
                    atoms = mol.atomtypes[atype]
                    charge = atoms[0].charge
                    outf.write("Charge=%.1f Atoms=%d\n"%(charge, len(atoms)))
                    for atom in atoms:
                        x, y, z = atom.coor[:]
                        outf.write("%s %10.6f %10.6f %10.6f\n"%(atype, x, y, z))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('foo')
    parser.add_argument('bar')
    args = parser.parse_args()
    # defalt
    kwargs = {"C":0.878, "H": 0.135, "O": 0.465}
    main(args.foo, args.bar, **kwargs) 
