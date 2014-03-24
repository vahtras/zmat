#!/usr/bin/env python
"""Script to convert from zmat format to other formats
 
    Supported output formats:
        xyz: standard format
        pot: as Dalton potential input file
        mol: as Dalton molecule input file
"""

import argparse
import os
import zmat

class UnsupportedMoltype(Exception): 
    """Basic exception for unregognized mol types"""
    pass
# make the class return error output?

FROMS = ['.zmat', ]
TOS = ['.xyz', '.mol', '.pot']
ALPHA_ALKENE = {"C":0.878, "H": 0.135, "O": 0.465}


def main(from_mol, to_mol, **kwargs):
    """Convert molecule file from_mol to to_mol

    from_mol -- filename of exisiting molecule 
    to_mol -- filename of generated file in new format

    The in- and out formats are extracted from extensions
    Currently these are
    in: zmat
    out: xyz, pot, mol
    """
    #
    # Determine type of from file by extension
    #
    _, from_ext  = os.path.splitext(from_mol)
    if from_ext not in FROMS: 
        raise(UnsupportedMoltype)

    #
    # Determine output file by extension
    #
    _, to_ext = os.path.splitext(to_mol)
    if to_ext not in TOS: 
        print to_ext, TOS
        raise(UnsupportedMoltype)


    #Process input file
    if from_ext == '.zmat':
        lines = [line.strip() for line in open(from_mol).readlines()]
        mol = zmat.Mol(lines)
        mol.update_cartesian()
        #
        # Filer out fake centers X in zmat input
        alist = [atom for atom in mol.atomlist if atom.label != "X"]

        #Open output flie
        with open(to_mol, 'w') as outf:
            if to_ext  == '.xyz':
                outf.write(
                    "%d\nGenerated from %s by molconvert\n" % (
                        len(alist), from_mol
                        )
                    )
                for atom in alist:
                    x, y, z = atom.coor[:]
                    outf.write(
                        "%s %10.6f %10.6f %10.6f\n" % (atom.label, x, y, z)
                        )
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
""" % (from_mol, len(mol.atomtypes))
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
    parser.add_argument('from_mol')
    parser.add_argument('to_mol')
    parser.add_argument('--alpha_C', default=ALPHA_ALKENE["C"])
    parser.add_argument('--alpha_H', default=ALPHA_ALKENE["H"])
    parser.add_argument('--alpha_O', default=ALPHA_ALKENE["O"])
    pargs = parser.parse_args()
    # defalt
    main(
        pargs.from_mol, pargs.to_mol, 
        C=pargs.alpha_C, H=pargs.alpha_H, O=pargs.alpha_O
        )
       
