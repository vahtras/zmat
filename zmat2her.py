import zmat
def main(zmatfile, herfile):
    lines = [line.strip() for line in open(zmatfile, 'r').readlines()]
    mol = zmat.Mol(lines)
    mol.update_cartesian()

    h = open(herfile, 'w')
    h.write("""BASIS
STO-3G


Atomtypes=%d
"""%len(mol.atomtypes))
    for atype in mol.atomtypes.keys():
        atoms = mol.atomtypes[atype]
        charge = atoms[0].charge
        h.write("Charge=%.1f Atoms=%d\n"%(charge, len(atoms)))
        for atom in atoms:
            x, y, z = atom.coor[:]
            h.write("%s %g %g %g\n"%(atype, x, y, z))
    h.close()
    
    

if __name__ == "__main__":
    import sys
    try:
        zmatfile = sys.argv[1]
        herfile = sys.argv[2]
    except(IndexError):
        print "Usage: %s zmatfile herfile"%sys.argv[0]
        raise SystemExit
    main(sys.argv[1], sys.argv[2])
    



