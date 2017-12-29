from indigo import Indigo


def fp_density(fp: str) -> float:
    ones = [c for c in fp if c == '1']
    return len(ones) / len(fp)


if __name__ == '__main__':
    smiles = [
            'CC(=O)OC1=CC=CC=C1C(=O)O',
            'CCCC1=NN(C2=C1NC(=NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C',
            'CC1=C(SC2=C1C(=NC(C3=NN=C(N32)C)CC(=O)OC(C)(C)C)C4=CC=C(C=C4)Cl)C',
            'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
            'CN1CCCC1C2=CN=CC=C2',
            'C1CN2C(=NN=C2C(F)(F)F)CN1C(=O)C[C@@H](CC3=CC(=C(C=C3F)F)F)N',
            'CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5',
            'C1CC(=O)NC2=C1C=CC(=C2)OCCCCN3CCN(CC3)C4=C(C(=CC=C4)Cl)Cl',
            'CNC1(CCCCC1=O)C2=CC=CC=C2Cl',
            'CN1C2CCC1C(C(C2)OC(=O)C3=CC=CC=C3)C(=O)OC',
            'CCCC1=NN(C2=C1NC(=NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C'
          ]

    for i, sml in enumerate(smiles):
        print("#%2d: %s" % (i, sml))

    df_descriptors = []
    indigo = Indigo()

    # indigo.setOption("similarity-type", "sim")
    # indigo.setOption("similarity-type", "chem")
    # indigo.setOption("similarity-type", "ECFP2")
    indigo.setOption("similarity-type", "ECFP4")
    # indigo.setOption("similarity-type", "ECFP6")
    # indigo.setOption("similarity-type", "ECFP8")

    fingerprints = []
    for sml in smiles:
        mol = indigo.loadMolecule(sml)
        fp = mol.fingerprint("sim")
        fingerprints += [fp]

    print("*** Fingerprints ***")

    for i, fp in enumerate(fingerprints):
        str = fp.toString()
        density = fp_density(str)
        print("#%2d: Density: %f ; FP: %s" % (i, density, str))

    print("*** Similarity matrix ***")

    for f0 in fingerprints:
        for f1 in fingerprints:
            similarity = indigo.similarity(f0, f1)
            print("   %.3f" % similarity, end="")
        print("\n")
