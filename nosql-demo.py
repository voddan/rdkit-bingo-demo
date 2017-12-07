from rdkit import Chem
from rdkit.Chem import Descriptors
from os import path


def get_custom_descriptors(df_in):
    """
    # USAGE: /given pd.DataFrame with columns/header ID, Structure <[PossibleOtherFileds]>/
    # df_fps = get_custom_fps( df_in )
    """
    import pandas as pd

    FP_names = [x[0] for x in Descriptors.descList]
    FP_funcs = [x[1] for x in Descriptors.descList]

    df1 = df_in  # df_in.copy()
    df2 = pd.DataFrame(index=range(df_in['Structure'].shape[0]), columns=FP_names)

    df_descriptors = pd.concat([df1, df2], axis=1, join='inner')

    print('Calculating Molecule Descriptors...')
    for idx, smile in enumerate(df_descriptors['Structure']):
        print("   Calculating descriptors...")
        m = Chem.MolFromSmiles(smile)
        fp4 = [f(m) for f in FP_funcs]
        df_descriptors.ix[idx, FP_names] = fp4

    print("Custom descriptors were successfully generated")
    return df_descriptors


def pack_normalized_descriptors_to_fingerprint(descriptors: [float], bit_size=512, density=0.3) -> int:
    """
    :param descriptors: list of numbers between 0.0 and 1.0 that characterise a molecule
    :param bit_size: size of the fingerprint in bits
    :param density: approximate density of '1's in the fingerprint
    :return: fingerprint as a number
    """
    length = len(descriptors)

    fingerprint = [False] * bit_size

    for id in range(length):
        # The magic formula. Works best for large bit_size. Will be adjusted later
        set_bits_num = int(descriptors[id] * (density * 10) * bit_size / length)

        hash = id
        for cnt in range(set_bits_num):
            hash = (hash * 0x8088405 + 1) % bit_size
            fingerprint[hash] = True

    fp_binary = 0
    for i in range(bit_size):
        fp_binary *= 2
        fp_binary += int(fingerprint[i])

    return fp_binary


def fp_density(fp: int) -> float:
    str = bin(fp)
    ones = [c for c in str if c == '1']
    return len(ones) / len(str)


# %% example of calling custom fingerprint routine:
if __name__ == '__main__':
    import pandas as pd

    df_in = pd.read_csv("smiles_input_table.csv")
    df_descriptors = get_custom_descriptors(df_in)

    print(df_descriptors)

    descriptor_lists = [list(df_descriptors.ix[i])[2:] for i in range(len(df_descriptors))]

    max_bounds = [20, 1, 20, 1, 1, 500, 500, 500, 200, 1, 1, -1, 1, 1, 10, 10, 10, 10, 1350, 40, 40, 40, 20, 20, 20, 10, 20, 10, 10, 10, 10, -1, 321337165, 40, 20, 10, 500, 20, 40, 20, 20, 10, 20, 20, 20, 20, 40, 40, 100, 100, 40, 40, 100, 1, 40, 40, 100, 100, 100, 1, 20, 20, 40, 10, 40, 100, 40, 40, 100, 100, 20, 40, 1, 200, 100, 40, 1, 100, 100, 100, 40, 40, 100, 100, 40, 1, 20, 1, 1, 1, 1, 1, 1, 100, 100, 1, 40, 10, 10, 1, 10, 10, 10, 10, 10, 10, 10, 20, 10, 1, 10, 10, 10, 10, 200]
    min_bounds = [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -10, 0]

    print("Add packed fingerprints to Indigo:")
    fp_size_bytes = 64  # or 512 bits, or 8 qwords

    from indigo.indigo import Indigo
    from indigo.bingo import Bingo
    indigo = Indigo()
    bingo = Bingo.createDatabaseFile(indigo, path.join('tempdb'), 'molecule', '')
    indigo.setOption("fp-sim-qwords", fp_size_bytes / 8)
    indigo.setOption("fp-ord-qwords", 0)
    indigo.setOption("fp-tau-qwords", 0)
    indigo.setOption("fp-any-qwords", 0)
    indigo.setOption("fp-ext-enabled", True)

    molecules = []
    for id in range(len(df_descriptors)):
        smiles = df_descriptors.ix[id, 'Structure']
        descriptors = list(df_descriptors.ix[id])[2:]
        normalized = [(f - min_bounds[i]) / (max_bounds[i] - min_bounds[i]) for i, f in enumerate(descriptors)]
        fp = pack_normalized_descriptors_to_fingerprint(normalized, density=0.3, bit_size=fp_size_bytes * 8)

        mol = indigo.loadMolecule(smiles)
        molecules += [mol]
        ext_fp = mol.fingerprintExt(hex(fp)[2:], fp_size_bytes)
        bingo.insertWithExtFP(mol, ext_fp)

    print("Similarity matrix:")
    for m1 in molecules:
        for m2 in molecules:
            similarity = indigo.similarity(m1, m2)
            print("%.4f  " % similarity, end="")
        print()
