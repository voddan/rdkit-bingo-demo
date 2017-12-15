from rdkit import Chem
from rdkit.Chem import Descriptors


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


def pack_normalized_descriptors_to_fingerprint(descriptors: [float], byte_size=64, density=0.3) -> str:
    """
    :param descriptors: list of numbers (roughly) between 0.0 and 1.0 that characterise a molecule
    :param byte_size: size of the fingerprint in bytes
    :param density: approximate density of '1's in the fingerprint
    :return: fingerprint as a hex string
    """
    length = len(descriptors)
    bit_size = byte_size * 8

    fingerprint = [False] * bit_size

    for id in range(length):
        # The magic formula. Works best for large bit_size. Will be adjusted later
        set_bits_num = int(descriptors[id] * (density * 10) * bit_size / length)

        hash = id
        for cnt in range(set_bits_num):
            hash = (hash * 0x8088405 + 1) % bit_size
            fingerprint[hash] = True

    str = ""
    for i in range(0, bit_size, 4):
        (a, b, c, d) = fingerprint[i: i + 4]
        digit = a + 2 * b + 4 * c + 8 * d
        str += hex(digit)[2:]

    return str


def fp_density(fp: str) -> float:
    ones = [c for c in fp if c == '1']
    return len(ones) / len(fp)


# %% example of calling custom fingerprint routine:
if __name__ == '__main__':
    import pandas as pd

    df_in = pd.read_csv("smiles_input_table.csv")
    df_descriptors = get_custom_descriptors(df_in)

    print(df_descriptors)

    descriptor_lists = [list(df_descriptors.ix[i])[2:] for i in range(len(df_descriptors))]

    print("Descriptors:")
    for list in descriptor_lists:
        print(list)

    max_bounds = [20, 1, 20, 1, 1, 500, 500, 500, 200, 1, 1, -1, 1, 1, 10, 10, 10, 10, 1350, 40, 40, 40, 20, 20, 20, 10, 20, 10, 10, 10, 10, -1, 321337165, 40, 20, 10, 500, 20, 40, 20, 20, 10, 20, 20, 20, 20, 40, 40, 100, 100, 40, 40, 100, 1, 40, 40, 100, 100, 100, 1, 20, 20, 40, 10, 40, 100, 40, 40, 100, 100, 20, 40, 1, 200, 100, 40, 1, 100, 100, 100, 40, 40, 100, 100, 40, 1, 20, 1, 1, 1, 1, 1, 1, 100, 100, 1, 40, 10, 10, 1, 10, 10, 10, 10, 10, 10, 10, 20, 10, 1, 10, 10, 10, 10, 200]
    min_bounds = [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -10, 0]

    print("Packed fingerprints:")
    for descrs in descriptor_lists:
        normalized = [(f - min_bounds[i]) / (max_bounds[i] - min_bounds[i]) for i, f in enumerate(descrs)]
        fp = pack_normalized_descriptors_to_fingerprint(normalized, density=0.3, byte_size=64)
        print("Density: %f ; FP: %s" % (fp_density(fp), fp))
