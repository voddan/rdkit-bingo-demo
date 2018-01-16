from indigo.bingo import *

FP_SIZE_BYTES = 64
MORGAN_RADIUS = 3
SIMILARITY_TYPE = "ecfp" + str(2 * MORGAN_RADIUS)

indigo = Indigo()

indigo.setOption("ignore-stereochemistry-errors", "1")
indigo.setOption("ignore-noncritical-query-features", "true")
indigo.setOption("ignore-bad-valence", "true")
indigo.setOption("ignore-stereochemistry-errors", "true")
indigo.setOption("deco-ignore-errors", "true")

indigo.setOption("similarity-type", SIMILARITY_TYPE)
indigo.setOption("fp-sim-qwords", FP_SIZE_BYTES / 8)
indigo.setOption("fp-ext-enabled", False)
indigo.setOption("fp-ord-qwords", 0)  # optional
indigo.setOption("fp-tau-qwords", 0)  # optional
indigo.setOption("fp-any-qwords", 0)  # optional


if __name__ == '__main__':
    path = os.path.join("..", "DATA", "chembl_23.sdf")

    all_zeros = "0" * (FP_SIZE_BYTES * 2)

    for i, molecule in enumerate(indigo.iterateSDFile(path)):
        fingerprint = molecule.fingerprint("sim")
        str = fingerprint.toString()

        if str == all_zeros:
            print("#%8d: %s" % (i, molecule.smiles()))

        if i % 10000 == 0:
            print("%d molecules have been scanned" % i)
