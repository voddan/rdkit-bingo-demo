import base64
from itertools import islice, starmap

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import Descriptors

import math
from indigo.indigo import *
from indigo.indigo_renderer import *
from indigo.bingo import *

FP_SIZE_BYTES = 64
THRESHOLD = 0.3
MORGAN_RADIUS = 3
SIMILARITY_TYPE = "ecfp" + str(2 * MORGAN_RADIUS)
IMPLEMENTATION = "rdkit"

indigo = Indigo()

indigo.setOption("ignore-stereochemistry-errors", "1")
indigo.setOption("ignore-noncritical-query-features", "true")
indigo.setOption("ignore-bad-valence", "true")
indigo.setOption("ignore-stereochemistry-errors", "true")
indigo.setOption("deco-ignore-errors", "true")

indigo.setOption("similarity-type", SIMILARITY_TYPE)
indigo.setOption("fp-sim-qwords", FP_SIZE_BYTES / 8)
indigo.setOption("fp-ext-enabled", True)
indigo.setOption("fp-ord-qwords", 0)  # optional
indigo.setOption("fp-tau-qwords", 0)  # optional
indigo.setOption("fp-any-qwords", 0)  # optional

root_mol_smiles = [
        ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        ("CETRIMONIUM", "CCCCCCCCCCCCCCCC[N+](C)(C)C"),
        ("Gleevec", "CN1CCN(Cc2ccc(cc2)C(=O)Nc3ccc(C)c(Nc4nccc(n4)c5cccnc5)c3)CC1.CS(=O)(=O)O"),
        ("Naphthalene", "c1ccc2ccccc2c1"),
        ("JQ1", "Cc1sc2c(C(=N[C@@H](CC(=O)OC(C)(C)C)c3nnc(C)n23)c4ccc(Cl)cc4)c1C"),
    ]


if __name__ == '__main__':
    db_path = os.path.join('DB', 'chembl_23.sdf', IMPLEMENTATION, SIMILARITY_TYPE)

    if not os.path.exists(db_path):
        os.makedirs(db_path)

    bingo = Bingo.createDatabaseFile(indigo, db_path, 'molecule', '')

    name, smiles = root_mol_smiles[0]
    molecule = indigo.loadMolecule(smiles)
    fingerprint = molecule.fingerprint("sim")

    path = os.path.join("..", "DATA", "chembl_23.sdf")

    for i, m in enumerate(indigo.iterateSDFile(path)):
        try:
            fingerprint_rdkit = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smiles), MORGAN_RADIUS, FP_SIZE_BYTES * 8)
            fingerprint_str = str(base64.b16encode(base64.b64decode(fingerprint_rdkit.ToBase64())))[2:-1].zfill(FP_SIZE_BYTES * 2)
            fingerprint = m.fingerprintExt(fingerprint_str, FP_SIZE_BYTES)
        except Exception as e:
            print(e)
            continue

        bingo.insert(m)

        if i % 10000 == 0:
            print("%d molecule loaded" % i)
