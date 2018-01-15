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

indigo = Indigo()
renderer = IndigoRenderer(indigo)
bingo = Bingo.createDatabaseFile(indigo, os.path.join('tempdb'), 'molecule', '')

FP_SIZE_BYTES = 64
THRESHOLD = 0.3
MORGAN_RADIUS = 3
SIMILARITY_TYPE = "ecfp" + str(2 * MORGAN_RADIUS)

indigo.setOption("ignore-stereochemistry-errors", "1")
indigo.setOption("ignore-noncritical-query-features", "true")
indigo.setOption("ignore-bad-valence", "true")

indigo.setOption("render-output-format", "svg")
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


def write_a_match(id: str, smiles: str, similarity: float, molecule, dir_path) -> str:
    name = str(id) + ".svg"
    path = os.path.join(dir_path, name)
    renderer.renderToFile(molecule, path)
    return """
          <tr>
            <th> <img width="100" height="100" src=%s> </th>
            <th>%.3f</th> 
            <th>%s</th>
          </tr>
        """ % (name, similarity, smiles)


# # RdKit FP
# if __name__ == '__main__' and USE_RDKIT:
#     name, smiles = root_mol_smiles[0]
#     molecule = indigo.loadMolecule(smiles)
#     fingerprint_rdkit = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smiles), MORGAN_RADIUS, FP_SIZE_BYTES * 8)
#     fingerprint_str = str(base64.b16encode(base64.b64decode(fingerprint_rdkit.ToBase64())))[2:-1].zfill(FP_SIZE_BYTES * 2)
#     fingerprint = molecule.fingerprintExt(fingerprint_str, FP_SIZE_BYTES)
#     rendering = str(renderer.renderToBuffer(molecule))
#
#     print(fingerprint_str)
#
#     path = os.path.join("..", "DATA", "pubchem_slice_100k.smiles")
#     file = open(path)
#     smiles_list = [l.rstrip() for l in file][:10000]
#     molecule_list = [Chem.MolFromSmiles(l) for l in smiles_list]
#
#     for sm in smiles_list:
#         m = indigo.loadMolecule(sm)
#         m_rdkit = Chem.MolFromSmiles(sm)
#         fp_rdkit = AllChem.GetMorganFingerprintAsBitVect(m_rdkit, MORGAN_RADIUS, FP_SIZE_BYTES * 8)
#         fp_str = str(base64.b16encode(base64.b64decode(fp_rdkit.ToBase64())))[2:-1].zfill(FP_SIZE_BYTES * 2)
#         fp_ext = m.fingerprintExt(fp_str, FP_SIZE_BYTES)
#         bingo.insertWithExtFP(m, fp_ext)
#
#     iterator = bingo.searchSimWithExtFP(molecule, THRESHOLD, 1.0, fingerprint, metric='tanimoto')
#     results = []
#
#     cur_mol = iterator.getIndigoObject()
#     while iterator.next():
#         sm = cur_mol.smiles()
#         mol = Chem.MolFromSmiles(sm)
#         r = str(renderer.renderToBuffer(cur_mol))
#         id = iterator.getCurrentId()
#         sim = iterator.getCurrentSimilarityValue()
#
#         results += [(id, sm, sim)]
#     iterator.close()
#
#     sorted_results = sorted(results, key=lambda result: result[2], reverse=True)
#
#     for res in sorted_results:
#         print("%.3f" % res[2])
#         print(res[1])
#

# Indigo FP
if __name__ == '__main__':
    name, smiles = root_mol_smiles[0]
    molecule = indigo.loadMolecule(smiles)
    fingerprint = molecule.fingerprint("sim")

    path = os.path.join("..", "DATA", "chembl_23.sdf")

    for i, m in enumerate(indigo.iterateSDFile(path)):
        bingo.insert(m)
        if i % 10000 == 0:
            print("%d molecule loaded" % i)

    iterator = bingo.searchSim(molecule, THRESHOLD, 1.0, metric='tanimoto')
    results = []

    cur_mol = iterator.getIndigoObject()
    while iterator.next():
        sm = cur_mol.smiles()
        mol = indigo.loadMolecule(sm)
        r = str(renderer.renderToBuffer(mol))
        sim = iterator.getCurrentSimilarityValue()
        id = iterator.getCurrentId()

        results += [(id, sm, sim)]
    iterator.close()

    sorted_results = sorted(results, key=lambda result: result[2], reverse=True)

    report_dir_path = os.path.join("reports", SIMILARITY_TYPE + '_' + name + '_' + str(THRESHOLD))

    if not os.path.exists(report_dir_path):
        os.makedirs(report_dir_path)

    report = open(os.path.join(report_dir_path, "report.html"), "w")

    report.write("""
    <!DOCTYPE html>
    <html>
    <head>
    <style>
        table, th, td {
            border: 1px solid black;
            border-collapse: collapse;
        }
    </style>
    </head>
    <body>
    
    <table style="width:100%">
      <tr>
        <th>Molecule</th>
        <th>Similarity</th> 
        <th>Smiles</th>
      </tr>
    """)

    report.write(write_a_match(name, smiles, 1.0, molecule, report_dir_path))
    for res in sorted_results:
        id, smiles, similarity = res
        print("%.3f" % similarity)
        print(smiles)
        report.write(write_a_match(id, smiles, similarity, bingo.getRecordById(id), report_dir_path))

    report.write("""
    </table>

    </body>
    </html>
    """)
    report.close()



