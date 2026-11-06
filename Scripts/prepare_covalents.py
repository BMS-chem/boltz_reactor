"""
prepare_covalents.py

Turns each covalent ligand into the "bound" form and works out the atom
attachment for Boltz.

A covalent ligand is given as its normal, unreacted SMILES plus a WARHEAD name.
This script applies that warhead's reaction, which opens one valence on the
attachment atom, and names that atom the way Boltz will name it.

NB: can't mix covalent and non-covalent ligands in the same experiment

Reads:  experiments/<experiment>/inputs/ligands.csv    (ID, SMILES, WARHEAD)
        experiments/<experiment>/inputs/proteins.csv   (ID, SEQUENCE, COFACTORS, COVALENT_SITE)
Writes: experiments/<experiment>/inputs/covalent_prepared.csv

Example:
    python Scripts/prepare_covalents.py --experiment covalent_test
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem

# RDKit prints its own warnings to stderr; the checks below report problems properly.
RDLogger.DisableLog("rdApp.*")


# ---------------------------------------------------------------------------
# 1. The warhead library
# ---------------------------------------------------------------------------
# Warheads are hardcoded here.
#
# The product keeps a filled valence in the SMARTS; the code below opens it by
# putting one radical electron on the mapped atom. That open-valence form is
# what Boltz takes, and the bond constraint joins it to the protein atom.
WARHEADS = {
    "acrylamide": {
        # Michael addition
        "smarts": "[CH2:1]=[CH:2][CX3:3](=[O:4])[NX3:5]>>[CH2:1][CH2:2][CX3:3](=[O:4])[NX3:5]",
        "residue": "C",
        "atom": "SG",
    },
    "chloroacetamide": {
        # SN2
        "smarts": "[Cl][CH2:1][CX3:2](=[O:3])[NX3:4]>>[CH2:1][CX3:2](=[O:3])[NX3:4]",
        "residue": "C",
        "atom": "SG",
    },
    "CF2_acrylamide": {
        # Michael addition
        "smarts": "[F][CH:1]([F])[CH:2]=[CH:3][CX3:4](=[O:5])[NX3:6]>>[F][CH:1]([F])[CH:2][CH2:3][CX3:4](=[O:5])[NX3:6]",
        "residue": "C",
        "atom": "SG",
    }
}


# ---------------------------------------------------------------------------
# 2. Read the command-line arguments
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Prepare covalent ligands for Boltz.")
parser.add_argument("--experiment", required=True, help="Name of the folder inside experiments/")
args = parser.parse_args()


# ---------------------------------------------------------------------------
# 3. Work out the folder paths
# ---------------------------------------------------------------------------
# This script lives in <project>/Scripts/, so the project folder is one level up.
project_dir = Path(__file__).resolve().parent.parent
experiment_dir = project_dir / "experiments" / args.experiment
inputs_dir = experiment_dir / "inputs"
prepared_file = inputs_dir / "covalent_prepared.csv"

if not experiment_dir.exists():
    sys.exit(f"ERROR: experiment folder not found: {experiment_dir}")


# ---------------------------------------------------------------------------
# 4. Read the two CSV files
# ---------------------------------------------------------------------------
# encoding="utf-8-sig" removes the invisible marker some text editors put at the start.
# dtype=str and keep_default_na=False keep every cell as text, with empty cells as "".
proteins = pd.read_csv(inputs_dir / "proteins.csv", encoding="utf-8-sig", dtype=str, keep_default_na=False)
ligands = pd.read_csv(inputs_dir / "ligands.csv", encoding="utf-8-sig", dtype=str, keep_default_na=False)

proteins.columns = proteins.columns.str.strip()
ligands.columns = ligands.columns.str.strip()
proteins = proteins.apply(lambda column: column.str.strip())
ligands = ligands.apply(lambda column: column.str.strip())


# ---------------------------------------------------------------------------
# 5. Check whether experiment is covalent
# ---------------------------------------------------------------------------
# (if any values are present in WARHEAD column)
if "WARHEAD" not in ligands.columns:
    warheads_named = []
else:
    warheads_named = [w for w in ligands["WARHEAD"] if w != ""]

if not warheads_named:
    print("Non-covalent experiment.")
    if prepared_file.exists():
        prepared_file.unlink()
        print(f"Removed a leftover {prepared_file.name} from an earlier run.")
    sys.exit(0)

if len(warheads_named) != len(ligands):
    sys.exit(f"ERROR: {len(warheads_named)} of {len(ligands)} ligands have a WARHEAD.\n"
             "  An experiment is either covalent or it is not; they cannot be mixed.\n"
             "  (Boltz cannot mix affinity and non-affinity jobs in one folder - see PLAN.md.)\n"
             "  Either give every ligand a WARHEAD, or split them into two experiments.")

print(f"Covalent experiment: {len(ligands)} ligands, {len(proteins)} proteins")


# ---------------------------------------------------------------------------
# 6. Check the warhead names and the protein attachment sites
# ---------------------------------------------------------------------------
# Every problem found is collected, so they can all be fixed in one go.
problems = []

# 6a. Check that warhead is encoded
for _, ligand in ligands.iterrows():
    if ligand["WARHEAD"] not in WARHEADS:
        problems.append(f"ligands.csv: '{ligand['ID']}' has an unknown WARHEAD "
                        f"'{ligand['WARHEAD']}'. Known: {sorted(WARHEADS)}")

if problems:
    sys.exit("ERROR:\n  " + "\n  ".join(problems))


# 6b. Check that protein site is correctly specified
if "COVALENT_SITE" not in proteins.columns:
    sys.exit("ERROR: this is a covalent experiment, but proteins.csv has no "
             "'COVALENT_SITE' column.\n"
             "  Add one, in the form <chain letter>:<residue letter><position>, e.g. A:C16.")

chain_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

for _, protein in proteins.iterrows():
    site = protein["COVALENT_SITE"]
    protein_id = protein["ID"]

    if site == "":
        problems.append(f"proteins.csv: '{protein_id}' has an empty COVALENT_SITE")
        continue

    # Check that the covalent site is in the expected format
    if ":" not in site:
        problems.append(f"proteins.csv: '{protein_id}' COVALENT_SITE '{site}' has no ':'. "
                        "Expected <chain>:<residue><position>, e.g. A:C16")
        continue

    chain_letter, residue_part = site.split(":", 1)
    chain_letter = chain_letter.upper()
    residue_part = residue_part.upper()

    if len(residue_part) < 2 or not residue_part[0].isalpha() or not residue_part[1:].isdigit():
        problems.append(f"proteins.csv: '{protein_id}' COVALENT_SITE '{site}' should be a "
                        "residue letter then a position, e.g. A:C16")
        continue

    residue_letter = residue_part[0]
    position = int(residue_part[1:])

    # Check the covalent attachment site is present and correct in the protein
    chains = protein["SEQUENCE"].split(":")
    if chain_letter not in chain_letters[:len(chains)]:
        problems.append(f"proteins.csv: '{protein_id}' COVALENT_SITE '{site}' names chain "
                        f"'{chain_letter}', but this protein has {len(chains)} chain(s) "
                        f"({chain_letters[:len(chains)]})")
        continue

    chain_sequence = chains[chain_letters.index(chain_letter)].upper()

    if not 1 <= position <= len(chain_sequence):
        problems.append(f"proteins.csv: '{protein_id}' COVALENT_SITE '{site}' position {position} "
                        f"is outside chain {chain_letter}, which is {len(chain_sequence)} residues long")
        continue

    # The provided residue has to match the residue in the sequence
    actual = chain_sequence[position - 1]
    if actual != residue_letter:
        problems.append(f"proteins.csv: '{protein_id}' COVALENT_SITE '{site}' says residue "
                        f"{position} is '{residue_letter}', but the sequence has '{actual}' there")
        continue

if problems:
    sys.exit("ERROR: please fix these problems in the input files:\n  " + "\n  ".join(problems))


# ---------------------------------------------------------------------------
# 7. React each ligand and name its attachment atom
# ---------------------------------------------------------------------------
# Function to generate Boltz atom names for a ligand molecule
def boltz_atom_names(mol):
    """This function is to figure out what Boltz will name each heavy atom in the ligand"""
    mol_h = Chem.AddHs(mol)
    ranks = list(Chem.CanonicalRankAtoms(mol_h))
    names = {}
    for atom in mol_h.GetAtoms():
        if atom.GetSymbol() != "H":
            names[atom.GetIdx()] = atom.GetSymbol() + str(ranks[atom.GetIdx()] + 1)
    return names


rows = []

for _, ligand in ligands.iterrows():
    ligand_id = ligand["ID"]
    warhead_name = ligand["WARHEAD"]
    spec = WARHEADS[warhead_name]

    mol = Chem.MolFromSmiles(ligand["SMILES"])
    if mol is None:
        problems.append(f"ligands.csv: '{ligand_id}' has a SMILES RDKit cannot read: {ligand['SMILES']}")
        continue

    reaction = AllChem.ReactionFromSmarts(spec["smarts"])
    products = reaction.RunReactants((mol,))

    if not products:
        problems.append(f"ligands.csv: '{ligand_id}' has no {warhead_name} group "
                        f"(SMILES {ligand['SMILES']})")
        continue

    # More than one match is an error:
    if len(products) > 1:
        problems.append(f"ligands.csv: '{ligand_id}' has {len(products)} {warhead_name} groups. "
                        "Pick one and write the ligand with a single warhead.")
        continue

    product = products[0][0]

    # Open the valence on the attachment atom
    attachment_idx = None
    for atom in product.GetAtoms(): # Iterate over all atoms to find the attachment atom
        if atom.GetPropsAsDict().get("old_mapno", atom.GetAtomMapNum()) == 1: # Attachment atom is mapped as :1
            atom.SetNoImplicit(True) # Stops RDKit from adding a hydrogen back again
            atom.SetNumRadicalElectrons(1) # Boltz expects open bond site to be a radical
            attachment_idx = atom.GetIdx() # Record the index of the attachment atom
        atom.SetAtomMapNum(0)

    if attachment_idx is None:
        problems.append(f"WARHEAD '{warhead_name}': its reaction SMARTS does not map an "
                        "attachment atom as :1. Fix the entry in WARHEADS.")
        continue

    # Sanitise to check valid molecule
    try:
        Chem.SanitizeMol(product)
    except Exception as error:
        problems.append(f"ligands.csv: '{ligand_id}' reacted, but the product will not "
                        f"sanitise: {error}")
        continue

    bound_smiles = Chem.MolToSmiles(product)

    # Name the atom from the SMILES as it will be WRITTEN by reading as Boltz would.
    reread = Chem.MolFromSmiles(bound_smiles)
    if reread is None:
        problems.append(f"ligands.csv: '{ligand_id}' produced a bound SMILES that cannot be "
                        f"read back: {bound_smiles}")
        continue

    open_valences = [a.GetIdx() for a in reread.GetAtoms() if a.GetNumRadicalElectrons() > 0]
    if len(open_valences) != 1:
        problems.append(f"ligands.csv: '{ligand_id}' bound form has {len(open_valences)} open "
                        f"valences, expected exactly 1: {bound_smiles}")
        continue

    # Check: Boltz will attempt to build a 3D conformer with ETKDG.
    embed_test = Chem.AddHs(Chem.MolFromSmiles(bound_smiles))
    if AllChem.EmbedMolecule(embed_test, randomSeed=0xC0FFEE) != 0:
        problems.append(f"ligands.csv: '{ligand_id}' bound form will not embed in 3D: {bound_smiles}")
        continue

    # Get the name of the attachment atom as Boltz would see it using the function defined earlier
    attachment_name = boltz_atom_names(reread)[open_valences[0]]

    rows.append({
        "ID": ligand_id,
        "WARHEAD": warhead_name,
        "ORIGINAL_SMILES": ligand["SMILES"],
        "BOUND_SMILES": bound_smiles,
        "ATTACHMENT_ATOM": attachment_name,
        "TARGET_RESIDUE": spec["residue"],
        "TARGET_ATOM": spec["atom"],
    })

if problems:
    sys.exit("ERROR: please fix these problems in the input files:\n  " + "\n  ".join(problems))


# ---------------------------------------------------------------------------
# 8. Write the table, and print it so it can be checked before anything runs
# ---------------------------------------------------------------------------
prepared = pd.DataFrame(rows)
prepared.to_csv(prepared_file, index=False)

print()
print(f"{'ID':<20} {'WARHEAD':<16} {'BOUND SMILES':<34} {'ATOM':<6} BONDS TO")
print("-" * 92)
for row in rows:
    print(f"{row['ID']:<20} {row['WARHEAD']:<16} {row['BOUND_SMILES']:<34} "
          f"{row['ATTACHMENT_ATOM']:<6} {row['TARGET_RESIDUE']} {row['TARGET_ATOM']}")

print()
print(f"Wrote {len(rows)} prepared ligand(s) to {prepared_file}")
print("Check the bound SMILES above before the jobs run: the open valence shown as")
print("[CH2] or similar is the atom that bonds to the protein.")
