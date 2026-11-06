# Reference data for the sigma-2 benchmark

Everything here was downloaded from public sources on 2026-11-05. Nothing in this
folder is produced by the pipeline, and nothing here is read by it: it is source
material and a record of where the numbers came from.

## The paper

Kim J, Correy GJ, Hall BW, Rachman MM, Mailhot O, Togo T, Gonciarz RL, Jaishankar P,
Neitz RJ, Hantz ER, Doruk YU, Stevens MGV, Diolaiti ME, Reid R, Gopalkrishnan S,
Krogan NJ, Renslo AR, Ashworth A, Shoichet BK, Fraser JS.
*Large-scale prospective evaluation of co-folding across hundreds of Mac1-ligand
complexes and three virtual screens.* eLife RP110475 (2026).
https://elifesciences.org/articles/110475

It compares AlphaFold3, Chai-1 and Boltz-2 against DOCK3.7 on the SARS-CoV-2 Mac1
macrodomain and on the tested hit lists of three large docking screens (AmpC
beta-lactamase, sigma-2 receptor, dopamine D4 receptor).

Code: https://github.com/jongbin99/cofold-em

## Files

| File | Source |
|---|---|
| `elife-110475-supp5.xlsx` | Supplementary file 5, "Data collection of all docked hits on AmpC, sigma2, D4". Downloaded from the eLife CDN. |
| `sigma2_compounds.csv` | The `Sigma2` sheet of that file, converted to CSV. |
| `ampc_compounds.csv`, `d4_compounds.csv` | The other two sheets, kept for completeness. Not used by this experiment. |
| `TMEM97_human_Q5BJF2.fasta` | Human TMEM97 (sigma-2 receptor), UniProt Q5BJF2, 176 aa. |
| `sigma2_bovine_7M95_construct.fasta` | The bovine sequence as crystallised in PDB 7M95, with and without its N-terminal GPGGSS tag. |
| `structures/*.cif` | mmCIF files from the RCSB PDB. |

### The compound set

`sigma2_compounds.csv` holds 506 compounds: 201 `True hits` and 305 `Non-hits`,
matching the counts in the paper. Columns include the SMILES, the ZINC ID, the
measured pKi, and the three scores the paper compared:

- `DOCK score (kcal/mol)` - DOCK3.7, more negative is better
- `L-pLDDT score` - AlphaFold3 ligand pLDDT, higher is better
- `Boltz2_pKi` - the paper's own Boltz-2 prediction, higher is better

The compounds come from the screen in Alon A, Lyu J, Braz JM, et al.,
*Structures of the sigma-2 receptor enable docking for bioactive ligand discovery*,
Nature 600:759-764 (2021), a 490-million-molecule docking screen.

How actives were defined (eLife Methods, "Benchmarking co-folding on
experimentally validated compounds"): more than 25% displacement of [3H]DTG,
relaxed from the original 50% to admit more molecules. Single-point 1 uM data were
converted to an apparent Ki with the Cheng-Prusoff equation, which the authors say
gives a non-binding cutoff of apparent Ki < 2.0 uM.

### Published results for sigma-2 (ROC AUC, actives vs inactives)

| | AUC |
|---|---|
| Boltz-2 pIC50 | 83.8 |
| DOCK3.7 | 78.8 |
| AF3 L-pLDDT | 56.1 |

Median Boltz-2 pIC50: 6.70 for actives, 5.38 for inactives.

Recomputing these from the `Boltz2_pKi`, `DOCK score` and `L-pLDDT score` columns
of `sigma2_compounds.csv`, as a plain ROC AUC over the 201/305 split, gives 83.9,
78.8 and 56.3, and medians of 6.70 and 5.38. So the published numbers are a plain
ROC AUC on these columns, with no extra weighting or filtering.

### Structures

All six are X-ray structures of the **bovine** sigma-2 receptor. There is no human
TMEM97 structure in the PDB (UniProt Q5BJF2 lists no PDB cross-references as of
2026-11-05). Human and bovine are 78% identical over 168 aligned residues.

| PDB | Resolution | Bound ligand |
|---|---|---|
| `7M93` | 2.94 A | PB28 |
| `7M94` | 2.71 A | Roluperidone |
| `7M95` | 2.41 A | Z1241145220 (a docking hit) |
| `7M96` | 2.41 A | Z4857158944 (a docking hit) |
| `7MFI` | 2.81 A | Cholesterol |
| `9ZIR` | 2.74 A | Zervimesine (CT1812). From a later paper, ACS Med Chem Lett 2026. |

`7M93`-`7M96` and `7MFI` are the structures from the Alon 2021 paper, so they are
what the original docking screen was built on.

## Which sequence to give Boltz, and why

Use **human TMEM97 (Q5BJF2), full length, no cofactors**.

The rule is to match the sequence to the protein the measurements were made on,
because those measurements are the labels being predicted. It is not to match
whatever was crystallised: a structure is needed for docking, which works on
coordinates, but co-folding takes only a sequence.

From the Methods of Alon et al. 2021 (Nature 600:759-764, PMC8867396):

- The affinities of the docking hits were measured on the **human** receptor:
  "The human sigma-2 receptor was cloned into pcDNA3.1 (Invitrogen) mammalian
  expression vector", using "membrane samples from Expi293 cells (Thermo Fisher
  Scientific) expressing wild-type or mutant sigma-2 receptor".
- The **bovine** receptor appears only in the crystallography: "The bovine sigma-2
  receptor was cloned into pVL1392" and "This receptor construct was expressed in
  Sf9 insect cells". The human receptor was expressed but was not tractable for
  structural work.
- The docking screen itself docked into the bovine cholesterol-bound structure:
  "The sigma-2 receptor bound to cholesterol (PDB ID: 7MFI) was used in the
  docking calculations."

So the `DOCK score` column was produced by docking into a bovine structure to
predict affinities measured on the human protein. Boltz has no such constraint,
since it needs no coordinates, and can use the human sequence that matches the
assay.

The bovine sequence is kept here anyway, as `sigma2_bovine_7M95_construct.fasta`,
in case it is worth running as a second protein row for comparison. The pipeline
pairs every protein with every ligand, so adding that row doubles the run.

## What is not known

The eLife paper does not say which sigma-2 sequence or construct it gave the
co-folding methods: no UniProt accession, no PDB ID, no sequence. It also does not
give the Boltz command line, the MSA settings, the number of diffusion samples,
recycling steps or seeds, or whether results were averaged over seeds. Its Boltz
input generator (`scripts/preprocessing/input_yaml_generator.py` in the repo)
builds a protein chain A plus a ligand chain B with `properties: affinity,
binder: B`, and reads the protein part from a template YAML that is not in the
repository.

It also does not state which field of Boltz's output `Boltz2_pKi` came from.
Boltz writes `affinity_pred_value` (log10 IC50 in uM) and
`affinity_probability_binary`. Since pIC50 = 6 - log10(IC50 in uM), their value is
most likely `6 - affinity_pred_value`. That transform is monotonic, so it does not
affect a rank-based measure such as AUC either way.

Note that Boltz-2's training cutoff is 1 June 2023, and the Alon 2021 structures
were deposited before that, so the sigma-2 receptor was available to it in
training.
