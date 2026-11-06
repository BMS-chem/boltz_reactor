# Sigma-2 benchmark

A check that the pipeline produces sensible numbers, by running it on a published
compound set and comparing against the result reported for the same set.

## The benchmark

Kim et al., *Large-scale prospective evaluation of co-folding across hundreds of
Mac1-ligand complexes and three virtual screens*, eLife RP110475 (2026)
([article](https://elifesciences.org/articles/110475)) compares AlphaFold3, Chai-1
and Boltz-2 against DOCK3.7 on three virtual screens. One of them is the sigma-2
receptor, with 506 compounds taken from a 490-million-molecule docking screen
(Alon et al., *Nature* 600:759-764, 2021) and tested experimentally: 201 actives
and 305 inactives, by the paper's definition of more than 25% displacement of
[3H]DTG.

Its supplementary file 5 gives, per compound, the SMILES, the measured activity,
and the scores from each method. That makes it possible both to run the compounds
through this pipeline and to recompute the paper's own numbers from its own
columns.

## What was run

506 ligands against one protein, affinity prediction on, otherwise Boltz-2 v2.2.0
defaults. One Slurm job on a single A100: 3 h 42 m end to end, of which 2 h 11 m
was structure prediction and the rest the MSA fetch and the affinity stage.

The protein is **human TMEM97** (UniProt Q5BJF2), full length, no cofactors.

That choice is worth stating, because it differs from what the docking study used.
The affinities being predicted were measured on the human receptor, in Expi293
cells (Alon et al. 2021, Methods). Every sigma-2 crystal structure in the PDB is
*bovine*, and the docking screen was run against one of them (7MFI), because
docking needs coordinates. Co-folding needs only a sequence, so it has no reason
to inherit that constraint, and the sequence was matched to the assay that
produced the labels instead. Human and bovine are 78% identical over 168 aligned
residues.

Cofactors were left empty. The paper's own Boltz input generator builds a protein
chain plus a ligand chain and nothing else, and the cholesterol present in 7MFI
would compete with the ligand being scored.

## Result

Ranking the 506 compounds by Boltz's predicted affinity (`AFFINITY_PRED_VALUE`,
log10 IC50 in uM, so negated to make higher mean tighter):

| Method | ROC AUC |
| --- | --- |
| **This pipeline, Boltz-2 affinity** | **0.848** |
| Boltz-2, as published | 0.838 |
| DOCK3.7, as published | 0.788 |
| AlphaFold3 ligand pLDDT, as published | 0.561 |

![ROC curve](../experiments/001_eLife_Sigma2_test/enrichment/roc_affinity_hTMEM97.png)

All 506 labelled compounds had a prediction, so nothing was dropped.

The published Boltz-2 figure is reproduced closely, slightly above it. The
difference is not interpretable as an improvement: the paper does not state its
receptor sequence or construct, its Boltz command line, MSA settings, number of
diffusion samples, recycling steps, or seeds, so this is not a controlled
comparison. It is a check that an independently built pipeline lands in the same
place.

## Checking the measurement itself

Before trusting the 0.848, the AUC code was checked against values it should
already know the answer to: the paper's own `Boltz2_pKi` and `DOCK score` columns
were passed through the same script, which returned **0.839** and **0.788**
against the published 83.8 and 78.8.

That also establishes what the published metric is: a plain ROC AUC over the
201/305 split, with no weighting or filtering. The medians come out at 6.70 for
actives and 5.38 for inactives, exactly as reported.

`pIC50 = 6 - affinity_pred_value` is a monotonic transform, so it makes no
difference to a rank-based measure which of the two the paper used.

## Caveats

- **The receptor was in Boltz-2's training data.** Its training cutoff is 1 June
  2023 and the Alon 2021 sigma-2 structures were deposited before that. This is a
  caveat on what the number means, not on whether the pipeline works.
- The comparison is uncontrolled, for the reasons above.
- One protein, one target class. It says nothing about how the pipeline behaves
  on a target Boltz has not seen.

## Reproducing it

The inputs, the labels and the collected results are in the repository; the
predicted structures (several hundred MB) are not. So the AUC can be recomputed
without a GPU:

```
python Scripts/calculate_enrichment.py --experiment 001_eLife_Sigma2_test
```

To re-run the predictions themselves, on a cluster with Boltz available:

```
./Scripts/submit_boltz.sh 001_eLife_Sigma2_test
```

`experiments/001_eLife_Sigma2_test/reference/README.md` records where every file
in that folder came from.
