# Cysteine map: UniProt P0A183 ↔ PDB 5CRL

Canonical sequence: UniProtKB reviewed entry **P0A183** (`MERR_PSEAEI` / gene `merR`), mercuric resistance operon regulatory protein from *Pseudomonas aeruginosa* (plasmid pVS1 / Tn501). Length **144 aa**. FASTA: [`data/P0A183.fasta`](P0A183.fasta).

Sources (re-fetched for this map):

- FASTA: https://rest.uniprot.org/uniprotkb/P0A183.fasta
- Features: https://rest.uniprot.org/uniprotkb/P0A183.json
- Structure: PDB **5CRL** (X-ray, 2.80 Å; UniProt cross-reference chains A/B = residues 1–134)

## Numbering alignment

PDB 5CRL maps auth chains A and B residues **1–134** to UniProt P0A183 **1–134** (`DBREF`). Residue numbers are identical (**offset = 0**). The C-terminal UniProt residues 135–144 are not modelled in 5CRL and contain no cysteines.

## All cysteines

| UniProt position | Residue | PDB 5CRL resseq | Chains in 5CRL | Role | UniProt evidence | PDB 5CRL evidence |
| ---: | :---: | ---: | :---: | --- | --- | --- |
| 82 | C | 82 | A, B | Hg(2+) ligand | Binding site for Hg(2+) (ChEBI:CHEBI:16793). Mutagenesis C→A abolishes transcriptional activation (PubMed:2551364). | SG–Hg donors ≤3.2 Å: A82 2.565 Å (HG A202), B82 2.547 Å (HG A201). |
| 115 | C | 115 | A, B | Non-ligand | Not annotated as an Hg(2+) binding site. Mutagenesis C→A: slight increase in transcriptional activation (PubMed:2551364). | Modelled on both chains; SG ~11.44 Å from the nearer Hg. Not a LINK donor. |
| 117 | C | 117 | A, B | Hg(2+) ligand | Binding site for Hg(2+) (ChEBI:CHEBI:16793). Mutagenesis C→A decreases transcriptional activation (PubMed:2551364). | SG–Hg donors ≤3.2 Å: A117 2.568 Å (HG A201), B117 2.553 Å (HG A202). |
| 126 | C | 126 | A, B | Hg(2+) ligand | Binding site for Hg(2+) (ChEBI:CHEBI:16793). Mutagenesis C→S abolishes transcriptional activation; loss of Hg binding (PubMed:2551364). | SG–Hg donors ≤3.2 Å: A126 2.562 Å (HG A201), B126 2.568 Å (HG A202). |

There are **four** cysteines in P0A183: C82, C115, C117, C126.

## Hg-binding sites C82 / C117 / C126

UniProt annotates **only** Cys82, Cys117 and Cys126 as Hg(2+) binding-site residues. PDB 5CRL supports the same three cysteines as the coordinating thiolates of the dimeric trigonal HgS3 site (two HETATM HG: A201 and A202):

- HG A201 donors: A117, A126, B82
- HG A202 donors: A82, B117, B126

Cys115 is present in both UniProt and the 5CRL models and is labelled **non-ligand** (no UniProt Hg-binding-site feature; SG far outside the 3.2 Å donor cutoff).

## Sequence context

```text
...LEDGTH C82 EEASSLAEHKLKDVREKMADLARMEAVLSELV C115 A C117 HARRGNVS C126 PLIASL...
```
