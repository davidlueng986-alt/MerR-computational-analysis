# MerR computational analysis

**Mechanistic modelling of Hg(II) speciation, ligand exchange, and transfer into the MerR Hg-binding site for an E. coli whole-cell mercury biosensor.**

> Project context: iGEM 2026 — soil mercury sensing using chloride-based extraction coupled to an **E. coli + MerR/pMer + dTomato** reporter circuit.

![MerR/pMer/dTomato circuit](assets/merR_pMer_dTomato_circuit.svg)

## Why this repository exists

The difficult part of a MerR whole-cell mercury sensor is not simply whether MerR can bind Hg(II). The real engineering problem is whether Hg extracted from soil in a chloride-rich chemical environment can be converted into a near-neutral biological assay environment while retaining enough **MerR-accessible Hg** to generate a strong and reproducible reporter response.

A chloride extract may contain Hg as a mixture of `HgCl+`, `HgCl2(aq)`, `HgCl3-`, `HgCl4^2-`, hydrolysed species such as `HgOH+` / `Hg(OH)2(aq)`, and complexes with dissolved organic matter or sulfur ligands. Therefore:

> **HgCl2 reagent added ≠ HgCl2 is the equilibrium species ≠ Hg is membrane-accessible ≠ Hg reaches MerR.**

This repository builds the computational layer needed to separate these steps rather than treating all dissolved Hg as equivalent.
