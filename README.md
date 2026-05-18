# KPC-204-FINAL
# KPC β-lactamase + Avibactam — Molecular Dynamics Simulations

> Comparative MD study of KPC-2 (wild-type) and KPC-204 (V204 variant) β-lactamases  
> in complex with avibactam, using crystallographic and homology-modeled structures.

---

## 📋 Table of Contents

- [Project Overview](#project-overview)
- [Repository Structure](#repository-structure)
- [Dependencies](#dependencies)
- [Workflow](#workflow)
  - [1. Ligand Preparation](#1-ligand-preparation)
  - [2. Molecular Docking](#2-molecular-docking)
  - [3. System Preparation](#3-system-preparation)
  - [4. MD Simulations](#4-md-simulations)
  - [5. Analysis](#5-analysis)
- [Systems](#systems)
- [Key Results](#key-results)
- [Known Issues & Notes](#known-issues--notes)
- [Citation](#citation)


## Project Overview

This project investigates the binding dynamics of **avibactam** (a non-β-lactam β-lactamase inhibitor)  
within the active site of two KPC variants:

| Protein | Structure Source | Ligand |
|---------|-----------------|--------|
| KPC-2 | X-ray crystallography | Avibactam (CID 9835049) |
| KPC-204 (V204 variant) | SwissModel homology model | Avibactam (CID 9835049) |

Each system underwent 100 ns of classical MD simulation using **GROMACS 2025.4**  
with the **AMBER99SB-ILDN** force field and **GAFF2** parameters for the ligand.


## Repository Structure
```
├── docking/
│   ├── avibactam_REAL.pdbqt          # Correct avibactam ligand (CID 9835049)
│   ├── KPC2_cristallo/
│   │   ├── KPC2_cristallo_receptor.pdbqt
│   │   ├── KPC2_cristallo_docked.pdbqt
│   │   ├── KPC2_cristallo_best_pose.pdb
│   │   └── vina_config.txt
│   └── KPC204_swissmodel/
│       ├── KPC204_swissmodel_receptor.pdbqt
│       ├── KPC204_swissmodel_docked.pdbqt
│       ├── KPC204_swissmodel_best_pose.pdb
│       └── vina_config.txt
│
├── ligand/
│   └── avibactam_REAL.acpype/
│       ├── avibactam_REAL_GMX.itp    # GAFF2 topology
│       ├── avibactam_REAL_GMX.gro    # Ligand coordinates
│       └── posre_avibactam_REAL.itp  # Position restraints
│
├── mdp/
│   ├── ions.mdp                      # Minimal MDP for genion
│   ├── em.mdp                        # Energy minimization
│   ├── nvt.mdp                       # NVT equilibration (300K, 100 ps)
│   ├── npt.mdp                       # NPT equilibration (1 bar, 100 ps)
│   └── md.mdp                        # Production MD (100 ns)
│
├── systems/
│   ├── KPC2_cristallo_v2/            # ← ACTIVE (correct ligand)
│   │   ├── protein_clean.pdb
│   │   ├── MOL.itp
│   │   ├── topol.top
│   │   ├── complex_ions.gro
│   │   ├── em.gro / nvt.gro / npt.gro
│   │   └── md.xtc / md.tpr
│   └── KPC204_swissmodel_v2/         # ← ACTIVE (correct ligand)
│       ├── protein_clean.pdb
│       ├── MOL.itp
│       ├── topol.top
│       └── ...
│
└── analysis/
├── plot_analysis.py              # RMSD / RMSF / Rg plots
├── rmsd.xvg
├── rmsf.xvg
├── gyrate.xvg
└── dist_ser70_C7.xvg             # Ser70–Avibactam C7 distance
# help from claude for the format beyong
```
