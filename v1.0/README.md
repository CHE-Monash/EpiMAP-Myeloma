# EpiMAP Myeloma v1.0

**Initial Multiple Myeloma Simulation Model**

## Publication

This model version supports:
> Irving A, Petrie D, Harris A, Fanning L, Wood EM, Moore E, Wellard C, Waters N, Huynh K, Augustson B, Cook G. Developing and validating a discrete-event simulation model of multiple myeloma disease outcomes and treatment pathways using a national clinical registry. Plos one. 2024 Aug 27;19(8):e0308812.

## What This Version Does

v1.0 provides a general-purpose simulation for multiple myeloma patients with:

- **30 risk equations** based on MRDR patient data
- **Treatment pathways** for up to 9 lines of therapy
- **Patient characteristics** including age, sex, ECOG, ISS
- **Hypothetical dataset** with 1,000 test patients

## Quick Start

```stata
cd v1.0
do "EpiMAP_Myeloma_v1.0.do"
