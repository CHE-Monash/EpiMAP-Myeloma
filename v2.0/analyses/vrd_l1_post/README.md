# VRd LoT 1 Post-Market Analysis

## Research Question

**"What is the clinical impact of VRd at LoT 1 compared to historical practice?"**

## Key Methodology

This analysis **excludes VRd from the risk equations** to simulate what would have happened to VRd patients under historical treatment patterns.

## Treatment Regimens in Risk Equations

### What's Excluded
- **VRd** (Bortezomib, Lenalidomide, Dexamethasone) ❌ **EXCLUDED**

### What's Included
- **VCd** (Bortezomib, Cyclophosphamide, Dexamethasone) ✅ 
- **VTd** (Bortezomib, Thalidomide, Dexamethasone) ✅
- **Other historical regimens** ✅

## How It Works

1. **Identify VRd patients**: Patients who would receive VRd based on their characteristics
2. **Simulate without VRd**: What regimen would they have gotten historically?
3. **Compare outcomes**: VRd effectiveness vs historical alternatives

## Analysis Scenarios

### Scenario 1: Historical Practice (SoC)
```stata
do "analyses/vrd_l1_post/EpiMAP_Myeloma_VRd_L1_Post.do" vrdpost SoC 1 VRd Predicted 1 4884 0
```
- VRd patients get alternative regimens (VCd, VTd, Other)
- Shows what outcomes would have been without VRd

### Scenario 2: VRd Available
```stata  
do "analyses/vrd_l1_post/EpiMAP_Myeloma_VRd_L1_Post.do" vrdpost VRd 1 VRd Predicted 1 4884 0
```
- VRd patients receive VRd
- Shows actual VRd impact

## Key Files

- **Coefficients**: `EpiMAP_Myeloma_Coefficients_VRd_L1_Post.dta` (VRd excluded)
- **Patient cohort**: `EpiMAP_Myeloma_Patients_VRd_L1_Post.dta` (VRd-eligible patients)
- **Bootstrap coefficients**: 500 samples in `bootstrap/` folder

## Expected Results

Compare survival, response rates, and treatment pathways between scenarios to quantify VRd's clinical benefit.

## Files Structure

```
analyses/vrd_l1_post/
├── README.md                                           # This file
├── EpiMAP_Myeloma_VRd_L1_Post.do                      # Main analysis script  
└── data/
    ├── coefficients/
    │   ├── EpiMAP_Myeloma_Coefficients_VRd_L1_Post.dta # VRd-excluded coefficients
    │   └── bootstrap/                                   # 500 bootstrap samples
    └── patients/
        └── EpiMAP_Myeloma_Patients_VRd_L1_Post.dta     # VRd-eligible cohort
```
