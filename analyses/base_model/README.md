# Base Model Analysis

## Overview

This is the base EpiMAP Myeloma model that includes **all available treatment regimens** in the risk equations. This represents the full treatment landscape as observed in the MRDR registry data.

## What's Different About This Model

The base model includes **ALL** treatment regimens in the risk equations:

### Line 1 (LoT 1) Regimens Included
- **VCd** (Bortezomib, Cyclophosphamide, Dexamethasone) - 58%
- **VRd** (Bortezomib, Lenalidomide, Dexamethasone) - 15% ✅ **Included**
- **VTd** (Bortezomib, Thalidomide, Dexamethasone) - 8% 
- **Other regimens** - 19%

### Line 2 (LoT 2) Regimens Included  
- **Rd** (Lenalidomide, Dexamethasone) - 16%
- **DVd** (Daratumumab, Bortezomib, Dexamethasone) - 11%
- **Other regimens** - 73%

## When to Use This Model

Use the base model for:
- Standard population projections with current treatment practices
- Comparing effectiveness across all available regimens
- Baseline estimates for health economic models

## Key Files

- **Coefficients**: `EpiMAP_Coefficients_v2.mmat` (includes all regimen effects)
- **Population**: Any population from `data/populations_2025-2030/`
- **Main script**: `EpiMAP_Myeloma_v2.0.do`

## Usage

```stata
// Basic run with all regimens included
do "EpiMAP_Myeloma_v2.0.do" base SoC 1 Base Population 1 4884 0
```

## Comparison with Other Analyses

| Analysis | VRd Status | Purpose |
|----------|------------|---------|
| **Base Model** | ✅ Included | Current practice projections |
| **VRd Post-Market** | ❌ Excluded | VRd impact assessment |

The VRd post-market analysis excludes VRd to estimate what outcomes would have been without VRd availability.
