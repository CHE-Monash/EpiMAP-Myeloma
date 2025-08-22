# EpiMAP Myeloma v2.0

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Stata](https://img.shields.io/badge/Stata-15.0%2B-red.svg)](https://www.stata.com/)
[![DOI](https://img.shields.io/badge/DOI-10.1371%2Fjournal.pone.0308812-blue.svg)](https://doi.org/10.1371/journal.pone.0308812)

**Epidemiological Modelling of Australian Patients with Myeloma**

A comprehensive discrete-event simulation model for multiple myeloma disease outcomes and treatment pathways, developed through collaboration between Monash University's Centre for Health Economics and Transfusion Research Unit.

## What's New in v2.0

- **Reorganised Architecture**: Clear version-specific folder structure (`v1.0/`, `v2.0/`)
- **Extended Treatment Options**: Additional pathways for later lines of therapy (LoT 3+)
- **Robust Testing Framework**: Complete model validation suite with benchmark results
- **Comprehensive Documentation**: Detailed user guides and technical specifications
- **Improved Accuracy**: Updated risk equations reflecting current clinical practice

## Model Overview

EpiMAP Myeloma simulates the complete treatment journey of multiple myeloma patients using **30 evidence-based risk equations** derived from the Australia and New Zealand Myeloma and Related Diseases Registry (MRDR).

### Key Features

- **Patient Characteristics**: Age, sex, ECOG performance score, R-ISS staging
- **Treatment Pathways**: Comprehensive modelling of up to 9 lines of therapy
- **Clinical Outcomes**: Best Clinical Response (BCR) and Overall Survival (OS)
- **ASCT Modelling**: Separate pathways for transplant-eligible patients
- **Maintenance Therapy**: Post-induction treatment modelling
- **Parametric Survival Models**: Time-to-event analysis for all outcomes

## Quick Start

### Prerequisites

- **Stata 15.0 or higher** (valid licence required)
- Windows, macOS, or Linux operating system

### Installation

```bash
git clone https://github.com/your-org/epimap-myeloma.git
cd epimap-myeloma/v2.0
```

### Basic Usage

```stata
cd "path/to/epimap-myeloma/v2.0"
do "EpiMAP_Myeloma_v2.0.do" [Analysis] [Intervention] [Line] [Coefficients] [Data] [MinID] [MaxID] [Bootstrap] [MinBS] [MaxBS]
```

#### Arguments Description

| Position | Argument | Description | Example Values |
|:--------:|:---------|:------------|:---------------|
| 1 | **Analysis** | Analysis identifier | `vrdpost` |
| 2 | **Intervention** | Treatment intervention | `VRd`, `SoC` |
| 3 | **Line** | Line of therapy (1-9) | `1`, `2`, `3` |
| 4 | **Coefficients** | Coefficient set to use | `VRd`, `SoC` |
| 5 | **Data** | Dataset type | `Predicted`, `Population` |
| 6 | **MinID** | Minimum patient ID | `1` |
| 7 | **MaxID** | Maximum patient ID | `4884`, `1000` |
| 8 | **Bootstrap** | Bootstrap flag (0/1) | `0` (no), `1` (yes) |
| 9 | **MinBS** | Minimum bootstrap sample | `1` |
| 10 | **MaxBS** | Maximum bootstrap sample | `5`, `100` |

Example Commands

The simulation will generate `EpiMAP_Simulated_v2.dta` containing comprehensive patient outcomes.

## Repository Structure

- **v2.0/**
  - `EpiMAP_Myeloma_v2.0.do` - Main simulation script
  - **functions/** - Core simulation functions
    - `SIM_BCR_L*.do` - Best Clinical Response models
    - `SIM_CR_L*.do` - Chemotherapy regimen selection
    - `SIM_OS_*.do` - Overall survival models
  - **data/**
    - **coefficients/** - Universal risk equation coefficients
      - `EpiMAP_Coefficients_v2.mmat`
    - **populations_2025-2030/** - Multiple MM population realisations for 2025-2030
      - `EpiMAP_Myeloma_Population_1.dta` 
  - **analyses/** - Specialised analysis scripts
    - **vrdpost/** - VRd post-market analysis
    - **dvdpre/** - DVd pre-market analysis
  - **documentation/** - Comprehensive user guides
    - `User_Guide.md`
    - `Technical_Specifications.md`
    - `Parameter_Reference.md`

## Clinical Applications

### Treatment Regimens Modelled

**Line 1 (LoT 1)**
- VCd (Bortezomib, Cyclophosphamide, Dexamethasone): 58%
- VRd (Bortezomib, Lenalidomide, Dexamethasone): 15%
- Other regimens: 26%

**Line 2 (LoT 2)**
- Rd (Lenalidomide, Dexamethasone): 16%
- DVd (Daratumumab, Bortezomib, Dexamethasone): 11%
- Other regimens: 73%

**Lines 3-9 (LoT 3+)**
- Averaged survival benefit approach for emerging therapies

### Risk Equations (30 Total)

| Category | Equations | Purpose |
|----------|-----------|---------|
| **Survival** | Overall Survival | Primary endpoint modelling |
| **Treatment Planning** | ASCT eligibility, regimen selection | Clinical decision support |
| **Response Prediction** | BCR for each LoT | Treatment effectiveness |
| **Time Intervals** | Chemotherapy duration, treatment-free periods | Disease progression modelling |

## Validation and Performance

- **Calibrated** against MRDR registry data (2008-2023)
- **Validated** survival curves match observed patient outcomes
- **Benchmarked** treatment distributions reflect real-world practice
- **Tested** across diverse patient populations and scenarios

## Advanced Usage

### Custom Patient Populations

Replace the hypothetical dataset with your institution's data:

```stata
// Edit line 24 in EpiMAP_Myeloma_v2.0.do
use "your_patient_data.dta", clear
```

**Required variables**: `age`, `male`, `ecog`, `iss`

### Scenario Analyses

Run specific treatment comparisons:

```stata
// VRd post-market analysis
do "analyses/vrdpost/EpiMAP_Myeloma_Analysis_VRd-Post.do"

// DVd pre-market analysis  
do "analyses/dvdpre/EpiMAP_Myeloma_Analysis_DVd-Pre.do"
```

### Bootstrap Validation

Enable uncertainty quantification:

```stata
global Bootstrap = 1
global BSMin = 1
global BSMax = 100
do "EpiMAP_Myeloma_v2.0.do"
```

## Output Variables

The simulation generates comprehensive outcomes for each patient:

- **Survival**: Overall survival time, mortality indicators
- **Treatment**: Regimen assignments, duration, response rates
- **Progression**: Treatment-free intervals, line progression
- **Clinical**: ASCT receipt, maintenance therapy, ECOG changes

## Contributing

We welcome collaborative research! Please see our [contribution guidelines](CONTRIBUTING.md) for:

- Bug reports and feature requests
- Model validation with external datasets
- Extension to new treatment regimens
- Integration with health economic models

## Documentation

- **[User Guide](documentation/User_Guide.md)**: Step-by-step instructions
- **[Technical Specifications](documentation/Technical_Specifications.md)**: Detailed model description
- **[Parameter Reference](documentation/Parameter_Reference.md)**: Complete coefficient listing
- **[Migration Guide](documentation/Migration_Guide.md)**: Upgrading from v1.0

## Citation

### Primary Publication
> Irving A, Petrie D, Harris A, et al. Developing and validating a discrete-event simulation model of multiple myeloma disease outcomes and treatment pathways using a national clinical registry. *PLOS ONE*. 2024;19(8):e0308812. [doi:10.1371/journal.pone.0308812](https://doi.org/10.1371/journal.pone.0308812)

### Software Citation
```bibtex
@software{epimap_myeloma_v2,
  title = {EpiMAP Myeloma: Epidemiological Modelling of Australian Patients with Myeloma},
  author = {Irving, Adam and Petrie, Dennis and Harris, Anthony and others},
  version = {2.0.0},
  year = {2025},
  url = {https://github.com/your-org/epimap-myeloma},
  institution = {Monash University}
}
```

## Research Team

**Health Economists**: Adam Irving, Dennis Petrie, Anthony Harris, Laura Fanning

**Clinical Experts**: Zoe K McQuilten, Erica M Wood, Bradley Augustson, Gordon Cook, Francesca Gay, Georgia McCaughan, Peter Mollee, Andrew Spencer

**Registry Team**: Elizabeth Moore, Cameron Wellard, Neil Waters

**Consumer Representative**: Andrew Marks

## Data Access

### MRDR Patient Data
For access to genuine patient data from the Australia and New Zealand Myeloma and Related Diseases Registry:
- Submit applications to the MRDR Steering Committee
- Visit: [mrdr.net.au](https://www.mrdr.net.au/)

## Licence

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Support

- **Model Questions**: [adam.irving@monash.edu](mailto:adam.irving@monash.edu)
- **Technical Issues**: [Create an issue](https://github.com/your-org/epimap-myeloma/issues)
- **Collaboration Enquiries**: Contact the research team

## Related Resources

- [Australia and New Zealand Myeloma and Related Diseases Registry](https://www.mrdr.net.au/)
- [Monash Centre for Health Economics](https://www.monash.edu/business/che)
- [Monash Transfusion Research Unit](https://www.monash.edu/medicine/sphpm/units/transfusion-research)

## Acknowledgements

The EpiMAP Myeloma project was supported by grant 1200706 to Prof Zoe K McQuilten. We thank patients, clinicians, and research staff at participating centres for their invaluable contributions to the MRDR.

---

**Important**: This model is designed for research purposes. Clinical decisions should always involve qualified healthcare professionals and consider individual patient circumstances.
