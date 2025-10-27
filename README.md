# EpiMAP Myeloma v2.1

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Stata](https://img.shields.io/badge/Stata-15.0%2B-red.svg)](https://www.stata.com/)
[![DOI](https://img.shields.io/badge/DOI-10.1371%2Fjournal.pone.0308812-blue.svg)](https://doi.org/10.1371/journal.pone.0308812)

**Epidemiological Modelling of Australian Patients with Myeloma**

A comprehensive discrete-event simulation model for multiple myeloma disease outcomes and treatment pathways, developed through collaboration between Monash University's Centre for Health Economics and Transfusion Research Unit.

## What's New in v2.1

- **Vectorised Implementation**: Complete rewrite using Mata's vector and matrix operations for dramatically improved performance
- **Enhanced Code Organisation**: `mata_setup.do` provides clean, maintainable vector-based architecture
- **Improved Efficiency**: Matrix-based computations replace patient-level loops for faster large-scale simulations
- **Modernised Repository**: Git-based versioning with clean structure (no version folders)

## Model Overview

EpiMAP Myeloma simulates the complete treatment journey of multiple myeloma patients using **30 evidence-based risk equations** derived from the Australia and New Zealand Myeloma and Related Diseases Registry (MRDR).

### Key Features

- **Patient Characteristics**: Age, sex, ECOG performance score, R-ISS staging & co-morbidity score
- **Treatment Pathways**: Comprehensive modelling of up to 9 lines of therapy
- **Clinical Outcomes**: Best Clinical Response (BCR) and Overall Survival (OS)
- **ASCT Modelling**: Separate pathways for transplant-eligible patients
- **Maintenance Therapy**: Post-induction treatment modelling
- **Parametric Survival Models**: Time-to-event analysis for all outcomes
- **High Performance**: Mata implementation for efficient large-scale simulations

## Quick Start

### Prerequisites

- **Stata 15.0 or higher** (valid licence required)
- Windows, macOS, or Linux operating system

### Installation

```bash
git clone https://github.com/CHE-Monash/EpiMAP-Myeloma.git
cd EpiMAP-Myeloma
```

### Basic Usage

```stata
cd "path/to/epimap-myeloma"
do "run.do"
```

The `run.do` file contains example parameters. Edit this file to customise your simulation with the following arguments:

#### Arguments Description

| Position | Argument | Description | Example Values |
|:--------:|:---------|:------------|:---------------|
| 1 | **Analysis** | Analysis identifier | `base_model`, `vrd_l1_post` |
| 2 | **Intervention** | Treatment intervention | `VRd`, `SoC`, `all` |
| 3 | **Line** | Line of therapy (1-9) | `1`, `2`, `3` |
| 4 | **Coefficients** | Coefficient set to use | `base_model`, `VRd`, `SoC` |
| 5 | **Data** | Dataset type | `Predicted`, `Population` |
| 6 | **MinID** | Minimum patient ID | `1` |
| 7 | **MaxID** | Maximum patient ID | `10`, `1000`, `4884` |
| 8 | **Bootstrap** | Bootstrap flag (0/1) | `0` (no), `1` (yes) |
| 9 | **MinBS** | Minimum bootstrap sample | `1` |
| 10 | **MaxBS** | Maximum bootstrap sample | `5`, `100` |
| 11| **Report** | Report flat | `0' (no), `1' (yes) |

### Example Commands

```stata
// Quick test with 10 patients
do "main.do" base_model all 0 base_model population 1 10 0 0 0 0

// Full population simulation
do "main.do" base_model all 0 base_model population 1 4884 0 0 0 1

// Bootstrap analysis
do "main.do" base_model VRd 1 base_model predicted 1 1000 1 1 100 1
```

The simulation will generate results in `analyses/[analysis_name]/data/simulated/`.

## Repository Structure

```
EpiMAP-Myeloma/
├── core/                  # Core simulation ine
│   ├── matrix_setup.do   # Matrix initialisation
│   ├ load_patients.do  # Patient data loading
│   ├ mata_functions.do # Mata utility functions
│   └── outcomes/         # Outcome simulation modules
├── analyses/             # Analysis-specific configurations
│   ├── base_model/       # Base model with all regimens
│   └── vrd_l1_post/      # VRd post-market analysis
├── patients/             # Patient cohort data
│   ├── population/       # Population projections
│   └── predicted/        # Predicted cohorts
├── tests/                # Validation test suite
│   ├── validate_vectors.do
│   └── test_*.do
├─EpiMAP_Myeloma.do    # Main simulation dispatcher
├─run.do               # Example wrapper script
├─README.md
├─CHANGELOG.md
└── LICENSE
```

## Performance Improvements in v2.1

The vectorised implementation provides significant performance benefits:

- **Faster Execution**: Matrix operations process all patients simultaneously
- **Better Memory Efficiency**: Reduced overhead from loop-based processing
- **Scalability**: Handles large cohorts (10,000+ patients) more efficiently
- **Maintainability**: Cleaner code structure with modular vector setup

## Validation

The model has been comprehensively validated:

- **Published Validation**: See Irving et al. (2024) in PLOS ONE
- **Out-of-Sample Testing**: 70/30 split validation with 100 bootstrap iterations
- **Survival Curve Accuracy**: No significant difference in 90% of 120 months post-diagnosis
- **Vectorisation Validation**: Comprehensive test suite confirms identical results to original implementation

## Version History

Access previous versions via Git tags:

- **v2.1**: Current version (vectorised implementation)
- **v2.0**: Reorganised architecture with extended treatment options
- **v1.0**: Initial public release (August 2024) - `git checkout v1.0`

## Citation

### Model Paper

Irving A, Petrie D, Harris A, Fanning L, Wood EM, Moore E, et al. Developing and validating a discrete-event simulation model of multiple myeloma disease outcomes and treatment pathways using a national clinical registry. *PLOS ONE*. 2024;19(8):e0308812. [doi:10.1371/journal.pone.0308812](https://doi.org/10.1371/journal.pone.0308812)

### Software Citation

```bibtex
@software{epimap_myeloma_v2_1,
  title = {EpiMAP Myeloma: Epidemiological Modelling of Australian Patients with Myeloma},
  author = {Irving, Adam and Petrie, Dennis and Harris, Anthony and Fanning, Laura and 
            Wood, Erica M and Moore, Elizabeth and Wellard, Cameron and Waters, Neil and
            Augustson, Bradley and Cook, Gordon and Gay, Francesca and McCaughan, Georgia and
            Mollee, Peter and Spencer, Andrew and McQuilten, Zoe K},
  version = {2.1},
  year = {2025},
  url = {https://github.com/CHE-Monash/EpiMAP-Myeloma},
  doi = {10.1371/journal.pone.0308812},
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

## Support

- **Model Questions**: [adam.irving@monash.edu](mailto:adam.irving@monash.edu)
- **Technical Issues**: [Create an issue](https://github.com/CHE-Monash/EpiMAP-Myeloma/issues)
- **Collaboration Enquiries**: Contact the research team

## Licence

This project is licenced under the GNU General Public Licence v3.0 - see the [LICENSE](LICENSE) file for details.

## Related Resources

- [Australia and New Zealand Myeloma and Related Diseases Registry](https://www.mrdr.net.au/)
- [Monash Centre for Health Economics](https://www.monash.edu/business/che)
- [Monash Transfusion Research Unit](https://www.monash.edu/medicine/sphpm/units/transfusion-research)

## Acknowledgements

The EpiMAP Myeloma project was supported by Medical Research Future Fund GNT1200706 and GNT2017480 and by National Health and Medical Research Council GNT1189490, GNT2024876 and GNT2036025. We thank patients, clinicians, and research staff at participating centres for their invaluable contributions to the MRDR.

---

**Important**: This model is designed for research purposes. Clinical decisions should always involve qualified healthcare professionals and consider individual patient circumstances.
