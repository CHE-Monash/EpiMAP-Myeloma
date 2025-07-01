# Changelog

All notable changes to the EpiMAP Myeloma project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Planned
- Enhanced validation framework for model parameters
- Additional treatment regimens for LoT 3+
- Integration with R for post-processing analysis

## [2.0.0] - 2025-01-XX

### Added
- Reorganised repository structure with clear version folders (v1.0/, v2.0/)
- Enhanced documentation with detailed user guide
- Improved parameter validation and error checking
- New treatment pathway options for later lines of therapy
- Comprehensive testing framework for model validation
- Better handling of edge cases in survival calculations

### Changed
- **BREAKING**: Repository structure reorganised - models now in version-specific folders
- **BREAKING**: Updated file naming conventions for better clarity
- Improved simulation performance through code optimisation
- Enhanced random number generation for better reproducibility
- Updated coefficient matrix structure for additional parameters
- Refined treatment-free interval calculations

### Fixed
- Corrected edge case in ASCT eligibility determination
- Fixed rare numerical precision issues in survival probability calculations
- Resolved memory allocation issues in large simulation runs
- Improved handling of missing data in patient characteristics

### Model Updates
- Updated risk equations based on latest MRDR data analysis
- Refined treatment regimen proportions to reflect current clinical practice
- Enhanced Best Clinical Response prediction accuracy
- Improved Overall Survival curve fitting

### Documentation
- Complete rewrite of user documentation
- Added model specification document with technical details
- Created parameter reference guide
- Added validation report with benchmark results

## [1.0.0] - 2024-XX-XX

### Added
- Initial public release of EpiMAP Myeloma simulation model
- Discrete-event simulation framework for multiple myeloma outcomes
- 30 risk equations based on MRDR patient-level data
- Support for up to 9 Lines of Therapy (LoTs)
- Hypothetical patient dataset with 1,000 patients for testing
- Complete documentation and usage instructions

### Model Features
- **Patient Characteristics**: Age, sex, ECOG performance score, ISS stage
- **Treatment Pathways**: Comprehensive modelling of treatment sequences
- **Survival Analysis**: Parametric survival models for overall survival
- **Clinical Response**: Ordered logit models for Best Clinical Response
- **ASCT Support**: Separate analysis paths for transplant-eligible patients
- **Maintenance Therapy**: Modelling of post-induction maintenance treatment

### Risk Equations
1. Overall Survival
2. Planned ASCT eligibility
3. Diagnosis to treatment interval
4. LoT 1 chemotherapy regimen selection
5-7. LoT 1 chemotherapy duration (with manual splines for ASCT patients)
8. LoT 1 chemotherapy duration (non-ASCT patients)
9. LoT 1 Best Clinical Response
10. Receipt of ASCT
11. ASCT Best Clinical Response
12. Receipt of maintenance therapy
13-14. LoT 1 to LoT 2 treatment-free intervals (ASCT vs non-ASCT)
15-17. LoT 2 treatment pathways and outcomes
18-30. LoTs 3-6+ treatment pathways and outcomes

### Chemotherapy Regimens
- **LoT 1**: VCd (58%), VRd (15%), Other (26%)
- **LoT 2**: Rd (16%), DVd (11%), Other (73%)
- **LoT 3+**: Averaged survival benefit approach

### Technical Specifications
- **Platform**: Stata 15.0 or higher required
- **Programming**: Stata with Mata matrix operations
- **Data Format**: Stata .dta files for patient data
- **Coefficients**: Mata .mmat matrix files for risk equations
- **Output**: Comprehensive simulated patient outcomes dataset

### Validation
- Model calibrated against MRDR registry data (20XX-20XX)
- Survival curves validated against observed patient outcomes
- Treatment pathway distributions match registry patterns
- Extensive sensitivity analysis conducted

## [0.1.0] - Internal Development

### Added
- Initial model framework development
- Risk equation estimation from MRDR data
- Prototype simulation algorithms
- Internal validation and testing

---

## Release Notes

Each release includes:
- **Complete Model**: All necessary files to run simulations
- **Documentation**: User guides, technical specifications, and examples
- **Validation Results**: Benchmark tests and model performance metrics
- **Example Data**: Hypothetical patient datasets for testing

## Migration Guides

### Upgrading from v1.0 to v2.0
- Repository structure has changed - see Documentation/Migration_Guide.md
- No changes to input data format required
- Enhanced output includes additional variables
- Improved parameter validation may flag previously undetected issues

## Contact

- **Model Questions**: adam.irving@monash.edu
- **Technical Issues**: Create an issue on GitHub
- **MRDR Data Access**: Visit mrdr.net.au
