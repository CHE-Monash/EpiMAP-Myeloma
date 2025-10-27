# Changelog

All notable changes to the EpiMAP Myeloma project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Planned
- Integration with R for post-processing analysis
- Extended documentation for health economic applications

## [2.1] - 2025-10-27

### Added
- **Vectorised Implementation**: Complete rewrite of core simulation engine using Mata vector and matrix operations
- **`core/vector_setup.do`**: New modular vector setup module for efficient data preparation
- **Comprehensive Test Suite**: New validation framework in `tests/` directory
  - `validate_vectors.do`: Validates vectorised implementation against original
  - Additional outcome-specific validation tests
- **Improved Error Handling**: Enhanced validation of vector dimensions and data consistency
- **Performance Metrics**: Built-in validation summaries with detailed reporting

### Changed
- **Core Engine**: Replaced patient-level loops with vectorised Mata operations throughout
- **Data Processing**: All patient characteristics now processed as vectors for simultaneous operations
- **Matrix Operations**: Optimised matrix algebra for risk equation calculations
- **Code Structure**: Modular architecture with clear separation between setup, computation, and validation
- **Repository Organisation**: Modernised Git-based versioning (removed version folders)
- **File Naming**: 
  - `EpiMAP_Start.do` → `run.do` (clearer purpose)
  - Simplified main dispatcher naming
- **Documentation**: Updated README with vectorisation details and performance notes

### Performance
- **Execution Speed**: Significantly faster for large cohorts (10,000+ patients)
- **Memory Efficiency**: Reduced memory overhead through bulk operations
- **Scalability**: Better handling of bootstrap iterations and large-scale simulations

### Fixed
- Memory allocation issues in large simulation runs (now handled via vectorisation)
- Numerical precision edge cases through consistent vector operations
- Random number generation consistency across parallel operations

### Validation
- All vectorised outcomes validated to produce identical results to v2.0
- Comprehensive testing confirms bit-for-bit equivalence with loop-based implementation
- Extended validation suite for ongoing quality assurance

### Technical Details
- Implementation uses Mata's native vector and matrix operations
- Pre-allocated vectors for all patient characteristics (age, sex, ECOG, R-ISS, comorbidities)
- Matrix-based risk equation calculations with element-wise operations
- Consistent random number seeding for reproducibility

## [2.0] - 2025-01

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

## [1.0] - 2024-08

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
- Model calibrated against MRDR registry data (2009-2023)
- Survival curves validated against observed patient outcomes
- Treatment pathway distributions match registry patterns
- Out-of-sample validation with 70/30 split (2,884/1,237 patients)
- 100 bootstrap iterations for robustness testing
- No significant difference in 90% of 120 months post-diagnosis

---

## Migration Guides

### Upgrading from v2.0 to v2.1

**Good News**: v2.1 is fully backward compatible with v2.0 inputs!

**What's Different**:
- Internal implementation is vectorised, but all inputs/outputs remain the same
- No changes to data formats, coefficient files, or simulation parameters
- Results are identical (validated bit-for-bit equivalence)
- Significantly faster performance, especially for large cohorts

**Action Required**:
- None for basic usage - existing scripts will work unchanged
- Optional: Review new test suite in `tests/` for validation examples
- Optional: Update file references if using old naming (`EpiMAP_Start.do` → `run.do`)

**Performance Benefits**:
- ~2-5x faster for typical simulations (depends on cohort size)
- Better scaling for bootstrap analyses
- More efficient memory usage

### Upgrading from v1.0 to v2.1

**Breaking Changes**:
- Repository structure has changed (no more version folders)
- File paths for analyses updated
- Coefficient file organisation modified

**Action Required**:
1. Review new repository structure
2. Update file paths in custom scripts
3. Verify coefficient file locations in `analyses/*/data/coefficients/`
4. Test simulations with small cohort first

**New Features Available**:
- Extended treatment regimens
- Vectorised performance
- Comprehensive validation tools

## Release Notes

Each release includes:
- **Complete Model**: All necessary files to run simulations
- **Documentation**: User guides, technical specifications, and examples
- **Validation Results**: Benchmark tests and model performance metrics
- **Test Data**: Example datasets for testing and validation

## Contact

- **Model Questions**: adam.irving@monash.edu
- **Technical Issues**: Create an issue on GitHub
- **MRDR Data Access**: Visit mrdr.net.au
- **Collaboration**: Contact the research team via email

## Version Naming Convention

This project uses semantic versioning (MAJOR.MINOR format):
- **Major version** (1.0 → 2.0): Significant changes, may break compatibility
- **Minor version** (2.0 → 2.1): New features, improvements, backward compatible

Examples:
- v1.0: Initial release
- v2.0: Reorganised structure (breaking changes)
- v2.1: Vectorised implementation (backward compatible)
