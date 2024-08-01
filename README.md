# EpiMAP Myeloma

**Epidemiological Modelling of Australian Patients with Myeloma**

A Monash University collaboration between the Centre for Health Economics & the Transfusion Research Unit.

The EpiMAP Myeloma model is a discrete-event simulation model of  multiple myeloma disease outcomes and treatment pathways. It is based on a series of risk equations estimated using patient-level data from Monash University's Australia and New Zealand Myeloma and Related Diseases Registry (MRDR). An open-access peer-reviewed manuscript has been accepted for publication which contains additional details regarding the structure and assumptions used to build the model, once published this repository will be updated to include the publication.

This repository houses the three items needed to simulate the disease outcomes and treatment pathways of patients with multiple myeloma. The code is written in Stata and a valid Stata license is required to execute the simulation.

1. Risk equation coefficients (EpiMAP Coefficients.mmat)
2. Hypoethetical patient data file (EpiMAP Hypothetical Patients.dta)
3. Simulation code (EpiMAP Simulation.do and the separate SIM do files located in the Sub folder)

**Running the Model**

To run the model, download the folder for the version you require (e.g., V1) ensuring not to alter the structure of the folder. Set your directory ('cd' command in Stata) to the downloaded folder and run the simulation code (EpiMAP Simulation.do) in Stata. The model uses Stata's in-built matrix language Mata to store simulated data, which is converted into a Stata data file at the end of the simulation (EpiMAP Simulated.dta).

The diagnostic characteristics of 1,000 hypothetical multiple myeloma patients are included in this repository (EpiMAP Hypothetical Patients.dta) so that the simulation code can be executed and verified. Do not draw inference from the results of simulating these patients as the relationships between the randomly assigned patient characeristics may not reflect reality. If you are a researcher interested in acquiring genuine patient data from the MRDR you may submit an application to the Steering Committee. Full details can be found at https://www.mrdr.net.au/

Alternatively, if you are an institution with access to individual multiple myeloma patient data at diagnosis including age, sex, ECOG and ISS, you may construct a patient data file with the same structure as the patient data file provided (EpiMAP Hypothetical Patients.dta) and edit Line 24 of the simulation code (EpiMAP Simulation.do) to use the new patient data file.

**Model Specification**

At diagnosis, MM patients are assigned four diagnostic characteristics - age, sex, Eastern Cooperative Oncology Group (ECOG) performance score and International Staging System (ISS) score. These diagnostic characteristics are included in all of the subsequent risk equations. The model uses Best Clinical Response (BCR) to treatment as a surrogate outcome for Overall Survival (OS).

![EpiMAP Myeloma - Model Framework](https://github.com/user-attachments/assets/3ae58966-9eb4-4a9a-b9b6-509fb8b1d952)

Blue health states represent treatments where BCR is predicted, either to chemotherapy or autologous stem cell transplant (ASCT). Grey health states represent treatment-free intervals or post-induction maintenance therapy. The EpiMAP Myeloma model (V1) considers up to nine LoTs and BCR is predicted after each treatment (chemotherapy or ASCT). The time between LoTs is defined as a treatment-free interval.

Given the wide variety of chemotherapy regimens used to treat myeloma in Australia, the EpiMAP Myeloma model includes coefficients for specific chemotherapy regimens at LoT 1 and LoT 2 only, alongside the 'other' category. The proportion of simulated patients that receive each chemotherapy regimen is based on observed data from the MRDR.

| LoT | Regimen | Proportion | 
| --- | ------- | ---------- |
| 1 | VCd (bortezomib, cyclophosphamide & dexamethasone) | 58% |
| 1 | VRd (bortezomib, lenalidomide & dexamethasone) | 15% |
| 1 | Other | 26% |
| 2 | Rd (lenalidomide & dexamethasone) | 16% |
| 2 | DVd (daratumumab, bortezomib & dexamethasone) | 11% |
| 2 | Other | 73% |

From LoT 3 onwards, all patients receive the 'other' category. It is not possible to simulate the impact of chemotherapy regimens not included in the model without restimating risk equations on individual patient-level data.

Please contact adam.irving@monash.edu with queries related to using the model, or if you find any bugs in the code.

**Risk Equations**

The EpiMAP Myeloma model contains 28 risk equations estimated using the Myeloma and Related Diseases Registry:

| No. | Outcome |
| --- | ------- | 
| 1 | Diagnosis Treatment-free Interval |
| 2 | LoT 1 Chemotherapy Regimen |
| 3 | LoT 1 Chemotherapy Duration ASCT Planned, Spline 1 (0-60 days) |
| 4 | LoT 1 Chemotherapy Duration ASCT Planned, Spline 2 (60-100 days) | 
| 5 | LoT 1 Chemotherapy Duration ASCT Planned, Spline 3 (≥100 days) |
| 6 | LoT 1 Chemotherapy Duration No ASCT Planned |
| 7 | LoT 1 Best Clinical Response |
| 8 | Receipt of ASCT |
| 9 | ASCT Best Clinical Response |
| 10 | Receipt of Maintenance Therapy | 
| 11 | LoT 1 Treatment-free Interval – ASCT Patients |
| 12 | LoT 1 Treatment-free Interval – No ASCT Patients |
| 13 | LoT 2 Chemotherapy Regimen |
| 14 | LoT 2 Chemotherapy Duration |
| 15 | LoT 2 Best Clinical Response |
| 16 | LoT 2 Treatment-free Interval | 
| 17 | LoT 3 Chemotherapy Duration |
| 18 | LoT 3 Best Clinical Response |
| 19 | LoT 3 Treatment-free Interval |
| 20 | LoT 4 Chemotherapy Duration |
| 21 | LoT 4 Best Clinical Response |
| 22 | LoT 4 Treatment-free Interval |
| 23 | LoT 5 Chemotherapy Duration |
| 24 | LoT 5 Best Clinical Response |
| 25 | LoT 5 Treatment-free Interval |
| 26 | LoT 6+ Chemotherapy Duration |
| 27 | LoT 6+ Best Clinical Response |
| 28 | LoT 6+ Treatment-free Interval |
