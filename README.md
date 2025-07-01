# EpiMAP Myeloma

**Epidemiological Modelling of Australian Patients with Myeloma**

_A Monash University collaboration between the Centre for Health Economics & Transfusion Research Unit_

Research Team: Adam Irving, Dennis Petrie, Anthony Harris, Laura Fanning, Erica M Wood, Elizabeth Moore, Cameron Wellard, Neil Waters, Bradley Augustson, Gordon Cook, Francesca Gay, Georgia McCaughan, Peter Mollee, Andrew Spencer, Zoe K McQuilten

The EpiMAP Myeloma model is a discrete-event simulation model of multiple myeloma disease outcomes and treatment pathways. It is based on a series of risk equations estimated using patient-level data from Monash University's Australia and New Zealand Myeloma and Related Diseases Registry (MRDR).

This repository houses the three items needed to simulate the disease outcomes and treatment pathways of patients with multiple myeloma. The code is written in Stata and a valid Stata license is required to execute the simulation.

1. Risk equation coefficients (EpiMAP Coefficients.mmat)
2. Synthetic patient data file (EpiMAP Hypothetical Patients.dta)
3. Simulation code (EpiMAP Simulation.do and the separate SIM do files located in the Sub folder)

**Running the Model**

To run the model, download the repository and extract the V1 folder, ensuring not to alter the structure of files inside the folder. Set your directory ('cd' command in Stata) to the downloaded V1 folder and run the simulation code (EpiMAP Simulation.do) in Stata. The model uses Stata's in-built matrix language Mata to store simulated data, which is converted into a Stata data file at the end of the simulation (EpiMAP Simulated.dta).

The diagnostic characteristics of 1,000 hypothetical multiple myeloma patients are included in this repository (EpiMAP Hypothetical Patients.dta) so that the simulation code can be executed and verified. Do not draw inference from the results of simulating these patients as the relationships between the randomly assigned patient characeristics may not reflect reality. If you are a researcher interested in acquiring genuine patient data from the MRDR you may submit an application to the Steering Committee. Full details can be found at https://www.mrdr.net.au/

Alternatively, if you are an institution with access to individual multiple myeloma patient data at diagnosis including age, sex, ECOG and ISS, you may construct a patient data file with the same structure as the patient data file provided (EpiMAP Hypothetical Patients.dta) and edit Line 24 of the simulation code (EpiMAP Simulation.do) to use the new patient data file.

**Model Specification**

At diagnosis, MM patients are assigned four diagnostic characteristics - age, sex, Eastern Cooperative Oncology Group (ECOG) performance score and International Staging System (ISS) score. These diagnostic characteristics are included in all of the subsequent risk equations. The model uses Best Clinical Response (BCR) to treatment as a surrogate outcome for Overall Survival (OS).

![EpiMAP Myeloma - Model Framework](https://github.com/user-attachments/assets/3ae58966-9eb4-4a9a-b9b6-509fb8b1d952)

Blue health states represent treatments where BCR is predicted, either to chemotherapy or autologous stem cell transplant (ASCT). Grey health states represent treatment-free intervals or post-induction maintenance therapy. The EpiMAP Myeloma model (V1) considers up to nine LoTs and BCR is predicted after each treatment (chemotherapy or ASCT). The time between LoTs is defined as a treatment-free interval. Post-induction maintenance therapy is not expected to improve BCR; therefore, receiving maintenance therapy only increases the time between induction chemotherapy (LoT 1) and 1st subsequent chemotherapy (LoT 2). For LoT 1 and LoT 2 we modelled those specific chemotherapy regimens that were used by at least 10% of MRDR patients, alongside ‘other’ chemotherapy. Given the wide variety of combinations of chemotherapeutic agents recorded in the MRDR, for LoT 3 onwards specific chemotherapy regimens were not modelled, with all patients instead receiving the survival benefit of the average chemotherapy regimen. The proportion of simulated patients that receive each chemotherapy regimen is based on observed data from the MRDR.

| LoT | Regimen | Proportion | 
| --- | ------- | ---------- |
| 1 | VCd (bortezomib, cyclophosphamide & dexamethasone) | 58% |
| 1 | VRd (bortezomib, lenalidomide & dexamethasone) | 15% |
| 1 | Other | 26% |
| 2 | Rd (lenalidomide & dexamethasone) | 16% |
| 2 | DVd (daratumumab, bortezomib & dexamethasone) | 11% |
| 2 | Other | 73% |

It is not possible to simulate the impact of chemotherapy regimens not included in the model without restimating risk equations on individual patient-level data. Additionally, beyond LoT 3 BCR was collapsed into a 3-item scale due to low numbers in some of the categories (1 - Complete Remission or Very Good Partial Response, 2 - Partial Response or Minimal Response, and 3 - Stable Disease or Progressive Disease). LoTs 6-9 were assumed to be equivalent, sharing the same OS coefficients and risk equations for chemotherapy duration, treatment-free interval, and BCR.

**Discrete-event Simulation**

Patient characteristics of the hypothetical MM patients to be simulated were based on the patients from the MRDR. For each decision point in the model framework (e.g., after induction chemotherapy does this patient receive an ASCT?), the coefficients from the respective risk equation were used to calculate the likelihood p of each patient experiencing the outcome. This is then compared to a random number r drawn from a uniform distribution between 0 and 1 to determine the outcome (e.g., if p < r this patient receives an ASCT). As a discrete-event simulation model, all time periods were modelled explicitly for each patient. OS was predicted at the start of each health state with random numbers drawn from a uniform distribution between 0 and each patient’s current position on their specific OS curve. When assessing the time to competing events such as death versus end of chemotherapy duration, the model selected whichever event happened first. An age limit was also used to curtail patients who survival was estimated beyond 100 years old. As the MRDR contained some missing chemotherapy end dates, chemotherapy duration and treatment-free intervals are curtailed at the maximum observed in the data. 

**Risk Equations**

Several types of risk equations are required to simulate the disease, treatment, and outcome trajectory of patients through the EpiMAP Myeloma model framework. Parametric survival analysis is used to model time-to-event data, including OS, chemotherapy duration, and treatment-free intervals. The parametric distribution for these outcomes were selected based on information criteria and visual inspection of the goodness-of-fit. Logit regression is used to model binary outcomes – planned ASCT, receipt of ASCT, and receipt of maintenance therapy. Multinomial logit regression is used to predict chemotherapy regimen and ordered logit for BCR. These risk equations were all estimated using the patient-level MRDR data. Covariates included in the risk equations were chosen based on availability of data and clinical plausibility via consultation with our expert advisory group. OS depended on patient characteristics including BCR to each LoT. BCR depended on patient characteristics, chemotherapy regimen, and BCR to the previous LoT. 

Given significant differences in outcomes, separate ASCT and non-ASCT regressions were performed for LoT 1 chemotherapy duration and LoT 1 to LoT 2 treatment-free interval. Furthermore, induction chemotherapy duration for patients with planned ASCT utilised three manual splines to capture rapidly changes hazard rates and patients whose BCR to LoT 1 chemotherapy was Stable Disease or Progressive Disease were ineligible for ASCT. In total the EpiMAP Myeloma model contains 30 risk equations.

| No. | Outcome |
| --- | ------- | 
| 1 | Overall Survival | 
| 2 | Planned ASCT |
| 3 | Diagnosis Treatment-free Interval |
| 4 | LoT 1 Chemotherapy Regimen |
| 5 | LoT 1 Chemotherapy Duration ASCT Planned, Spline 1 (0-60 days) |
| 6 | LoT 1 Chemotherapy Duration ASCT Planned, Spline 2 (60-100 days) | 
| 7 | LoT 1 Chemotherapy Duration ASCT Planned, Spline 3 (≥100 days) |
| 8 | LoT 1 Chemotherapy Duration No ASCT Planned |
| 9 | LoT 1 Best Clinical Response |
| 10 | Receipt of ASCT |
| 11 | ASCT Best Clinical Response |
| 12 | Receipt of Maintenance Therapy | 
| 13 | LoT 1 Treatment-free Interval – ASCT Patients |
| 14 | LoT 1 Treatment-free Interval – No ASCT Patients |
| 15 | LoT 2 Chemotherapy Regimen |
| 16 | LoT 2 Chemotherapy Duration |
| 17 | LoT 2 Best Clinical Response |
| 18 | LoT 2 Treatment-free Interval | 
| 19 | LoT 3 Chemotherapy Duration |
| 20 | LoT 3 Best Clinical Response |
| 21 | LoT 3 Treatment-free Interval |
| 22 | LoT 4 Chemotherapy Duration |
| 23 | LoT 4 Best Clinical Response |
| 24 | LoT 4 Treatment-free Interval |
| 25 | LoT 5 Chemotherapy Duration |
| 26 | LoT 5 Best Clinical Response |
| 27 | LoT 5 Treatment-free Interval |
| 28 | LoT 6+ Chemotherapy Duration |
| 29 | LoT 6+ Best Clinical Response |
| 30 | LoT 6+ Treatment-free Interval |

By making this model open source we hope to begin collaborative research that can help improve outcomes for patients with multiple myeloma. Please contact adam.irving@monash.edu with queries related to using the model, or if you find any bugs in the code.

_The EpiMAP Myeloma project was directly supported by grant 1200706 to investigator Prof Zoe K. McQuilten. The EpiMAP Myeloma team thanks Andrew Marks for his involvement in the Expert Advisory Group as a consumer representative. The MRDR thanks patients, clinicians, and research staff at participating centres for their support._
