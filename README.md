# EpiMAP Myeloma
Epidemiological Modelling of Australian Patients with Myeloma.
A Monash University collaboration between the Centre for Health Economics & the Transfusion Research Unit.

The EpiMAP Myeloma model is a discrete-event simulation model of the multiple myeloma disease and treatment pathway. It is based on a series of risk equations estimated using patient-level data from Monash University's Myeloma and Related Diseases Registry. 

This repository houses the three items needed to simulate the disease outcomes and treatment pathways of patients with multiple myeloma. The code is written in Stata and a valid Stata license is required to execute the simulation.

1. Risk equation coefficients (EpiMAP Coefficients.mmat)
2. Hypoethetical patient data file (EpiMAP Hypothetical Patients.dta)
3. Simulation code (EpiMAP Simulation.do and the separate SIM do files located in the Sub folder)

To run the model, download the folder for the version you require (e.g., V1) ensuring not to alter the structure of the folder, set your directory ('cd' command in Stata) to the downloaded folder and run the simulation code (EpiMAP Simulation.do) in Stata.

The diagnostic characteristics of 1,000 hypothetical multiple myeloma patients are included in this repository (EpiMAP Hypothetical Patients.dta) so that the simulation code can be executed and verified. Do not draw inference from the results of simulating these patients as the relationships between the randomly assigned patient characeristics may not reflect reality. If you are a researcher interested in acquiring genuine patient data from the Myeloma and Related Disease Registry you may submit an application to the Steering Committee. Full details can be found at https://www.mrdr.net.au/

Alternatively, if you are an institution with access to individual multiple myeloma patient data at diagnosis, you may construct a patient data file with the same structure as the patient data file provided (EpiMAP Hypothetical Patients.dta) and edit Line 24 of the simulation code (EpiMAP Simulation.do) to use the new patient data file.

Note that the model uses Stata's in-built matrix language Mata to store simulated data, which is converted into a Stata data file at the end of the simulation (EpiMAP Simulated.dta)

Please contact adam.irving@monash.edu with queries related to using the model

