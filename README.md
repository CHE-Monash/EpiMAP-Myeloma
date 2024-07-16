# EpiMAP Myeloma
Epidemiological Modelling of Australian Patients with Myeloma. A Monash University collaboration between the Centre for Health Economics & the Transfusion Research Unit.

This repository houses the three items needed to simulate the disease outcomes and treatment pathways of Australian patients with multiple myeloma.

1. Risk equation coefficients generated using the Australia & New Zealand Myeloma and Related Diseases Registy (EpiMAP Coefficients.mmat)
2. Multiply imputed patient characteristics at diagnosis data (EpiMAP MI Patients.dta)
3. Simulation Stata code (EpiMAP Simulation.do and the separate do files located in the Sub folder)

To run the model, download the folder for the version you require (e.g., V1) ensuring not to alter the structure of the folder and run "EpiMAP Simulation.do" in Stata.

Note that the model uses Stata's in-built matrix language Mata to store simulated data, which is converted into a Stata data file at the end of the simulation (EpiMAP Simulated.dta)
