// Test 1: Simple population simulation
cd "/Users/adami/Documents/Monash/Research/Blood Disorders/Myeloma/EpiMAP/Github/"

// Define settings
local analysis      "base_model"    // Analysis folder name
local intervention  "all"           // Treatment intervention
local line          "0"             // Start from diagnosis
local coefficients  "base_model"    // Coefficient set to use
local data_type     "population"    // Use population data
local min_id        "1"             // First patient ID
local max_id        "1000"            // Last patient ID  
local bootstrap     "0"             // No bootstrap
local min_bootstrap ""             // (Not used when bootstrap=0)
local max_bootstrap ""             // (Not used when bootstrap=0)

// Execute simulation
do "EpiMAP_Myeloma_v2.0.do" ///
    `analysis' `intervention' `line' `coefficients' `data_type' ///
    `min_id' `max_id' `bootstrap' `min_bootstrap' `max_bootstrap'
