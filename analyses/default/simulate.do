**********
* Monash Myeloma Model - Simulate (default dispatcher)
*
* Purpose: the reference analysis. One dispatcher, two modes selected by $scenario:
*            $scenario ""          -> PROJECTION: full-registry fit (coeffs=full) applied to the
*                                     synthetic incidence population (data=synthetic); costed, reported.
*            $scenario "outsample" -> OUT-OF-SAMPLE VALIDATION: 70%-train fit (coeffs=train) applied to
*                                     the held-out real 30% (data=test); compared to observed targets by
*                                     validate_outsample.do.
* Usage:   orchestrated by run.do; on the HPC it is sbatch'd directly (never sources run.do). Point
*          estimate: $boot 0. Bootstrap: $boot 1 with $min_bs/$max_bs over the coefficient resamples.
*          Optional positional args (run.do / HPC arrays): boot min_bs max_bs [scenario].
**********

* Optional positional args, read into locals BEFORE clear all:
local a_boot  `1'
local a_minbs `2'
local a_maxbs `3'
local a_scen  `4'   // scenario override ("" projection / "outsample"); default set below

clear all
set more off

if "$repo_path" != "" cd "$repo_path"   // cd to repo root only if config.do set it
capture run "config.do"     // machine-specific paths (git-ignored; see config.example.do)

**********
* Configuration
**********

// Common settings (shared by both scenarios)
global analysis     "default"           // Analysis name
global int          "all"               // Intervention
global line         "0"                 // Line being assessed (0-9)
global min_year     "1995"              // Patients diagnosed from (>= 1995)
global max_year     "2040"              // Patients diagnosed until (<= 2040)
global min_id       "1"                 // First patient ID (>= 1)
global max_id       "101212"            // Last patient ID
global drate        "0.05"              // Annual discount rate (PBAC = 5%)

// Scenario select (positional arg 4 overrides; default "" = projection)
global scenario     ""
if "`a_scen'" != "" global scenario "`a_scen'"

// Scenario presets: coeffs (fit) x data (cohort) + cost year + report
if ("$scenario" == "outsample") {
    global coeffs       "train"                                     // fit on the 70% training fold
    global data         "test"                                      // simulate the held-out real 30%
    global cohort_file  "analyses/$analysis/patients/patients_test.dta"
    global cost_year    "2025"                                      // Cost year (AUD)
    global report       "0"                                         // no PDF report for validation
}
else {
    global coeffs       "full"                                      // fit on the full (100%) registry
    global data         "synthetic"                                 // simulate the synthetic incidence population
    global cost_year    "2025"                                      // Cost year (AUD); falls back to latest treatment_costs_*.csv
    global report       ""                                         // costed PDF projection report
}

// Bootstrap settings
global boot         "0"                 // Bootstrap flag (0/1)
if "`a_boot'"  != "" global boot   "`a_boot'"
if "`a_minbs'" != "" global min_bs "`a_minbs'"
if "`a_maxbs'" != "" global max_bs "`a_maxbs'"

**********
* Set Paths
**********

global coefficients_path    "analyses/$analysis/coefficients"
global outcomes_path        "analyses/$analysis/outcomes"
global patients_path        "analyses/$analysis/patients"
global simulated_path       "analyses/$analysis/simulated"

// Output partition for the simulated .dta: scenario is an optional subfolder
global sim_out = cond("$scenario" == "", "$simulated_path", "$simulated_path/$scenario")

**********
* Load Programs
**********

run "core/export_results.do"
run "core/run_pipeline.do"

**********
* Execute Simulation
**********

if ("$boot" == "0") {

// No Bootstrapping

    // Load coefficients
    qui mata: mata matuse "$coefficients_path/coefficients_$coeffs"

    // Execute pipeline
    run_pipeline

    // Save results
    capture mkdir "$sim_out"
    save "$sim_out/${int}_${line}_${data}.dta", replace

    // Validate results
    run "core/validation.do"

    // Generate report (projection only)
    if ("$report" == "1") qui do "core/generate_report.do"

}
else {

    // Bootstrapping: one simulated dataset per coefficient resample
    forvalues b = $min_bs/$max_bs {
    global BSIteration "`b'"
    mata: mata clear

        // Load coefficients
        qui mata: mata matuse "$coefficients_path/bootstrap/coefficients_${coeffs}_B`b'"

        // Execute pipeline
        run_pipeline

        // Save results
        capture mkdir "$sim_out"
        capture mkdir "$sim_out/bootstrap"
        save "$sim_out/bootstrap/${int}_${line}_${data}_B`b'.dta", replace

        di as text "Iteration `b' completed"
    }
}
