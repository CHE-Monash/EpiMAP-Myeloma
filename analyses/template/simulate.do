**********
* Monash Myeloma Model - Simulate (template dispatcher)
*
* Purpose: the standard simulation dispatcher -- a configuration block of globals, then the shared engine
*          pass (run_pipeline) run once per arm. Supports a two-arm comparison (intervention vs comparator)
*          or a single arm. COPY THIS FOLDER to analyses/<your_analysis>/, adapt, and fill in the
*          <...> placeholders in the Configuration block.
* Usage:   orchestrated by run.do; on the HPC it is sbatch'd directly (never sources run.do). Point
*          estimate: $boot 0. Bootstrap: $boot 1 with $min_bs/$max_bs over coefficient resamples.
*          Optional positional args (run.do / HPC arrays): boot min_bs max_bs [scenario].
* Notes:   Two arms are set by $int1 (intervention) and $int0 (comparator); leave $int0 "" for a single
*          arm. Both arms simulate the SAME cohort with the SAME coefficients -- they differ only via the
*          outcomes/sim_*_override.do files, which branch on $int (e.g. `if ("$int" == "<intervention>")`). A
*          comparator label with no override branch therefore gets the standard model (e.g. the natural
*          regimen blend). Worked examples -- analyses/default/ (single arm; start here);
*          analyses/transport_dvd/ (two arms + scenarios + outcome overrides).
**********

* Optional positional args, read into locals BEFORE clear all:
local a_boot  `1'
local a_minbs `2'
local a_maxbs `3'
local a_scen  `4'   // scenario override; delete if your analysis has no scenarios

clear all
set more off

if "$repo_path" != "" cd "$repo_path"   // cd to repo root only if config.do set it
capture run "config.do"     // machine-specific paths (git-ignored; see config.example.do)

**********
* Configuration  -- EDIT THESE
**********

// Analysis settings
global analysis     "template"          // <- your analysis name (this folder under analyses/)
global int1         "<intervention>"    // Intervention arm label (drives override branches; names output)
global int0         ""                  // Comparator arm label ("" = single arm; e.g. a standard-care blend)
global line         "0"                 // Line assessed (0 = all lines; 1-9 = a single line, e.g. for an override)
global coeffs       "template"          // Coefficient set -> coefficients_<coeffs>.mmat + outcomes/txr_<coeffs>.do
global data         "synthetic"         // Patient cohort (synthetic / train / test / predicted)
global min_year     "1995"              // Patients diagnosed from (>= 1995)
global max_year     "2040"              // Patients diagnosed until (<= 2040)
global min_id       "1"                 // First patient ID
global max_id       "999999"            // Last patient ID
global cost_year    "2025"              // Cost year (AUD)
global drate        "0.05"              // Annual discount rate (PBAC = 5%)
global report       "0"                 // Generate PDF report (0/1; single arm only)
global scenario     ""                  // Scenario label (optional; partitions outputs + selects overrides)

local single_arm = ("$int0" == "")

// Bootstrap settings
global boot         "0"                 // Bootstrap flag (0/1)
if "`a_boot'"  != "" global boot     "`a_boot'"
if "`a_minbs'" != "" global min_bs   "`a_minbs'"
if "`a_maxbs'" != "" global max_bs   "`a_maxbs'"
if "`a_scen'"  != "" global scenario "`a_scen'"   // run.do / HPC arrays loop scenarios via this arg

**********
* Set Paths
**********

global coefficients_path    "analyses/$analysis/coefficients"
global outcomes_path		"analyses/$analysis/outcomes"   // core/outcomes/sim_bcr.do & sim_txr.do
global patients_path        "analyses/$analysis/patients"   //   auto-run sim_bcr_override.do /
global simulated_path       "analyses/$analysis/simulated"  //   sim_txr_override.do from here (if present)

// Output partition for the simulated .dta: scenario is an optional subfolder
global sim_out = cond("$scenario" == "", "$simulated_path", "$simulated_path/$scenario")

// Cohort file -- set explicitly EVERY run: clear all does NOT clear globals, so a prior analysis's
// $cohort_file (e.g. default's outsample -> patients_test.dta) would otherwise persist and be loaded
// instead. Empty "" -> load_patients uses the convention path patients_${analysis}_${line}.dta.
global cohort_file  ""    // or point at a non-standard cohort: "$patients_path/<your_cohort>.dta"

**********
* Load Programs
**********

run "core/export_results.do"
run "core/run_pipeline.do"

**********
* Execute Simulation
**********

if ("$boot" == "0") {

// No Bootstrapping -- one pipeline pass per arm

    // Fail loudly if the coefficient set is missing -- the qui matuse below would otherwise
    // swallow the "file not found" error and the do-file would abort silently.
    capture confirm file "$coefficients_path/coefficients_${coeffs}.mmat"
    if _rc {
        di as error "Coefficients not found: $coefficients_path/coefficients_${coeffs}.mmat"
        di as error "  -> fit them (prep/risk_equations.do) or set \$coeffs to an existing set."
        exit 601
    }

    capture mkdir "$simulated_path"
    capture mkdir "$sim_out"

    // ---- Intervention arm ----
    global int "$int1"

        // Load coefficients
        qui mata: mata matuse "$coefficients_path/coefficients_$coeffs"

        // Execute pipeline
        run_pipeline

        // Save results
        save "$sim_out/${int}_${line}_${data}.dta", replace

        // Validate + report the single arm (a two-arm run is compared below, not validated to targets)
        if (`single_arm') {
            run "core/validation.do"
            if ("$report" == "1") qui do "core/generate_report.do"
        }

    // ---- Comparator arm ----
    if (!`single_arm') {

        global int "$int0"
        qui mata: mata clear              // reset Mata, then reload coefficients for this arm

        // Load coefficients
        qui mata: mata matuse "$coefficients_path/coefficients_$coeffs"

        // Execute pipeline
        run_pipeline

        // Save results
        save "$sim_out/${int}_${line}_${data}.dta", replace

        // Incremental summary (point estimate only; bootstrap gives the CI). Comparator (arm 0) is in
        // memory; append the saved intervention arm (arm 1) and compare mean discounted outcomes.
        gen arm = 0
        preserve
            use "$sim_out/${int1}_${line}_${data}.dta", clear
            gen arm = 1
            tempfile intervention
            qui save `intervention'
        restore
        append using `intervention'

        collapse (mean) cost=cost_total_d qaly=qaly_total_d ly=OC_TIME, by(arm)
        // OC_TIME is in MONTHS and undiscounted; cost_total_d / qaly_total_d are discounted years.
        // Divide LY by 12 for a years readout (it stays undiscounted -- descriptive, the ICER uses QALYs).
        local dcost = cost[2] - cost[1]
        local dqaly = qaly[2] - qaly[1]
        local dly   = (ly[2]  - ly[1]) / 12
        di as text "Inc Cost (disc):  $" as result %9.0fc `dcost'
        di as text "Inc LY (undisc):  " as result %6.2f `dly'
        di as text "Inc QALY (disc):  " as result %6.2f `dqaly'
        di as text "ICER (disc):      $" as result %9.0fc (`dcost' / `dqaly')

        // Two-arm landscape comparison report (outcomes side-by-side)
        if ("$report" == "1") {
            global report_twoarm "1"
            qui do "core/generate_report.do"
            global report_twoarm ""
        }
    }
}
else {

// Bootstrapping -- one simulated dataset per coefficient resample, per arm

    capture mkdir "$simulated_path"
    capture mkdir "$sim_out"
    capture mkdir "$sim_out/bootstrap"

    forvalues b = $min_bs/$max_bs {
    global b "`b'"

        // ---- Intervention arm ----
        global int "$int1"
        mata: mata clear

            // Load coefficients
            qui mata: mata matuse "$coefficients_path/bootstrap/coefficients_${coeffs}_B`b'"

            // Execute pipeline
            run_pipeline

            // Save results
            save "$sim_out/bootstrap/${int}_${line}_${data}_B`b'.dta", replace

        // ---- Comparator arm ----
        if (!`single_arm') {
            global int "$int0"
            mata: mata clear

                // Load coefficients
                qui mata: mata matuse "$coefficients_path/bootstrap/coefficients_${coeffs}_B`b'"

                // Execute pipeline
                run_pipeline

                // Save results
                save "$sim_out/bootstrap/${int}_${line}_${data}_B`b'.dta", replace
        }

        di as text "Iteration `b' completed"
    }
}
