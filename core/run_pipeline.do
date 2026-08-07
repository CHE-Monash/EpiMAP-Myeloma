**********
* Monash Myeloma Model - Run Pipeline (shared lean engine pass)
*
* Purpose: The single definition of the core simulation pass (load_patients -> mata_setup
*          -> simulation -> process_data) used by every orchestrator, WITHOUT CSV export
*          (callers run core/export_results.do separately).
* Notes:   This file defines the four stage programs as well as run_pipeline, so a caller needs
*          `run "core/run_pipeline.do"' alone. The stage files are run HERE, outside the program
*          body, so they are parsed once rather than on every run_pipeline call.
*          Caller must still load coefficients (mata matuse), set the usual globals ($data, $int,
*          $line, $scenario, $coeffs, $boot, paths), and run core/export_results.do separately.
**********

run "core/load_patients.do"
run "core/mata_setup.do"
run "core/simulation_engine.do"
run "core/process_data.do"

capture program drop run_pipeline
program define run_pipeline
	// Compile the persistent Mata utility functions once per (cleared) Mata
	// state. The definition files error ("... already exists") if re-run while
	// their functions are still compiled, so load them only when absent: true on
	// a fresh session and again after the bootstrap loop's `mata clear`, but a
	// no-op if run_pipeline is called repeatedly within one live Mata state.
	capture mata: _rp_probe = &get_txr_coef()
	if (_rc) {
		run "core/mata_functions.do"
		run "core/rng_slots.do"
	}
	load_patients
	mata_setup
	simulation
	process_data
end
