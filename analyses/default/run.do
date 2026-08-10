**********
* Monash Myeloma Model - default runbook
*
* Purpose: The reference analysis, in two modes (simulate.do $scenario select):
*            PROJECTION  ($scenario "")          - full-registry fit (coeffs=full) on the synthetic
*                                                  incidence population (data=synthetic); costed, reported.
*            OUT-OF-SAMPLE ($scenario "outsample")- 70%-train fit (coeffs=train) on the held-out real 30%
*                                                  (data=test); validated against observed targets.
* Notes:   DETERMINISTIC (point estimate) tracks are light and run locally, top to bottom. BOOTSTRAP
*          (prediction intervals) is heavy and runs on the HPC -- this file is not run there; the
*          bootstrap sections are the canonical record (commented), with the sbatch + rsync/submit/pull
*          plumbing in the /* ... */ shell block at the end. The 70/30 split and per-fold imputation need
*          the restricted MRDR drive ($data_path via config.do); the rest uses the imputed outputs.
**********

if "$repo_path" != "" cd "$repo_path"
capture run "config.do"

**********
* PROJECTION -- full fit x synthetic population (the primary use)
**********

* P0. Multiple imputation on the FULL registry -> ${data_path}/MRDR Long MI.dta  (needs MRDR drive).
*     Re-run whenever the extraction changes a kept variable (e.g. MND_L1, refr_len_*); the risk
*     equations and in-sample benchmarks below read this file. args: imp boot min_bs max_bs sample
*     (empty 5th arg = full cohort; also clears any $sample left over from an OOS run).
do "prep/multiple_imputation.do" 10 0 . .

* P1. Risk equations on the FULL registry (100%) -> analyses/default/coefficients/coefficients_full.mmat
*     args: analysis coeffs min_year max_year boot min_bs max_bs   (loads outcomes/txr_full.do)
do "prep/risk_equations.do" default full 1995 2040 0

* P2. Synthetic incidence population cohort(s) -> patients/synthetic_1995_2040_*.dta  (needs MRDR drive)
*do "prep/synthetic_1995_2040.do"

* P3. Simulate the synthetic population (projection; costed PDF report) -> simulated/all_0_synthetic.dta
do "analyses/default/simulate.do"

* P4. In-sample benchmarks (full-registry observed outcomes) -> scratch/benchmarks/  (needs MRDR drive)
do "prep/generate_benchmarks.do"

* P5. In-sample FACE VALIDITY: synthetic projection vs the in-sample benchmarks (NOT the mainstay
*     validation -- synthetic cohort vs whole-registry, not held out; see validate_insample.do header)
do "analyses/default/validate_insample.do"


**********
* OUT-OF-SAMPLE VALIDATION -- 70%-train fit x held-out 30% (the mainstay validation)
* DETERMINISTIC (point estimate) -- runs locally, top to bottom
**********

* O0. Split patients 70/30 -> ${data_path}/oos/oos_split.dta   (run ONCE, fixed seed; needs MRDR drive)
do "analyses/default/prep/split.do"

* O1. Multiple imputation, each fold imputed SEPARATELY (split BEFORE imputation = no leakage)
*     args: imp boot min_bs max_bs sample
do "prep/multiple_imputation.do" 2 0 . . train         // train (70%) -> coefficients
do "prep/multiple_imputation.do" 2 0 . . test          // test  (30%) -> targets + cohort

* O2. Risk equations on the TRAIN fold -> analyses/default/coefficients/coefficients_train.mmat
*     (loads outcomes/txr_train.do, which sources the canonical txr_full.do)
do "prep/risk_equations.do" default train 1995 2040 0 . . train

* O3. Held-out 30% validation targets (observed outcomes) -> analyses/default/targets/ (18 csv)
do "analyses/default/prep/test_targets.do"

* O4. Held-out 30% simulation cohort -> analyses/default/patients/patients_test.dta
do "analyses/default/prep/test_cohort.do"

* O5. Simulate the 30% with the 70%-trained coefficients ($boot 0, $scenario outsample)
*     -> analyses/default/simulated/outsample/all_0_test.dta
do "analyses/default/simulate.do" 0 . . outsample

* O6. Validate the point estimate vs the observed targets (fixed tolerances)
do "analyses/default/validate_outsample.do"


**********
* BOOTSTRAP (prediction intervals) -- HEAVY; runs on the HPC. The `do ...` lines are the LOCAL (serial,
* slow) equivalent; the actual HPC plumbing is the /* ... */ shell block at the END of this file (a Stata
* block comment, inert here) -- SELECT-AND-RUN those lines in a VS Code TERMINAL (bash). The slurm scripts
* (hpc/*.script) are generic array jobs 1-500; per-analysis args ride on the `sbatch --export=` lines.
**********

*  Projection PIs (optional): bootstrap the full-fit coefficients, then 500 projection sims.
* do "prep/risk_equations.do" default full 1995 2040 1 1 500
* do "analyses/default/simulate.do" 1 1 500

*  Out-of-sample PIs (the headline validation metric):
*  (a) TRAIN bootstrap: MI + risk equations -> coefficients/bootstrap/coefficients_train_B1..500
* do "prep/multiple_imputation.do" 2 1 1 500 train
* do "prep/risk_equations.do" default train 1995 2040 1 1 500 train
*  (b) 500 bootstrap sims of the held-out 30% -> simulated/outsample/bootstrap/
* do "analyses/default/simulate.do" 1 1 500 outsample
*  (c) Percentile-method PI coverage vs the held-out observed -> analyses/default/results/
* do "analyses/default/bootstrap_validation.do" 500


/* ================================================================================================
   HPC plumbing (out-of-sample bootstrap) -- RUN BY HAND IN A VS CODE TERMINAL (bash), NOT in Stata.
   (Stata block comment, ignored on `do`; select the lines you need and run them in a terminal.)
   Prereq (on the Mac, needs the drive): do analyses/default/prep/split.do -> ${data_path}/oos/oos_split.dta

   # FIRST, load the machine-specific paths (git-ignored; sets user/repo_path/hpc/drive_path):
   source env.sh
   analysis=default

   # Directories on hpc (mirror the repo tree + the data tree)
   ssh $hpc "mkdir -p em76/$user/data/$data/oos/bootstrap"
   ssh $hpc "mkdir -p em76/$user/prep/inputs"
   ssh $hpc "mkdir -p em76/$user/hpc"
   ssh $hpc "mkdir -p em76/$user/analyses/$analysis/outcomes"
   ssh $hpc "mkdir -p em76/$user/analyses/$analysis/coefficients/bootstrap"
   ssh $hpc "mkdir -p em76/$user/analyses/$analysis/patients"
   ssh $hpc "mkdir -p em76/$user/analyses/$analysis/simulated/outsample/bootstrap"
   ssh $hpc "mkdir -p em76/$user/analyses/$analysis/results"

   # Config + prep/engine code + the OOS regimen definition risk_equations loads
   rsync -auvzce ssh "$repo_path"/hpc/config.do $hpc:~/em76/$user/config.do
   rsync -auvzce ssh "$repo_path"/prep/multiple_imputation.do $hpc:~/em76/$user/prep/
   rsync -auvzce ssh "$repo_path"/prep/risk_equations.do $hpc:~/em76/$user/prep/
   rsync -auvzce ssh "$repo_path"/analyses/$analysis/outcomes/ $hpc:~/em76/$user/analyses/$analysis/outcomes/

   # Data: MRDR Long + the 70/30 split crosswalk
   rsync -auvzce ssh "$drive_path/data/$data/MRDR Long.dta" $hpc:~/em76/$user/data/$data/
   rsync -auvzce ssh "$drive_path/data/$data/oos/oos_split.dta" $hpc:~/em76/$user/data/$data/oos/

   # Slurm scripts -> hpc/ (generic array jobs 1-500; per-analysis args via --export at submit)
   rsync -auvzce ssh "$repo_path"/hpc/multiple_imputation.script $hpc:~/em76/$user/hpc/
   rsync -auvzce ssh "$repo_path"/hpc/risk_equations.script $hpc:~/em76/$user/hpc/
   rsync -auvzce ssh "$repo_path"/hpc/simulate.script $hpc:~/em76/$user/hpc/

   # == Step (a): TRAIN-fold MI, then TRAIN-fold risk equations. Both are 500-job arrays and the
   #    second reads the first's B<b> datasets, so the order is fixed. Submitted SEPARATELY below so
   #    either can be run on its own; the chained one-liner that does both is at the end of the step.
   #    NOTE this is the OOS TRAIN fold (IMP=2, SAMPLE=train). The transport_dvd / car_t bootstraps
   #    need the MAIN-MODEL MI instead -- full data, no fold, IMP=10 -- which is step (a) in
   #    analyses/transport_dvd/run.do. Running this one will not satisfy those.

   # (a1) TRAIN-fold bootstrap MI -- 500-job array. --parsable prints just the job id; note it as
   #      MIJOB if you want to queue (a2) behind it instead of waiting.
   ssh $hpc "cd em76/$user ; sbatch --parsable --mail-user=$hpc_email --export=ALL,IMP=2,SAMPLE=train hpc/multiple_imputation.script"

   #      Watch it finish before submitting (a2):
   ssh $hpc "squeue -u $user"

   # (a2) TRAIN-fold risk equations -- 500-job array. Run ONLY once (a1) has finished.
   ssh $hpc "cd em76/$user ; sbatch --mail-user=$hpc_email --export=ALL,ANALYSIS=default,COEFFS=train,MINYR=1995,MAXYR=2040,SAMPLE=train hpc/risk_equations.script"

   #      Or, to queue (a2) behind (a1) without waiting, using the id (a1) printed:
   # ssh $hpc "cd em76/$user ; sbatch --mail-user=$hpc_email --dependency=afterok:<MIJOB> --export=ALL,ANALYSIS=default,COEFFS=train,MINYR=1995,MAXYR=2040,SAMPLE=train hpc/risk_equations.script"

   #      Or both in one go, chained (the original single-command form):
   # ssh $hpc "cd em76/$user ; mi=\$(sbatch --parsable --mail-user=$hpc_email --export=ALL,IMP=2,SAMPLE=train hpc/multiple_imputation.script) ; echo MI job \$mi ; sbatch --mail-user=$hpc_email --dependency=afterok:\$mi --export=ALL,ANALYSIS=default,COEFFS=train,MINYR=1995,MAXYR=2040,SAMPLE=train hpc/risk_equations.script"

   # Pull the bootstrap coefficients back (optional; already on hpc for step b)
   rsync -auvzce ssh $hpc:~/em76/$user/analyses/$analysis/coefficients/bootstrap/ "$repo_path"/analyses/$analysis/coefficients/bootstrap/

   # == Step (b): 500 bootstrap simulations of the held-out 30% (after step a's coefficients are on hpc) ==
   # Needs core/ + simulate.do + cost inputs + the test cohort, which the MI/risk phase did not send:
   rsync -auvzce ssh "$repo_path"/core/ $hpc:~/em76/$user/core/
   rsync -auvzce ssh "$repo_path"/prep/inputs/ $hpc:~/em76/$user/prep/inputs/
   rsync -auvzce ssh "$repo_path"/analyses/$analysis/simulate.do $hpc:~/em76/$user/analyses/$analysis/
   rsync -auvzce ssh "$repo_path"/analyses/$analysis/patients/patients_test.dta $hpc:~/em76/$user/analyses/$analysis/patients/

   # Submit the simulation array (simulate.do 1 <task> <task> outsample; --export picks the analysis):
   ssh $hpc "cd em76/$user ; sbatch --mail-user=$hpc_email --export=ALL,ANALYSIS=default,SCENARIO=outsample hpc/simulate.script"

   # == Step (c): PI coverage. Aggregate the 500 sims on hpc; pull only results/. ==
   rsync -auvzce ssh "$repo_path"/analyses/$analysis/targets/ $hpc:~/em76/$user/analyses/$analysis/targets/
   rsync -auvzce ssh "$repo_path"/analyses/$analysis/bootstrap_validation.do $hpc:~/em76/$user/analyses/$analysis/
   ssh $hpc "cd em76/$user ; sbatch --mail-user=$hpc_email --time=0-00:30:00 --wrap=\"module load stata/16u2021 && stata-mp analyses/$analysis/bootstrap_validation.do 500\""
   rsync -auvzce ssh $hpc:~/em76/$user/analyses/$analysis/results/ "$repo_path"/analyses/$analysis/results/
   ================================================================================================ */
