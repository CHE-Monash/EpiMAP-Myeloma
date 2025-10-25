**********
	*EpiMAP Myeloma - Overall Survival Results
**********

capture confirm variable OC_TIME
if _rc {
	di as error "Run process_data.do first"
	exit 111
}

capture confirm existence $simulated_path
if _rc {
	global simulated_path "analyses/$Analysis/data/simulated"
}
capture mkdir "$simulated_path/results"
capture mkdir "$simulated_path/results/figures"


// Summary statistics
di _n "{bf:OVERALL SURVIVAL STATISTICS}"
di "{hline 60}"
quietly summarize OC_TIME, detail
di "N: " %9.0fc r(N)
di "Mean (SD): " %6.2f r(mean) " (" %5.2f r(sd) ")"
di "Median [IQR]: " %6.2f r(p50) " [" %6.2f r(p25) "-" %6.2f r(p75) "]"

// Survival at key time points
di _n "{bf:SURVIVAL AT KEY TIME POINTS}"
di "{hline 60}"
foreach year in 1 2 3 5 10 {
	quietly count if OC_TIME/12 >= `year'
	local pct = (r(N) / _N) * 100
	local se = sqrt(`pct' * (100 - `pct') / _N)
	local ci_l = max(0, `pct' - 1.96 * `se')
	local ci_u = min(100, `pct' + 1.96 * `se')
	di "`year'-year: " %5.1f `pct' "% (95% CI: " %5.1f `ci_l' "%-" %5.1f `ci_u' "%)"
}

// Overall K-M curve
preserve
stset OC_TIME if OC_TIME < 240, failure(OC_MORT)

sts graph, ///
	title("Overall Survival") ///
	xtitle("Months") ///
	ytitle("Probability") ///
	ylabel(0(0.2)1, angle(0) format(%3.1f)) ///
	xlabel(0(24)240) ///
	ci risktable legend(off) ///
	graphregion(color(white)) ///
	name(os_overall, replace)

graph export "$simulated_path/results/figures/os_overall.png", replace width(1200)
restore


// By BCR_L1 / ASCT
capture confirm variable BCR_L1
if !_rc {
	preserve
	gen bcr_group = .
	replace bcr_group = 1 if BCR_L1 == 1 | BCR_SCT == 1
	replace bcr_group = 2 if BCR_L1 == 2 | BCR_SCT == 2
	replace bcr_group = 3 if BCR_L1 == 3 | BCR_SCT == 3
	replace bcr_group = 4 if BCR_L1 == 4 | BCR_SCT == 4
	replace bcr_group = 5 if BCR_L1 == 5
	replace bcr_group = 6 if BCR_L1 == 6
	label define bcr_lbl 1 "≥VGPR" 2 "PR" 3 "≤MR"
	label values bcr_group bcr_lbl
	
	stset OC_TIME if OC_TIME < 240, failure(OC_MORT)
	
	di _n "{bf:SURVIVAL BY RESPONSE}"
	table bcr_group, stat(mean OC_TIME) stat(count OC_TIME)
	sts test bcr_group
	
	sts graph, by(bcr_group) ///
		title("OS by BCR L1 / ASCT") ///
		xtitle("Months") ///
		ytitle("Probability") ///
		ylabel(0(0.2)1, angle(0)) ///
		xlabel(0(24)240) ///
		graphregion(color(white)) ///
		legend(label(1 "CR") label(2 "VGPR") label(3 "PR") label(4 " MR") label(5 "SD") label(6 "PD") rows(2)) ///
		name(os_bcr, replace)
	
	graph export "$simulated_path/results/figures/os_bcr.png", replace width(1200)
	restore
}


// By ASCT
capture confirm variable SCT_L1
if !_rc {
	preserve
	gen asct = SCT_L1
	label define asct_lbl 0 "No ASCT" 1 "ASCT"
	label values asct asct_lbl
	
	stset OC_TIME if OC_TIME < 240, failure(OC_MORT)
	
	di _n "{bf:SURVIVAL BY ASCT}"
	table asct, stat(mean OC_TIME) stat(count OC_TIME)
	sts test asct
	
	sts graph, by(asct) ///
		title("OS by ASCT Status") ///
		xtitle("Months") ///
		ytitle("Probability") ///
		ylabel(0(0.2)1, angle(0)) ///
		xlabel(0(24)240) ///
		graphregion(color(white)) ///		
		risktable(, rowtitle("No ASCT") group(#1)) ///
		risktable(, rowtitle("ASCT") group(#2)) ///
		legend(label(1 "No ASCT") label(2 "ASCT")) ///
		name(os_asct, replace)
	
	graph export "$simulated_path/results/figures/os_asct.png", replace width(1200)
	restore
}


// By Age
preserve
gen age_group = .
replace age_group = 1 if Age < 65
replace age_group = 2 if Age >= 65 & Age < 75
replace age_group = 3 if Age >= 75 & Age < .
label define age_lbl 1 "<65" 2 "65-74" 3 "≥75"
label values age_group age_lbl

stset OC_TIME if OC_TIME < 240, failure(OC_MORT)

di _n "{bf:SURVIVAL BY AGE}"
table age_group, stat(mean OC_TIME) stat(count OC_TIME)
sts test age_group

sts graph, by(age_group) ///
	title("OS by Age at Diagnosis") ///
	xtitle("Months") ///
	ytitle("Probability") ///
	ylabel(0(0.2)1, angle(0)) ///
	xlabel(0(24)240) ///
	graphregion(color(white)) ///			
	risktable(, rowtitle("<65") group(#1)) ///
	risktable(, rowtitle("66 - 75") group(#2)) ///
	risktable(, rowtitle("75+") group(#3)) ///	
	legend(label(1 "<65") label(2 "66 - 75") label(3 "75+") rows(1)) ///	
	name(os_age, replace)

graph export "$simulated_path/results/figures/os_age.png", replace width(1200)
restore


// By R-ISS
capture confirm variable RISS
if !_rc {
	preserve
	label define riss_lbl 1 "R-ISS I" 2 "R-ISS II" 3 "R-ISS III"
	label values RISS riss_lbl
	
	stset OC_TIME if OC_TIME < 240, failure(OC_MORT)
	
	di _n "{bf:SURVIVAL BY R-ISS}"
	table RISS, stat(mean OC_TIME) stat(count OC_TIME)
	sts test RISS
	
	sts graph, by(RISS) ///
		title("OS by R-ISS Stage") ///
		xtitle("Months") ///
		ytitle("Probability") ///
		ylabel(0(0.2)1, angle(0)) ///
		xlabel(0(24)240) ///
		graphregion(color(white)) ///			
		risktable(, rowtitle("RISS-1") group(#1)) ///
		risktable(, rowtitle("RISS-2") group(#2)) ///
		risktable(, rowtitle("RISS-3") group(#3)) ///	
		legend(label(1 "RISS-1") label(2 "RISS-2") label(3 "RISS-3") rows(1)) ///	
		name(os_riss, replace)
	
	graph export "$simulated_path/results/figures/os_riss.png", replace width(1200)
	restore
}


// By ECOG
capture confirm variable ECOGcc
if !_rc {
	preserve
	label define ecog_lbl 1 "ECOG 0" 2 "ECOG 1" 3 "ECOG 2+"
	label values ECOGcc ecog_lbl
	
	stset OC_TIME if OC_TIME < 240, failure(OC_MORT)
	
	di _n "{bf:SURVIVAL BY ECOG}"
	table ECOGcc, stat(mean OC_TIME) stat(count OC_TIME)
	sts test ECOGcc
	
	sts graph, by(ECOGcc) ///
		title("OS by ECOG Score") ///
		xtitle("Months") ///
		ytitle("Probability") ///
		ylabel(0(0.2)1, angle(0)) ///
		xlabel(0(24)240) ///
		graphregion(color(white)) ///			
		risktable(, rowtitle("ECOG 0") group(#1)) ///
		risktable(, rowtitle("ECOG 1") group(#2)) ///
		risktable(, rowtitle("ECOG 2+") group(#3)) ///	
		legend(label(1 "ECOG 0") label(2 "ECOG 1") label(3 "ECOG 2+") rows(1)) ///	
		name(os_ecog, replace)
	
	graph export "$simulated_path/results/figures/os_ecog.png", replace width(1200)
	restore
}

di _n "Complete"
