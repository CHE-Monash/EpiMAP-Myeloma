****************************************************************************
* Generate Baseline Characteristics Report
* 
* Purpose: Create PDF report of baseline characteristics using putpdf
*          Integrates with EpiMAP v2.0 dispatcher system
* 
* Requirements:
*   - Stata 15+ (putpdf built-in)
*   - Global variables from dispatcher:
*       $Analysis, $Int, $Line, $Data, $MinID, $MaxID, $simulated_path
* 
* Output:
*   - baseline_characteristics_[Analysis]_[Int]_[Data].pdf
*
* Author: EpiMAP Research Team
* Date: October 2025
****************************************************************************

di as text _n(2) "{hline 78}"
di as text "{bf:EpiMAP Myeloma v2.0 - Simulation Report}"
di as text "{hline 78}"

****************************************************************************
* VALIDATE REQUIRED GLOBALS
****************************************************************************

capture confirm existence $Analysis
if _rc {
    di as error "Error: Global variables not set"
    di as error "This script must be called after running EpiMAP_Myeloma_v2.0.do"
    exit 198
}

di as text _n "Simulation Parameters:"
di as text "  Analysis:     " as result "$Analysis"
di as text "  Intervention: " as result "$Int"
di as text "  Line:         " as result "$Line"  
di as text "  Data:         " as result "$Data"
di as text "  Patient IDs:  " as result "$MinID to $MaxID"

****************************************************************************
* CREATE OUTPUT DIRECTORY
****************************************************************************

capture mkdir "$simulated_path/reports"
local report_dir "$simulated_path/reports"

di as text "  Output:       " as result "`report_dir'"

****************************************************************************
* LOAD DATA
****************************************************************************

local datafile "$simulated_path/$Int $Line $Data $MinID $MaxID.dta"

capture confirm file "`datafile'"
if _rc {
    di as error _n "Error: Simulated data file not found:"
    di as error "  `datafile'"
    exit 601
}

quietly use "`datafile'", clear
di as result "✓ Data loaded: " as text _N " patients"

****************************************************************************
* Start PDF
****************************************************************************

capture putpdf clear
set graphics off
putpdf begin

// Title page
putpdf paragraph, halign(center)
putpdf text ("EpiMAP Myeloma v2.0"), bold font(,18)

putpdf paragraph
putpdf text ("Simulation Report"), bold font(,16)

putpdf table settings = (7, 2), border(all)
putpdf table settings(1,1) = ("Setting"), bold
putpdf table settings(1,2) = ("Value"), bold
putpdf table settings(2,1) = ("Analysis")
putpdf table settings(2,2) = ("$Analysis")
putpdf table settings(3,1) = ("Intervention")
putpdf table settings(3,2) = ("$Int")
putpdf table settings(4,1) = ("Line")
putpdf table settings(4,2) = ("$Line")
putpdf table settings(5,1) = ("Data")
putpdf table settings(5,2) = ("$Data")
putpdf table settings(6,1) = ("Patient IDs")
putpdf table settings(6,2) = ("$MinID to $MaxID")
putpdf table settings(7,1) = ("Report Date")
putpdf table settings(7,2) = ("`c(current_date)'")

****************
* Demographics
****************

putpdf paragraph
putpdf text ("Demographics"), bold font(,16)

// Sample size
quietly count
local total_n = string(r(N), "%9.0fc")

putpdf paragraph
putpdf text ("Total Sample Size: `total_n'"), bold

// Age
quietly summarize Age
local mean_age = string(r(mean), "%4.1f")
local sd_age = string(r(sd), "%4.1f")
local min_age = string(r(min), "%4.0f")
local max_age = string(r(max), "%4.0f")

quietly summarize Age, detail
local p25_age = string(r(p25), "%4.1f")
local p50_age = string(r(p50), "%4.1f")
local p75_age = string(r(p75), "%4.1f")

putpdf paragraph
putpdf text ("Age Distribution"), bold

putpdf table age_tbl = (7, 2), border(all)
putpdf table age_tbl(1,1) = ("Statistic"), bold
putpdf table age_tbl(1,2) = ("Value (years)"), bold
putpdf table age_tbl(2,1) = ("Mean (SD)")
putpdf table age_tbl(2,2) = ("`mean_age' (`sd_age')")
putpdf table age_tbl(3,1) = ("Min")
putpdf table age_tbl(3,2) = ("`min_age'")
putpdf table age_tbl(4,1) = ("25th percentile")
putpdf table age_tbl(4,2) = ("`p25_age'")
putpdf table age_tbl(5,1) = ("Median")
putpdf table age_tbl(5,2) = ("`p50_age'")
putpdf table age_tbl(6,1) = ("75th percentile")
putpdf table age_tbl(6,2) = ("`p75_age'")
putpdf table age_tbl(7,1) = ("Max")
putpdf table age_tbl(7,2) = ("`max_age'")

// Sex
quietly count if Male == 1
local male_n = string(r(N), "%9.0fc")
local male_pct = string(100*r(N)/_N, "%4.1f")
quietly count if Male == 0
local female_n = string(r(N), "%9.0fc")
local female_pct = string(100*r(N)/_N, "%4.1f")

putpdf paragraph
putpdf text ("Sex Distribution"), bold

putpdf table sex_tbl = (3, 3), border(all)
putpdf table sex_tbl(1,1) = ("Sex"), bold
putpdf table sex_tbl(1,2) = ("N"), bold
putpdf table sex_tbl(1,3) = ("%"), bold
putpdf table sex_tbl(2,1) = ("Female")
putpdf table sex_tbl(2,2) = ("`female_n'")
putpdf table sex_tbl(2,3) = ("`female_pct'")
putpdf table sex_tbl(3,1) = ("Male")
putpdf table sex_tbl(3,2) = ("`male_n'")
putpdf table sex_tbl(3,3) = ("`male_pct'")

****************************************************************************
* CLINICAL CHARACTERISTICS
****************************************************************************

putpdf paragraph
putpdf text ("Clinical Characteristics at Diagnosis"), bold font(,16)

// ECOG
quietly count if ECOGcc == 0
local ecog0_n = string(r(N), "%9.0fc")
local ecog0_pct = string(100*r(N)/_N, "%4.1f")
quietly count if ECOGcc == 1
local ecog1_n = string(r(N), "%9.0fc")
local ecog1_pct = string(100*r(N)/_N, "%4.1f")
quietly count if ECOGcc == 2
local ecog2_n = string(r(N), "%9.0fc")
local ecog2_pct = string(100*r(N)/_N, "%4.1f")

putpdf paragraph
putpdf text ("ECOG Performance Status"), bold

putpdf table ecog_tbl = (4, 3), border(all)
putpdf table ecog_tbl(1,1) = ("ECOG Score"), bold
putpdf table ecog_tbl(1,2) = ("N"), bold
putpdf table ecog_tbl(1,3) = ("%"), bold
putpdf table ecog_tbl(2,1) = ("0")
putpdf table ecog_tbl(2,2) = ("`ecog0_n'")
putpdf table ecog_tbl(2,3) = ("`ecog0_pct'")
putpdf table ecog_tbl(3,1) = ("1")
putpdf table ecog_tbl(3,2) = ("`ecog1_n'")
putpdf table ecog_tbl(3,3) = ("`ecog1_pct'")
putpdf table ecog_tbl(4,1) = ("2+")
putpdf table ecog_tbl(4,2) = ("`ecog2_n'")
putpdf table ecog_tbl(4,3) = ("`ecog2_pct'")

putpdf pagebreak

// R-ISS
quietly count if RISS == 1
local riss1_n = string(r(N), "%9.0fc")
local riss1_pct = string(100*r(N)/_N, "%4.1f")
quietly count if RISS == 2
local riss2_n = string(r(N), "%9.0fc")
local riss2_pct = string(100*r(N)/_N, "%4.1f")
quietly count if RISS == 3
local riss3_n = string(r(N), "%9.0fc")
local riss3_pct = string(100*r(N)/_N, "%4.1f")

putpdf paragraph
putpdf text ("Revised International Staging System (R-ISS)"), bold

putpdf table riss_tbl = (4, 3), border(all)
putpdf table riss_tbl(1,1) = ("R-ISS Stage"), bold
putpdf table riss_tbl(1,2) = ("N"), bold
putpdf table riss_tbl(1,3) = ("%"), bold
putpdf table riss_tbl(2,1) = ("I")
putpdf table riss_tbl(2,2) = ("`riss1_n'")
putpdf table riss_tbl(2,3) = ("`riss1_pct'")
putpdf table riss_tbl(3,1) = ("II")
putpdf table riss_tbl(3,2) = ("`riss2_n'")
putpdf table riss_tbl(3,3) = ("`riss2_pct'")
putpdf table riss_tbl(4,1) = ("III")
putpdf table riss_tbl(4,2) = ("`riss3_n'")
putpdf table riss_tbl(4,3) = ("`riss3_pct'")

// Autologous Stem Cell Transplant
capture confirm variable SCT_L1
if _rc == 0 {
    quietly count if SCT_L1 == 1
    local sct_n = string(r(N), "%9.0fc")
    local sct_pct = string(100*r(N)/_N, "%4.1f")
    quietly count if SCT_L1 == 0
    local nosct_n = string(r(N), "%9.0fc")
    local nosct_pct = string(100*r(N)/_N, "%4.1f")
    
    putpdf paragraph
    putpdf text ("Receipt of ASCT"), bold
    
    putpdf table sct_tbl = (3, 3), border(all)
    putpdf table sct_tbl(1,1) = ("Received ASCT"), bold
    putpdf table sct_tbl(1,2) = ("N"), bold
    putpdf table sct_tbl(1,3) = ("%"), bold
    putpdf table sct_tbl(2,1) = ("No")
    putpdf table sct_tbl(2,2) = ("`nosct_n'")
    putpdf table sct_tbl(2,3) = ("`nosct_pct'")
    putpdf table sct_tbl(3,1) = ("Yes")
    putpdf table sct_tbl(3,2) = ("`sct_n'")
    putpdf table sct_tbl(3,3) = ("`sct_pct'")
}

******************
* Overall Survival
******************

putpdf paragraph
putpdf text ("Overall Survival Results"), bold font(,16)

// Summary statistics
quietly summarize OC_TIME, detail
local n = string(r(N), "%9.0fc")
local mean = string(r(mean), "%6.2f")
local sd = string(r(sd), "%5.2f")
local median = string(r(p50), "%6.2f")
local p25 = string(r(p25), "%6.2f")
local p75 = string(r(p75), "%6.2f")

putpdf paragraph
putpdf text ("Summary Statistics"), bold

putpdf table os_sum = (4, 2), border(all)
putpdf table os_sum(1,1) = ("Statistic"), bold
putpdf table os_sum(1,2) = ("Value (months)"), bold
putpdf table os_sum(2,1) = ("N")
putpdf table os_sum(2,2) = ("`n'")
putpdf table os_sum(3,1) = ("Mean (SD)")
putpdf table os_sum(3,2) = ("`mean' (`sd')")
putpdf table os_sum(4,1) = ("Median [IQR]")
putpdf table os_sum(4,2) = ("`median' [`p25'-`p75']")

// Survival at key time points
putpdf paragraph
putpdf text ("Survival at Key Time Points"), bold

putpdf table surv_time = (6, 2), border(all)
putpdf table surv_time(1,1) = ("Time Point"), bold
putpdf table surv_time(1,2) = ("Survival % (95% CI)"), bold

local row = 2
foreach year in 1 2 3 5 10 {
    quietly count if OC_TIME/12 >= `year'
    local pct = (r(N) / _N) * 100
    local se = sqrt(`pct' * (100 - `pct') / _N)
    local ci_l = max(0, `pct' - 1.96 * `se')
    local ci_u = min(100, `pct' + 1.96 * `se')
    
    local pct_str = string(`pct', "%5.1f")
    local ci_l_str = string(`ci_l', "%5.1f")
    local ci_u_str = string(`ci_u', "%5.1f")
    
    putpdf table surv_time(`row',1) = ("`year'-year")
    putpdf table surv_time(`row',2) = ("`pct_str'% (`ci_l_str'%-`ci_u_str'%)")
    local row = `row' + 1
}

// Generate and insert KM curves
preserve
capture mkdir "$simulated_path/results"
capture mkdir "$simulated_path/results/figures"

// Overall KM
stset OC_TIME if OC_TIME < 240, failure(OC_MORT)
sts graph, ///
    title("Overall Survival") ///
    xtitle("Months") ytitle("Probability") ///
    ylabel(0(0.2)1, angle(0) format(%3.1f)) ///
    xlabel(0(24)240) ci risktable legend(off) ///
    graphregion(color(white)) name(os_overall, replace)
graph export "$simulated_path/results/figures/os_overall.png", replace width(1200)

putpdf paragraph
putpdf image "$simulated_path/results/figures/os_overall.png", width(6)

// By ASCT
capture confirm variable SCT_L1
if !_rc {
    gen asct = SCT_L1
    label define asct_lbl 0 "No ASCT" 1 "ASCT"
    label values asct asct_lbl
    
    sts graph, by(asct) ///
        title("OS by ASCT Status") ///
        xtitle("Months") ytitle("Probability") ///
        ylabel(0(0.2)1, angle(0)) xlabel(0(24)240) ///
        graphregion(color(white)) ///
        legend(label(1 "No ASCT") label(2 "ASCT")) ///
        name(os_asct, replace)
    graph export "$simulated_path/results/figures/os_asct.png", replace width(1200)
    
    putpdf paragraph
    putpdf text ("Overall Survival by ASCT Status"), bold
    putpdf image "$simulated_path/results/figures/os_asct.png", width(6)
}

restore

// By BCR
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
    
    stset OC_TIME if OC_TIME < 240, failure(OC_MORT)
    
    sts graph, by(bcr_group) ///
        title("OS by Best Clinical Response") ///
        xtitle("Months") ytitle("Probability") ///
        ylabel(0(0.2)1, angle(0)) xlabel(0(24)240) ///
        graphregion(color(white)) ///
        legend(label(1 "CR") label(2 "VGPR") label(3 "PR") label(4 "MR") label(5 "SD") label(6 "PD") rows(2)) ///
        name(os_bcr, replace)
    graph export "$simulated_path/results/figures/os_bcr.png", replace width(1200)
    
    putpdf paragraph
    putpdf text ("Overall Survival by Best Clinical Response"), bold
    putpdf image "$simulated_path/results/figures/os_bcr.png", width(6)
    restore
}

// By Age
preserve
gen age_group = .
replace age_group = 1 if Age < 65
replace age_group = 2 if Age >= 65 & Age < 75
replace age_group = 3 if Age >= 75 & Age < .

stset OC_TIME if OC_TIME < 240, failure(OC_MORT)

sts graph, by(age_group) ///
    title("OS by Age at Diagnosis") ///
    xtitle("Months") ytitle("Probability") ///
    ylabel(0(0.2)1, angle(0)) xlabel(0(24)240) ///
    graphregion(color(white)) ///
    legend(label(1 "<65") label(2 "65-74") label(3 "≥75") rows(1)) ///
    name(os_age, replace)
graph export "$simulated_path/results/figures/os_age.png", replace width(1200)

putpdf paragraph
putpdf text ("Overall Survival by Age Group"), bold
putpdf image "$simulated_path/results/figures/os_age.png", width(6)
restore

// By R-ISS
capture confirm variable RISS
if !_rc {
    preserve
    stset OC_TIME if OC_TIME < 240, failure(OC_MORT)
    
    sts graph, by(RISS) ///
        title("OS by R-ISS Stage") ///
        xtitle("Months") ytitle("Probability") ///
        ylabel(0(0.2)1, angle(0)) xlabel(0(24)240) ///
        graphregion(color(white)) ///
        legend(label(1 "Stage I") label(2 "Stage II") label(3 "Stage III") rows(1)) ///
        name(os_riss, replace)
    graph export "$simulated_path/results/figures/os_riss.png", replace width(1200)
    
    putpdf paragraph
    putpdf text ("Overall Survival by R-ISS Stage"), bold
    putpdf image "$simulated_path/results/figures/os_riss.png", width(6)
    restore
}

// By ECOG
capture confirm variable ECOGcc
if !_rc {
    preserve
    stset OC_TIME if OC_TIME < 240, failure(OC_MORT)
    
    sts graph, by(ECOGcc) ///
        title("OS by ECOG Performance Status") ///
        xtitle("Months") ytitle("Probability") ///
        ylabel(0(0.2)1, angle(0)) xlabel(0(24)240) ///
        graphregion(color(white)) ///
        legend(label(1 "ECOG 0") label(2 "ECOG 1") label(3 "ECOG 2+") rows(1)) ///
        name(os_ecog, replace)
    graph export "$simulated_path/results/figures/os_ecog.png", replace width(1200)
    
    putpdf paragraph
    putpdf text ("Overall Survival by ECOG Performance Status"), bold
    putpdf image "$simulated_path/results/figures/os_ecog.png", width(6)
    restore
}

****************************************************************************
* SAVE PDF
****************************************************************************

set graphics on
local output_file "`report_dir'/report_${Analysis}_${Int}_${Data}.pdf"
putpdf save "`output_file'", replace

****************************************************************************
* SUMMARY
****************************************************************************

di as text _n(2) "{hline 78}"
di as result "{bf:Report Generation Complete}"
di as text "{hline 78}"
di as result "✓ " as text "PDF created: `output_file'"
di as text "{hline 78}"

****************************************************************************
* End of File
****************************************************************************
