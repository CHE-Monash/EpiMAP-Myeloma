**********
* Monash Myeloma Model - Treatment cost calculator (PBS first-principles DPMQ method)
*
* Purpose: Build the per-cycle regimen drug costs and the (inflated) non-drug costs the economic
*          model uses, for a chosen price year, from the committed PBS inputs. Replaces the old
*          "DPMA x vials" spreadsheet port: costs are now built from first principles as the PBS
*          Dispensed Price for Maximum Quantity (DPMQ), reproducible from a dated PBS extract
*          (see prep/extract_pbs_costs.do).
* Usage:   do "prep/treatment_costs.do" [target_year]      (default 2026)
* Inputs:  prep/inputs/treatment_regimens.csv - dosing spec (dose, basis, schedule, cycles)
*          prep/inputs/pbs_prices_<year>.csv   - AEMP per (drug, strength, program, pack), cheapest brand
*          prep/inputs/pbs_fees_<year>.csv     - EFC preparation / dispensing fees per program
*          prep/inputs/pbs_markups_<year>.csv  - General-Schedule wholesale + pharmacy mark-up bands
*          prep/inputs/pbs_copayments_<year>.csv - general / concessional co-payment (perspective switch)
*          (<year> = target price year; the matching schedule's extract, from extract_pbs_costs.do)
*          prep/inputs/other_costs.csv         - non-drug costs (ASCT, hospital, community, emergency)
*          prep/inputs/cost_index.csv          - year x ABS price index (inflates the non-drug costs)
* Output:  prep/inputs/treatment_costs_<year>.csv - the c* per-cycle drug costs + inflated non-drug costs
* Method:  DPMQ per program (perspective = Australian health system: full DPMQ, co-payment NOT netted):
*            EFC injectable (IN): fewest vials (min count, then min AEMP) to cover the dose + EFC
*                                 preparation fee, per administration.
*            S100 HSD oral (HB) : pack AEMP + dispensing fee.
*            General oral   (GE) : pack AEMP + wholesale + pharmacy mark-up + dispensing fee.
*          Oral pack policy (`oral_policy'): "wholepack" (default) costs whole packs per cycle (ceil to
*            the pack boundary, so dispensing wastage is included) — this is the PBAC costing convention;
*            "prorata" costs only units consumed. Injectables are always fewest-vials (a vial cannot be
*            split; its wastage is inherent). Set `net_copay' = 1 for a cost-to-government view.
*          Per-drug per-cycle = per-admin (injectable) or per-cycle (oral) cost x admins-per-cycle.
*          Regimen cost = sum over its drugs. "Other" = mean(VTd,TCd,Td,Vd); maintenance = MRDR
*          usage-weighted blend of lenalidomide (R) and thalidomide (T). Body size is a population
*          average (BSA / weight from MRDR); per-patient dosing is a future extension.
* Notes:   Does NOT modify core/process_data.do (which still carries hardcoded locals). This script
*          regenerates and validates the c* values; wiring process_data to read the output CSV is a
*          separate change.
**********

clear all
set more off
if "$repo_path" != "" cd "$repo_path"
capture run "config.do"

* ---- Config (args: 1=year 2=oral_policy 3=net_copay; all optional) ----
local target = "`1'"
if "`target'" == "" local target 2026
local oral_policy = "`2'"
if "`oral_policy'" == "" local oral_policy "wholepack"  // wholepack (default, PBAC: waste costed) | prorata
local net_copay = "`3'"
if "`net_copay'" == "" local net_copay 0               // 0 = full DPMQ (health-system); 1 = net co-pay (govt)
local BSA         = 1.93          // population-average body surface area (m^2), MRDR
local WEIGHT      = 81.10597      // population-average weight (kg), MRDR
local setting     "public"        // public path: EFC IN, S100 HSD HB, General Schedule GE
local copay_class "concessional"  // concessional (default) | general  (only used if net_copay==1)
local IN          "prep/inputs"

di as text _n "=== Treatment costs, price year `target' (`setting', orals=`oral_policy', net_copay=`net_copay') ===" _n

* ---- Load the PBS inputs into frames for the Mata core ----
* year-stamped PBS inputs: the extract writes pbs_*_<schedule year>.csv, so multiple price years coexist
import delimited "`IN'/pbs_prices_`target'.csv", varnames(1) case(preserve) clear
keep if setting == "`setting'" | (kind == "oral" & setting == "community")
frame copy default fprices, replace

* fees, markups, co-payments
import delimited "`IN'/pbs_fees_`target'.csv",       varnames(1) case(preserve) clear
frame copy default ffees, replace
import delimited "`IN'/pbs_markups_`target'.csv",    varnames(1) case(preserve) clear
frame copy default fmark, replace
import delimited "`IN'/pbs_copayments_`target'.csv", varnames(1) case(preserve) clear
quietly summarize copay if type == "`copay_class'", meanonly
scalar copay_amt = cond(`net_copay', r(mean), 0)

* regimens (some are phased: DVd load/mid/tail, Kd load/maint - front-loaded schedules)
import delimited "`IN'/treatment_regimens.csv", varnames(1) case(preserve) clear
capture confirm variable phase
if _rc gen phase = ""
frame copy default fregs, replace

* cost keys = regimen, or regimen_phase for phased regimens (DVd/Kd phases p1/p2/p3)
local reglist "VCd VRd Rd Kd_p1 Kd_p2 DVd_p1 DVd_p2 DVd_p3 Pd Vd VTd TCd Td R T VPd DPd Dd PCd D KCd KRd"

* scalars the Mata block reads (locals are not visible inside mata:)
scalar BSA_    = `BSA'
scalar WEIGHT_ = `WEIGHT'

* ===========================================================================
* Mata core: build DPMQ per drug/dose, roll up to per-cycle regimen costs
* ===========================================================================
mata:
mata clear

// ---- pull inputs from frames ----
st_framecurrent("fprices")
p_drug = st_sdata(., "drug"); p_prog = st_sdata(., "program"); p_kind = st_sdata(., "kind")
p_str  = st_data(., "strength_mg"); p_pack = st_data(., "pack_size"); p_aemp = st_data(., "aemp")

st_framecurrent("ffees")
f_prog = st_sdata(., "program")
f_disp = st_data(., "disp_fee"); f_prep = st_data(., "prep_fee")

st_framecurrent("fmark")
m_type = st_sdata(., "markup"); m_lo = st_data(., "lo_limit")
m_pct  = st_data(., "pct"); m_off = st_data(., "off"); m_fix = st_data(., "fixed_amt")

copay = st_numscalar("copay_amt")

// ---- fee lookups ----
real scalar fee_disp(string scalar prog, string colvector fp, real colvector fd) {
    for (i=1;i<=rows(fp);i++) if (fp[i]==prog) return(fd[i]==. ? 0 : fd[i])
    return(0)
}
real scalar fee_prep(string scalar prog, string colvector fp, real colvector fr) {
    for (i=1;i<=rows(fp);i++) if (fp[i]==prog) return(fr[i]==. ? 0 : fr[i])
    return(0)
}

// ---- General-Schedule banded mark-up (wholesale or pharmacy) ----
real scalar markup(real scalar aemp, string scalar kind, string colvector mt,
                   real colvector ml, real colvector mp, real colvector mo, real colvector mf) {
    real scalar bestlo, pct, off, fix
    bestlo = -1
    for (i=1;i<=rows(mt);i++) {
        if (mt[i]==kind & ml[i] <= aemp & ml[i] > bestlo) {
            bestlo = ml[i]; pct = mp[i]; off = mo[i]; fix = mf[i]
        }
    }
    if (bestlo < 0) return(0)
    return(pct/100*(aemp + off) + fix)
}

// ---- injectable: fewest vials (min count, then min AEMP) covering dose ----
// returns sum of chosen vial AEMPs (prep fee added by caller)
real scalar inj_vialcost(real scalar D, real colvector S, real colvector A) {
    real scalar L, n, k, i, j
    real colvector sums, costs, ns, nc, us, idx
    L = max(S)
    n = max((1, ceil(D/L - 1e-9)))
    sums = S; costs = A                      // 1-vial reachable (sum, min cost)
    for (k=2; k<=n; k++) {
        ns = J(0,1,.); nc = J(0,1,.)
        for (i=1;i<=rows(sums);i++) {
            for (j=1;j<=rows(S);j++) {
                ns = ns \ (sums[i]+S[j])
                nc = nc \ (costs[i]+A[j])
            }
        }
        us = uniqrows(ns)                    // collapse to min cost per unique sum
        sums = us; costs = J(rows(us),1,.)
        for (i=1;i<=rows(us);i++) costs[i] = min(nc[selectindex(ns:==us[i])])
    }
    idx = selectindex(sums :>= D - 1e-9)
    return(min(costs[idx]))
}

// ---- oral pack DPMQ (per pack) ----
real scalar oral_packdpmq(real scalar aemp, string scalar prog,
                          string colvector mt, real colvector ml, real colvector mp,
                          real colvector mo, real colvector mf,
                          string colvector fp, real colvector fd) {
    real scalar d
    d = aemp + fee_disp(prog, fp, fd)
    if (prog=="GE") d = d + markup(aemp,"wholesale",mt,ml,mp,mo,mf) + markup(aemp,"pharmacy",mt,ml,mp,mo,mf)
    return(d)
}

// ---- loop over regimen rows ----
st_framecurrent("fregs")
r_reg = st_sdata(., "regimen"); r_drug = st_sdata(., "drug"); r_basis = st_sdata(., "basis")
r_str = st_data(., "strength_mg"); r_dose = st_data(., "dose")
r_freq = st_data(., "freq_per_cycle"); r_phase = st_sdata(., "phase")
BSA = st_numscalar("BSA_"); WT = st_numscalar("WEIGHT_")
policy = st_local("oral_policy")

regs = tokens(st_local("reglist"))
totals = J(1, cols(regs), 0)

for (i=1;i<=rows(r_reg);i++) {
    // body size and dose; per-cycle admins = freq_per_cycle (within this phase, if phased)
    bs = (r_basis[i]=="m2" ? BSA : (r_basis[i]=="kg" ? WT : 1))
    dose_mg = r_dose[i]*bs
    admpc = r_freq[i]

    // injectable if this drug appears with kind=="injectable" in the price table
    isinj = anyof(select(p_kind, p_drug:==r_drug[i]), "injectable")

    if (isinj) {
        sel = selectindex(p_drug:==r_drug[i] :& p_kind:=="injectable")
        S = p_str[sel]; A = p_aemp[sel]; prog = p_prog[sel[1]]
        vc = inj_vialcost(dose_mg, S, A)
        admin = vc + fee_prep(prog, f_prog, f_prep) - copay
        pc = admin*admpc
    }
    else {
        sel = selectindex(p_drug:==r_drug[i] :& abs(p_str:-r_str[i]):<1e-6)
        packs = p_pack[sel]; aemps = p_aemp[sel]; progs = p_prog[sel]
        units_admin = round(dose_mg/r_str[i])
        units_cyc = units_admin*admpc
        best = .
        for (j=1;j<=rows(sel);j++) {
            pd = oral_packdpmq(aemps[j], progs[j], m_type,m_lo,m_pct,m_off,m_fix, f_prog,f_disp) - copay
            if (policy=="wholepack") cost = ceil(units_cyc/packs[j] - 1e-9)*pd
            else                     cost = units_cyc*(pd/packs[j])
            if (best==. | cost<best) best = cost
        }
        pc = best
    }
    // accumulate into the matching cost key (regimen, or regimen_phase for phased regimens)
    key = (r_phase[i]=="" ? r_reg[i] : r_reg[i]+"_"+r_phase[i])
    for (g=1; g<=cols(regs); g++) if (regs[g]==key) totals[g] = totals[g] + pc
}

// export per-regimen totals to Stata scalars c_<reg>
for (g=1; g<=cols(regs); g++) st_numscalar("c_"+regs[g], totals[g])
end

* ---- Derived blends: pooled "Other" (code 0) and usage-weighted maintenance ----
* c_Other is the fallback: a fixed, unweighted basket of four older regimens. It is superseded by the
* per-line rates below wherever prep/inputs/other_mix.csv is available.
scalar c_Other = (c_VTd + c_TCd + c_Td + c_Vd)/4

* ---- Per-line "Other" blends, usage-weighted from the registry ----
* scratch/regimen_freq.do (MRDR-gated) writes prep/inputs/other_mix.csv: for each line, the share of
* that line's PRICEABLE 'other' bucket taken by each regimen. Weighting by those shares gives a rate
* that reflects what is actually prescribed outside the modelled regimens, by line, instead of one
* fixed basket everywhere. Lines are L1..L5 with L6+ collapsed (line 6 = L6 onwards).
* ASSUMPTION: the registry's own 'other'/'Unknown' rows, which have no drug composition and so cannot
* be priced, are assumed to resemble the priceable mix. That is what weighting within the priceable
* subset implies, and it must be disclosed.
forval g = 1/6 {
    scalar c_Other_L`g' = .
}
capture confirm file "`IN'/other_mix.csv"
if (_rc == 0) {
    preserve
    import delimited "`IN'/other_mix.csv", varnames(1) case(preserve) clear
    forval g = 1/6 {
        qui count if lineno == `g'
        if (r(N) > 0) {
            scalar c_Other_L`g' = 0
            qui levelsof costkey if lineno == `g', local(_ks) clean
            foreach k of local _ks {
                qui summarize weight if lineno == `g' & costkey == "`k'", meanonly
                local w = r(mean)
                capture local cv = scalar(c_`k')
                if (_rc == 0 & "`cv'" != "" & !missing(scalar(c_`k'))) {
                    scalar c_Other_L`g' = scalar(c_Other_L`g') + `w' * scalar(c_`k')
                }
            }
        }
    }
    restore
    di as text _n "Per-line 'Other' blends (usage-weighted from other_mix.csv):"
    forval g = 1/6 {
        local nm = cond(`g' == 6, "L6+", "L`g'")
        di as text "  c_Other_`nm'" _col(18) %12.2f scalar(c_Other_L`g')
    }
}
else {
    di as text _n "(prep/inputs/other_mix.csv not found - per-line Other blends left missing;"
    di as text " process_data.do will fall back to the single c_Other. Run scratch/regimen_freq.do.)"
}
scalar c_MNT   = (1002/1504)*c_R + (502/1504)*c_T          // MRDR R/T usage weights

di as text _n "Per-cycle drug costs (full DPMQ, `oral_policy' orals):"
foreach rg in VCd VRd Rd Kd_p1 Kd_p2 DVd_p1 DVd_p2 DVd_p3 Pd Vd Other MNT {
    di as text "  c`rg' " _col(16) %12.2f scalar(c_`rg')
}

* ===========================================================================
* Non-drug costs, inflated source-year -> target via the ABS index
* ===========================================================================
frame change default
import delimited "`IN'/cost_index.csv", varnames(1) case(preserve) clear
* index for the target year, carrying the latest year <= target forward if that year is missing
quietly summarize year if year <= `target', meanonly
scalar idx_yr = r(max)
quietly summarize index if year == idx_yr, meanonly
scalar idx_t = r(mean)
if missing(idx_t) {
    di as text "  (no ABS index at or before `target'; non-drug costs left un-inflated)"
    scalar idx_t = .
}
else if idx_yr < `target' {
    di as text "  (no ABS index for `target'; carrying `=idx_yr' forward for non-drug costs)"
}
tempfile idx
keep year index
rename (year index) (source_year idx_s)
save `idx'

import delimited "`IN'/other_costs.csv", varnames(1) case(preserve) clear
merge m:1 source_year using `idx', keep(master match) nogen
gen double factor = cond(missing(idx_t) | missing(idx_s), 1, idx_t/idx_s)
gen double cost_infl = cost*factor

* ASCT one-off (no phase). NEP-based (AR-DRG price weight x National Efficient Price), set administratively
* each year -> NOT CPI-inflated: pick the cASCT row for the latest source_year <= target and use its raw cost.
quietly summarize source_year if component == "cASCT" & source_year <= `target', meanonly
scalar asct_yr = r(max)
quietly summarize cost if component == "cASCT" & source_year == asct_yr, meanonly
scalar cASCT = r(mean)
if asct_yr < `target' di as text "  (no NEP for `target'; carrying `=asct_yr' ASCT cost forward)"

* Phase-based non-treatment (Yap 2025), kept SPLIT by component so a breakdown can be reported:
* cHosp (APDC admitted-hospital), cMBS (out-of-hospital Medicare), cEmer (emergency). PBS is already
* excluded (not in the input; drugs are costed separately). The transplant admission sits in Yap's
* INITIAL-phase hospital, so net the cohort-average ASCT out of cHosp_initial to avoid double-counting
* cASCT (applied patient-specifically in process_data): subtract asct_prev x cASCT, asct_prev = Yap's
* ASCT prevalence (STable10: 125/520). process_data sums the components per patient after phase allocation.
scalar asct_prev = 0.24
foreach comp in cHosp cMBS cEmer {
    foreach ph in initial continuing terminal {
        quietly summarize cost_infl if component == "`comp'" & phase == "`ph'"
        scalar `comp'_`ph' = r(sum)
    }
}
scalar cHosp_initial = cHosp_initial - asct_prev*cASCT

di as text _n "Non-treatment costs/yr by component (Yap 2025, `target'; PBS & ASCT excluded):"
di as text "  phase" _col(16) "cHosp" _col(30) "cMBS" _col(44) "cEmer"
foreach ph in initial continuing terminal {
    di as text "  " %-12s "`ph'" _col(14) %12.2f cHosp_`ph' _col(28) %12.2f cMBS_`ph' _col(42) %12.2f cEmer_`ph'
}
di as text "  cASCT (one-off)" _col(14) %12.2f cASCT

* ===========================================================================
* Validation: injectable admin DPMQ against known 2026 PBS figures
* ===========================================================================
if `net_copay' == 0 {
di as text _n "Validation (full DPMQ; the reference DPMQs below are verified 2026 PBS figures):"
local nfail 0
* Externally-verifiable published PBS DPMQs (policy-independent; verified against the 2026 schedule -
* asserted only when target==2026, since other price years legitimately differ):
*   Carfilzomib 108.08 mg (56 mg/m2) admin = 2 x 60 mg + prep      = 2,503.54
*   Bortezomib  2.51 mg   (1.3 mg/m2) admin = 1 x 3.5 mg + prep    = 113.03
*   Daratumumab 1297.7 mg (16 mg/kg)  admin = 3x400 + 1x100 + prep = 7,307.30
*   Lenalidomide 25 mg x21 pack (S100 HSD) = 769.49 + 9.24 disp    = 778.73
*   Kd per cycle (dexa pack is exact, so policy-independent)        = 12,629.70
* (reuses the persisted price vectors p_* and helper functions from the core Mata block above)
mata:
prepc = fee_prep("IN", f_prog, f_prep)
selc = selectindex(p_drug:=="Carfilzomib"  :& p_kind:=="injectable")
selb = selectindex(p_drug:=="Bortezomib"   :& p_kind:=="injectable")
seld = selectindex(p_drug:=="Daratumumab"  :& p_kind:=="injectable")
sell = selectindex(p_drug:=="Lenalidomide" :& abs(p_str:-25):<1e-6 :& p_pack:==21)
st_numscalar("v_c56", inj_vialcost(56*1.93,      p_str[selc], p_aemp[selc]) + prepc)
st_numscalar("v_cb",  inj_vialcost(1.3*1.93,     p_str[selb], p_aemp[selb]) + prepc)
st_numscalar("v_cd",  inj_vialcost(16*81.10597,  p_str[seld], p_aemp[seld]) + prepc)
st_numscalar("v_lena", oral_packdpmq(p_aemp[sell[1]], p_prog[sell[1]], m_type,m_lo,m_pct,m_off,m_fix, f_prog,f_disp))
end
foreach pair in "carf56 v_c56 2503.54" "bort v_cb 113.03" "dara v_cd 7307.30" "lena25 v_lena 778.73" "Kd_p1 c_Kd_p1 12629.70" {
    tokenize "`pair'"
    local got = scalar(`2')
    local d = abs(`got' - `3')
    if `target' == 2026 {
        local ok = cond(`d' < 0.5, "OK", "FAIL")
        if "`ok'" == "FAIL" local ++nfail
        di as text "  `1' " _col(10) %12.2f `got' "  expect " %12.2f `3' "   `ok'"
    }
    else di as text "  `1' " _col(10) %12.2f `got' "  (2026 ref " %12.2f `3' ")"
}
if `target' == 2026 {
    if `nfail' == 0 di as result "  All checks pass."
    else            di as error  "  `nfail' check(s) FAILED - investigate before use."
}
else di as text "  (reference figures are 2026 PBS; `target' prices differ by schedule - shown for context, not asserted)"
}
else di as text _n "(validation skipped: reference DPMQs are full-price; net_copay=1)"

* ===========================================================================
* Write the output cost set for this year
* ===========================================================================
clear
set obs 40
gen str14 parameter = ""
gen double value = .
gen str8 unit = ""
gen str48 note = ""
local i 0
foreach p in cVCd cVRd cRd cKd_p1 cKd_p2 cDVd_p1 cDVd_p2 cDVd_p3 cPd cVd cOther cMNT cVTd cTCd cTd cR cT cVPd cDPd cDd cPCd cD cKCd cKRd ///
             cOther_L1 cOther_L2 cOther_L3 cOther_L4 cOther_L5 cOther_L6 {
    local ++i
    local nm = subinstr("`p'", "c", "", 1)
    replace parameter = "`p'" in `i'
    replace value = scalar(c_`nm') in `i'
    replace unit = "AUD/cyc" in `i'
    replace note = "PBS DPMQ `target' (`oral_policy')" in `i'
}
local ++i
replace parameter = "cASCT" in `i'
replace value = scalar(cASCT) in `i'
replace unit = "AUD" in `i'
replace note = "AR-DRG R06 x NEP" in `i'
foreach comp in cHosp cMBS cEmer {
    foreach ph in initial continuing terminal {
        local ++i
        replace parameter = "`comp'_`ph'" in `i'
        replace value = scalar(`comp'_`ph') in `i'
        replace unit = "AUD/yr" in `i'
        replace note = "Yap 2025 `comp' `ph' phase (to `target')" in `i'
    }
}
replace value = round(value, 0.01)
format value %12.2f
export delimited "`IN'/treatment_costs_`target'.csv", replace datafmt
di as result _n "Wrote `IN'/treatment_costs_`target'.csv"
list parameter value unit, noobs sep(0)
