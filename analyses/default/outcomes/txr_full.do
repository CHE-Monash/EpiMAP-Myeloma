**********
* Monash Myeloma Model - regimen lists (default analysis)
*
* Purpose: declare every regimen list this analysis models - per-line treatment (TXR_L1..L9), L1
*          maintenance (MNR_L1), and the len-refractory gate (LENREFR_regimens). gen_txr in
*          prep/risk_equations.do builds TXR_L1..L9 from these; any regimen not listed for a line falls
*          into 0 = 'other'. This is the CANONICAL regimen list for the default analysis; the train fit
*          (txr_train.do) sources it, so the in-sample/out-of-sample validation uses the same regimens.
* Notes:   for per-line regimen frequencies before choosing a list, see scratch/regimen_freq.do.
**********

* Regimen codes:
*   2  Thal/Cycl/Dexa    4  Bort/Cycl/Dexa    7  Lena/Dexa       9  Bort/Thal/Dexa
*  31  Bort/Lena/Dexa   49  Carf/Dexa        56  Poma/Dexa      80  Dara/Bort/Dexa

* Rd (7) was trialled in L1 and L4 and reverted: TXD is fitted on i.TXR_L{line}, so adding a
* continuous-until-progression regimen to a line moves that line's duration equation. L4 TXD then
* over-predicted time on treatment by 24-39pp out of sample. Re-test it on its own branch, against
* the L1/L4 TXD blocks, before re-adding.
global TXR_L1 "4 31"
global TXR_L2 "7 80"
global TXR_L3 "7 49"
global TXR_L4 "49 56"

* Lenalidomide-containing regimen codes, for the refractory gate. The FIT conditions on the true
* drug binary (Lenalidomide == 1); the engine has no drug field, so sim_lenrefr.do counts a line as
* lenalidomide iff the DRAWN regimen is one of these. 'other' (code 0) is treated as non-len, which
* under-counts: some of that bucket contains lenalidomide. That asymmetry is deliberate and is one
* of the two causes of the treatment-arm under-generation in docs/refractory.md 5(6).
* Declared per analysis alongside the TXR lists so a change of modelled regimens carries the gate.
global LENREFR_regimens "7 31"
* L5-L9 unset => all 'other'

**********
* MAINTENANCE regimens (MNR_L1). gen_mnr builds MNR_L1 from this list; any drug not listed falls to
* 0 = 'other'. Drug codes, per docs/refractory.md 2:
*   0 none/other   1 lenalidomide   2 daratumumab   3 carfilzomib   4 bortezomib   5 thalidomide
*
* ONE list serves both eras, so there is no per-analysis switch. Historically thalidomide was the
* MAJORITY regimen until 2020 (51.6% of starts in 2019, none from 2021), so len + thal covers ~85%
* of the OOS window explicitly; in a modern window thalidomide simply empties out and the r(r) == 1
* guard in risk_equations.do assigns everyone lenalidomide. Mix by year: scratch/refractory/mnr_recency.do.
*
* SIMPLE-FIRST: the fits restrict to inlist(MNR_L1, 1, 5), so bortezomib, daratumumab and
* carfilzomib maintenance are excluded from estimation and the engine never draws an 'other'
* maintenance regimen. Bortezomib is 5-11% of starts but has no PBS maintenance DPMQ to price
* separately (MSAG guideline, June 2022); dara/carf maintenance is a couple of dozen mostly
* later-line patients. Rationale: docs/refractory.md 7.4.
global MNR_L1 "1 5"
