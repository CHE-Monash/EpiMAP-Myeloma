**********
* Monash Myeloma Model - TXR regimens per line (default analysis)
*
* Purpose: declare the per-line regimen code lists (MRDR Regimen codes). gen_txr in
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

