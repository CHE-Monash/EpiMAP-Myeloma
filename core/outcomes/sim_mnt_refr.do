**********
* Monash Myeloma Model - Sim maintenance len-refractory (L1)
*
* Purpose: set the LATCHED lenalidomide-refractory state for patients who progress on, or within 60
*          days of, their L1 lenalidomide MAINTENANCE. Writes into the same vLenRefr_in the
*          treatment arm uses - one flag, not two (see below).
*
* Notes:   NO FITTED MODEL. This is arithmetic on quantities the engine has already drawn:
*
*              tail = TFI_L1 - MND_L1        refractory if tail < 60 days
*
*          TFI_L1 is the L1-end-to-L2 gap and MND_L1 the maintenance duration, so the tail is the
*          time from maintenance ending to the next line starting. Under IMWG that IS the 60-day
*          clock (docs/refractory.md 1.1), because in the engine L2 starting is the progression.
*
*          THIS ONLY WORKS BECAUSE OF THE DURATION REWORK. Before it, MND and TFI were drawn
*          independently and process_data clipped the overshoot, so the tail was an artefact of the
*          clip rather than a quantity. Now sim_tfi_l1.do draws the gap truncated below at the
*          maintenance, so the tail is meaningful. It lands at 25.9% of simulated lenalidomide
*          maintenance patients under 60 days, against a registry 27.6% - which is why the fitted
*          logit this file used to carry (bL1_MNTREFR) is no longer needed.
*
*          CONSUMES NO RANDOMNESS, so it needs no CRN slot and cannot perturb the common-random-
*          number layout. Deterministic given TFI_L1 and MND_L1.
*
*          ONE FLAG. The treatment-dose and maintenance-dose flags were kept separate while it was
*          open whether they carried different prognostic weight. os_lenrefr_check.do settled that
*          at p = 0.97: they collapse, so both arms write to vLenRefr_in and the OS and TXR
*          equations read the union.
*
*          KNOWN CHARACTERISTIC, recorded rather than hidden: 87.7% of the simulated sub-60-day
*          tails are exactly zero, because the truncated TFI draw piles mass at the bound, against
*          roughly three in ten in the registry. The flag is binary and nothing downstream reads the
*          tail's value, so this does not affect any consumer - but it does mean which patients are
*          flagged depends on where the truncation binds. Whether that tracks prognosis is the open
*          question the OS-by-status check exists to answer (docs/refractory.md 5(6)).
*
* ORDER:   AFTER sim_tfi_l1.do (needs the drawn gap) and sim_mnd.do (needs the duration), at L1E.
*          The flag applies from L2 onward, which is satisfied: L2's equations run at a later OMC.
* Reads:   mTFI[.,2] (gap, months), vMND (duration, months), vMNT, vMNR.
* Writes:  vLenRefr_in (latched 0 -> 1).
**********

mata {
	// Only lenalidomide maintenance can make a patient lenalidomide-refractory. sim_mnr never draws
	// 'other', so this is len against thal.
	idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC) ///
	                  :& (vMNT :== 1) :& (vMNR :== 1) :& (vLenRefr_in :== 0))

	if (rows(idx) > 0) {
		// 60 days on the engine's clock. Both inputs are months.
		wk60 = 60 / 30.4375

		vTail = mTFI[idx, 2] :- vMND[idx]

		// A missing duration means no maintenance was costed (no MND model, or an arm that produced
		// none), and a missing gap means no L2 was drawn. Neither can be called refractory.
		vAcq = (vTail :< wk60) :& (vTail :< .) :& (mTFI[idx, 2] :< .) :& (vMND[idx] :< .)

		vLenRefr_in[idx] = vLenRefr_in[idx] :| vAcq
	}
}

// Check for override file, execute if it exists
local override_file "${outcomes_path}/sim_mnt_refr_override.do"
capture confirm file "`override_file'"
if _rc == 0 {
	qui do `override_file'
}
