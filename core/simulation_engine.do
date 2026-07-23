**********
* Monash Myeloma Model - Simulation Engine
*
* Purpose: Drive the simulation pass, dispatching each pathway point (DN, L1S..L9E) to
*          its outcome module in core/outcomes/ (ASCT, TFI, OS, mortality, BCR, TXR, TXD).
**********

cap program drop simulation
program define simulation

di "Running simulation"

**********
*Diagnosis (DN)
	mata: OMC = 1
	mata: Line = 0
			
	di "DN - SCT"	
		qui do "core/outcomes/sim_asct_dn.do"
		*mata: _matrix_list(bDN_SCT, rbDN_SCT, cbDN_SCT)
		*mata: _matrix_list(vSCT_DN)		
		
	di "DN - Treatment-free Interval"		
		qui do "core/outcomes/sim_tfi_dn.do"
		*mata: _matrix_list(bDN_TFI, rbDN_TFI, cbDN_TFI)
		*mata: _matrix_list(mTFI, rmTFI, cmTFI)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)	

	di "DN - Overall Survival" 
		qui do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)
				
	di "DN - Mortality"
		qui do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
		
**********
*Line 1 Start (L1S)	
	mata: OMC = 2
	mata: Line = 1
	
	di "L1S - Age"
		qui do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)
		
	di "L1S - Treatment Regimen"
		qui do "core/outcomes/sim_txr.do"
		*mata: _matrix_list(bL1_TXR, rbL1_TXR, cbL1_TXR)
		*mata: _matrix_list(mTXR, rmTXR, cmTXR)
	di "L1S - Best Clinical Response"
		qui do "core/outcomes/sim_bcr.do"
		*mata: _matrix_list(bL1_BCR, rbL1_BCR, cbL1_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)

	di "L1S - Treatment Duration"
		qui do "core/outcomes/sim_txd_l1.do"
		*mata: _matrix_list(bL1_TXD_ASCT_S1, rbL1_TXD_ASCT_S1, cbL1_TXD_ASCT_S1)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L1S - Overall Survival"
		qui do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)
				
	di "L1S - Mortality"
		qui do "core/outcomes/sim_mort.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
		
**********
*Line 1 End (L1E)
	mata: OMC = 3

	di "L1E - Age"
		qui do "core/outcomes/sim_age.do"
		*mata: _matrix_list(mAge, rmAge, cmAge)	

	di "L1E - SCT"
		qui do "core/outcomes/sim_asct_l1.do"			
		*mata: _matrix_list(bL1_SCT, rbL1_SCT, cbL1_SCT)
		*mata: _matrix_list(vSCT_L1)

	di "L1E - SCT Best Clinical Response"
		qui do "core/outcomes/sim_bcr_asct.do"			
		*mata: _matrix_list(bSCT_BCR, rbSCT_BCR, cbSCT_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)

	di "L1E - MNT"
		qui do "core/outcomes/sim_mnt.do"
		*mata: _matrix_list(bMNT, rbMNT, cbMNT)
		*mata: _matrix_list(vMNT)

	di "L1E - MNT Regimen"
		qui do "core/outcomes/sim_mnr.do"
		*mata: _matrix_list(bL1_MNR, rbL1_MNR, cbL1_MNR)
		*mata: _matrix_list(vMNR)

	// ORDER REVERSED, July 2026: MND now comes BEFORE TFI. The maintenance duration is drawn from a
	// fit that no longer conditions on the gap - which is what frees it from the complete-gap
	// selection - and sim_tfi_l1 then draws the gap TRUNCATED below at that duration, so the gap
	// always contains the maintenance it has to hold. Reversing these two is the whole design; the
	// previous order existed only so MND could read the drawn gap as a covariate.
	// No CRN consequence: both files consume the same slots they always did.
	di "L1E - MNT Duration (drawn first; the gap is then truncated to contain it)"
		qui do "core/outcomes/sim_mnd.do"
		*mata: _matrix_list(bL1_MND_LEN, rbL1_MND_LEN, cbL1_MND_LEN)
		*mata: _matrix_list(vMND)

	di "L1E - Treatment-free Interval (truncated below at vMND for maintenance patients)"
		qui do "core/outcomes/sim_tfi_l1.do"			
		*mata: _matrix_list(bL1_TFI, rbL1_TFI, cbL1_TFI)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)	
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)	


	di "L1E - Overall Survival"
		qui do "core/outcomes/sim_os.do"

	// LenRefr (treatment lines): fires AFTER this line's OS, so every equation reads the state from
	// STRICTLY PRIOR lines - which is how the covariate was fitted. Latched 0 -> 1.
	di "L1E - LenRefr (treatment)"
		qui do "core/outcomes/sim_lenrefr.do"

	// Maintenance len-refractory. The tail (TFI_L1 - MND_L1) needs both draws, which happened
	// earlier at this OMC, but this runs AFTER L1E's OS for the same reason the treatment arm does:
	// the refractory event occurs at the END of the gap, when L2 starts, so it must not influence
	// survival ACROSS that gap. Read from L2 on. Arithmetic, not a fitted model, so no CRN slot.
	di "L1E - Maintenance len-refractory (tail < 60 days)"
		qui do "core/outcomes/sim_mnt_refr.do"
		*mata: _matrix_list(vLenRefr_in)
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)


	di "L1E - Mortality"
		qui do "core/outcomes/sim_mort.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
		
**********
*Line 2 Start (L2S)	
	mata: OMC = 4
	mata: Line = 2
	
	di "L2S - Age"
		qui do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
		
	di "L2S - Treatment Regimen"
		qui do "core/outcomes/sim_txr.do"
		*mata: _matrix_list(bL2_TXR, rbL2_TXR, cbL2_TXR)
		*mata: _matrix_list(mTXR, rmTXR, cmTXR)
	di "L2S - Best Clinical Response"
		qui do "core/outcomes/sim_bcr.do"
		*mata: _matrix_list(bL2_BCR, rbL2_BCR, cbL2_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)	

	di "L2S - Treatment Duration"
		qui do "core/outcomes/sim_txd.do"
		*mata: _matrix_list(bL2_TXD, rbL2_TXD, cbL2_TXD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)

	di "L2S - Overall Survival"
		qui do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)
			
	di "L2S - Mortality"
		qui do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
			
**********
*Line 2 End (L2E)	
	mata: OMC = 5

	di "L2E - Age"
		qui do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	

	di "L2E - Treatment-free Interval"
		qui do "core/outcomes/sim_tfi.do"
		*mata: _matrix_list(bL2_TFI, rbL2_TFI, cbL2_TFI)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)	

	di "L2E - Overall Survival"
		qui do "core/outcomes/sim_os.do"

	// LenRefr (treatment lines): fires AFTER this line's OS, so every equation reads the state from
	// STRICTLY PRIOR lines - which is how the covariate was fitted. Latched 0 -> 1.
	di "L2E - LenRefr (treatment)"
		qui do "core/outcomes/sim_lenrefr.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS))


	di "L2E - Mortality"
		qui do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Line 3 Start (L3S)
	mata: OMC = 6
	mata: Line = 3
		
	di "L3S - Age"
		qui do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)

	di "L3S - Treatment Regimen"
		qui do "core/outcomes/sim_txr.do"
		*mata: _matrix_list(bL3_TXR, rbL3_TXR, cbL3_TXR)
		*mata: _matrix_list(mTXR, rmTXR, cmTXR)
	di "L3S - Best Clinical Response"
		qui do "core/outcomes/sim_bcr.do"
		*mata: _matrix_list(bL3_BCR, rbL3_BCR, cbL3_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)				
				
	di "L3S - Treatment Duration"
		qui do "core/outcomes/sim_txd.do"	
		*mata: _matrix_list(bL3_TXD, rbL3_TXD, cbL3_TXD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L3S - Overall Survival"
		qui do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)
		
	di "L3S - Mortality"
		qui do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Line 3 End (L3E)		
	mata: OMC = 7
	
	di "L3E - Age"
		qui do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)

	di "L3E - Treatment-free Interval"
		qui do "core/outcomes/sim_tfi.do"		
		*mata: _matrix_list(bL3_TFI, rbL3_TFI, cbL3_TFI)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)

	di "L3E - Overall Survival"
		qui do "core/outcomes/sim_os.do"

	// LenRefr (treatment lines): fires AFTER this line's OS, so every equation reads the state from
	// STRICTLY PRIOR lines - which is how the covariate was fitted. Latched 0 -> 1.
	di "L3E - LenRefr (treatment)"
		qui do "core/outcomes/sim_lenrefr.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)


	di "L3E - Mortality"
		qui do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Line 4 Start (L4S)	
	mata: OMC = 8
	mata: Line = 4
		
	di "L4S - Age"
		qui do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	

	di "L4S - Treatment Regimen"
		qui do "core/outcomes/sim_txr.do"
		*mata: _matrix_list(bL4_TXR, rbL4_TXR, cbL4_TXR)
		*mata: _matrix_list(mTXR, rmTXR, cmTXR)
	di "L4S - Best Clinical Response"
		qui do "core/outcomes/sim_bcr.do"
		*mata: _matrix_list(bL4_BCR, rbL4_BCR, cbL4_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)	
		
	di "L4S - Treatment Duration"
		qui do "core/outcomes/sim_txd.do"
		*mata: _matrix_list(bL4_TXD, rbL4_TXD, cbL4_TXD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L4S - Overall Survival"
		qui do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)
		
	di "L4S - Mortality"
		qui do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Line 4 End (L4E)
	mata: OMC = 9
		
	di "L4E - Age"
		qui do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)
		
	di "L4E - Treatment-free Interval"
		qui do "core/outcomes/sim_tfi.do"	
		*mata: _matrix_list(bL4_TFI, rbL4_TFI, cbL4_TFI)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)

	di "L4E - Overall Survival"
		qui do "core/outcomes/sim_os.do"

	// LenRefr (treatment lines): fires AFTER this line's OS, so every equation reads the state from
	// STRICTLY PRIOR lines - which is how the covariate was fitted. Latched 0 -> 1.
	di "L4E - LenRefr (treatment)"
		qui do "core/outcomes/sim_lenrefr.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)


	di "L4E - Mortality"
		qui do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Line 5 Start (L5S)	
	mata: OMC = 10
	mata: Line = 5
		
	di "L5S - Age" 
		qui do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)
		
	di "L5S - Treatment Regimen"
		qui do "core/outcomes/sim_txr.do"
		*mata: _matrix_list(mTXR, rmTXR, cmTXR)

	di "L5S - Best Clinical Response"
		qui do "core/outcomes/sim_bcr.do"
		*mata: _matrix_list(bL5_BCR, rbL5_BCR, cbL5_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)		

	di "L5S - Treatment Duration"
		qui do "core/outcomes/sim_txd.do"		
		*mata: _matrix_list(bL5_TXD, rbL5_TXD, cbL5_TXD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L5S - Overall Survival"
		qui do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)	
	
	di "L5S - Mortality"
		qui do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Line 5 End (L5E)	
	mata: OMC = 11
		
	di "L5E - Age"
		qui do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)

	di "L5E - Treatment-free Interval"
		qui do "core/outcomes/sim_tfi.do"				
		*mata: _matrix_list(bL5_TFI, rbL5_TFI, cbL5_TFI)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)

	di "L5E - Overall Survival"
		qui do "core/outcomes/sim_os.do"

	// LenRefr (treatment lines): fires AFTER this line's OS, so every equation reads the state from
	// STRICTLY PRIOR lines - which is how the covariate was fitted. Latched 0 -> 1.
	di "L5E - LenRefr (treatment)"
		qui do "core/outcomes/sim_lenrefr.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)


	di "L5E - Mortality"
		qui do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
		
**********
*Line 6 Start (L6S) 		
	mata: OMC = 12
	mata: Line = 6
		
	di "L6S - Age" 
		qui do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)
		
	di "L6S - Treatment Regimen"
		qui do "core/outcomes/sim_txr.do"
		*mata: _matrix_list(mTXR, rmTXR, cmTXR)
		
	di "L6S - Best Clinical Response"
		qui do "core/outcomes/sim_bcr.do"
		*mata: _matrix_list(bL6_BCR, rbL6_BCR, cbL6_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)		

	di "L6S - Treatment Duration"
		qui do "core/outcomes/sim_txd.do"	
		*mata: _matrix_list(bL6_TXD, rbL6_TXD, cbL6_TXD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L6S - Overall Survival"
		qui do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)	
	
	di "L6S - Mortality"
		qui do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Line 6 End (L6E)
	mata: OMC = 13
	
	di "L6E - Age"
		qui do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)
		
	di "L6E - Treatment-free Interval"
		qui do "core/outcomes/sim_tfi.do"				
		*mata: _matrix_list(bL6_TFI, rbL6_TFI, cbL6_TFI)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)	

	di "L6E - Overall Survival"
		qui do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L6E - Mortality" 
		qui do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
					
**********
*Line 7 Start (L7S) 	
	mata: OMC = 14
	mata: Line = 7
		
	di "L7S - Age" 
		qui do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)
		
	di "L7S - Treatment Regimen"
		qui do "core/outcomes/sim_txr.do"
		*mata: _matrix_list(mTXR, rmTXR, cmTXR)
		
	di "L7E - Best Clinical Response" 
		qui do "core/outcomes/sim_bcr.do"
		*mata: _matrix_list(bL7_BCR, rbL7_BCR, cbL7_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)	

	di "L7S - Treatment Duration"
		qui do "core/outcomes/sim_txd.do"	
		*mata: _matrix_list(bL7_TXD, rbL7_TXD, cbL7_TXD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L7S - Overall Survival"
		qui do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L7S - Mortality"
		qui do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Line 7 End (L7E)
	mata: OMC = 15
		
	di "L7E - Age"
		qui do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)
		
	di "L7E - Treatment-free Interval"
		qui do "core/outcomes/sim_tfi.do"				
		*mata: _matrix_list(bL7_TFI, rbL7_TFI, cbL7_TFI)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)		

	di "L7E - Overall Survival"
		qui do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L7E - Mortality" 
		qui do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
		
**********
*Line 8 Start (L8S) 
	mata: OMC = 16
	mata: Line = 8
	
	di "L8S - Age" 
		qui do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)
		
	di "L8S - Treatment Regimen"
		qui do "core/outcomes/sim_txr.do"
		*mata: _matrix_list(mTXR, rmTXR, cmTXR)
		
	di "L8S - Best Clinical Response"
		qui do "core/outcomes/sim_bcr.do"
		*mata: _matrix_list(bL8_BCR, rbL8_BCR, cbL8_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)

	di "L8S - Treatment Duration"
		qui do "core/outcomes/sim_txd.do"	
		*mata: _matrix_list(bL8_TXD, rbL8_TXD, cbL8_TXD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L8S - Overall Survival"
		qui do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)		
	
	di "L8S - Mortality"
		qui do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
		
**********
*Line 8 End (L8E)
	mata: OMC = 17
		
	di "L8E - Age"
		qui do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)

	di "L8E - Treatment-free Interval"
		qui do "core/outcomes/sim_tfi.do"				
		*mata: _matrix_list(bL8_TFI, rbL8_TFI, cbL8_TFI)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)	

	di "L8E - Overall Survival"
		qui do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)		

	di "L8E - Mortality" 
		qui do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Line 9 Start (L9S) 
	mata: OMC = 18
	mata: Line = 9
	
	di "L9S - Age" 
		qui do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)
		
	di "L9S - Treatment Regimen"
		qui do "core/outcomes/sim_txr.do"
		*mata: _matrix_list(mTXR, rmTXR, cmTXR)
		
	di "L9E - Best Clinical Response"
		qui do "core/outcomes/sim_bcr.do"
		*mata: _matrix_list(bL9_BCR, rbL9_BCR, cbL9_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)

	di "L9S - Treatment Duration"
		qui do "core/outcomes/sim_txd.do"	
		*mata: _matrix_list(bL9_TXD, rbL9_TXD, cbL9_TXD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L9S - Overall Survival"
		qui do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)	
	
	di "L9S - Mortality"
		qui do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
		
**********
*Line 9 End (L9E)
	mata: OMC = 19
	
	di "L9E - Age"
		qui do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
		
	di "L9E - Overall Survival"
		qui do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L9E - Mortality" 
		*Everyone still alive dies at predicted OS or Limit
		mata {
			// Identify patients still alive
			idxAlive = selectindex(mMOR[., OMC-1] :== 0)
			
			if (rows(idxAlive) > 0) {
				// Cap OS at age limit for those who would exceed it
				vExceedsLimit = ((mAge[idxAlive, 1] :+ (mOS[idxAlive, OMC] :/ 12)) :> Limit)
				mOS[idxAlive, OMC] = vExceedsLimit :* ((Limit :- mAge[idxAlive, 1]) :* 12) :+ 
									  (!vExceedsLimit) :* mOS[idxAlive, OMC]
				
				// Mark everyone as dead
				mMOR[idxAlive, OMC] = J(rows(idxAlive), 1, 1)
				
				// Set outcome time and mortality indicator
				mOC[idxAlive, 1] = mOS[idxAlive, OMC]
				mOC[idxAlive, 2] = J(rows(idxAlive), 1, 1)
			}
		}
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
		*mata: _matrix_list(mOC, rmOC, cmOC)
		
end
