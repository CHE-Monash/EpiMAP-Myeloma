**********
	*EpiMAP Myeloma - Simulation Engine
**********

capture program drop simulation
program define simulation

di "Running simulation"

**********
*Diagnosis (DN)
	mata: OMC = 1
	mata: Line = 0
	mata: LX = 0
			
	di "DN - SCT"	
		quietly do "core/outcomes/sim_asct_dn.do"
		*mata: _matrix_list(bDN_SCT, rbDN_SCT, cbDN_SCT)
		*mata: _matrix_list(mSCT, rmSCT, cmSCT)		
		
	di "DN - Treatment-free Interval"		
		quietly do "core/outcomes/sim_tfi_dn.do"
		*mata: _matrix_list(bDN_TFI, rbDN_TFI, cbDN_TFI)
		*mata: _matrix_list(mTFI, rmTFI, cmTFI)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)	

	di "DN - Overall Survival" 
		quietly do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "DN - Mortality"
		quietly do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Line 1 Start (L1S)	
	mata: OMC = 2

	di "L1S - Age"
		quietly do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)

	di "L1S - Treatment Regimen"
		quietly do "core/outcomes/sim_txr.do"
		*mata: _matrix_list(bL1_TXR, rbL1_TXR, cbL1_TXR)
		*mata: _matrix_list(mTXR, rmTXR, cmTXR)
		
		if ("$Data" == "Population" & $Line == 1) {
			exit
		}

	di "L1S - Treatment Duration"
		quietly do "core/outcomes/sim_txd_l1.do"
		*mata: _matrix_list(bL1_TXD, rbL1_TXD, cbL1_TXD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L1S - Overall Survival"
		quietly do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L1S - Mortality"
		quietly do "core/outcomes/sim_mort.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Line 1 End (L1E)
	mata: OMC = 3
	mata: Line = 1
	mata: LX = 1

	di "L1E - Age"
		quietly do "core/outcomes/sim_age.do"
		*mata: _matrix_list(mAge, rmAge, cmAge)	

	di "L1E - Best Clinical Response"
		quietly do "core/outcomes/sim_bcr.do"
		*mata: _matrix_list(bL1_BCR, rbL1_BCR, cbL1_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)

	di "L1E - SCT"
		quietly do "core/outcomes/sim_asct_l1.do"			
		*mata: _matrix_list(bL1_SCT, rbL1_SCT, cbL1_SCT)
		*mata: _matrix_list(mSCT, rmSCT, cmSCT)

	di "L1E - SCT Best Clinical Response"
		quietly do "core/outcomes/sim_bcr_asct.do"			
		*mata: _matrix_list(bSCT_BCR, rbSCT_BCR, cbSCT_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)

	di "L1E - MNT"
		quietly do "core/outcomes/sim_mnt.do"		
		*mata: _matrix_list(bMNT, rbMNT, cbMNT)
		*mata: _matrix_list(vMNT)

	di "L1E - Treatment-free Interval"
		quietly do "core/outcomes/sim_tfi_l1.do"			
		*mata: _matrix_list(bL1_TFI, rbL1_TFI, cbL1_TFI)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)	
	
	di "L1E - Overall Survival"
		quietly do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L1E - Mortality"
		quietly do "core/outcomes/sim_mort.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Line 2 Start (L2S)	
	mata: OMC = 4
		
	di "L2S - Age"
		quietly do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
		
	di "L2S - Treatment Regimen"
		quietly do "core/outcomes/sim_txr.do"
		*mata: _matrix_list(bL2_TXR, rbL2_TXR, cbL2_TXR)
		*mata: _matrix_list(mTXR, rmTXR, cmTXR)
		
		if ("$Data" == "Population" & $Line == 2) {
			exit
		}

	di "L2S - Treatment Duration"
		quietly do "core/outcomes/sim_txd.do"
		*mata: _matrix_list(bL2_TXD, rbL2_TXD, cbL2_TXD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)

	di "L2S - Overall Survival"
		quietly do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L2S - Mortality"
		quietly do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
	
**********
*Line 2 End (L2E)	
	mata: OMC = 5
	mata: Line = 2
	mata: LX = 2

	di "L2E - Age"
		quietly do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
	
	di "L2E - Best Clinical Response"
		quietly do "core/outcomes/sim_bcr.do"
		*mata: _matrix_list(bL2_BCR, rbL2_BCR, cbL2_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)

	di "L2E - Treatment-free Interval"
		quietly do "core/outcomes/sim_tfi.do"
		*mata: _matrix_list(bL2_TFI, rbL2_TFI, cbL2_TFI)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)

	di "L2E - Overall Survival"
		quietly do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS))

	di "L2E - Mortality"
		quietly do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Line 3 Start (L3S)
	mata: OMC = 6
		
	di "L3S - Age"
		quietly do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)
		
	di "L3S - Treatment Regimen"
		quietly do "core/outcomes/sim_txr.do"
		*mata: _matrix_list(bL3_TXR, rbL3_TXR, cbL3_TXR)
		*mata: _matrix_list(mTXR, rmTXR, cmTXR)
		
		if ("$Data" == "Population" & $Line == 3) {
			exit
		}
		
/*		mata: oL3_TXR = 0
		scalar o = "oL3_TXR"
		forval i = 1/`=Obs' {
			mata {
				if (mMOR[`i',`=OMC'-1] == 0) mCore[`i',cCR] = 0
				if (mMOR[`i',`=OMC'-1] == 0) mCR[`i',`=LX'+2] = 0
				if (mMOR[`i',`=OMC'-1] != 0) mCore[`i',cCR] = .					
			}
		}
*/
		*mata: _matrix_list(mL3_TXR, rmL3_TXR, cmL3_TXR)
				
	di "L3S - Treatment Duration"
		quietly do "core/outcomes/sim_txd.do"	
		*mata: _matrix_list(bL3_TXD, rbL3_TXD, cbL3_TXD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L3S - Overall Survival"
		quietly do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L3S - Mortality"
		quietly do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Line 3 End (L3E)		
	mata: OMC = 7
	mata: Line = 3
	mata: LX = 3
	
	di "L3E - Age"
		quietly do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)
	
	di "L3E - Best Clinical Response"
		quietly do "core/outcomes/sim_bcr.do"
		*mata: _matrix_list(bL3_BCR, rbL3_BCR, cbL3_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		
	di "L3E - Treatment-free Interval"
		quietly do "core/outcomes/sim_tfi.do"		
		*mata: _matrix_list(bL3_TFI, rbL3_TFI, cbL3_TFI)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)

	di "L3E - Overall Survival"
		quietly do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L3E - Mortality"
		quietly do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Line 4 Start (L4S)	
	mata: OMC = 8
		
	di "L4S - Age"
		quietly do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)

	di "L4S - Treatment Regimen"
		quietly do "core/outcomes/sim_txr.do"
		*mata: _matrix_list(bL4_TXR, rbL4_TXR, cbL4_TXR)
		*mata: _matrix_list(mTXR, rmTXR, cmTXR)
		
		if ("$Data" == "Population" & $Line == 4) {
			exit
		}

/*		mata {
			// Identify alive and dead patients
			idxAlive = selectindex(mMOR[., OMC-1] :== 0)
			idxDead = selectindex(mMOR[., OMC-1] :!= 0)
			
			// Set treatment regimen to 0 for alive patients
			if (rows(idxAlive) > 0) {
				mTXR[idxAlive, LX+1] = J(rows(idxAlive), 1, 0)
			}
			
			// Set to missing for dead patients
			if (rows(idxDead) > 0) {
				mTXR[idxDead, LX+1] = J(rows(idxDead), 1, .)
			}
		}
*/		
		*mata: _matrix_list(mL4_TXR, rmL4_TXR, cmL4_TXR)
			
	di "L4S - Treatment Duration"
		quietly do "core/outcomes/sim_txd.do"
		*mata: _matrix_list(bL4_TXD, rbL4_TXD, cbL4_TXD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L4S - Overall Survival"
		quietly do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L4S - Mortality"
		quietly do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Line 4 End (L4E)
	mata: OMC = 9
	mata: Line = 4
	mata: LX = 4
		
	di "L4E - Age"
		quietly do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)
	
	di "L4E - Best Clinical Response"
		quietly do "core/outcomes/sim_bcr.do"
		*mata: _matrix_list(bL4_BCR, rbL4_BCR, cbL4_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		
	di "L4E - Treatment-free Interval"

		quietly do "core/outcomes/sim_tfi.do"	
		*mata: _matrix_list(bL4_TFI, rbL4_TFI, cbL4_TFI)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)

	di "L4E - Overall Survival"
		quietly do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L4E - Mortality"
		quietly do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
			
**********
*Line 5 Start (L5S)	
	mata: OMC = 10
		
	di "L5S - Age" 
		quietly do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)
		
	di "L5S - Treatment Regimen"
	mata {
		// Identify alive and dead patients
		idxAlive = selectindex(mMOR[., OMC-1] :== 0)
		idxDead = selectindex(mMOR[., OMC-1] :!= 0)
		
		// Set treatment regimen to 0 for alive patients
		if (rows(idxAlive) > 0) {
			mTXR[idxAlive, LX+1] = J(rows(idxAlive), 1, 0)
		}
		
		// Set to missing for dead patients
		if (rows(idxDead) > 0) {
			mTXR[idxDead, LX+1] = J(rows(idxDead), 1, .)
		}
	}

	di "L5S - Treatment Duration"
		scalar m = "mL5_TXD"
		scalar b = "bLX_TXD"
		scalar c = "cSU_"
		quietly do "core/outcomes/sim_txd.do"		
		*mata: _matrix_list(bL5_TXD, rbL5_TXD, cbL5_TXD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L5S - Overall Survival"
		quietly do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L5S - Mortality"
		quietly do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Line 5 End (L5E)	
	mata: OMC = 11
	mata: Line = 5
	mata: LX = 5
		
	di "L5E - Age"
		quietly do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)

	di "L5E - Best Clinical Response"
		quietly do "core/outcomes/sim_bcr.do"
		*mata: _matrix_list(bL5_BCR, rbL5_BCR, cbL5_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
	
	di "L5E - Treatment-free Interval"
		quietly do "core/outcomes/sim_tfi.do"				
		*mata: _matrix_list(bL5_TFI, rbL5_TFI, cbL5_TFI)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)

	di "L5E - Overall Survival"
		quietly do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L5E - Mortality"
		quietly do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
		
**********
*Line 6 Start (L6S) 		
	mata: OMC = 12
	
	di "L6S - Age" 
		quietly do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)
		
	di "L6S - Treatment Regimen"
		mata {
			// Identify alive and dead patients
			idxAlive = selectindex(mMOR[., OMC-1] :== 0)
			idxDead = selectindex(mMOR[., OMC-1] :!= 0)
			
			// Set treatment regimen to 0 for alive patients
			if (rows(idxAlive) > 0) {
				mTXR[idxAlive, LX+1] = J(rows(idxAlive), 1, 0)
			}
			
			// Set to missing for dead patients
			if (rows(idxDead) > 0) {
				mTXR[idxDead, LX+1] = J(rows(idxDead), 1, .)
			}
		}

	di "L6S - Treatment Duration"
		scalar m = "mL6_TXD"
		scalar b = "bLX_TXD" 		
		scalar c = "cSU_"
		quietly do "core/outcomes/sim_txd.do"	
		*mata: _matrix_list(bL6_TXD, rbL6_TXD, cbL6_TXD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L6S - Overall Survival"
		quietly do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L6S - Mortality"
		quietly do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Line 6 End (L6E)
	mata: OMC = 13
	mata: Line = 6
	mata: LX = 6
		
	di "L6E - Age"
		quietly do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)
	
	di "L6E - Best Clinical Response"
		quietly do "core/outcomes/sim_bcr.do"
		*mata: _matrix_list(bL6_BCR, rbL6_BCR, cbL6_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		
	di "L6E - Treatment-free Interval"
		quietly do "core/outcomes/sim_tfi.do"				
		*mata: _matrix_list(bL6_TFI, rbL6_TFI, cbL6_TFI)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)	

	di "L6E - Overall Survival"
		quietly do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L6E - Mortality" 
		quietly do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
					
**********
*Line 7 Start (L7S) 	
	mata: OMC = 14
		
	di "L7S - Age" 
		quietly do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)
		
	di "L7S - Treatment Regimen"
		mata {
			// Identify alive and dead patients
			idxAlive = selectindex(mMOR[., OMC-1] :== 0)
			idxDead = selectindex(mMOR[., OMC-1] :!= 0)
			
			// Set treatment regimen to 0 for alive patients
			if (rows(idxAlive) > 0) {
				mTXR[idxAlive, LX+1] = J(rows(idxAlive), 1, 0)
			}
			
			// Set to missing for dead patients
			if (rows(idxDead) > 0) {
				mTXR[idxDead, LX+1] = J(rows(idxDead), 1, .)
			}
		}

	di "L7S - Chemo Duration"
		quietly do "core/outcomes/sim_txd.do"	
		*mata: _matrix_list(bL7_TXD, rbL7_TXD, cbL7_TXD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L7S - Overall Survival"
		quietly do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)
		
	di "L7S - Mortality"
		quietly do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Line 7 End (L7E)
	mata: OMC = 15
	mata: LX = 7
		
	di "L7E - Age"
		quietly do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)
	
	di "L7E - Best Clinical Response" 
		quietly do "core/outcomes/sim_bcr.do"
		*mata: _matrix_list(bL7_BCR, rbL7_BCR, cbL7_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR))
		
	di "L7E - Treatment-free Interval"
		quietly do "core/outcomes/sim_tfi.do"				
		*mata: _matrix_list(bL7_TFI, rbL7_TFI, cbL7_TFI)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)		

	di "L7E - Overall Survival"
		quietly do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L7E - Mortality" 
		quietly do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
		
**********
*Line 8 Start (L8S) 
	mata: OMC = 16
		
	di "L8S - Age" 
		quietly do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)
		
	di "L8S - Treatment Regimen"
		mata {
			// Identify alive and dead patients
			idxAlive = selectindex(mMOR[., OMC-1] :== 0)
			idxDead = selectindex(mMOR[., OMC-1] :!= 0)
			
			// Set treatment regimen to 0 for alive patients
			if (rows(idxAlive) > 0) {
				mTXR[idxAlive, LX+1] = J(rows(idxAlive), 1, 0)
			}
			
			// Set to missing for dead patients
			if (rows(idxDead) > 0) {
				mTXR[idxDead, LX+1] = J(rows(idxDead), 1, .)
			}
		}

	di "L8S - Treatment Duration"
		quietly do "core/outcomes/sim_txd.do"	
		*mata: _matrix_list(bL8_TXD, rbL8_TXD, cbL8_TXD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L8S - Overall Survival"
		quietly do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L8S - Mortality"
		quietly do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Line 8 End (L8E)
	mata: OMC = 17
	mata: LX = 8
		
	di "L8E - Age"
		quietly do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)
	
	di "L8E - Best Clinical Response"
		quietly do "core/outcomes/sim_bcr.do"
		*mata: _matrix_list(bL8_BCR, rbL8_BCR, cbL8_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		
	di "L8E - Treatment-free Interval"
		quietly do "core/outcomes/sim_tfi.do"				
		*mata: _matrix_list(bL8_TFI, rbL8_TFI, cbL8_TFI)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)	

	di "L8E - Overall Survival"
		quietly do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L8E - Mortality" 
		quietly do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Line 9 Start (L9S) 
	mata: OMC = 18
		
	di "L9S - Age" 
		quietly do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)
		
	di "L9S - Treatment Regimen"
		mata {
			// Identify alive and dead patients
			idxAlive = selectindex(mMOR[., OMC-1] :== 0)
			idxDead = selectindex(mMOR[., OMC-1] :!= 0)
			
			// Set treatment regimen to 0 for alive patients
			if (rows(idxAlive) > 0) {
				mTXR[idxAlive, LX+1] = J(rows(idxAlive), 1, 0)
			}
			
			// Set to missing for dead patients
			if (rows(idxDead) > 0) {
				mTXR[idxDead, LX+1] = J(rows(idxDead), 1, .)
			}
		}

	di "L9S - Treatment Duration"
		quietly do "core/outcomes/sim_txd.do"	
		*mata: _matrix_list(bL9_TXD, rbL9_TXD, cbL9_TXD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L9S - Overall Survival"
		quietly do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L9S - Mortality"
		quietly do "core/outcomes/sim_mort.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Line 9 End (L9E)
	mata: OMC = 19
	mata: LX = 9
		
	di "L9E - Age"
		quietly do "core/outcomes/sim_age.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
	
	di "L9E - Best Clinical Response"
		quietly do "core/outcomes/sim_bcr.do"
		*mata: _matrix_list(bL9_BCR, rbL9_BCR, cbL9_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		
	di "L9E - Overall Survival"
		quietly do "core/outcomes/sim_os.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L9E - Mortality" 
		*Everyone still alive dies at predicted OS or Limit
		mata {
			// Identify patients still alive
			idxAlive = selectindex(mMOR[., OMC-1] :== 0)
			
			if (rows(idxAlive) > 0) {
				// Cap OS at age limit for those who would exceed it
				vExceedsLimit = ((mAge[idxAlive, 1] :+ mOS[idxAlive, OMC]) :> Limit)
				mOS[idxAlive, OMC] = vExceedsLimit :* (Limit :- mAge[idxAlive, 1]) :+ 
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
