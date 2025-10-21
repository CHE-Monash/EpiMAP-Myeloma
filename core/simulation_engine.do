**********
	*EpiMAP Myeloma - Simulation Engine
**********

capture program drop simulation
program define simulation

di "Running simulation"

**********
*Diagnosis (DN)
	scalar OMC = 1
	scalar Line = 0
	scalar LX = 0
	
	mata: OMC = 1
	mata: Line = 0
	mata: LX = 0
			
	di "DN - SCT"	
		quietly do "core/outcomes/SIM SCT DN Vector.do"
		*mata: _matrix_list(bDN_SCT, rbDN_SCT, cbDN_SCT)
		*mata: _matrix_list(mSCT, rmSCT, cmSCT)		
		
	di "DN - Treatment-free Interval"		
		quietly do "core/outcomes/SIM TFI DN Vector.do"
		*mata: _matrix_list(bDN_CI, rbDN_CI, cbDN_CI) // Needs renaming to bDN_TFI / bTFI_DN
		*mata: _matrix_list(mTFI, rmTFI, cmTFI)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)	

	di "DN - Overall Survival" 
		quietly do "core/outcomes/SIM OS Vector.do"
		*mata: _matrix_list(bOS, rbOS, cbOS)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "DN - Mortality"
		quietly do "core/outcomes/SIM MORT Vector.do"
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 1 Start (L1S)
	scalar OMC = 2
	
	mata: OMC = 2

	di "L1S - Age"
		quietly do "core/outcomes/SIM AGE Vector.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)

	di "L1S - Treatment Regimen"
		quietly do "core/outcomes/SIM CR L1 Vector.do"
		*mata: _matrix_list(bL1_CR, rbL1_CR, cbL1_CR)
		*mata: _matrix_list(mTXD, rmTXD, cmTXD)
		
		if ("$Data" == "Population" & $Line == 1) {
			exit
		}
		
	di "L1S - Chemo Duration"
		scalar m = "mL1_CD"
		scalar c = "cSU_"
		quietly do "core/outcomes/SIM CD L1 Vector.do"
		*mata: _matrix_list(bL1_CD, rbL1_CD, cbL1_CD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L1S - Overall Survival"
		scalar m = "mOS_L1S"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "core/outcomes/SIM OS Vector.do"
		*mata: _matrix_list(mOS_L1S, rmOS_L1S, cmOS_L1S)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L1S - Mortality"
		quietly do "core/outcomes/SIM MORT Vector.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 1 End (L1E)
	scalar OMC = 3
	scalar Line = 1
	scalar LX = 1

	di "L1E - Age"
		quietly do "core/outcomes/SIM AGE.do"
		*mata: _matrix_list(mAge, rmAge, cmAge)	

	di "L1E - Best Clinical Response"
		scalar m = "mL1_BCR"
		scalar b = "bL1_BCR"
		scalar c = "cOL_"
		quietly do "core/outcomes/SIM BCR L1.do"		
		*mata: _matrix_list(mL1_BCR, rmL1_BCR, cmL1_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)

	di "L1E - SCT"
		scalar m = "mL1_SCT"
		scalar b = "bL1_SCT"
		scalar c = "cLO_"
		quietly do "core/outcomes/SIM SCT L1.do"			
		*mata: _matrix_list(mL1_SCT, rmL1_SCT, cmL1_SCT)
		*mata: _matrix_list(mSCT, rmSCT, cmSCT)

	di "L1E - SCT Best Clinical Response"
		scalar m = "mSCT_BCR"
		scalar b = "bSCT_BCR"
		scalar c = "cOL_"
		quietly do "core/outcomes/SIM BCR SCT.do"			
		*mata: _matrix_list(mSCT_BCR, rmSCT_BCR, cmSCT_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)

	di "L1E - MNT"
		scalar m = "mMNT"
		scalar b = "bMNT"
		scalar c = "cLO_"
		quietly do "core/outcomes/SIM MNT.do"		
		*mata: _matrix_list(mMNT, rmMNT, cmMNT)
		*mata: _matrix_list(mCore, rmCore, cmCore)

	di "L1E - Chemo Interval"
		scalar m = "mL1_CI"
		scalar c = "cSU_"
		quietly do "core/outcomes/SIM CI L1.do"			
		*mata: _matrix_list(mL1_CI, rmL1_CI, cmL1_CI)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)		
	
	di "L1E - Overall Survival"
		scalar m = "mOS_L1E"
		scalar b = "bOS"	
		scalar c = "cLO_"
		quietly do "core/outcomes/SIM OS.do"	
		*mata: _matrix_list(mOS_L1E, rmOS_L1E, cmOS_L1E)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L1E - Mortality"
		quietly do "core/outcomes/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 2 Start (L2S)
	scalar OMC = 4
		
	di "L2S - Age"
		quietly do "core/outcomes/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)		
		
	di "L2S - Chemo Regimen"
		scalar m = "mL2_CR"
		scalar b = "bL2_CR"
		scalar o = "oL2_CR"
		scalar c = "cML_"
		quietly do "core/outcomes/SIM CR L2.do"
		*mata: _matrix_list(mL2_CR, rmL2_CR, cmL2_CR)
		
		if ("$Data" == "Population" & $Line == 2) {
			exit
		}

	di "L2S - Chemo Duration"
		scalar m = "mL2_CD"
		scalar b = "bL2_CD"
		scalar c = "cSU_"
		quietly do "core/outcomes/SIM CD L2.do"
		*mata: _matrix_list(mL2_CD, rmL2_CD, cmL2_CD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)

	di "L2S - Overall Survival"
		scalar m = "mOS_L2S"
		scalar b = "bOS"	
		scalar c = "cLO_"
		quietly do "core/outcomes/SIM OS.do"	
		*mata: _matrix_list(mOS_L2S, rmOS_L2S, cmOS_L2S)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L2S - Mortality"
		quietly do "core/outcomes/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
	
**********
*Chemo Line 2 End (L2E)
	scalar OMC = 5
	scalar Line = 2
	scalar LX = 2

	di "L2E - Age"
		quietly do "core/outcomes/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)		
	
	di "L2E - Best Clinical Response"
		scalar m = "mL2_BCR"
		scalar b = "bL2_BCR"
		scalar c = "cOL_"
		quietly do "core/outcomes/SIM BCR L2.do"
		*mata: _matrix_list(mL2_BCR, rmL2_BCR, cmL2_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)

	di "L2E - Chemo Interval"
		scalar m = "mL2_CI"
		scalar b = "bL2_CI"
		scalar c = "cSU_"
		quietly do "core/outcomes/SIM CI L2.do"
		*mata: _matrix_list(mL2_CI, rmL2_CI, cmL2_CI)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)

	di "L2E - Overall Survival"
		scalar m = "mOS_L2E"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "core/outcomes/SIM OS.do"	
		*mata: _matrix_list(mOS_L2E, rmOS_L2E, cmOS_L2E)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L2E - Mortality"
		quietly do "core/outcomes/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 3 Start (L3S)
	scalar OMC = 6
		
	di "L3S - Age"
		quietly do "core/outcomes/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
		
	di "L3S - Chemo Regimen"
		scalar m = "mL3_CR"
		scalar b = "bL3_CR"
		scalar o = "oL3_CR"
		scalar c = "cML_"
		quietly do "core/outcomes/SIM CR L3.do"
		
		if ("$Data" == "Population" & $Line == 3) {
			exit
		}
		
/*		mata: oL3_CR = 0
		scalar o = "oL3_CR"
		forval i = 1/`=Obs' {
			mata {
				if (mMOR[`i',`=OMC'-1] == 0) mCore[`i',cCR] = 0
				if (mMOR[`i',`=OMC'-1] == 0) mCR[`i',`=LX'+2] = 0
				if (mMOR[`i',`=OMC'-1] != 0) mCore[`i',cCR] = .					
			}
		}
*/
		*mata: _matrix_list(mL3_CR, rmL3_CR, cmL3_CR)
				
	di "L3S - Chemo Duration"
		scalar m = "mL3_CD"
		scalar b = "bL3_CD"
		scalar c = "cSU_"
		quietly do "core/outcomes/SIM CD L3.do"	
		*mata: _matrix_list(mL3_CD, rmL3_CD, cmL3_CD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L3S - Overall Survival"
		scalar m = "mOS_L3S"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "core/outcomes/SIM OS.do"
		*mata: _matrix_list(mOS_L3S, rmOS_L3S, cmOS_L3S)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L3S - Mortality"
		quietly do "core/outcomes/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)	

**********
*Chemo Line 3 End (L3E)
	scalar OMC = 7
	scalar Line = 3
	scalar LX = 3
		
	di "L3E - Age"
		quietly do "core/outcomes/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
	
	di "L3E - Best Clinical Response"
		scalar m = "mL3_BCR"
		scalar b = "bL3_BCR"
		scalar c = "cOL_"
		quietly do "core/outcomes/SIM BCR L3.do"		
		*mata: _matrix_list(mL3_BCR, rmL3_BCR, cmL3_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		
	di "L3E - Chemo Interval"
		scalar m = "mL3_CI"
		scalar b = "bL3_CI"
		scalar c = "cSU_"
		quietly do "core/outcomes/SIM CI L3 L4.do"		
		*mata: _matrix_list(mL3_CI, rmL3_CI, cmL3_CI)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)

	di "L3E - Overall Survival"
		scalar m = "mOS_L3E"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "core/outcomes/SIM OS.do"
		*mata: _matrix_list(mOS_L3E, rmOS_L3E, cmOS_L3E)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L3E - Mortality"
		quietly do "core/outcomes/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 4 Start (L4S)
	scalar OMC = 8
		
	di "L4S - Age"
		quietly do "core/outcomes/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)		

	di "L4S - Chemo Regimen"
		scalar m = "mL4_CR"
		scalar b = "bL4_CR"
		scalar o = "oL4_CR"
		scalar c = "cML_"
		quietly do "core/outcomes/SIM CR L4.do"
		
		if ("$Data" == "Population" & $Line == 4) {
			exit
		}

/*		mata: oL4_CR = 0
		scalar o = "oL4_CR"
		forval i = 1/`=Obs' {
			mata {
				if 	(mMOR[`i',`=OMC'-1] == 0) mCore[`i',cCR] = 0
				if 	(mMOR[`i',`=OMC'-1] != 0) mCore[`i',cCR] = .
				if 	(mMOR[`i',`=OMC'-1] == 0) mCR[`i',`=LX'+2] = 0
			}
		}
*/		
		*mata: _matrix_list(mL4_CR, rmL4_CR, cmL4_CR)
			
	di "L4S - Chemo Duration"
		scalar m = "mL4_CD"
		scalar b = "bL4_CD"
		scalar c = "cSU_"
		quietly do "core/outcomes/SIM CD L4.do"		
		*mata: _matrix_list(mL4_CD, rmL4_CD, cmL4_CD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L4S - Overall Survival"
		scalar m = "mOS_L4S"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "core/outcomes/SIM OS.do"
		*mata: _matrix_list(mOS_L4S, rmOS_L4S, cmOS_L4S)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L4S - Mortality"
		quietly do "core/outcomes/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 4 End (L4E)
	scalar OMC = 9
	scalar Line = 4
	scalar LX = 4
		
	di "L4E - Age"
		quietly do "core/outcomes/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
	
	di "L4E - Best Clinical Response"
		scalar m = "mL4_BCR"
		scalar b = "bL4_BCR"
		scalar c = "cOL_"
		quietly do "core/outcomes/SIM BCR L4.do"	
		*mata: _matrix_list(mL4_BCR, rmL4_BCR, cmL4_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		
	di "L4E - Chemo Interval"
		scalar m = "mL4_CI"
		scalar b = "bL4_CI"
		scalar c = "cSU_"
		quietly do "core/outcomes/SIM CI L3 L4.do"	
		*mata: _matrix_list(mL4_CI, rmL4_CI, cmL4_CI)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)

	di "L4E - Overall Survival"
		scalar m = "mOS_L4E"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "core/outcomes/SIM OS.do"	
		*mata: _matrix_list(mOS_L4E, rmOS_L4E, cmOS_L4E)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L4E - Mortality"
		quietly do "core/outcomes/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
			
**********
*Chemo Line 5 Start (L5S)
	scalar OMC = 10
		
	di "L5S - Age" 
		quietly do "core/outcomes/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
		
	di "L5S - Chemo Regimen"
		forval i = 1/`=Obs' {
			mata {
				if 	(mMOR[`i',`=OMC'-1] == 0) mCore[`i',cCR] = 0
				if 	(mMOR[`i',`=OMC'-1] != 0) mCore[`i',cCR] = .
				if 	(mMOR[`i',`=OMC'-1] == 0) mCR[`i',`=LX'+1] = 0
			}
		}

	di "L5S - Chemo Duration"
		scalar m = "mL5_CD"
		scalar b = "bLX_CD"
		scalar c = "cSU_"
		quietly do "core/outcomes/SIM CD LX (5 - 9).do"		
		*mata: _matrix_list(mL5_CD, rmL5_CD, cmL5_CD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L5S - Overall Survival"
		scalar m = "mOS_L5S"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "core/outcomes/SIM OS.do"
		*mata: _matrix_list(mOS_L5S, rmOS_L5S, cmOS_L5S)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L5S - Mortality"
		quietly do "core/outcomes/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 5 End (L5E)	
	scalar OMC = 11
	scalar Line = 5
	scalar LX = 5
		
	di "L5E - Age"
		quietly do "core/outcomes/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	

	di "L5E - Best Clinical Response"
		scalar m = "mL5_BCR"
		scalar b = "bLX_BCR"
		scalar c = "cOL_"
		quietly do "core/outcomes/SIM BCR LX (5 - 9).do"
		*mata: _matrix_list(mL5_BCR, rmL5_BCR, cmL5_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
	
	di "L5E - Chemo Interval"
		scalar m = "mL5_CI"
		scalar b = "bLX_CI"
		scalar c = "cSU_"
		quietly do "core/outcomes/SIM CI LX (5 - 8).do"				
		*mata: _matrix_list(mL5_CI, rmL5_CI, cmL5_CI)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)

	di "L5E - Overall Survival"
		scalar m = "mOS_L5E"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "core/outcomes/SIM OS.do"
		*mata: _matrix_list(mOS_L5E, rmOS_L5E, cmOS_L5E)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L5E - Mortality"
		quietly do "core/outcomes/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
		
**********
*Chemo Line 6 Start (L6S) 
	scalar OMC = 12
		
	di "L6S - Age" 
		quietly do "core/outcomes/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
		
	di "L6S - Chemo Regimen"
		forvalues i = 1/`=Obs' {
			mata {
				if 	(mMOR[`i',`=OMC'-1] == 0) mCore[`i',cCR] = 0
				if 	(mMOR[`i',`=OMC'-1] != 0) mCore[`i',cCR] = .
				if 	(mMOR[`i',`=OMC'-1] == 0) mCR[`i',`=LX'+1] = 0
			}
		}

	di "L6S - Chemo Duration"
		scalar m = "mL6_CD"
		scalar b = "bLX_CD" 		
		scalar c = "cSU_"
		quietly do "core/outcomes/SIM CD LX (5 - 9).do"	
		*mata: _matrix_list(mL6_CD, rmL6_CD, cmL6_CD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L6S - Overall Survival"
		scalar m = "mOS_L6S"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "core/outcomes/SIM OS.do"
		*mata: _matrix_list(mOS_L6S, rmOS_L6S, cmOS_L6S)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L6S - Mortality"
		quietly do "core/outcomes/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 6 End (L6E)
	scalar OMC = 13
	scalar Line = 6
	scalar LX = 6
		
	di "L6E - Age"
		quietly do "core/outcomes/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
	
	di "L6E - Best Clinical Response"
		scalar m = "mL6_BCR"
		scalar b = "bLX_BCR"
		scalar c = "cOL_"
		quietly do "core/outcomes/SIM BCR LX (5 - 9).do"	 	
		*mata: _matrix_list(mL6_BCR, rmL6_BCR, cmL6_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		
	di "L6E - Chemo Interval"
		scalar m = "mL6_CI"
		scalar b = "bLX_CI"
		scalar c = "cSU_"
		quietly do "core/outcomes/SIM CI LX (5 - 8).do"				
		*mata: _matrix_list(mL6_CI, rmL6_CI, cmL6_CI)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)	

	di "L6E - Overall Survival"
		scalar m = "mOS_L6E"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "core/outcomes/SIM OS.do"
		*mata: _matrix_list(mOS_L6E, rmOS_L6E, cmOS_L6E)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L6E - Mortality" 
		quietly do "core/outcomes/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
					
**********
*Chemo Line 7 Start (L7S) 
	scalar OMC = 14
		
	di "L7S - Age" 
		quietly do "core/outcomes/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
		
	di "L7S - Chemo Regimen"
		forvalues i = 1/`=Obs' {
			mata {
				if 	(mMOR[`i',`=OMC'-1] == 0) mCore[`i',cCR] = 0
				if 	(mMOR[`i',`=OMC'-1] != 0) mCore[`i',cCR] = .
				if 	(mMOR[`i',`=OMC'-1] == 0) mCR[`i',`=LX'+1] = 0
			}
		}

	di "L7S - Chemo Duration"
		scalar m = "mL7_CD"
		scalar b = "bLX_CD" 		
		scalar c = "cSU_"
		quietly do "core/outcomes/SIM CD LX (5 - 9).do"	
		*mata: _matrix_list(mL7_CD, rmL7_CD, cmL7_CD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L7S - Overall Survival"
		scalar m = "mOS_L7S"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "core/outcomes/SIM OS.do"
		*mata: _matrix_list(mOS_L7S, rmOS_L7S, cmOS_L7S)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L7S - Mortality"
		quietly do "core/outcomes/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 7 End (L7E)
	scalar OMC = 15
	scalar LX = 7
		
	di "L7E - Age"
		quietly do "core/outcomes/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
	
	di "L7E - Best Clinical Response" 
		scalar m = "mL7_BCR"
		scalar b = "bLX_BCR"
		scalar c = "cOL_"
		quietly do "core/outcomes/SIM BCR LX (5 - 9).do"	 	
		*mata: _matrix_list(mL7_BCR, rmL7_BCR, cmL7_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		
	di "L7E - Chemo Interval"
		scalar m = "mL7_CI"
		scalar b = "bLX_CI"
		scalar c = "cSU_"
		quietly do "core/outcomes/SIM CI LX (5 - 8).do"				
		*mata: _matrix_list(mL7_CI, rmL7_CI, cmL7_CI)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)		

	di "L7E - Overall Survival"
		scalar m = "mOS_L7E"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "core/outcomes/SIM OS.do"
		*mata: _matrix_list(mOS_L7E, rmOS_L7E, cmOS_L7E)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L7E - Mortality" 
		quietly do "core/outcomes/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
		
**********
*Chemo Line 8 Start (L8S) 
	scalar OMC = 16
		
	di "L8S - Age" 
		quietly do "core/outcomes/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
		
	di "L8S - Chemo Regimen"
		forvalues i = 1/`=Obs' {
			mata {
				if 	(mMOR[`i',`=OMC'-1] == 0) mCore[`i',cCR] = 0
				if 	(mMOR[`i',`=OMC'-1] != 0) mCore[`i',cCR] = .
				if 	(mMOR[`i',`=OMC'-1] == 0) mCR[`i',`=LX'+1] = 0
			}
		}

	di "L8S - Chemo Duration"
		scalar m = "mL8_CD"
		scalar b = "bLX_CD" 		
		scalar c = "cSU_"
		quietly do "core/outcomes/SIM CD LX (5 - 9).do"	
		*mata: _matrix_list(mL8_CD, rmL8_CD, cmL8_CD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L8S - Overall Survival"
		scalar m = "mOS_L8S"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "core/outcomes/SIM OS.do"
		*mata: _matrix_list(mOS_L8S, rmOS_L8S, cmOS_L8S)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L8S - Mortality"
		quietly do "core/outcomes/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 8 End (L8E)
	scalar OMC = 17
	scalar LX = 8
		
	di "L8E - Age"
		quietly do "core/outcomes/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
	
	di "L8E - Best Clinical Response"
		scalar m = "mL8_BCR"
		scalar b = "bLX_BCR"
		scalar c = "cOL_"
		quietly do "core/outcomes/SIM BCR LX (5 - 9).do"	
		*mata: _matrix_list(mL8_BCR, rmL8_BCR, cmL8_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		
	di "L8E - Chemo Interval"
		scalar m = "mL8_CI"
		scalar b = "bLX_CI"
		scalar c = "cSU_"
		quietly do "core/outcomes/SIM CI LX (5 - 8).do"				
		*mata: _matrix_list(mL8_CI, rmL8_CI, cmL8_CI)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)	

	di "L8E - Overall Survival"
		scalar m = "mOS_L8E"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "core/outcomes/SIM OS.do"
		*mata: _matrix_list(mOS_L8E, rmOS_L8E, cmOS_L8E)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L8E - Mortality" 
		quietly do "core/outcomes/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)	

**********
*Chemo Line 9 Start (L9S) 
	scalar OMC = 18
		
	di "L9S - Age" 
		quietly do "core/outcomes/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
		
	di "L9S - Chemo Regimen"
		forvalues i = 1/`=Obs' {
			mata {
				if 	(mMOR[`i',`=OMC'-1] == 0) mCore[`i',cCR] = 0
				if 	(mMOR[`i',`=OMC'-1] != 0) mCore[`i',cCR] = .
				if 	(mMOR[`i',`=OMC'-1] == 0) mCR[`i',`=LX'+1] = 0
			}
		}

	di "L9S - Chemo Duration"
		scalar m = "mL9_CD"
		scalar b = "bLX_CD" 		
		scalar c = "cSU_"
		quietly do "core/outcomes/SIM CD LX (5 - 9).do"	
		*mata: _matrix_list(mL9_CD, rmL9_CD, cmL9_CD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L9S - Overall Survival"
		scalar m = "mOS_L9S"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "core/outcomes/SIM OS.do"
		*mata: _matrix_list(mOS_L9S, rmOS_L9S, cmOS_L9S)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L9S - Mortality"
		quietly do "core/outcomes/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 9 End (L9E)
	scalar OMC = 19
	scalar LX = 9
		
	di "L9E - Age"
		quietly do "core/outcomes/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
	
	di "L9E - Best Clinical Response"
		scalar m = "mL9_BCR"
		scalar b = "bLX_BCR"
		scalar c = "cOL_"
		quietly do "core/outcomes/SIM BCR LX (5 - 9).do"
		*mata: _matrix_list(mL9_BCR, rmL9_BCR, cmL9_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		
	di "L9E - Overall Survival"
		scalar m = "mOS_L9E"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "core/outcomes/SIM OS.do"
		*mata: _matrix_list(mOS_L9E, rmOS_L9E, cmOS_L9E)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L9E - Mortality" 
		*Everyone still alive dies at predicted OS or Limit
		forvalues i = 1/`=Obs'{
			mata {
				if (mMOR[`i',`=OMC'-1] == 0) { // Alive filter
					if	((mAge[`i',1] + mOS[`i',`=OMC']) > `=Limit') mOS[`i',`=OMC'] = `=Limit' - mAge[`i',1] // Set mOS to max if Age > Limit
					mMOR[`i',`=OMC'] = 1 // Patient dies
					mOC[`i',1] = mOS[`i',`=OMC'] // Set OC Time
					mOC[`i',2] = 1 // Set OC Outcome
				}		
			}	
		}
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
		*mata: _matrix_list(mOC, rmOC, cmOC)		

end
