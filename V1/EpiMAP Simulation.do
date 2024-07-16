**********
	*BCR MM - Simulation
**********

	clear
	clear matrix
	set more off
	scalar drop _all

**********
	*Copied below from Core - Simulation 230419
	*Turned off CR for L3+ 
	*Using BCRc in L1_CI, LX_CI, L3_BCR+, L4_CD+
	*BCR SCT - now looks to see if there are 18 or 19 coefficients (very few BCR ==6)
**********

forvalues bs = 1/1 {
	clear mata
	scalar drop _all
	
	*Load risk equation coefficients
		mata: mata matuse "Data/Coefficients/EpiMAP Coefficients `bs'"

	*Number of patients
		scalar Obs = 10000
		
	*Age limit
		scalar Limit = 100
		
	*Draw sample of patients from MRDR, save in Stata matrices
{	
		use "Data/EpiMAP V Wide MI.dta"
		preserve
		bsample `=Obs'
		*stID
			mkmat ID, matrix(stID)
		*stAge
			mkmat Age, matrix(stAge)
		*stPts
			mkmat ID Age Male ECOGcc ISS, matrix(stPts)
		*stCom	
			mkmat Age70 Age75 CMCard CMPulm CMDiab CMLive CMPNeu CMMali, matrix(stCom)
		restore
	
	*mCore - Core patient matrix
{
		mata: mCore = st_matrix("stPts")
		mata: rmCore = st_matrixrowstripe("stPts")
		mata: cmCore = st_matrixcolstripe("stPts")
		
		*Add columns
			*SCT
				mata: mCore = mCore , J(rows(mCore),1,0)
			*MNT
				mata: mCore = mCore , J(rows(mCore),1,0)
			*CR & CD
				mata: mCore = mCore , J(rows(mCore),2,.)
			*BCR (set to 5 during Line 0)
				mata: mCore = mCore , J(rows(mCore),1,5)
	
		*Add rows
			mata: cmCore = cmCore \ J(5,cols(cmCore),"")
			*Label
				mata: cmCore[6,2] = "SCT"
				mata: cmCore[7,2] = "MNT"
				mata: cmCore[8,2] = "CR"
				mata: cmCore[9,2] = "CD"
				mata: cmCore[10,2] = "BCR"
				
		*Add mata column references
				mata: cAge = 2
				mata: cMale = 3
				mata: cECOGc = 4
				mata: cISS = 5
				mata: cSCT = 6
				mata: cMNT = 7
				mata: cCR = 8
				mata: cCD = 9
				mata: cBCR = 10
}	
		*mata: _matrix_list(mCore, rmCore, cmCore)
		
	*mCom - Comorbidity matrix (for SCT only)
{
		mata: mCom = st_matrix("stCom")
		mata: rmCom = st_matrixrowstripe("stCom")
		mata: cmCom = st_matrixcolstripe("stCom")
}
		*mata: _matrix_list(mCom, rmCom, cmCom)
	
	*mAge - Age (at event)
{	
		mata: mAge = st_matrix("stID")
		mata: mAge = mAge , st_matrix("stAge")
		mata: rmAge = st_matrixrowstripe("stID")
		mata: cmAge = st_matrixcolstripe("stID")
		mata: cmAge = cmAge \ st_matrixcolstripe("stAge")
		
		*Add columns
			mata: mAge = mAge , J(rows(mAge),18,.)
				
		*Add rows
			mata: cmAge = cmAge \ J(18,cols(cmAge),"")
			*Label
				mata: cmAge[2,2] = "Age_DN"
				mata: cmAge[3,2] = "Age_L1S"
				mata: cmAge[4,2] = "Age_L1E"
				mata: cmAge[5,2] = "Age_L2S"
				mata: cmAge[6,2] = "Age_L2E"
				mata: cmAge[7,2] = "Age_L3S"
				mata: cmAge[8,2] = "Age_L3E"
				mata: cmAge[9,2] = "Age_L4S"
				mata: cmAge[10,2] = "Age_L4E"
				mata: cmAge[11,2] = "Age_L5S"
				mata: cmAge[12,2] = "Age_L5E"
				mata: cmAge[13,2] = "Age_L6S"
				mata: cmAge[14,2] = "Age_L6E"
				mata: cmAge[15,2] = "Age_L7S"
				mata: cmAge[16,2] = "Age_L7E"
				mata: cmAge[17,2] = "Age_L8S"
				mata: cmAge[18,2] = "Age_L8E"
				mata: cmAge[19,2] = "Age_L9S"
				mata: cmAge[20,2] = "Age_L9E"
}		
		*mata: _matrix_list(mAge, rmAge, cmAge)		

	*mOS - Overall Survival (from DN)
{	
		mata: mOS = st_matrix("stID")
		mata: rmOS = st_matrixrowstripe("stID")
		mata: cmOS = st_matrixcolstripe("stID")
		
		*Add columns
			mata: mOS = mOS , J(rows(mOS),19,.)
				
		*Add rows
			mata: cmOS = cmOS \ J(19,cols(cmOS),"")
			*Label
				mata: cmOS[2,2] = "OS_DN"
				mata: cmOS[3,2] = "OS_L1S"
				mata: cmOS[4,2] = "OS_L1E"
				mata: cmOS[5,2] = "OS_L2S"
				mata: cmOS[6,2] = "OS_L2E"
				mata: cmOS[7,2] = "OS_L3S"
				mata: cmOS[8,2] = "OS_L3E"
				mata: cmOS[9,2] = "OS_L4S"
				mata: cmOS[10,2] = "OS_L4E"
				mata: cmOS[11,2] = "OS_L5S"
				mata: cmOS[12,2] = "OS_L5E"
				mata: cmOS[13,2] = "OS_L6S"
				mata: cmOS[14,2] = "OS_L6E"
				mata: cmOS[15,2] = "OS_L7S"
				mata: cmOS[16,2] = "OS_L7E"
				mata: cmOS[17,2] = "OS_L8S"
				mata: cmOS[18,2] = "OS_L8E"
				mata: cmOS[19,2] = "OS_L9S"
				mata: cmOS[20,2] = "OS_L9E"
}			
		*mata: _matrix_list(mOS, rmOS, cmOS)		

	*mTNE - Time to Next Event
{	
		mata: mTNE = st_matrix("stID")
		mata: rmTNE = st_matrixrowstripe("stID")
		mata: cmTNE = st_matrixcolstripe("stID")
		
		*Add columns
			mata: mTNE = mTNE , J(rows(mTNE),19,.)
				
		*Add rows
			mata: cmTNE = cmTNE \ J(19,cols(cmTNE),"")
			*Label
				mata: cmTNE[2,2] = "TNE_DN"
				mata: cmTNE[3,2] = "TNE_L1S"
				mata: cmTNE[4,2] = "TNE_L1E"
				mata: cmTNE[5,2] = "TNE_L2S"
				mata: cmTNE[6,2] = "TNE_L2E"
				mata: cmTNE[7,2] = "TNE_L3S"
				mata: cmTNE[8,2] = "TNE_L3E"
				mata: cmTNE[9,2] = "TNE_L4S"
				mata: cmTNE[10,2] = "TNE_L4E"
				mata: cmTNE[11,2] = "TNE_L5S"
				mata: cmTNE[12,2] = "TNE_L5E"
				mata: cmTNE[13,2] = "TNE_L6S"
				mata: cmTNE[14,2] = "TNE_L6E"
				mata: cmTNE[15,2] = "TNE_L7S"
				mata: cmTNE[16,2] = "TNE_L7E"
				mata: cmTNE[17,2] = "TNE_L8S"
				mata: cmTNE[18,2] = "TNE_L8E"
				mata: cmTNE[19,2] = "TNE_L9S"
				mata: cmTNE[20,2] = "TNE_L9E"
}			
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)	
		
	*mTSD - Time Since Diagnosis
{	
		mata: mTSD = st_matrix("stID")
		mata: rmTSD = st_matrixrowstripe("stID")
		mata: cmTSD = st_matrixcolstripe("stID")
		
		*Add columns
			mata: mTSD = mTSD , J(rows(mTSD),19,.)
				
		*Add rows
			mata: cmTSD = cmTSD \ J(19,cols(cmTSD),"")
			*Label
				mata: cmTSD[2,2] = "TSD_DN"
				mata: cmTSD[3,2] = "TSD_L1S"
				mata: cmTSD[4,2] = "TSD_L1E"
				mata: cmTSD[5,2] = "TSD_L2S"
				mata: cmTSD[6,2] = "TSD_L2E"
				mata: cmTSD[7,2] = "TSD_L3S"
				mata: cmTSD[8,2] = "TSD_L3E"
				mata: cmTSD[9,2] = "TSD_L4S"
				mata: cmTSD[10,2] = "TSD_L4E"
				mata: cmTSD[11,2] = "TSD_L5S"
				mata: cmTSD[12,2] = "TSD_L5E"
				mata: cmTSD[13,2] = "TSD_L6S"
				mata: cmTSD[14,2] = "TSD_L6E"
				mata: cmTSD[15,2] = "TSD_L7S"
				mata: cmTSD[16,2] = "TSD_L7E"
				mata: cmTSD[17,2] = "TSD_L8S"
				mata: cmTSD[18,2] = "TSD_L8E"
				mata: cmTSD[19,2] = "TSD_L9S"
				mata: cmTSD[20,2] = "TSD_L9E"
}			
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)		
		
	*mMOR - Mortality (=1 if patient dies before next event) 
{	
		mata: mMOR= st_matrix("stID")
		mata: rmMOR= st_matrixrowstripe("stID")
		mata: cmMOR= st_matrixcolstripe("stID")
		
		*Add columns
			mata: mMOR= mMOR , J(rows(mMOR),19,.)
				
		*Add rows
			mata: cmMOR= cmMOR \ J(19,cols(cmMOR),"")
			*Label
				mata: cmMOR[2,2] = "MOR_DN"
				mata: cmMOR[3,2] = "MOR_L1S"
				mata: cmMOR[4,2] = "MOR_L1E"
				mata: cmMOR[5,2] = "MOR_L2S"
				mata: cmMOR[6,2] = "MOR_L2E"
				mata: cmMOR[7,2] = "MOR_L3S"
				mata: cmMOR[8,2] = "MOR_L3E"
				mata: cmMOR[9,2] = "MOR_L4S"
				mata: cmMOR[10,2] = "MOR_L4E"
				mata: cmMOR[11,2] = "MOR_L5S"
				mata: cmMOR[12,2] = "MOR_L5E"
				mata: cmMOR[13,2] = "MOR_L6S"
				mata: cmMOR[14,2] = "MOR_L6E"
				mata: cmMOR[15,2] = "MOR_L7S"
				mata: cmMOR[16,2] = "MOR_L7E"
				mata: cmMOR[17,2] = "MOR_L8S"
				mata: cmMOR[18,2] = "MOR_L8E"
				mata: cmMOR[19,2] = "MOR_L9S"
				mata: cmMOR[20,2] = "MOR_L9E"
}			
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

	*mOC - Outcome
{	
		mata: mOC = st_matrix("stID")
		mata: rmOC = st_matrixrowstripe("stID")
		mata: cmOC = st_matrixcolstripe("stID")
			
		*Add columns
			mata: mOC = mOC , J(rows(mOC),2,.)
				
		*Add rows
			mata: cmOC= cmOC \ J(2,cols(cmOC),"")
			*Label
				mata: cmOC[2,2] = "Time"
				mata: cmOC[3,2] = "Mort"
}		
		*mata: _matrix_list(mOC, rmOC, cmOC)	
		
	*mCR - Chemotherapy Regimen 
{		
		mata: mCR = st_matrix("stID")
		mata: rmCR = st_matrixrowstripe("stID")
		mata: cmCR = st_matrixcolstripe("stID")
		
		*Add columns
			forvalues i = 1/9{
				mata: mCR = mCR , J(rows(mCR),1,.)
			}
	
		*Add rows
			mata: cmCR = cmCR \ J(9,cols(cmCR),"")
			*Label
				mata: cmCR[2,2] = "CR_L1"
				mata: cmCR[3,2] = "CR_L2"
				mata: cmCR[4,2] = "CR_L3"
				mata: cmCR[5,2] = "CR_L4"
				mata: cmCR[6,2] = "CR_L5"
				mata: cmCR[7,2] = "CR_L6"
				mata: cmCR[8,2] = "CR_L7"
				mata: cmCR[9,2] = "CR_L8"
				mata: cmCR[10,2] = "CR_L9"
}		
		*mata: _matrix_list(mCR, rmCR, cmCR)
		
	*mCD - Chemotherapy Duration
{		
		mata: mCD = st_matrix("stID")
		mata: rmCD = st_matrixrowstripe("stID")
		mata: cmCD = st_matrixcolstripe("stID")

		*Add columns
			forvalues i = 1/9{
				mata: mCD = mCD , J(rows(mCD),1,.)
			}
	
		*Add rows
			mata: cmCD = cmCD \ J(9,cols(cmCD),"")
			*Label
				mata: cmCD[2,2] = "CD_L1"
				mata: cmCD[3,2] = "CD_L2"
				mata: cmCD[4,2] = "CD_L3"
				mata: cmCD[5,2] = "CD_L4"
				mata: cmCD[6,2] = "CD_L5"
				mata: cmCD[7,2] = "CD_L6"
				mata: cmCD[8,2] = "CD_L7"
				mata: cmCD[9,2] = "CD_L8"
				mata: cmCD[10,2] = "CD_L9"
}		
		*mata: _matrix_list(mCD, rmCD, cmCD)	
		
	*mBCR - Best Clinical Response 
{	
		mata: mBCR = st_matrix("stID")
		mata: rmBCR = st_matrixrowstripe("stID")
		mata: cmBCR = st_matrixcolstripe("stID")
		
		*Add columns
			mata: mBCR = mBCR , J(rows(mBCR),9,.)
				
		*Add rows
			mata: cmBCR = cmBCR \ J(9,cols(cmBCR),"")
			*Label
				mata: cmBCR[2,2] = "BCR_L1"
				mata: cmBCR[3,2] = "BCR_L2"
				mata: cmBCR[4,2] = "BCR_L3"
				mata: cmBCR[5,2] = "BCR_L4"
				mata: cmBCR[6,2] = "BCR_L5"
				mata: cmBCR[7,2] = "BCR_L6"
				mata: cmBCR[8,2] = "BCR_L7"
				mata: cmBCR[9,2] = "BCR_L8"
				mata: cmBCR[10,2] = "BCR_L9"
}			
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		
	*mCI - Chemotherapy Interval 
{		
		mata: mCI = st_matrix("stID")
		mata: rmCI = st_matrixrowstripe("stID")
		mata: cmCI = st_matrixcolstripe("stID")
		
		*Add columns
			forvalues i = 1/9{
				mata: mCI = mCI , J(rows(mCI),1,.)
			}
	
		*Add rows
			mata: cmCI = cmCI \ J(9,cols(cmCI),"")
			*Label
				mata: cmCI[2,2] = "CI_L1"
				mata: cmCI[3,2] = "CI_L2"
				mata: cmCI[4,2] = "CI_L3"
				mata: cmCI[5,2] = "CI_L4"
				mata: cmCI[6,2] = "CI_L5"
				mata: cmCI[7,2] = "CI_L6"
				mata: cmCI[8,2] = "CI_L7"
				mata: cmCI[9,2] = "CI_L8"
				mata: cmCI[10,2] = "CI_L9"
}		
		*mata: _matrix_list(mCI, rmCI, cmCI)
		
	*mNFT - No Further Treatment
{		
		mata: mNFT = st_matrix("stID")
		mata: rmNFT = st_matrixrowstripe("stID")
		mata: cmNFT = st_matrixcolstripe("stID")
		
		*Add columns
			forvalues i = 1/9{
				mata: mNFT = mNFT , J(rows(mNFT),1,.)
			}
	
		*Add rows
			mata: cmNFT = cmNFT \ J(9,cols(cmCI),"")
			*Label
				mata: cmNFT[2,2] = "NFT_DN"
				mata: cmNFT[3,2] = "NFT_L1"
				mata: cmNFT[4,2] = "NFT_L2"
				mata: cmNFT[5,2] = "NFT_L3"
				mata: cmNFT[6,2] = "NFT_L4"
				mata: cmNFT[7,2] = "NFT_L5"
				mata: cmNFT[8,2] = "NFT_L6"
				mata: cmNFT[9,2] = "NFT_L7"
				mata: cmNFT[10,2] = "NFT_L8"
}		
		*mata: _matrix_list(mNFT, rmNFT, cmNFT)
		
	*mSU - Survival
{
		mata: mSU = J(`=Obs',1,0)
		mata: mSU = mSU , J(`=Obs',2,.)
		mata: cmSU = J(3,2,"")
		
		*Label
			mata: cmSU[1,2] = "SU_XB"
			mata: cSU_XB = cols(mCore) + 1
			mata: cmSU[2,2] = "SU_RN"
			mata: cSU_RN = cols(mCore) + 2
			mata: cmSU[3,2] = "SU_OC"
			mata: cSU_OC = cols(mCore) + 3
}		
		*mata: _matrix_list(mSU, rmCore, cmSU)
	
	*mLO - Logit
{
		mata: mLO = J(`=Obs',1,0)
		mata: mLO = mLO , J(`=Obs',3,.)
		mata: cmLO = J(4,2,"")
		
		*Label
			mata: cmLO[1,2] = "LO_XB"
			mata: cLO_XB = cols(mCore) + 1
			mata: cmLO[2,2] = "LO_PR"
			mata: cLO_PR = cols(mCore) + 2
			mata: cmLO[3,2] = "LO_RN"
			mata: cLO_RN = cols(mCore) + 3
			mata: cmLO[4,2] = "LO_OC"
			mata: cLO_OC = cols(mCore) + 4
}			
		*mata: _matrix_list(mLO, rmCore, cmLO)
	
	*mML - Multinomial Logit
{
		mata: mML = J(`=Obs',4,0)
		mata: mML = mML , J(`=Obs',6,.)
		mata: cmML = J(10,2,"")
		
		*Label
			mata: cmML[1,2] = "ML_XB1"
			mata: cML_XB1 = cols(mCore) + 1
			mata: cmML[2,2] = "ML_XB2"
			mata: cML_XB2 = cols(mCore) + 2
			mata: cmML[3,2] = "ML_XB3"
			mata: cML_XB3 = cols(mCore) + 3
			mata: cmML[4,2] = "ML_XB4"
			mata: cML_XB4 = cols(mCore) + 4
			mata: cmML[5,2] = "ML_PR1"
			mata: cML_PR1 = cols(mCore) + 5
			mata: cmML[6,2] = "ML_PR2"
			mata: cML_PR2 = cols(mCore) + 6
			mata: cmML[7,2] = "ML_PR3"
			mata: cML_PR3 = cols(mCore) + 7
			mata: cmML[8,2] = "ML_PR4"
			mata: cML_PR4 = cols(mCore) + 8
			mata: cmML[9,2] = "ML_RN"
			mata: cML_RN = cols(mCore) + 9
			mata: cmML[10,2] = "ML_OC"
			mata: cML_OC = cols(mCore) + 10
}			
		*mata: _matrix_list(mML, rmCore, cmML)
	
	*mOL - Ordered Logit
{
		mata: mOL = J(`=Obs',1,0)
		mata: mOL = mOL , J(`=Obs',8,.)
		mata: cmOL = J(9,2,"")
		
		*Label
			mata: cmOL[1,2] = "OL_XB"
			mata: cOL_XB = cols(mCore) + 1
			mata: cmOL[2,2] = "OL_PR1"
			mata: cOL_PR1 = cols(mCore) + 2
			mata: cmOL[3,2] = "OL_PR2"
			mata: cOL_PR2 = cols(mCore) + 3
			mata: cmOL[4,2] = "OL_PR3"
			mata: cOL_PR3 = cols(mCore) + 4
			mata: cmOL[5,2] = "OL_PR4"
			mata: cOL_PR4 = cols(mCore) + 5
			mata: cmOL[6,2] = "OL_PR5"
			mata: cOL_PR5 = cols(mCore) + 6
			mata: cmOL[7,2] = "OL_PR6"
			mata: cOL_PR6 = cols(mCore) + 7
			mata: cmOL[8,2] = "OL_RN"
			mata: cOL_RN = cols(mCore) + 8
			mata: cmOL[9,2] = "OL_OC"
			mata: cOL_OC = cols(mCore) + 9
}			
		*mata: _matrix_list(mOL, rmCore, cmOL)
}
		*mata: _matrix_list(mCore, rmCore, cmCore)
		*mata: _matrix_list(mCom, rmCom, cmCom)
		*mata: _matrix_list(mAge, rmAge, cmAge)
		*mata: _matrix_list(mOS, rmOS, cmOS)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)	
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
		*mata: _matrix_list(mOC, rmOC, cmOC)
		*mata: _matrix_list(mCR, rmCR, cmCR)
		*mata: _matrix_list(mCD, rmCD, cmCD)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		*mata: _matrix_list(mCI, rmCI, cmCI)
		*mata: _matrix_list(mNFT, rmNFT, cmNFT)

**********
*Diagnosis (DN)
	scalar OMC = 2
	scalar Line = 0
	scalar LX = 0
	scalar NFT = 2
			
	di "DN - SCT"
		scalar m = "mDN_SCT"
		scalar b = "bDN_SCT"
		scalar c = "cLO_"
		quietly do "Sub/SIM/SIM SCT DN.do"			
		*mata: _matrix_list(mDN_SCT, rmDN_SCT, cmDN_SCT)

	di "DN - Chemo Interval"
		scalar m = "mDN_CI"
		scalar b = "bDN_CI"
		scalar c = "cSU_"
		quietly do "Sub/SIM/SIM CI DN.do"
		*mata: _matrix_list(mDN_CI, rmDN_CI, cmDN_CI)
		*mata: _matrix_list(mNFT, rmNFT, cmNFT)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)

	di "DN - Overall Survival"
		scalar m = "mOS_DN"
		scalar b = "bOS"
		scalar c = "cSU_"	
		quietly do "Sub/SIM/SIM OS DN.do"	
		*mata: _matrix_list(mDN_CI, rmDN_CI, cmDN_CI)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "DN - Mortality"
		forvalues i = 1/`=Obs' {
			mata {
				if 	(mTSD[`i',`=OMC'+1] > mOS[`i',`=OMC']) 	mMOR[`i',`=OMC'] = 1
				if 	(mTSD[`i',`=OMC'+1] <= mOS[`i',`=OMC']) mMOR[`i',`=OMC'] = 0
				if 	(mMOR[`i',`=OMC'] == 1)	{
					if	((mAge[`i',2] + mOS[`i',`=OMC']) > `=Limit') mOS[`i',`=OMC'] = `=Limit' - mAge[`i',2] // Set mOS to max at Limit
					mOC[`i',2] = mOS[`i',`=OMC']
					mOC[`i',3] = 1
					mTSD[`i',`=OMC'+1] = . 
					mTNE[`i',`=OMC'] = .
					mCI[`i', `=OMC'/2] = mOC[`i',2] - mTSD[`i',`=OMC']
				}
			}
		}
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
		
**********
*Chemo Line 1 Start (L1S)
	scalar OMC = 3

	di "L1S - Age"
		quietly do "Sub/SIM/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)

	di "L1S - Chemo Regimen"
		scalar m = "mL1_CR"
		scalar b = "bL1_CR"
		scalar o = "oL1_CR"
		scalar c = "cML_"
		quietly do "Sub/SIM/SIM CR L1.do"	
		*mata: _matrix_list(mL1_CR, rmL1_CR, cmL1_CR)
	
	di "L1S - Chemo Duration"
		scalar m = "mL1_CD"
		scalar c = "cSU_"
		quietly do "Sub/SIM/SIM CD L1.do"	
		*mata: _matrix_list(mL1_CD, rmL1_CD, cmL1_CD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L1S - Overall Survival"
		scalar m = "mOS_L1S"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM/SIM OS.do"
		*mata: _matrix_list(mOS_L1S, rmOS_L1S, cmOS_L1S)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L1S - Mortality"
		quietly do "Sub/SIM/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 1 End (L1E)
	scalar OMC = 4
	scalar Line = 1
	scalar LX = 1
	scalar NFT = 3

	di "L1E - Age"
		quietly do "Sub/SIM/SIM AGE.do"
		*mata: _matrix_list(mAge, rmAge, cmAge)	
	
	di "L1E - Best Clinical Response"
		scalar m = "mL1_BCR"
		scalar b = "bL1_BCR"
		scalar c = "cOL_"
		quietly do "Sub/SIM/SIM BCR L1.do"		
		*mata: _matrix_list(mL1_BCR, rmL1_BCR, cmL1_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)

	di "L1E - SCT"
		scalar m = "mL1_SCT"
		scalar b = "bL1_SCT"
		scalar c = "cLO_"
		quietly do "Sub/SIM/SIM SCT L1.do"			
		*mata: _matrix_list(mL1_SCT, rmL1_SCT, cmL1_SCT)	

	di "L1E - SCT Best Clinical Response"
		scalar m = "mSCT_BCR"
		scalar b = "bSCT_BCR"
		scalar c = "cOL_"
		quietly do "Sub/SIM/SIM BCR SCT.do"			
		*mata: _matrix_list(mSCT_BCR, rmSCT_BCR, cmSCT_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)

	di "L1E - MNT"
		scalar m = "mMNT"
		scalar b = "bMNT"
		scalar c = "cLO_"
		quietly do "Sub/SIM/SIM MNT.do"		
		*mata: _matrix_list(mMNT, rmMNT, cmMNT)
		*mata: _matrix_list(mCore, rmCore, cmCore)
		
	di "L1E - Chemo Interval"
		scalar m = "mL1_CI"
		scalar c = "cSU_"
		quietly do "Sub/SIM/SIM CI L1.do"			
		*mata: _matrix_list(mL1_CI, rmL1_CI, cmL1_CI)
		*mata: _matrix_list(mNFT, rmNFT, cmNFT)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)		
	
	di "L1E - Overall Survival"
		scalar m = "mOS_L1E"
		scalar b = "bOS"	
		scalar c = "cLO_"
		quietly do "Sub/SIM/SIM OS.do"	
		*mata: _matrix_list(mOS_L1E, rmOS_L1E, cmOS_L1E)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L1E - Mortality"
		quietly do "Sub/SIM/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 2 Start (L2S)
	scalar OMC = 5
		
	di "L2S - Age"
		quietly do "Sub/SIM/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)		
		
	di "L2S - Chemo Regimen"
		scalar m = "mL2_CR"
		scalar b = "bL2_CR"
		scalar o = "oL2_CR"
		scalar c = "cML_"
		quietly do "Sub/SIM/SIM CR L2.do"
		*mata: _matrix_list(mL2_CR, rmL2_CR, cmL2_CR)

	di "L2S - Chemo Duration"
		scalar m = "mL2_CD"
		scalar b = "bL2_CD"
		scalar c = "cSU_"
		quietly do "Sub/SIM/SIM CD.do"
		*mata: _matrix_list(mL2_CD, rmL2_CD, cmL2_CD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L2S - Overall Survival"
		scalar m = "mOS_L2S"
		scalar b = "bOS"	
		scalar c = "cLO_"
		quietly do "Sub/SIM/SIM OS.do"	
		*mata: _matrix_list(mOS_L2S, rmOS_L2S, cmOS_L2S)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L2S - Mortality"
		quietly do "Sub/SIM/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
	
**********
*Chemo Line 2 End (L2E)
	scalar OMC = 6
	scalar Line = 2
	scalar LX = 2
	scalar NFT = 4

	di "L2E - Age"
		quietly do "Sub/SIM/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)		
	
	di "L2E - Best Clinical Response"
		scalar m = "mL2_BCR"
		scalar b = "bL2_BCR"
		scalar c = "cOL_"
		quietly do "Sub/SIM/SIM BCR L2.do"
		*mata: _matrix_list(mL2_BCR, rmL2_BCR, cmL2_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		
	di "L2E - Chemo Interval"
		scalar m = "mL2_CI"
		scalar b = "bL2_CI"
		scalar c = "cSU_"
		quietly do "Sub/SIM/SIM CI.do"
		*mata: _matrix_list(mL2_CI, rmL2_CI, cmL2_CI)
		*mata: _matrix_list(mNFT, rmNFT, cmNFT)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)

	di "L2E - Overall Survival"
		scalar m = "mOS_L2E"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM/SIM OS.do"	
		*mata: _matrix_list(mOS_L2E, rmOS_L2E, cmOS_L2E)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L2E - Mortality"
		quietly do "Sub/SIM/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
	
**********
*Chemo Line 3 Start (L3S)
	scalar OMC = 7
		
	di "L3S - Age"
		quietly do "Sub/SIM/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
		
	di "L3S - Chemo Regimen"
		forvalues i = 1/`=Obs'{
			mata {
				if 	(mMOR[`i',`=OMC'-1] == 0) mCore[`i',cCR] = 0
				if 	(mMOR[`i',`=OMC'-1] != 0) mCore[`i',cCR] = .
				if 	(mMOR[`i',`=OMC'-1] == 0) mCR[`i',`=Line'+2] = 0
			}
		}
				
	di "L3S - Chemo Duration"
		scalar m = "mL3_CD"
		scalar b = "bL3_CD"
		scalar c = "cSU_"
		quietly do "Sub/SIM/SIM CD L3.do"	
		*mata: _matrix_list(mL3_CD, rmL3_CD, cmL3_CD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L3S - Overall Survival"
		scalar m = "mOS_L3S"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM/SIM OS.do"
		*mata: _matrix_list(mOS_L3S, rmOS_L3S, cmOS_L3S)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L3S - Mortality"
		quietly do "Sub/SIM/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)	

**********
*Chemo Line 3 End (L3E)
	scalar OMC = 8
	scalar Line = 3
	scalar LX = 3
	scalar NFT = 5
		
	di "L3E - Age"
		quietly do "Sub/SIM/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
	
	di "L3E - Best Clinical Response"
		scalar m = "mL3_BCR"
		scalar b = "bL3_BCR"
		scalar c = "cOL_"
		quietly do "Sub/SIM/SIM BCR L3.do"		
		*mata: _matrix_list(mL3_BCR, rmL3_BCR, cmL3_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		
	di "L3E - Chemo Interval"
		scalar m = "mL3_CI"
		scalar b = "bL3_CI"
		scalar c = "cSU_"
		quietly do "Sub/SIM/SIM CI LX.do"		
		*mata: _matrix_list(mL3_CI, rmL3_CI, cmL3_CI)
		*mata: _matrix_list(mNFT, rmNFT, cmNFT)	
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)

	di "L3E - Overall Survival"
		scalar m = "mOS_L3E"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM/SIM OS.do"
		*mata: _matrix_list(mOS_L3E, rmOS_L3E, cmOS_L3E)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L3E - Mortality"
		quietly do "Sub/SIM/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 4 Start (L4S)
	scalar OMC = 9
		
	di "L4S - Age"
		quietly do "Sub/SIM/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)		

	di "L4S - Chemo Regimen"
		forvalues i = 1/`=Obs'{
			mata {
				if 	(mMOR[`i',`=OMC'-1] == 0) mCore[`i',cCR] = 0
				if 	(mMOR[`i',`=OMC'-1] != 0) mCore[`i',cCR] = .
				if 	(mMOR[`i',`=OMC'-1] == 0) mCR[`i',`=Line'+2] = 0
			}
		}
			
	di "L4S - Chemo Duration"
		scalar m = "mL4_CD"
		scalar b = "bL4_CD"
		scalar c = "cSU_"
		quietly do "Sub/SIM/SIM CD LX.do"		
		*mata: _matrix_list(mL4_CD, rmL4_CD, cmL4_CD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L4S - Overall Survival"
		scalar m = "mOS_L4S"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM/SIM OS.do"
		*mata: _matrix_list(mOS_L4S, rmOS_L4S, cmOS_L4S)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L4S - Mortality"
		quietly do "Sub/SIM/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 4 End (L4E)
	scalar OMC = 10
	scalar Line = 4
	scalar LX = 4
	scalar NFT = 6
		
	di "L4E - Age"
		quietly do "Sub/SIM/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
	
	di "L4E - Best Clinical Response"
		scalar m = "mL4_BCR"
		scalar b = "bL4_BCR"
		scalar c = "cOL_"
		quietly do "Sub/SIM/SIM BCR LX.do"	
		*mata: _matrix_list(mL4_BCR, rmL4_BCR, cmL4_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		
	di "L4E - Chemo Interval"
		scalar m = "mL4_CI"
		scalar b = "bL4_CI"
		scalar c = "cSU_"
		quietly do "Sub/SIM/SIM CI LX.do"	
		*mata: _matrix_list(mL4_CI, rmL4_CI, cmL4_CI)
		*mata: _matrix_list(mNFT, rmNFT, cmNFT)	
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)

	di "L4E - Overall Survival"
		scalar m = "mOS_L4E"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM/SIM OS.do"	
		*mata: _matrix_list(mOS_L4E, rmOS_L4E, cmOS_L4E)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L4E - Mortality"
		quietly do "Sub/SIM/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
			
**********
*Chemo Line 5 Start (L5S)
	scalar OMC = 11
		
	di "L5S - Age" 
		quietly do "Sub/SIM/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
		
	di "L5S - Chemo Regimen"
		forvalues i = 1/`=Obs'{
			mata {
				if 	(mMOR[`i',`=OMC'-1] == 0) mCore[`i',cCR] = 0
				if 	(mMOR[`i',`=OMC'-1] != 0) mCore[`i',cCR] = .
				if 	(mMOR[`i',`=OMC'-1] == 0) mCR[`i',`=Line'+2] = 0
			}
		}

	di "L5S - Chemo Duration"
		scalar m = "mL5_CD"
		scalar b = "bL5_CD"
		scalar c = "cSU_"
		quietly do "Sub/SIM/SIM CD LX.do"		
		*mata: _matrix_list(mL5_CD, rmL5_CD, cmL5_CD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L5S - Overall Survival"
		scalar m = "mOS_L5S"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM/SIM OS.do"
		*mata: _matrix_list(mOS_L5S, rmOS_L5S, cmOS_L5S)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L5S - Mortality"
		quietly do "Sub/SIM/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 5 End (L5E)
	scalar OMC = 12
	scalar Line = 5
	scalar LX = 5
	scalar NFT = 7
		
	di "L5E - Age"
		quietly do "Sub/SIM/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	

	di "L5E - Best Clinical Response"
		scalar m = "mL5_BCR"
		scalar b = "bL5_BCR"
		scalar c = "cOL_"
		quietly do "Sub/SIM/SIM BCR LX.do"
		*mata: _matrix_list(mL5_BCR, rmL5_BCR, cmL5_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
	
	di "L5E - Chemo Interval"
		scalar m = "mL5_CI"
		scalar b = "bL5_CI"
		scalar c = "cSU_"
		quietly do "Sub/SIM/SIM CI LX.do"				
		*mata: _matrix_list(mL5_CI, rmL5_CI, cmL5_CI)
		*mata: _matrix_list(mNFT, rmNFT, cmNFT)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)

	di "L5E - Overall Survival"
		scalar m = "mOS_L5E"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM/SIM OS.do"
		*mata: _matrix_list(mOS_L5E, rmOS_L5E, cmOS_L5E)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L5E - Mortality"
		quietly do "Sub/SIM/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
		
**********
*Chemo Line 6 Start (L6S) 
	scalar OMC = 13
		
	di "L6S - Age" 
		quietly do "Sub/SIM/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
		
	di "L6S - Chemo Regimen"
		forvalues i = 1/`=Obs'{
			mata {
				if 	(mMOR[`i',`=OMC'-1] == 0) mCore[`i',cCR] = 0
				if 	(mMOR[`i',`=OMC'-1] != 0) mCore[`i',cCR] = .
				if 	(mMOR[`i',`=OMC'-1] == 0) mCR[`i',`=Line'+2] = 0
			}
		}

	di "L6S - Chemo Duration"
		scalar m = "mL6_CD"
		scalar b = "bLX_CD" 		
		scalar c = "cSU_"
		quietly do "Sub/SIM/SIM CD LXX.do"	
		*mata: _matrix_list(mL6_CD, rmL6_CD, cmL6_CD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L6S - Overall Survival"
		scalar m = "mOS_L6S"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM/SIM OS.do"
		*mata: _matrix_list(mOS_L6S, rmOS_L6S, cmOS_L6S)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L6S - Mortality"
		quietly do "Sub/SIM/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 6 End (L6E)
	scalar OMC = 14
	scalar Line = 6
	scalar LX = 6
	scalar NFT = 8
		
	di "L6E - Age"
		quietly do "Sub/SIM/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
	
	di "L6E - Best Clinical Response"
		scalar m = "mL6_BCR"
		scalar b = "bLX_BCR"
		scalar c = "cOL_"
		quietly do "Sub/SIM/SIM BCR LXX.do"	 	
		*mata: _matrix_list(mL6_BCR, rmL6_BCR, cmL6_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		
	di "L6E - Chemo Interval"
		scalar m = "mL6_CI"
		scalar b = "bLX_CI"
		scalar c = "cSU_"
		quietly do "Sub/SIM/SIM CI LXX.do"				
		*mata: _matrix_list(mL6_CI, rmL6_CI, cmL6_CI)
		*mata: _matrix_list(mNFT, rmNFT, cmNFT)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)	

	di "L6E - Overall Survival"
		scalar m = "mOS_L6E"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM/SIM OS.do"
		*mata: _matrix_list(mOS_L6E, rmOS_L6E, cmOS_L6E)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L6E - Mortality" 
		quietly do "Sub/SIM/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
					
**********
*Chemo Line 7 Start (L7S) 
	scalar OMC = 15
		
	di "L7S - Age" 
		quietly do "Sub/SIM/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
		
	di "L7S - Chemo Regimen"
		forvalues i = 1/`=Obs'{
			mata {
				if 	(mMOR[`i',`=OMC'-1] == 0) mCore[`i',cCR] = 0
				if 	(mMOR[`i',`=OMC'-1] != 0) mCore[`i',cCR] = .
				if 	(mMOR[`i',`=OMC'-1] == 0) mCR[`i',`=Line'+2] = 0
			}
		}

	di "L7S - Chemo Duration"
		scalar m = "mL7_CD"
		scalar b = "bLX_CD" 		
		scalar c = "cSU_"
		quietly do "Sub/SIM/SIM CD LXX.do"	
		*mata: _matrix_list(mL7_CD, rmL7_CD, cmL7_CD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L7S - Overall Survival"
		scalar m = "mOS_L7S"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM/SIM OS.do"
		*mata: _matrix_list(mOS_L7S, rmOS_L7S, cmOS_L7S)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L7S - Mortality"
		quietly do "Sub/SIM/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 7 End (L7E)
	scalar OMC = 16
	scalar LX = 7
	scalar NFT = 9
		
	di "L7E - Age"
		quietly do "Sub/SIM/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
	
	di "L7E - Best Clinical Response" 
		scalar m = "mL7_BCR"
		scalar b = "bLX_BCR"
		scalar c = "cOL_"
		quietly do "Sub/SIM/SIM BCR LXX.do"	 	
		*mata: _matrix_list(mL7_BCR, rmL7_BCR, cmL7_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		
	di "L7E - Chemo Interval"
		scalar m = "mL7_CI"
		scalar b = "bLX_CI"
		scalar c = "cSU_"
		quietly do "Sub/SIM/SIM CI LXX.do"				
		*mata: _matrix_list(mL7_CI, rmL7_CI, cmL7_CI)
		*mata: _matrix_list(mNFT, rmNFT, cmNFT)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)		

	di "L7E - Overall Survival"
		scalar m = "mOS_L7E"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM/SIM OS.do"
		*mata: _matrix_list(mOS_L7E, rmOS_L7E, cmOS_L7E)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L7E - Mortality" 
		quietly do "Sub/SIM/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
		
**********
*Chemo Line 8 Start (L8S) 
	scalar OMC = 17
		
	di "L8S - Age" 
		quietly do "Sub/SIM/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
		
	di "L8S - Chemo Regimen"
		forvalues i = 1/`=Obs'{
			mata {
				if 	(mMOR[`i',`=OMC'-1] == 0) mCore[`i',cCR] = 0
				if 	(mMOR[`i',`=OMC'-1] != 0) mCore[`i',cCR] = .
				if 	(mMOR[`i',`=OMC'-1] == 0) mCR[`i',`=Line'+2] = 0
			}
		}

	di "L8S - Chemo Duration"
		scalar m = "mL8_CD"
		scalar b = "bLX_CD" 		
		scalar c = "cSU_"
		quietly do "Sub/SIM/SIM CD LXX.do"	
		*mata: _matrix_list(mL8_CD, rmL8_CD, cmL8_CD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L8S - Overall Survival"
		scalar m = "mOS_L8S"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM/SIM OS.do"
		*mata: _matrix_list(mOS_L8S, rmOS_L8S, cmOS_L8S)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L8S - Mortality"
		quietly do "Sub/SIM/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 8 End (L8E)
	scalar OMC = 18
	scalar LX = 8
	scalar NFT = 10
		
	di "L8E - Age"
		quietly do "Sub/SIM/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
	
	di "L8E - Best Clinical Response"
		scalar m = "mL8_BCR"
		scalar b = "bLX_BCR"
		scalar c = "cOL_"
		quietly do "Sub/SIM/SIM BCR LXX.do"	
		*mata: _matrix_list(mL8_BCR, rmL8_BCR, cmL8_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		
	di "L8E - Chemo Interval"
		scalar m = "mL8_CI"
		scalar b = "bLX_CI"
		scalar c = "cSU_"
		quietly do "Sub/SIM/SIM CI LXX.do"				
		*mata: _matrix_list(mL8_CI, rmL8_CI, cmL8_CI)
		*mata: _matrix_list(mNFT, rmNFT, cmNFT)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)	

	di "L8E - Overall Survival"
		scalar m = "mOS_L8E"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM/SIM OS.do"
		*mata: _matrix_list(mOS_L8E, rmOS_L8E, cmOS_L8E)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L8E - Mortality" 
		quietly do "Sub/SIM/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)	

**********
*Chemo Line 9 Start (L9S) 
	scalar OMC = 19
		
	di "L9S - Age" 
		quietly do "Sub/SIM/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
		
	di "L9S - Chemo Regimen"
		forvalues i = 1/`=Obs'{
			mata {
				if 	(mMOR[`i',`=OMC'-1] == 0) mCore[`i',cCR] = 0
				if 	(mMOR[`i',`=OMC'-1] != 0) mCore[`i',cCR] = .
				if 	(mMOR[`i',`=OMC'-1] == 0) mCR[`i',`=Line'+2] = 0
			}
		}

	di "L9S - Chemo Duration"
		scalar m = "mL9_CD"
		scalar b = "bLX_CD" 		
		scalar c = "cSU_"
		quietly do "Sub/SIM/SIM CD LXX.do"	
		*mata: _matrix_list(mL9_CD, rmL9_CD, cmL9_CD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L9S - Overall Survival"
		scalar m = "mOS_L9S"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM/SIM OS.do"
		*mata: _matrix_list(mOS_L9S, rmOS_L9S, cmOS_L9S)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L9S - Mortality"
		quietly do "Sub/SIM/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 9 End (L9E)
	scalar OMC = 20
	scalar LX = 9
	scalar NFT = 11
		
	di "L9E - Age"
		quietly do "Sub/SIM/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
	
	di "L9E - Best Clinical Response"
		scalar m = "mL9_BCR"
		scalar b = "bLX_BCR"
		scalar c = "cOL_"
		quietly do "Sub/SIM/SIM BCR LXX.do"
		*mata: _matrix_list(mL9_BCR, rmL9_BCR, cmL9_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		
	di "L9E - Overall Survival"
		scalar m = "mOS_L9E"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM/SIM OS.do"
		*mata: _matrix_list(mOS_L9E, rmOS_L9E, cmOS_L9E)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L9E - Mortality" 
		*Everyone still alive dies at predicted OS or Limit
		forvalues i = 1/`=Obs'{
			mata {
				if	(mMOR[`i',`=OMC'-1] == 0) {
					*If mOS would take them over Limit
						if	((mAge[`i',2] + mOS[`i',`=OMC']) > `=Limit') mOS[`i',`=OMC'] = `=Limit' - mAge[`i',2] // Set mOS to max if Age > Limit
					mMOR[`i',`=OMC'] = 1
					mOC[`i',2] = mOS[`i',`=OMC']
					mOC[`i',3] = 1
				}		
			}	
		}
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
		*mata: _matrix_list(mOC, rmOC, cmOC)		
		
**********
*Summary Matrix

	*Create mSum in Mata 
		mata: mSum = mCore , mAge , mOS , mTNE , mTSD , mMOR , mOC , mCR , mCD , mBCR , mSCT_BCR[.,19] , mCI , mNFT
	
	*Convert mSum to stSum
		mata: st_matrix("stSum", mSum)
		drop _all
		
	*Convert stSum to variables
		svmat double stSum
	
	*Name variables
		local varnames ID Age Male ECOGc ISS SCT MNT CR CD BCR ///
			Age_ID Age_DN Age_L1S Age_L1E Age_L2S Age_L2E Age_L3S Age_L3E Age_L4S Age_L4E Age_L5S Age_L5E Age_L6S Age_L6E Age_L7S Age_L7E Age_L8S Age_L8E Age_L9S Age_L9E ///
			OS_ID OS_DN OS_L1S OS_L1E OS_L2S OS_L2E OS_L3S OS_L3E OS_L4S OS_L4E OS_L5S OS_L5E OS_L6S OS_L6E OS_L7S OS_L7E OS_L8S OS_L8E OS_L9S OS_L9E ///
			TNE_ID TNE_DN TNE_L1S TNE_L1E TNE_L2S TNE_L2E TNE_L3S TNE_L3E TNE_L4S TNE_L4E TNE_L5S TNE_L5E TNE_L6S TNE_L6E TNE_L7S TNE_L7E TNE_L8S TNE_L8E TNE_L9S TNE_L9E ///
			TSD_ID TSD_DN TSD_L1S TSD_L1E TSD_L2S TSD_L2E TSD_L3S TSD_L3E TSD_L4S TSD_L4E TSD_L5S TSD_L5E TSD_L6S TSD_L6E TSD_L7S TSD_L7E TSD_L8S TSD_L8E TSD_L9S TSD_L9E ///
			MOR_ID MOR_DN MOR_L1S MOR_L1E MOR_L2S MOR_L2E MOR_L3S MOR_L3E MOR_L4S MOR_L4E MOR_L5S MOR_L5E MOR_L6S MOR_L6E MOR_L7S MOR_L7E MOR_L8S MOR_L8E MOR_L9S MOR_L9E ///
			OC_ID OC_TIME OC_MORT ///
			CR_ID CR_L1 CR_L2 CR_L3 CR_L4 CR_L5 CR_L6 CR_L7 CR_L8 CR_L9 ///
			CD_ID CD_L1 CD_L2 CD_L3 CD_L4 CD_L5 CD_L6 CD_L7 CD_L8 CD_L9 ///
			BCR_ID BCR_L1 BCR_L2 BCR_L3 BCR_L4 BCR_L5 BCR_L6 BCR_L7 BCR_L8 BCR_L9 BCR_SCT ///
			CI_ID CI_L1 CI_L2 CI_L3 CI_L4 CI_L5 CI_L6 CI_L7 CI_L8 CI_L9 /// 
			NFT_ID NFT_DN NFT_L1 NFT_L2 NFT_L3 NFT_L4 NFT_L5 NFT_L6 NFT_L7 NFT_L8 ///
		
		local varlength : word count `varnames'
		
		forvalues i = 1/`varlength'{
			local currentvar : word `i' of `varnames'
			rename stSum`i' `currentvar'
		}		
	
	*Drop unnecessary variables
		drop Age_ID OS_ID TNE_ID TSD_ID MOR_ID OC_ID CR_ID CD_ID BCR_ID CI_ID NFT_ID CR CD BCR TSD_DN TNE_L9E CI_L9
		order ID Age Male ECOGc ISS SCT MNT
		
	*Change ID to reflect bootstrap sample
		replace ID = ID + `bs'0000000

	*Save EpiMAP Simulated Dataset
		save "Data/Simulated/EpiMAP V Wide Simulated `bs'.dta", replace
}		
/*
	
	*Code to clear matrix xb and outcome
		forvalues i = 1/`=Obs'{
			mata: `=m'[`i',`=c'XB] = 0
			mata: `=m'[`i',`=c'OC] = .
		}
		mata: _matrix_list(`=b', r`=b', c`=b')
