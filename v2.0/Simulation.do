**********
	*EpiMAP - Simulation
**********
	
	clear
	clear mata
	clear matrix
	set more off
	
**********
*Capture arguments 
	global Analysis `1'
	global Int `2'
	global Line `3'
	global Coeffs `4'
	global Data `5'
	global MinID `6'
	global MaxID `7'
	global Boot `8' 
	global MinBS `9' 
	global MaxBS `10'

**********
*Simulation function
	capture program drop simulation
	program define simulation
	
	*Additional variables needed for mCore2 to transfer to matrix multiplication for XB
		gen ECOGcc0 = (ECOGcc == 0)
		gen ECOGcc1 = (ECOGcc == 1)
		gen ECOGcc2 = (ECOGcc == 2)
		gen RISS1 = (RISS == 1)
		gen RISS2 = (RISS == 2)
		gen RISS3 = (RISS == 3)
		gen Cons = 1
	
	*Set up matrices
{	

	*mState - State matrix
{
		mata {
			mState = st_data(.,"ID State DateDN")
			rmState = J(rows(mState),2,"")
			rmState[.,2] = strofreal(mState[.,1])
			cmState = J(3,2,"")
			cmState[1,1] = "ID" 
			cmState[2,1] = "State"
			cmState[3,1] = "DateDN"
		}	
}	
		*_matrix_list(mState, rmState, cmState)
	
	*mCore - Core patient matrix
{
		mata {
			mCore = st_data(., "ID Age Male ECOGcc RISS SCT MNT CR CD BCR")
			rmCore = J(rows(mCore), 2, "")
			rmCore[.,2] = strofreal(mCore[.,1])
			cmCore = J(10, 2, "")
			cmCore[1,1] = "ID" 
			cmCore[2,1] = "Age"
			cmCore[3,1] = "Male" 
			cmCore[4,1] = "ECOGcc" 
			cmCore[5,1] = "RISS" 
			cmCore[6,1] = "SCT" 
			cmCore[7,1] = "MNT" 
			cmCore[8,1] = "CR" 
			cmCore[9,1] = "CD" 
			cmCore[10,1] = "BCR" 
		
		// Add mata column references
			cAge = 2
			cMale = 3
			cECOG = 4
			cRISS = 5
			cSCT = 6
			cMNT = 7
			cCR = 8
			cCD = 9
			cBCR = 10
		}
}	
		*mata: _matrix_list(mCore, rmCore, cmCore)
		
	*mCore2 - Core patient matrix for matrix multiplication
{
		mata {
			mCore2 = st_data(., "Age Age2 Male ECOGcc0 ECOGcc1 ECOGcc2 RISS1 RISS2 RISS3")
			rmCore2 = J(rows(mCore2),2,"")
			rmCore2[.,2] = strofreal(1::rows(mCore2)) 
			cmCore2 = J(9,2,"")
			cmCore2[1,1] = "Age" 
			cmCore2[2,1] = "Age2"
			cmCore2[3,1] = "Male" 
			cmCore2[4,1] = "ECOGcc0" 
			cmCore2[5,1] = "ECOGcc1" 
			cmCore2[6,1] = "ECOGcc2" 
			cmCore2[7,1] = "RISS1" 
			cmCore2[8,1] = "RISS2" 
			cmCore2[9,1] = "RISS3" 
		}
}
		*mata: _matrix_list(mCore2, rmCore2, cmCore2)
	
	*mCom - Comorbidity matrix
{
		mata {
			mCom = st_data(.,"Age70 Age75 CMScore")
			rmCom = J(rows(mCom),2,"")
			rmCom[.,2] = strofreal(1::rows(mCom))
			cmCom = J(3,2,"")
			cmCom[1,1] = "Age70" 
			cmCom[2,1] = "Age75"
			cmCom[3,1] = "CMScore" 
		}
}
		*mata: _matrix_list(mCom, rmCom, cmCom)
		
	*mSCT
{
		mata {
			mSCT = st_data(., "SCT SCT") // Currently using same variable twice
			rmSCT = J(rows(mSCT),2,"")
			rmSCT[.,2] = strofreal(1::rows(mSCT))
			cmSCT = J(2,2,"")
			cmSCT[1,1] = "SCT_DN"
			cmSCT[2,1] = "SCT_L1"
		}
}		
		*mata: _matrix_list(mSCT, rmSCT, cmSCT)
		
	*mMNT
{
		mata {
			mMNT = st_data(.,"MNT")
			rmMNT = J(rows(mMNT),2,"")
			rmMNT[.,2] = strofreal(1::rows(mMNT))
			cmMNT = J(1,2,"")
			cmMNT[1,1] = "MNT"
		}
}		
		*mata: _matrix_list(mMNT, rmMNT, cmMNT)	
		
	*mCons
{
		mata {
			mCons = st_data(.,"Cons")
			rmCons = J(rows(mCons), 2,"")
			rmCons[.,2] = strofreal(1::rows(mCons))
			cmCons = J(1,2,"")
			cmCons[1,1] = "Cons"
		}
}		
		*mata: _matrix_list(mCons, rmCons, cmCons)	
	
	*mAge - Age (at event)
{	
		mata {
			mAge = st_data(., "ID Age_DN Age_L1S Age_L1E Age_L2S Age_L2E Age_L3S Age_L3E Age_L4S Age_L4E Age_L5S Age_L5E Age_L6S Age_L6E Age_L7S Age_L7E Age_L8S Age_L8E Age_L9S Age_L9E")
			rmAge = J(rows(mAge), 2, "")
			rmAge[., 2] = strofreal(1::rows(mAge))
			cmAge = J(20, 2, "")
			cmAge[1,1] = "ID"
			cmAge[2,1] = "Age_DN"
			cmAge[3,1] = "Age_L1S"
			cmAge[4,1] = "Age_L1E"
			cmAge[5,1] = "Age_L2S"
			cmAge[6,1] = "Age_L2E"
			cmAge[7,1] = "Age_L3S"
			cmAge[8,1] = "Age_L3E"
			cmAge[9,1] = "Age_L4S"
			cmAge[10,1] = "Age_L4E"
			cmAge[11,1] = "Age_L5S"
			cmAge[12,1] = "Age_L5E"
			cmAge[13,1] = "Age_L6S"
			cmAge[14,1] = "Age_L6E"
			cmAge[15,1] = "Age_L7S"
			cmAge[16,1] = "Age_L7E"
			cmAge[17,1] = "Age_L8S"
			cmAge[18,1] = "Age_L8E"
			cmAge[19,1] = "Age_L9S"
			cmAge[20,1] = "Age_L9E"
		}
}		
		*mata: _matrix_list(mAge, rmAge, cmAge)		

	*mOS - Overall Survival (from DN)
{	
		mata {
			mOS = st_data(., "ID")
			mOS = mOS , J(rows(mOS),19,.)
			rmOS = J(rows(mOS), 2, "")
			rmOS[., 2] = strofreal(1::rows(mOS))
			cmOS = J(20, 2, "")
			cmOS[1,2] = "ID"
			cmOS[2,2] = "OS_DN"
			cmOS[3,2] = "OS_L1S"
			cmOS[4,2] = "OS_L1E"
			cmOS[5,2] = "OS_L2S"
			cmOS[6,2] = "OS_L2E"
			cmOS[7,2] = "OS_L3S"
			cmOS[8,2] = "OS_L3E"
			cmOS[9,2] = "OS_L4S"
			cmOS[10,2] = "OS_L4E"
			cmOS[11,2] = "OS_L5S"
			cmOS[12,2] = "OS_L5E"
			cmOS[13,2] = "OS_L6S"
			cmOS[14,2] = "OS_L6E"
			cmOS[15,2] = "OS_L7S"
			cmOS[16,2] = "OS_L7E"
			cmOS[17,2] = "OS_L8S"
			cmOS[18,2] = "OS_L8E"
			cmOS[19,2] = "OS_L9S"
			cmOS[20,2] = "OS_L9E"
		}
}			
		*mata: _matrix_list(mOS, rmOS, cmOS)		

	*mTNE - Time to Next Event
{	
		mata {
			mTNE = st_data(., "ID TNE_DN TNE_L1S TNE_L1E TNE_L2S TNE_L2E TNE_L3S TNE_L3E TNE_L4S TNE_L4E TNE_L5S TNE_L5E TNE_L6S TNE_L6E TNE_L7S TNE_L7E TNE_L8S TNE_L8E TNE_L9S TNE_L9E")
			rmTNE = J(rows(mTNE), 2, "")
			rmTNE[., 2] = strofreal(1::rows(mTNE))
			cmTNE = J(20, 2, "")
			cmTNE[1,1] = "ID"
			cmTNE[2,1] = "TNE_DN"
			cmTNE[3,1] = "TNE_L1S"
			cmTNE[4,1] = "TNE_L1E"
			cmTNE[5,1] = "TNE_L2S"
			cmTNE[6,1] = "TNE_L2E"
			cmTNE[7,1] = "TNE_L3S"
			cmTNE[8,1] = "TNE_L3E"
			cmTNE[9,1] = "TNE_L4S"
			cmTNE[10,1] = "TNE_L4E"
			cmTNE[11,1] = "TNE_L5S"
			cmTNE[12,1] = "TNE_L5E"
			cmTNE[13,1] = "TNE_L6S"
			cmTNE[14,1] = "TNE_L6E"
			cmTNE[15,1] = "TNE_L7S"
			cmTNE[16,1] = "TNE_L7E"
			cmTNE[17,1] = "TNE_L8S"
			cmTNE[18,1] = "TNE_L8E"
			cmTNE[19,1] = "TNE_L9S"
			cmTNE[20,1] = "TNE_L9E"
		}
}			
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)	
		
	*mTSD - Time Since Diagnosis
{	
		mata {
			mTSD = st_data(., "ID TSD_DN TSD_L1S TSD_L1E TSD_L2S TSD_L2E TSD_L3S TSD_L3E TSD_L4S TSD_L4E TSD_L5S TSD_L5E TSD_L6S TSD_L6E TSD_L7S TSD_L7E TSD_L8S TSD_L8E TSD_L9S TSD_L9E")
			rmTSD = J(rows(mTSD), 2, "")
			rmTSD[., 2] = strofreal(1::rows(mTSD))
			cmTSD = J(20, 2, "")
			cmTSD[1,1] = "ID"
			cmTSD[2,1] = "TSD_DN"
			cmTSD[3,1] = "TSD_L1S"
			cmTSD[4,1] = "TSD_L1E"
			cmTSD[5,1] = "TSD_L2S"
			cmTSD[6,1] = "TSD_L2E"
			cmTSD[7,1] = "TSD_L3S"
			cmTSD[8,1] = "TSD_L3E"
			cmTSD[9,1] = "TSD_L4S"
			cmTSD[10,1] = "TSD_L4E"
			cmTSD[11,1] = "TSD_L5S"
			cmTSD[12,1] = "TSD_L5E"
			cmTSD[13,1] = "TSD_L6S"
			cmTSD[14,1] = "TSD_L6E"
			cmTSD[15,1] = "TSD_L7S"
			cmTSD[16,1] = "TSD_L7E"
			cmTSD[17,1] = "TSD_L8S"
			cmTSD[18,1] = "TSD_L8E"
			cmTSD[19,1] = "TSD_L9S"
			cmTSD[20,1] = "TSD_L9E"
		}
}			
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)		
		
	*mMOR - Mortality (=1 if patient dies before next event) 
{	
		mata {
			mMOR = st_data(., "ID MOR_DN MOR_L1S MOR_L1E MOR_L2S MOR_L2E MOR_L3S MOR_L3E MOR_L4S MOR_L4E MOR_L5S MOR_L5E MOR_L6S MOR_L6E MOR_L7S MOR_L7E MOR_L8S MOR_L8E MOR_L9S MOR_L9E")
			rmMOR = J(rows(mMOR), 2, "")
			rmMOR[., 2] = strofreal(1::rows(mMOR))
			cmMOR = J(20, 2, "")
			cmMOR[1,1] = "ID"
			cmMOR[2,1] = "MOR_DN"
			cmMOR[3,1] = "MOR_L1S"
			cmMOR[4,1] = "MOR_L1E"
			cmMOR[5,1] = "MOR_L2S"
			cmMOR[6,1] = "MOR_L2E"
			cmMOR[7,1] = "MOR_L3S"
			cmMOR[8,1] = "MOR_L3E"
			cmMOR[9,1] = "MOR_L4S"
			cmMOR[10,1] = "MOR_L4E"
			cmMOR[11,1] = "MOR_L5S"
			cmMOR[12,1] = "MOR_L5E"
			cmMOR[13,1] = "MOR_L6S"
			cmMOR[14,1] = "MOR_L6E"
			cmMOR[15,1] = "MOR_L7S"
			cmMOR[16,1] = "MOR_L7E"
			cmMOR[17,1] = "MOR_L8S"
			cmMOR[18,1] = "MOR_L8E"
			cmMOR[19,1] = "MOR_L9S"
			cmMOR[20,1] = "MOR_L9E"
		}
}			
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

	*mOC - Outcome
{	
		mata {
			mOC = st_data(., "ID")
			mOC = mOC , J(rows(mOC),2,.)				
			rmOC = J(rows(mOC), 2, "")
			rmOC[., 2] = strofreal(1::rows(mOC))
			cmOC = J(3, 2, "")
			cmOC[1,2] = "ID"
			cmOC[2,2] = "Time"
			cmOC[3,2] = "Mort"
		}
}		
		*mata: _matrix_list(mOC, rmOC, cmOC)	
		
	*mCR - Chemotherapy Regimen 
{		
		mata {
			mCR = st_data(., "ID CR_L1 CR_L2 CR_L3 CR_L4 CR_L5 CR_L6 CR_L7 CR_L8 CR_L9")		
			rmCR = J(rows(mCR), 2, "")
			rmCR[., 2] = strofreal(1::rows(mCR))
			cmCR = J(10, 2, "")
			cmCR[1,2] = "ID"
			cmCR[2,2] = "CR_L1"
			cmCR[3,2] = "CR_L2"
			cmCR[4,2] = "CR_L3"
			cmCR[5,2] = "CR_L4"
			cmCR[6,2] = "CR_L5"
			cmCR[7,2] = "CR_L6"
			cmCR[8,2] = "CR_L7"
			cmCR[9,2] = "CR_L8"
			cmCR[10,2] = "CR_L9"
		}
}		
		*mata: _matrix_list(mCR, rmCR, cmCR)
		
	*mCD - Chemotherapy Duration
{		
		mata {
			mCD = st_data(., "ID CD_L1 CD_L2 CD_L3 CD_L4 CD_L5 CD_L6 CD_L7 CD_L8 CD_L9")		
			rmCD = J(rows(mCD), 2, "")
			rmCD[., 2] = strofreal(1::rows(mCD))
			cmCD = J(10, 2, "")
			cmCD[1,2] = "ID"
			cmCD[2,2] = "CD_L1"
			cmCD[3,2] = "CD_L2"
			cmCD[4,2] = "CD_L3"
			cmCD[5,2] = "CD_L4"
			cmCD[6,2] = "CD_L5"
			cmCD[7,2] = "CD_L6"
			cmCD[8,2] = "CD_L7"
			cmCD[9,2] = "CD_L8"
			cmCD[10,2] = "CD_L9"
		}
}		
		*mata: _matrix_list(mCD, rmCD, cmCD)	
		
	*mBCR - Best Clinical Response 
{	
		mata {
			mBCR = st_data(., "ID BCR_L1 BCR_L2 BCR_L3 BCR_L4 BCR_L5 BCR_L6 BCR_L7 BCR_L8 BCR_L9 BCR_SCT")		
			rmBCR = J(rows(mBCR), 2, "")
			rmBCR[., 2] = strofreal(1::rows(mBCR))
			cmBCR = J(11, 2, "")
			cmBCR[1,2] = "ID"
			cmBCR[2,2] = "BCR_L1"
			cmBCR[3,2] = "BCR_L2"
			cmBCR[4,2] = "BCR_L3"
			cmBCR[5,2] = "BCR_L4"
			cmBCR[6,2] = "BCR_L5"
			cmBCR[7,2] = "BCR_L6"
			cmBCR[8,2] = "BCR_L7"
			cmBCR[9,2] = "BCR_L8"
			cmBCR[10,2] = "BCR_L9"
			cmBCR[11,2] = "BCR_SCT"
		}
}			
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		
	*mCI - Chemotherapy Interval 
{		
		mata {
			mCI = st_data(., "ID CI_L1 CI_L2 CI_L3 CI_L4 CI_L5 CI_L6 CI_L7 CI_L8 CI_L9")		
			rmCI = J(rows(mCI), 2, "")
			rmCI[., 2] = strofreal(1::rows(mCI))
			cmCI = J(10, 2, "")
			cmCI[1,2] = "ID"
			cmCI[2,2] = "CI_L1"
			cmCI[3,2] = "CI_L2"
			cmCI[4,2] = "CI_L3"
			cmCI[5,2] = "CI_L4"
			cmCI[6,2] = "CI_L5"
			cmCI[7,2] = "CI_L6"
			cmCI[8,2] = "CI_L7"
			cmCI[9,2] = "CI_L8"
			cmCI[10,2] = "CI_L9"
		}
}		
		*mata: _matrix_list(mCI, rmCI, cmCI)
		
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
		*mata: _matrix_list(mState, rmState, cmState)
		*mata: _matrix_list(mCore, rmCore, cmCore)
		*mata: _matrix_list(mCore2, rmCore2, cmCore2)		
		*mata: _matrix_list(mCom, rmCom, cmCom)
		*mata: _matrix_list(mSCT, rmSCT, cmSCT)
		*mata: _matrix_list(mMNT, rmMNT, cmMNT)	
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
		
		*OS/TSD/TNE in Years
		*CD/CI in Days

**********
*Diagnosis (DN)
	scalar OMC = 2
	scalar Line = 0
	scalar LX = 0
	scalar NFT = 2
			
	di "DN - SCT"
		scalar m = "DN_SCT"
		scalar b = "bDN_SCT"
		scalar c = "cLO_"
		quietly do "Sub/New/SIM SCT DN.do"			
		*mata: _matrix_list(pDN_SCT, rpDN_SCT, cpDN_SCT)
		*mata: _matrix_list(oDN_SCT, roDN_SCT, coDN_SCT)

	di "DN - Chemo Interval"
		scalar m = "DN_CI"
		scalar b = "bDN_CI"
		scalar c = "cSU_"
		quietly do "Sub/New/SIM CI DN.do"
		*mata: _matrix_list(mDN_CI, rmDN_CI, cmDN_CI)
		*mata: _matrix_list(mNFT, rmNFT, cmNFT)	
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)

	di "DN - Overall Survival" 
		scalar m = "OS_DN"
		scalar b = "bOS"
		scalar c = "cSU_"	
		quietly do "Sub/New/SIM OS DN.do"	
		*mata: _matrix_list(mOS_DN, rmOS_DN, cmOS_DN)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "DN - Mortality"
		forvalues i = 1/`=Obs' {
			mata {
				if (mState[`i',2] <= `=OMC') { // State filter
					if (mTSD[`i',`=OMC'+1] > mOS[`i',`=OMC']) { // If TSD > OS...
						mMOR[`i',`=OMC'] = 1 // Patient dies
						if ((mAge[`i',2] + mOS[`i',`=OMC']) > `=Limit') mOS[`i',`=OMC'] = `=Limit' - mAge[`i',2] // Set mOS to max if Age > Limit
						mOC[`i',2] = mOS[`i',`=OMC'] // Set OC Time
						mOC[`i',3] = 1 // Set OC Outcome
						mTSD[`i',`=OMC'+1] = . // Clear TSD
						mTNE[`i',`=OMC'] = mOS[`i',`=OMC'] - mTSD[`i',`=OMC'] // Truncate TNE to OC_TIME
						mCI[`i', `=OMC'/2] = (mOC[`i',2] - mTSD[`i',`=OMC'])*365.25 // Overwrite CI
						mCore[`i', cSCT] = 0 // Overwrite SCT
					}
					if (mTSD[`i',`=OMC'+1] <= mOS[`i',`=OMC']) mMOR[`i',`=OMC'] = 0 // If TSD < OS, Patient alive
				}
			}
		}
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 1 Start (L1S)
	scalar OMC = 3

	di "L1S - Age"
		quietly do "Sub/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)

	di "L1S - Chemo Regimen"
		scalar m = "mL1_CR"
		scalar b = "bL1_CR"
		scalar o = "oL1_CR"
		scalar c = "cML_"
		quietly do "Sub/SIM CR L1.do"
		*mata: _matrix_list(mL1_CR, rmL1_CR, cmL1_CR)
		*mata: _matrix_list(mCR, rmCR, cmCR)
		
		if ("$Data" == "Population" & $Line == 1) {
			exit
		}
		
	di "L1S - Chemo Duration"
		scalar m = "mL1_CD"
		scalar c = "cSU_"
		quietly do "Sub/SIM CD L1.do"
		*mata: _matrix_list(mL1_CD, rmL1_CD, cmL1_CD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)

	di "L1S - Overall Survival"
		scalar m = "mOS_L1S"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM OS.do"
		*mata: _matrix_list(mOS_L1S, rmOS_L1S, cmOS_L1S)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L1S - Mortality"
		quietly do "Sub/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 1 End (L1E)
	scalar OMC = 4
	scalar Line = 1
	scalar LX = 1
	scalar NFT = 3

	di "L1E - Age"
		quietly do "Sub/SIM AGE.do"
		*mata: _matrix_list(mAge, rmAge, cmAge)	

	di "L1E - Best Clinical Response"
		scalar m = "mL1_BCR"
		scalar b = "bL1_BCR"
		scalar c = "cOL_"
		quietly do "Sub/SIM BCR L1.do"		
		*mata: _matrix_list(mL1_BCR, rmL1_BCR, cmL1_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)

	di "L1E - SCT"
		scalar m = "mL1_SCT"
		scalar b = "bL1_SCT"
		scalar c = "cLO_"
		quietly do "Sub/SIM SCT L1.do"			
		*mata: _matrix_list(mL1_SCT, rmL1_SCT, cmL1_SCT)	

	di "L1E - SCT Best Clinical Response"
		scalar m = "mSCT_BCR"
		scalar b = "bSCT_BCR"
		scalar c = "cOL_"
		quietly do "Sub/SIM BCR SCT.do"			
		*mata: _matrix_list(mSCT_BCR, rmSCT_BCR, cmSCT_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)

	di "L1E - MNT"
		scalar m = "mMNT"
		scalar b = "bMNT"
		scalar c = "cLO_"
		quietly do "Sub/SIM MNT.do"		
		*mata: _matrix_list(mMNT, rmMNT, cmMNT)
		*mata: _matrix_list(mCore, rmCore, cmCore)

	di "L1E - Chemo Interval"
		scalar m = "mL1_CI"
		scalar c = "cSU_"
		quietly do "Sub/SIM CI L1.do"			
		*mata: _matrix_list(mL1_CI, rmL1_CI, cmL1_CI)
		*mata: _matrix_list(mNFT, rmNFT, cmNFT)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)		
	
	di "L1E - Overall Survival"
		scalar m = "mOS_L1E"
		scalar b = "bOS"	
		scalar c = "cLO_"
		quietly do "Sub/SIM OS.do"	
		*mata: _matrix_list(mOS_L1E, rmOS_L1E, cmOS_L1E)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L1E - Mortality"
		quietly do "Sub/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 2 Start (L2S)
	scalar OMC = 5
		
	di "L2S - Age"
		quietly do "Sub/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)		
		
	di "L2S - Chemo Regimen"
		scalar m = "mL2_CR"
		scalar b = "bL2_CR"
		scalar o = "oL2_CR"
		scalar c = "cML_"
		quietly do "Sub/SIM CR L2.do"
		*mata: _matrix_list(mL2_CR, rmL2_CR, cmL2_CR)
		
		if ("$Analysis" == "DVd-Post" & "$Int" == "All" & "$Data" == "Population" & $Line == 2) {
			exit
		}

	di "L2S - Chemo Duration"
		scalar m = "mL2_CD"
		scalar b = "bL2_CD"
		scalar c = "cSU_"
		quietly do "Sub/SIM CD L2.do"
		*mata: _matrix_list(mL2_CD, rmL2_CD, cmL2_CD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)

	di "L2S - Overall Survival"
		scalar m = "mOS_L2S"
		scalar b = "bOS"	
		scalar c = "cLO_"
		quietly do "Sub/SIM OS.do"	
		*mata: _matrix_list(mOS_L2S, rmOS_L2S, cmOS_L2S)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L2S - Mortality"
		quietly do "Sub/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
	
**********
*Chemo Line 2 End (L2E)
	scalar OMC = 6
	scalar Line = 2
	scalar LX = 2
	scalar NFT = 4

	di "L2E - Age"
		quietly do "Sub/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)		
	
	di "L2E - Best Clinical Response"
		scalar m = "mL2_BCR"
		scalar b = "bL2_BCR"
		scalar c = "cOL_"
		quietly do "Sub/SIM BCR L2.do"
		*mata: _matrix_list(mL2_BCR, rmL2_BCR, cmL2_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)

	di "L2E - Chemo Interval"
		scalar m = "mL2_CI"
		scalar b = "bL2_CI"
		scalar c = "cSU_"
		quietly do "Sub/SIM CI L2.do"
		*mata: _matrix_list(mL2_CI, rmL2_CI, cmL2_CI)
		*mata: _matrix_list(mNFT, rmNFT, cmNFT)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)

	di "L2E - Overall Survival"
		scalar m = "mOS_L2E"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM OS.do"	
		*mata: _matrix_list(mOS_L2E, rmOS_L2E, cmOS_L2E)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L2E - Mortality"
		quietly do "Sub/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 3 Start (L3S)
	scalar OMC = 7
		
	di "L3S - Age"
		quietly do "Sub/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
		
	di "L3S - Chemo Regimen"
		scalar m = "mL3_CR"
		scalar b = "bL3_CR"
		scalar o = "oL3_CR"
		scalar c = "cML_"
		quietly do "Sub/SIM CR L3.do"
		
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
		quietly do "Sub/SIM CD L3.do"	
		*mata: _matrix_list(mL3_CD, rmL3_CD, cmL3_CD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L3S - Overall Survival"
		scalar m = "mOS_L3S"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM OS.do"
		*mata: _matrix_list(mOS_L3S, rmOS_L3S, cmOS_L3S)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L3S - Mortality"
		quietly do "Sub/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)	

**********
*Chemo Line 3 End (L3E)
	scalar OMC = 8
	scalar Line = 3
	scalar LX = 3
	scalar NFT = 5
		
	di "L3E - Age"
		quietly do "Sub/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
	
	di "L3E - Best Clinical Response"
		scalar m = "mL3_BCR"
		scalar b = "bL3_BCR"
		scalar c = "cOL_"
		quietly do "Sub/SIM BCR L3.do"		
		*mata: _matrix_list(mL3_BCR, rmL3_BCR, cmL3_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		
	di "L3E - Chemo Interval"
		scalar m = "mL3_CI"
		scalar b = "bL3_CI"
		scalar c = "cSU_"
		quietly do "Sub/SIM CI L3 L4.do"		
		*mata: _matrix_list(mL3_CI, rmL3_CI, cmL3_CI)
		*mata: _matrix_list(mNFT, rmNFT, cmNFT)	
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)

	di "L3E - Overall Survival"
		scalar m = "mOS_L3E"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM OS.do"
		*mata: _matrix_list(mOS_L3E, rmOS_L3E, cmOS_L3E)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L3E - Mortality"
		quietly do "Sub/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 4 Start (L4S)
	scalar OMC = 9
		
	di "L4S - Age"
		quietly do "Sub/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)		

	di "L4S - Chemo Regimen"
		scalar m = "mL4_CR"
		scalar b = "bL4_CR"
		scalar o = "oL4_CR"
		scalar c = "cML_"
		quietly do "Sub/SIM CR L4.do"

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
		quietly do "Sub/SIM CD L4.do"		
		*mata: _matrix_list(mL4_CD, rmL4_CD, cmL4_CD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L4S - Overall Survival"
		scalar m = "mOS_L4S"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM OS.do"
		*mata: _matrix_list(mOS_L4S, rmOS_L4S, cmOS_L4S)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L4S - Mortality"
		quietly do "Sub/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 4 End (L4E)
	scalar OMC = 10
	scalar Line = 4
	scalar LX = 4
	scalar NFT = 6
		
	di "L4E - Age"
		quietly do "Sub/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
	
	di "L4E - Best Clinical Response"
		scalar m = "mL4_BCR"
		scalar b = "bL4_BCR"
		scalar c = "cOL_"
		quietly do "Sub/SIM BCR L4.do"	
		*mata: _matrix_list(mL4_BCR, rmL4_BCR, cmL4_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		
	di "L4E - Chemo Interval"
		scalar m = "mL4_CI"
		scalar b = "bL4_CI"
		scalar c = "cSU_"
		quietly do "Sub/SIM CI L3 L4.do"	
		*mata: _matrix_list(mL4_CI, rmL4_CI, cmL4_CI)
		*mata: _matrix_list(mNFT, rmNFT, cmNFT)	
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)

	di "L4E - Overall Survival"
		scalar m = "mOS_L4E"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM OS.do"	
		*mata: _matrix_list(mOS_L4E, rmOS_L4E, cmOS_L4E)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L4E - Mortality"
		quietly do "Sub/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
			
**********
*Chemo Line 5 Start (L5S)
	scalar OMC = 11
		
	di "L5S - Age" 
		quietly do "Sub/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
		
	di "L5S - Chemo Regimen"
		forval i = 1/`=Obs' {
			mata {
				if 	(mMOR[`i',`=OMC'-1] == 0) mCore[`i',cCR] = 0
				if 	(mMOR[`i',`=OMC'-1] != 0) mCore[`i',cCR] = .
				if 	(mMOR[`i',`=OMC'-1] == 0) mCR[`i',`=LX'+2] = 0
			}
		}

	di "L5S - Chemo Duration"
		scalar m = "mL5_CD"
		scalar b = "bLX_CD"
		scalar c = "cSU_"
		quietly do "Sub/SIM CD LX (5 - 9).do"		
		*mata: _matrix_list(mL5_CD, rmL5_CD, cmL5_CD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L5S - Overall Survival"
		scalar m = "mOS_L5S"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM OS.do"
		*mata: _matrix_list(mOS_L5S, rmOS_L5S, cmOS_L5S)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L5S - Mortality"
		quietly do "Sub/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 5 End (L5E)
	scalar OMC = 12
	scalar Line = 5
	scalar LX = 5
	scalar NFT = 7
		
	di "L5E - Age"
		quietly do "Sub/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	

	di "L5E - Best Clinical Response"
		scalar m = "mL5_BCR"
		scalar b = "bLX_BCR"
		scalar c = "cOL_"
		quietly do "Sub/SIM BCR LX (5 - 9).do"
		*mata: _matrix_list(mL5_BCR, rmL5_BCR, cmL5_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
	
	di "L5E - Chemo Interval"
		scalar m = "mL5_CI"
		scalar b = "bLX_CI"
		scalar c = "cSU_"
		quietly do "Sub/SIM CI LX (5 - 8).do"				
		*mata: _matrix_list(mL5_CI, rmL5_CI, cmL5_CI)
		*mata: _matrix_list(mNFT, rmNFT, cmNFT)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)

	di "L5E - Overall Survival"
		scalar m = "mOS_L5E"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM OS.do"
		*mata: _matrix_list(mOS_L5E, rmOS_L5E, cmOS_L5E)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L5E - Mortality"
		quietly do "Sub/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
		
**********
*Chemo Line 6 Start (L6S) 
	scalar OMC = 13
		
	di "L6S - Age" 
		quietly do "Sub/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
		
	di "L6S - Chemo Regimen"
		forvalues i = 1/`=Obs' {
			mata {
				if 	(mMOR[`i',`=OMC'-1] == 0) mCore[`i',cCR] = 0
				if 	(mMOR[`i',`=OMC'-1] != 0) mCore[`i',cCR] = .
				if 	(mMOR[`i',`=OMC'-1] == 0) mCR[`i',`=LX'+2] = 0
			}
		}

	di "L6S - Chemo Duration"
		scalar m = "mL6_CD"
		scalar b = "bLX_CD" 		
		scalar c = "cSU_"
		quietly do "Sub/SIM CD LX (5 - 9).do"	
		*mata: _matrix_list(mL6_CD, rmL6_CD, cmL6_CD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L6S - Overall Survival"
		scalar m = "mOS_L6S"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM OS.do"
		*mata: _matrix_list(mOS_L6S, rmOS_L6S, cmOS_L6S)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L6S - Mortality"
		quietly do "Sub/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 6 End (L6E)
	scalar OMC = 14
	scalar Line = 6
	scalar LX = 6
	scalar NFT = 8
		
	di "L6E - Age"
		quietly do "Sub/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
	
	di "L6E - Best Clinical Response"
		scalar m = "mL6_BCR"
		scalar b = "bLX_BCR"
		scalar c = "cOL_"
		quietly do "Sub/SIM BCR LX (5 - 9).do"	 	
		*mata: _matrix_list(mL6_BCR, rmL6_BCR, cmL6_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		
	di "L6E - Chemo Interval"
		scalar m = "mL6_CI"
		scalar b = "bLX_CI"
		scalar c = "cSU_"
		quietly do "Sub/SIM CI LX (5 - 8).do"				
		*mata: _matrix_list(mL6_CI, rmL6_CI, cmL6_CI)
		*mata: _matrix_list(mNFT, rmNFT, cmNFT)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)	

	di "L6E - Overall Survival"
		scalar m = "mOS_L6E"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM OS.do"
		*mata: _matrix_list(mOS_L6E, rmOS_L6E, cmOS_L6E)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L6E - Mortality" 
		quietly do "Sub/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
					
**********
*Chemo Line 7 Start (L7S) 
	scalar OMC = 15
		
	di "L7S - Age" 
		quietly do "Sub/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
		
	di "L7S - Chemo Regimen"
		forvalues i = 1/`=Obs' {
			mata {
				if 	(mMOR[`i',`=OMC'-1] == 0) mCore[`i',cCR] = 0
				if 	(mMOR[`i',`=OMC'-1] != 0) mCore[`i',cCR] = .
				if 	(mMOR[`i',`=OMC'-1] == 0) mCR[`i',`=LX'+2] = 0
			}
		}

	di "L7S - Chemo Duration"
		scalar m = "mL7_CD"
		scalar b = "bLX_CD" 		
		scalar c = "cSU_"
		quietly do "Sub/SIM CD LX (5 - 9).do"	
		*mata: _matrix_list(mL7_CD, rmL7_CD, cmL7_CD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L7S - Overall Survival"
		scalar m = "mOS_L7S"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM OS.do"
		*mata: _matrix_list(mOS_L7S, rmOS_L7S, cmOS_L7S)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L7S - Mortality"
		quietly do "Sub/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 7 End (L7E)
	scalar OMC = 16
	scalar LX = 7
	scalar NFT = 9
		
	di "L7E - Age"
		quietly do "Sub/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
	
	di "L7E - Best Clinical Response" 
		scalar m = "mL7_BCR"
		scalar b = "bLX_BCR"
		scalar c = "cOL_"
		quietly do "Sub/SIM BCR LX (5 - 9).do"	 	
		*mata: _matrix_list(mL7_BCR, rmL7_BCR, cmL7_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		
	di "L7E - Chemo Interval"
		scalar m = "mL7_CI"
		scalar b = "bLX_CI"
		scalar c = "cSU_"
		quietly do "Sub/SIM CI LX (5 - 8).do"				
		*mata: _matrix_list(mL7_CI, rmL7_CI, cmL7_CI)
		*mata: _matrix_list(mNFT, rmNFT, cmNFT)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)		

	di "L7E - Overall Survival"
		scalar m = "mOS_L7E"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM OS.do"
		*mata: _matrix_list(mOS_L7E, rmOS_L7E, cmOS_L7E)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L7E - Mortality" 
		quietly do "Sub/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
		
**********
*Chemo Line 8 Start (L8S) 
	scalar OMC = 17
		
	di "L8S - Age" 
		quietly do "Sub/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
		
	di "L8S - Chemo Regimen"
		forvalues i = 1/`=Obs' {
			mata {
				if 	(mMOR[`i',`=OMC'-1] == 0) mCore[`i',cCR] = 0
				if 	(mMOR[`i',`=OMC'-1] != 0) mCore[`i',cCR] = .
				if 	(mMOR[`i',`=OMC'-1] == 0) mCR[`i',`=LX'+2] = 0
			}
		}

	di "L8S - Chemo Duration"
		scalar m = "mL8_CD"
		scalar b = "bLX_CD" 		
		scalar c = "cSU_"
		quietly do "Sub/SIM CD LX (5 - 9).do"	
		*mata: _matrix_list(mL8_CD, rmL8_CD, cmL8_CD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L8S - Overall Survival"
		scalar m = "mOS_L8S"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM OS.do"
		*mata: _matrix_list(mOS_L8S, rmOS_L8S, cmOS_L8S)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L8S - Mortality"
		quietly do "Sub/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 8 End (L8E)
	scalar OMC = 18
	scalar LX = 8
	scalar NFT = 10
		
	di "L8E - Age"
		quietly do "Sub/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
	
	di "L8E - Best Clinical Response"
		scalar m = "mL8_BCR"
		scalar b = "bLX_BCR"
		scalar c = "cOL_"
		quietly do "Sub/SIM BCR LX (5 - 9).do"	
		*mata: _matrix_list(mL8_BCR, rmL8_BCR, cmL8_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		
	di "L8E - Chemo Interval"
		scalar m = "mL8_CI"
		scalar b = "bLX_CI"
		scalar c = "cSU_"
		quietly do "Sub/SIM CI LX (5 - 8).do"				
		*mata: _matrix_list(mL8_CI, rmL8_CI, cmL8_CI)
		*mata: _matrix_list(mNFT, rmNFT, cmNFT)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)	

	di "L8E - Overall Survival"
		scalar m = "mOS_L8E"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM OS.do"
		*mata: _matrix_list(mOS_L8E, rmOS_L8E, cmOS_L8E)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L8E - Mortality" 
		quietly do "Sub/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)	

**********
*Chemo Line 9 Start (L9S) 
	scalar OMC = 19
		
	di "L9S - Age" 
		quietly do "Sub/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
		
	di "L9S - Chemo Regimen"
		forvalues i = 1/`=Obs' {
			mata {
				if 	(mMOR[`i',`=OMC'-1] == 0) mCore[`i',cCR] = 0
				if 	(mMOR[`i',`=OMC'-1] != 0) mCore[`i',cCR] = .
				if 	(mMOR[`i',`=OMC'-1] == 0) mCR[`i',`=LX'+2] = 0
			}
		}

	di "L9S - Chemo Duration"
		scalar m = "mL9_CD"
		scalar b = "bLX_CD" 		
		scalar c = "cSU_"
		quietly do "Sub/SIM CD LX (5 - 9).do"	
		*mata: _matrix_list(mL9_CD, rmL9_CD, cmL9_CD)
		*mata: _matrix_list(mTNE, rmTNE, cmTNE)
		*mata: _matrix_list(mTSD, rmTSD, cmTSD)
		
	di "L9S - Overall Survival"
		scalar m = "mOS_L9S"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM OS.do"
		*mata: _matrix_list(mOS_L9S, rmOS_L9S, cmOS_L9S)
		*mata: _matrix_list(mOS, rmOS, cmOS)
	
	di "L9S - Mortality"
		quietly do "Sub/SIM MORT.do"	
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)

**********
*Chemo Line 9 End (L9E)
	scalar OMC = 20
	scalar LX = 9
	scalar NFT = 11
		
	di "L9E - Age"
		quietly do "Sub/SIM AGE.do"	
		*mata: _matrix_list(mAge, rmAge, cmAge)	
	
	di "L9E - Best Clinical Response"
		scalar m = "mL9_BCR"
		scalar b = "bLX_BCR"
		scalar c = "cOL_"
		quietly do "Sub/SIM BCR LX (5 - 9).do"
		*mata: _matrix_list(mL9_BCR, rmL9_BCR, cmL9_BCR)
		*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		
	di "L9E - Overall Survival"
		scalar m = "mOS_L9E"
		scalar b = "bOS"
		scalar c = "cLO_"
		quietly do "Sub/SIM OS.do"
		*mata: _matrix_list(mOS_L9E, rmOS_L9E, cmOS_L9E)
		*mata: _matrix_list(mOS, rmOS, cmOS)

	di "L9E - Mortality" 
		*Everyone still alive dies at predicted OS or Limit
		forvalues i = 1/`=Obs'{
			mata {
				if (mMOR[`i',`=OMC'-1] == 0) { // Alive filter
						if	((mAge[`i',2] + mOS[`i',`=OMC']) > `=Limit') mOS[`i',`=OMC'] = `=Limit' - mAge[`i',2] // Set mOS to max if Age > Limit
					mMOR[`i',`=OMC'] = 1 // Patient dies
					mOC[`i',2] = mOS[`i',`=OMC'] // Set OC Time
					mOC[`i',3] = 1 // Set OC Outcome
				}		
			}	
		}
		*mata: _matrix_list(mMOR, rmMOR, cmMOR)
		*mata: _matrix_list(mOC, rmOC, cmOC)		

	end
	
**********
*Process function
	capture program drop process
	program define process

	*Create mSum in Mata 
		mata: mSum = mCore , mAge , mOS , mTNE , mTSD , mMOR , mOC , mCR , mCD , mBCR , mCI , mState, mSCT
	
	*Convert mSum to stSum
		mata: st_matrix("stSum", mSum)
		drop _all
		
	*Convert stSum to variables
		svmat double stSum
	
	*Name variables
		local varnames ID Age Male ECOGcc RISS SCT MNT CR CD BCR ///
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
			State_ID State DateDN ///
			SCT_DN SCT_L1
		
		local varlength : word count `varnames'
		
		forvalues i = 1/`varlength'{
			local currentvar : word `i' of `varnames'
			rename stSum`i' `currentvar'
		}		
	
		format DateDN %td
	
	*Drop unnecessary variables
		drop Age_ID OS_ID TNE_ID TSD_ID MOR_ID OC_ID CR_ID CD_ID BCR_ID CI_ID State_ID CR CD BCR
		order ID Age Male ECOGcc RISS SCT MNT
		
	*Label
		label values State State_lbl
	
	*Generate Dates
		gen DateL1S = DateDN + (TNE_DN*365.25)
		gen DateL1E = DateL1S + (TNE_L1S*365.25)
		gen DateL2S = DateL1E + (TNE_L1E*365.25)
		gen DateL2E = DateL2S + (TNE_L2S*365.25)
		gen DateL3S = DateL2E + (TNE_L2E*365.25)
		gen DateL3E = DateL3S + (TNE_L3S*365.25)
		gen DateL4S = DateL3E + (TNE_L3E*365.25)
		gen DateL4E = DateL4S + (TNE_L4S*365.25)
		gen DateL5S = DateL4E + (TNE_L4E*365.25)
		gen DateL5E = DateL5S + (TNE_L5S*365.25)
		gen DateL6S = DateL5E + (TNE_L5E*365.25)
		gen DateL6E = DateL6S + (TNE_L6S*365.25)
		gen DateL7S = DateL6E + (TNE_L6E*365.25)
		gen DateL7E = DateL7S + (TNE_L7S*365.25)
		gen DateL8S = DateL7E + (TNE_L7E*365.25)
		gen DateL8E = DateL8S + (TNE_L8S*365.25)
		gen DateL9S = DateL8E + (TNE_L8E*365.25)
		gen DateL9E = DateL9S + (TNE_L9S*365.25)
		gen DateSCT = DateL1E + 1 if(SCT == 1) // Fix DateSCT 1 day after DateL1E
		gen DateMOR = DateDN + (OC_TIME*365.25)
		format Date* %td
	
	*Generate Years
		gen YearDN = yofd(DateDN)
		gen YearL1 = yofd(DateL1S)
		gen YearL2 = yofd(DateL2S)
		gen YearL3 = yofd(DateL3S)
		gen YearL4 = yofd(DateL4S)
		gen YearL5 = yofd(DateL5S)
		gen YearL6 = yofd(DateL6S)
		gen YearL7 = yofd(DateL7S)
		gen YearL8 = yofd(DateL8S)
		gen YearL9 = yofd(DateL9S)
		gen YearSCT = yofd(DateSCT)
		gen YearMOR = yofd(DateMOR)
	
	end
	
**********
*Execute based on arguments	

	if("$Boot" == "0") {
		
		*Load coefficient
			mata: mata matuse "Analysis/$Analysis/Coefficients/$Analysis $Coeffs Coefficients"

		*Load data
			use "Analysis/$Analysis/$Analysis $Data.dta", replace
			
		*State filter
			keep if State <= ($Line * 2) + 1
			replace ID = _n

		*ID filter
			keep if ID >= $MinID & ID <= $MaxID
			quietly sum ID
			scalar Obs = r(N)
	
		*Age limit
			scalar Limit = 100
			
		*Execute functions
			simulation
			process
		
		*Save Simulated Dataset
			save "Analysis/$Analysis/Simulated/$Int $Line $Data $MinID $MaxID.dta", replace
	}
	
	else if("$Boot" == "1") {
		
		forvalues b = $MinBS / $MaxBS {
			
			*Load coefficients
				clear mata
				mata: mata matuse "Analysis/$Analysis/Coefficients/$Analysis $Coeffs Coefficients B`b'"
				
			*Load data
				use "Analysis/$Analysis/$Analysis $Data.dta", replace
						
			*State filter
				keep if State <= ($Line * 2) + 1
				replace ID = _n
				
			*ID filter
				keep if ID >= $MinID & ID <= $MaxID
				quietly sum ID
				scalar Obs = r(N)
				
			*Age limit
				scalar Limit = 100	
			
			*Execute functions
				simulation
				process
				
			*Change ID to reflect bootstrap sample
				replace ID = ID + `b'0000000	
				
			*Save Simulated Dataset
				save "Analysis/$Analysis/Simulated/$Int $Line $Data B`b'.dta", replace
		}
	}
	

