**********
	*EpiMAP Myeloma - Matrix Setup
**********

capture program drop matrix_setup
program define matrix_setup
	
	di "Setting up matrices"
	
	*Additional variables needed for mCore2 to transfer to matrix multiplication for XB
		capture gen ECOGcc0 = (ECOGcc == 0)
		capture gen ECOGcc1 = (ECOGcc == 1)
		capture gen ECOGcc2 = (ECOGcc == 2)
		capture gen RISS1 = (RISS == 1)
		capture gen RISS2 = (RISS == 2)
		capture gen RISS3 = (RISS == 3)
		capture gen Cons = 1

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
	
	*mCom - Comorbidity matrix
{
		mata {
			mCom = st_data(.,"Age70 Age75 CMc")
			rmCom = J(rows(mCom),2,"")
			rmCom[.,2] = strofreal(1::rows(mCom))
			cmCom = J(3,2,"")
			cmCom[1,1] = "Age70" 
			cmCom[2,1] = "Age75"
			cmCom[3,1] = "CMc" 
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

end
