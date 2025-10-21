**********
	*EpiMAP Myeloma - Matrix Setup
**********

capture program drop matrix_setup
program define matrix_setup
	
	di "Setting up matrices"

	*mState - State matrix
{
	mata {
		mState = st_data(.,"State DateDN")
		rmState = J(rows(mState),2,"")
		rmState[.,2] = strofreal(mState[.,1])
		cmState = J(2,2,"")
		cmState[1,1] = "State"
		cmState[2,1] = "DateDN"
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
		mAge = st_data(., "Age_DN Age_L1S Age_L1E Age_L2S Age_L2E Age_L3S Age_L3E Age_L4S Age_L4E Age_L5S Age_L5E Age_L6S Age_L6E Age_L7S Age_L7E Age_L8S Age_L8E Age_L9S Age_L9E")
		rmAge = J(rows(mAge), 2, "")
		rmAge[., 2] = strofreal(1::rows(mAge))
		cmAge = J(19, 2, "")
		cmAge[1,1] = "Age_DN"
		cmAge[2,1] = "Age_L1S"
		cmAge[3,1] = "Age_L1E"
		cmAge[4,1] = "Age_L2S"
		cmAge[5,1] = "Age_L2E"
		cmAge[6,1] = "Age_L3S"
		cmAge[7,1] = "Age_L3E"
		cmAge[8,1] = "Age_L4S"
		cmAge[9,1] = "Age_L4E"
		cmAge[10,1] = "Age_L5S"
		cmAge[11,1] = "Age_L5E"
		cmAge[12,1] = "Age_L6S"
		cmAge[13,1] = "Age_L6E"
		cmAge[14,1] = "Age_L7S"
		cmAge[15,1] = "Age_L7E"
		cmAge[16,1] = "Age_L8S"
		cmAge[17,1] = "Age_L8E"
		cmAge[18,1] = "Age_L9S"
		cmAge[19,1] = "Age_L9E"
	}
}		
	*mata: _matrix_list(mAge, rmAge, cmAge)		

	*mOS - Overall Survival (from DN)
{	
	mata {
		mOS = J(rows(mCore),19,.)
		rmOS = J(rows(mOS), 2, "")
		rmOS[., 2] = strofreal(1::rows(mOS))
		cmOS = J(19, 2, "")
		cmOS[1,1] = "OS_DN"
		cmOS[2,1] = "OS_L1S"
		cmOS[3,1] = "OS_L1E"
		cmOS[4,1] = "OS_L2S"
		cmOS[5,1] = "OS_L2E"
		cmOS[6,1] = "OS_L3S"
		cmOS[7,1] = "OS_L3E"
		cmOS[8,1] = "OS_L4S"
		cmOS[9,1] = "OS_L4E"
		cmOS[10,1] = "OS_L5S"
		cmOS[11,1] = "OS_L5E"
		cmOS[12,1] = "OS_L6S"
		cmOS[13,1] = "OS_L6E"
		cmOS[14,1] = "OS_L7S"
		cmOS[15,1] = "OS_L7E"
		cmOS[16,1] = "OS_L8S"
		cmOS[17,1] = "OS_L8E"
		cmOS[18,1] = "OS_L9S"
		cmOS[19,1] = "OS_L9E"
	}
}			
	*mata: _matrix_list(mOS, rmOS, cmOS)		

	*mTNE - Time to Next Event
{	
	mata {
		mTNE = st_data(., "TNE_DN TNE_L1S TNE_L1E TNE_L2S TNE_L2E TNE_L3S TNE_L3E TNE_L4S TNE_L4E TNE_L5S TNE_L5E TNE_L6S TNE_L6E TNE_L7S TNE_L7E TNE_L8S TNE_L8E TNE_L9S TNE_L9E")
		rmTNE = J(rows(mTNE), 2, "")
		rmTNE[., 2] = strofreal(1::rows(mTNE))
		cmTNE = J(19, 2, "")
		cmTNE[1,1] = "TNE_DN"
		cmTNE[2,1] = "TNE_L1S"
		cmTNE[3,1] = "TNE_L1E"
		cmTNE[4,1] = "TNE_L2S"
		cmTNE[5,1] = "TNE_L2E"
		cmTNE[6,1] = "TNE_L3S"
		cmTNE[7,1] = "TNE_L3E"
		cmTNE[8,1] = "TNE_L4S"
		cmTNE[9,1] = "TNE_L4E"
		cmTNE[10,1] = "TNE_L5S"
		cmTNE[11,1] = "TNE_L5E"
		cmTNE[12,1] = "TNE_L6S"
		cmTNE[13,1] = "TNE_L6E"
		cmTNE[14,1] = "TNE_L7S"
		cmTNE[15,1] = "TNE_L7E"
		cmTNE[16,1] = "TNE_L8S"
		cmTNE[17,1] = "TNE_L8E"
		cmTNE[18,1] = "TNE_L9S"
		cmTNE[19,1] = "TNE_L9E"
	}
}			
	*mata: _matrix_list(mTNE, rmTNE, cmTNE)	
		
	*mTSD - Time Since Diagnosis
{	
	mata {
		mTSD = st_data(., "TSD_DN TSD_L1S TSD_L1E TSD_L2S TSD_L2E TSD_L3S TSD_L3E TSD_L4S TSD_L4E TSD_L5S TSD_L5E TSD_L6S TSD_L6E TSD_L7S TSD_L7E TSD_L8S TSD_L8E TSD_L9S TSD_L9E")
		mTSD[., 1] = J(rows(mTSD), 1, 0) // Set to TSD_DN to 0
		rmTSD = J(rows(mTSD), 2, "")
		rmTSD[., 2] = strofreal(1::rows(mTSD))
		cmTSD = J(19, 2, "")
		cmTSD[1,1] = "TSD_DN"
		cmTSD[2,1] = "TSD_L1S"
		cmTSD[3,1] = "TSD_L1E"
		cmTSD[4,1] = "TSD_L2S"
		cmTSD[5,1] = "TSD_L2E"
		cmTSD[6,1] = "TSD_L3S"
		cmTSD[7,1] = "TSD_L3E"
		cmTSD[8,1] = "TSD_L4S"
		cmTSD[9,1] = "TSD_L4E"
		cmTSD[10,1] = "TSD_L5S"
		cmTSD[11,1] = "TSD_L5E"
		cmTSD[12,1] = "TSD_L6S"
		cmTSD[13,1] = "TSD_L6E"
		cmTSD[14,1] = "TSD_L7S"
		cmTSD[15,1] = "TSD_L7E"
		cmTSD[16,1] = "TSD_L8S"
		cmTSD[17,1] = "TSD_L8E"
		cmTSD[18,1] = "TSD_L9S"
		cmTSD[19,1] = "TSD_L9E"
	}
}			
	*mata: _matrix_list(mTSD, rmTSD, cmTSD)		
		
	*mMOR - Mortality (=1 if patient dies before next event) 
{	
	mata {
		mMOR = st_data(., "MOR_DN MOR_L1S MOR_L1E MOR_L2S MOR_L2E MOR_L3S MOR_L3E MOR_L4S MOR_L4E MOR_L5S MOR_L5E MOR_L6S MOR_L6E MOR_L7S MOR_L7E MOR_L8S MOR_L8E MOR_L9S MOR_L9E")
		rmMOR = J(rows(mMOR), 2, "")
		rmMOR[., 2] = strofreal(1::rows(mMOR))
		cmMOR = J(19, 2, "")
		cmMOR[1,1] = "MOR_DN"
		cmMOR[2,1] = "MOR_L1S"
		cmMOR[3,1] = "MOR_L1E"
		cmMOR[4,1] = "MOR_L2S"
		cmMOR[5,1] = "MOR_L2E"
		cmMOR[6,1] = "MOR_L3S"
		cmMOR[7,1] = "MOR_L3E"
		cmMOR[8,1] = "MOR_L4S"
		cmMOR[9,1] = "MOR_L4E"
		cmMOR[10,1] = "MOR_L5S"
		cmMOR[11,1] = "MOR_L5E"
		cmMOR[12,1] = "MOR_L6S"
		cmMOR[13,1] = "MOR_L6E"
		cmMOR[14,1] = "MOR_L7S"
		cmMOR[15,1] = "MOR_L7E"
		cmMOR[16,1] = "MOR_L8S"
		cmMOR[17,1] = "MOR_L8E"
		cmMOR[18,1] = "MOR_L9S"
		cmMOR[19,1] = "MOR_L9E"
	}
}			
	*mata: _matrix_list(mMOR, rmMOR, cmMOR)

	*mOC - Outcome
{	
	mata {
		mOC = J(rows(mCore),2,.)				
		rmOC = J(rows(mOC), 2, "")
		rmOC[., 2] = strofreal(1::rows(mOC))
		cmOC = J(2, 2, "")
		cmOC[1,1] = "Time"
		cmOC[2,1] = "Mort"
	}
}		
	*mata: _matrix_list(mOC, rmOC, cmOC)	
		
	*mTXR - Treatment Regimen 
{		
	mata {
		mTXR = st_data(., "CR_L1 CR_L2 CR_L3 CR_L4 CR_L5 CR_L6 CR_L7 CR_L8 CR_L9")		
		rmTXR = J(rows(mTXR), 2, "")
		rmTXR[., 2] = strofreal(1::rows(mTXR))
		cmTXR = J(9, 2, "")
		cmTXR[1,1] = "TXR_L1"
		cmTXR[2,1] = "TXR_L2"
		cmTXR[3,1] = "TXR_L3"
		cmTXR[4,1] = "TXR_L4"
		cmTXR[5,1] = "TXR_L5"
		cmTXR[6,1] = "TXR_L6"
		cmTXR[7,1] = "TXR_L7"
		cmTXR[8,1] = "TXR_L8"
		cmTXR[9,1] = "TXR_L9"
	}
}		
	*mata: _matrix_list(mTXR, rmTXR, cmTXR)
		
	*mTXD - Treatment Duration
{		
	mata {
		mTXD = st_data(., "CD_L1 CD_L2 CD_L3 CD_L4 CD_L5 CD_L6 CD_L7 CD_L8 CD_L9")		
		rmTXD = J(rows(mTXD), 2, "")
		rmTXD[., 2] = strofreal(1::rows(mTXD))
		cmTXD = J(9, 2, "")
		cmTXD[1,1] = "TXD_L1"
		cmTXD[2,1] = "TXD_L2"
		cmTXD[3,1] = "TXD_L3"
		cmTXD[4,1] = "TXD_L4"
		cmTXD[5,1] = "TXD_L5"
		cmTXD[6,1] = "TXD_L6"
		cmTXD[7,1] = "TXD_L7"
		cmTXD[8,1] = "TXD_L8"
		cmTXD[9,1] = "TXD_L9"
	}
}		
	*mata: _matrix_list(mTXD, rmTXD, cmTXD)	
		
	*mBCR - Best Clinical Response 
{	
	mata {
		mBCR = st_data(., "BCR_L1 BCR_L2 BCR_L3 BCR_L4 BCR_L5 BCR_L6 BCR_L7 BCR_L8 BCR_L9 BCR_SCT")		
		rmBCR = J(rows(mBCR), 2, "")
		rmBCR[., 2] = strofreal(1::rows(mBCR))
		cmBCR = J(10, 2, "")
		cmBCR[1,1] = "BCR_L1"
		cmBCR[2,1] = "BCR_L2"
		cmBCR[3,1] = "BCR_L3"
		cmBCR[4,1] = "BCR_L4"
		cmBCR[5,1] = "BCR_L5"
		cmBCR[6,1] = "BCR_L6"
		cmBCR[7,1] = "BCR_L7"
		cmBCR[8,1] = "BCR_L8"
		cmBCR[9,1] = "BCR_L9"
		cmBCR[10,1] = "BCR_SCT"
	}
}			
	*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		
	*mTFI - Treatment-free Interval 
{		
	mata {
		mTFI = st_data(., "CI_L1 CI_L2 CI_L3 CI_L4 CI_L5 CI_L6 CI_L7 CI_L8 CI_L9")		
		rmTFI = J(rows(mTFI), 2, "")
		rmTFI[., 2] = strofreal(1::rows(mTFI))
		cmTFI = J(9, 2, "")
		cmTFI[1,1] = "TFI_DN"
		cmTFI[2,1] = "TFI_L1"
		cmTFI[3,1] = "TFI_L2"
		cmTFI[4,1] = "TFI_L3"
		cmTFI[5,1] = "TFI_L4"
		cmTFI[6,1] = "TFI_L5"
		cmTFI[7,1] = "TFI_L6"
		cmTFI[8,1] = "TFI_L7"
		cmTFI[9,1] = "TFI_L8"
	}
}		
	*mata: _matrix_list(mCI, rmCI, cmCI)
		
	*mSU - Survival
{
		mata: mSU = J(`=Obs',1,0)
		mata: mSU = mSU , J(`=Obs',2,.)
		mata: cmSU = J(3,2,"")
		
		*Label
			mata: cmSU[1,1] = "SU_XB"
			mata: cSU_XB = cols(mCore) + 1
			mata: cmSU[2,1] = "SU_RN"
			mata: cSU_RN = cols(mCore) + 2
			mata: cmSU[3,1] = "SU_OC"
			mata: cSU_OC = cols(mCore) + 3
}		
		*mata: _matrix_list(mSU, rmCore, cmSU)
	
	*mLO - Logit
{
		mata: mLO = J(`=Obs',1,0)
		mata: mLO = mLO , J(`=Obs',3,.)
		mata: cmLO = J(4,2,"")
		
		*Label
			mata: cmLO[1,1] = "LO_XB"
			mata: cLO_XB = cols(mCore) + 1
			mata: cmLO[2,1] = "LO_PR"
			mata: cLO_PR = cols(mCore) + 2
			mata: cmLO[3,1] = "LO_RN"
			mata: cLO_RN = cols(mCore) + 3
			mata: cmLO[4,1] = "LO_OC"
			mata: cLO_OC = cols(mCore) + 4
}			
		*mata: _matrix_list(mLO, rmCore, cmLO)
	
	*mML - Multinomial Logit
{
		mata: mML = J(`=Obs',4,0)
		mata: mML = mML , J(`=Obs',6,.)
		mata: cmML = J(10,2,"")
		
		*Label
			mata: cmML[1,1] = "ML_XB1"
			mata: cML_XB1 = cols(mCore) + 1
			mata: cmML[2,1] = "ML_XB2"
			mata: cML_XB2 = cols(mCore) + 2
			mata: cmML[3,1] = "ML_XB3"
			mata: cML_XB3 = cols(mCore) + 3
			mata: cmML[4,1] = "ML_XB4"
			mata: cML_XB4 = cols(mCore) + 4
			mata: cmML[5,1] = "ML_PR1"
			mata: cML_PR1 = cols(mCore) + 5
			mata: cmML[6,1] = "ML_PR2"
			mata: cML_PR2 = cols(mCore) + 6
			mata: cmML[7,1] = "ML_PR3"
			mata: cML_PR3 = cols(mCore) + 7
			mata: cmML[8,1] = "ML_PR4"
			mata: cML_PR4 = cols(mCore) + 8
			mata: cmML[9,1] = "ML_RN"
			mata: cML_RN = cols(mCore) + 9
			mata: cmML[10,1] = "ML_OC"
			mata: cML_OC = cols(mCore) + 10
}			
		*mata: _matrix_list(mML, rmCore, cmML)
	
	*mOL - Ordered Logit
{
		mata: mOL = J(`=Obs',1,0)
		mata: mOL = mOL , J(`=Obs',8,.)
		mata: cmOL = J(9,2,"")
		
		*Label
			mata: cmOL[1,1] = "OL_XB"
			mata: cOL_XB = cols(mCore) + 1
			mata: cmOL[2,1] = "OL_PR1"
			mata: cOL_PR1 = cols(mCore) + 2
			mata: cmOL[3,1] = "OL_PR2"
			mata: cOL_PR2 = cols(mCore) + 3
			mata: cmOL[4,1] = "OL_PR3"
			mata: cOL_PR3 = cols(mCore) + 4
			mata: cmOL[5,1] = "OL_PR4"
			mata: cOL_PR4 = cols(mCore) + 5
			mata: cmOL[6,1] = "OL_PR5"
			mata: cOL_PR5 = cols(mCore) + 6
			mata: cmOL[7,1] = "OL_PR6"
			mata: cOL_PR6 = cols(mCore) + 7
			mata: cmOL[8,1] = "OL_RN"
			mata: cOL_RN = cols(mCore) + 8
			mata: cmOL[9,1] = "OL_OC"
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
