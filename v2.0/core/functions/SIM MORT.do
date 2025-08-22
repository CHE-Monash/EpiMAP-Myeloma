**********
	*SIM MORT
**********

	forval i = 1/`=Obs' {
		mata {
			if (mState[`i',2] <= `=OMC') { // State filter
				if (mMOR[`i',`=OMC'-1] == 0 & mTSD[`i',`=OMC'+1] > mOS[`i',`=OMC']) { // If TSD > OS
					mMOR[`i',`=OMC'] = 1 // Patient dies
					if ((mAge[`i',2] + mOS[`i',`=OMC']) > `=Limit') mOS[`i',`=OMC'] = `=Limit' - mAge[`i',2] // Set mOS to max if Age > Limit
					mOC[`i',2] = mOS[`i',`=OMC'] // Set OC Time
					mOC[`i',3] = 1 // Set OC Outcome
					mTSD[`i',`=OMC'+1] = . // Clear TSD
					mTNE[`i',`=OMC'] = mOS[`i',`=OMC'] - mTSD[`i',`=OMC'] // Truncate TNE to OC_TIME
					mCore[`i',cCD] = . // No CD
					if (mod(`=OMC', 2) == 0) mCI[`i', `=OMC'/2] = (mOS[`i',`=OMC'] - mTSD[`i',`=OMC'])*365.25 // Overwrite CI (even OMCs)
					if (mod(`=OMC', 2) != 0) mCD[`i', (`=OMC'+1)/2] = (mOS[`i',`=OMC'] - mTSD[`i',`=OMC'])*365.25 // Overwrite CD (odd OMCs)
					if (`=OMC' <= 3) mCore[`i', cSCT] = 0 // Overwrite SCT if death before L1E
				}
				if (mMOR[`i',`=OMC'-1] == 0 & mTSD[`i',`=OMC'+1] <= mOS[`i',`=OMC']) mMOR[`i',`=OMC'] = 0 // If TSD < OS, Patient alive
			}
		} 
	}
