**********
	*SIM MORT
**********

	forvalues i = 1/`=Obs' {
		mata {
		*For patients with NFT == 0
			if	(mMOR[`i',`=OMC'-1] == 0 & mTSD[`i',`=OMC'+1] >= mOS[`i',`=OMC']) 	mMOR[`i',`=OMC'] = 1 
			if  (mMOR[`i',`=OMC'-1] == 0 & mTSD[`i',`=OMC'+1] < mOS[`i',`=OMC']) 	mMOR[`i',`=OMC'] = 0
		*For patients with NFT == 1
			*if (mMOR[`i',`=OMC'-1] == 0 & mNFT[`i',`=NFT'] == 1) 				mMOR[`i',`=OMC'] = 1 // Removed cure models
		*For patients who die
			if 	(mMOR[`i',`=OMC'] == 1) {
				*If mOS would take them over Limit
					if	((mAge[`i',2] + mOS[`i',`=OMC']) > `=Limit') {
						mOS[`i',`=OMC'] = `=Limit' - mAge[`i',2] // Set mOS to max if Age > Limit
						if	(mod(`=OMC', 2) == 0)	mCI[`i', `=OMC'/2] = (mOS[`i',`=OMC'] - mTSD[`i',`=OMC'])*365.25 // Overwrite CI (even OMCs)
						if	(mod(`=OMC', 2) != 0)	mCD[`i', (`=OMC'+1)/2] = (mOS[`i',`=OMC'] - mTSD[`i',`=OMC'])*365.25 // Overwrite CD (odd OMCs)
					}
				mOC[`i',2] = mOS[`i',`=OMC']
				mOC[`i',3] = 1
				mTSD[`i',`=OMC'+1] = .
				mTNE[`i',`=OMC'] = .
			}										 
		}
	}
