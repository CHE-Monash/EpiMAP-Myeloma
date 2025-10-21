**********
	*SIM AGE
**********

	forvalues i = 1/`=Obs' {
		mata {
		*For patients alive
			if (mMOR[`i',`=OMC'-1] == 0) {
				mAge[`i',`=OMC'] = mAge[`i',`=OMC'-1] + mTNE[`i',`=OMC'-1] // Update mAge
				mCore[`i',cAge] = mAge[`i',`=OMC'] // Update mCore	
			*For patients over Limit
				if	(mAge[`i',`=OMC'] > `=Limit' & mAge[`i',`=OMC'] != .) {
					mAge[`i',`=OMC'] = `=Limit' // Set Age to Limit
					mMOR[`i',`=OMC'-1] = 1 // Overwrite mMOR
					mCore[`i',cAge] = . // Overwrite mCore		
					mOC[`i',2] = `=Limit' - mAge[`i',2] // Set OC_TIME so patient dies at Limit
					mOC[`i',3] = 1
					*Code below similar to SIM MORT but OMC has now moved on by 1
					if	(mod((`=OMC'-1), 2) == 0) mCI[`i', (`=OMC'-1)/2] = (mOC[`i',2] - mTSD[`i',`=OMC'-1])*365.25 // Overwrite CI (even OMCs)
					if	(mod((`=OMC'-1), 2) != 0) mCD[`i', `=OMC'/2] = (mOC[`i',2] - mTSD[`i',`=OMC'-1])*365.25 // Overwrite CD (odd OMCs)
					mTSD[`i',`=OMC'] = .
					mTNE[`i',`=OMC'-1] = .
				}
			}
				
		} 
	}
