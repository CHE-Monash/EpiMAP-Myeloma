**********
	*SIM SCT DN
**********

	*REPLACE CORE2 WITH CORE ONCE UPDATE COMPLETE

	*Patient matrix
		mata: p`=m' = mCore2 , mCom , mCons
		mata: rp`=m' = rmCore2
		mata: cp`=m' = cmCore2 \ cmCom \ cmCons
		mata: _matrix_list(p`=m', rp`=m', cp`=m')
		
	*Outcome matrix
		mata: o`=m' = mLO
		mata: ro`=m' = rmCore2
		mata: co`=m' = cmLO
		
		mata: o`=m'[.,1] = p`=m' * `=b'' // XB
		mata: o`=m'[.,2] = 1:/(1:+exp(-o`=m'[.,1]))	// PR		
		mata: o`=m'[.,3] = runiform(`=Obs',1) // RN

		*Loop to determine outcome
			forval i = 1/`=Obs'{
				mata {
					if (mState[`i',2] <= `=OMC') { // State filter

					*Compare PR with RN 
						if (o`=m'[`i',2] > o`=m'[`i',3]) o`=m'[`i',4] = 1
						if (o`=m'[`i',2] < o`=m'[`i',3]) o`=m'[`i',4] = 0
							
					*Update mSCT
						mSCT[`i',1] = o`=m'[`i',4]
						
					*Update mCore REMOVE ONCE UPDATE COMPLETE
						mCore[`i',cSCT] = o`=m'[`i',4]
					}	
				}		
			}
		
		mata: _matrix_list(o`=m', ro`=m', co`=m')

