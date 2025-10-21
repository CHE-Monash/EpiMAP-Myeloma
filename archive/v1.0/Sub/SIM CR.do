**********
	*SIM CR
**********

	*Patient matrix	
		mata: `=m' = mCore , mML
		mata: r`=m' = rmCore
		mata: c`=m' = cmCore \ cmML
		mata: `=m'[.,`=c'RN] = runiform(`=Obs',1)	
		mata: _matrix_list(`=m', r`=m', c`=m')
		
		*Determine outcome
			forvalues i=1(1)`=Obs'{
				*Set e(XB1) to 1
					mata: `=m'[`i',`=c'XB1] = 1	
				
				*Calculate e(XB2)
					*Age
						mata: `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + (`=m'[`i',cAge] * `=b'[1,11])
					*Male
						mata: `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + (`=m'[`i',cMale] * `=b'[1,12])
					*ECOG 
						mata {
							if 		(`=m'[`i',cECOGc] == 1) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,14] 
							else if (`=m'[`i',cECOGc] == 2) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,15] 
							else if (`=m'[`i',cECOGc] == 3) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,16] 
						}
					*ISS 
						mata {
							if 		(`=m'[`i',cISS] == 2) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,18] 
							else if (`=m'[`i',cISS] == 3) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,19]
						}
					*cons
						mata: `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,20]
					*Exponent
						mata: `=m'[`i',`=c'XB2] = exp(`=m'[`i',`=c'XB2])
				
				*Calculate e(XB3)
					*Age
						mata: `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + (`=m'[`i',cAge] * `=b'[1,21])
					*Male
						mata: `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + (`=m'[`i',cMale] * `=b'[1,22])
					*ECOG 
						mata {
							if 		(`=m'[`i',cECOGc] == 1) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,24] 
							else if (`=m'[`i',cECOGc] == 2) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,25] 
							else if (`=m'[`i',cECOGc] == 3) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,26] 
						}
					*ISS 
						mata {
							if 		(`=m'[`i',cISS] == 2) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,28] 
							else if (`=m'[`i',cISS] == 3) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,29]
						}
					*cons
						mata: `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,30]
					*Exponent
						mata: `=m'[`i',`=c'XB3] = exp(`=m'[`i',`=c'XB3])
				
				*Calculate e(XB4)
					*Age
						mata: `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + (`=m'[`i',cAge] * `=b'[1,31])
					*Male
						mata: `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + (`=m'[`i',cMale] * `=b'[1,32])
					*ECOG 
						mata {
							if 		(`=m'[`i',cECOGc] == 1) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,34] 
							else if (`=m'[`i',cECOGc] == 2) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,35] 
							else if (`=m'[`i',cECOGc] == 3) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,36] 
						}
					*ISS 
						mata {
							if 		(`=m'[`i',cISS] == 2) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,38] 
							else if (`=m'[`i',cISS] == 3) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,39]
						}
					*cons
						mata: `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,40]
					*Exponent
						mata: `=m'[`i',`=c'XB4] = exp(`=m'[`i',`=c'XB4])
						
				*Calculate probabilities
					mata: `=m'[`i',`=c'PR1] = `=m'[`i',`=c'XB1]/(`=m'[`i',`=c'XB1] + `=m'[`i',`=c'XB2] + `=m'[`i',`=c'XB3] + `=m'[`i',`=c'XB4])
					mata: `=m'[`i',`=c'PR2] = `=m'[`i',`=c'PR1] + `=m'[`i',`=c'XB2]/(`=m'[`i',`=c'XB1] + `=m'[`i',`=c'XB2] + `=m'[`i',`=c'XB3] + `=m'[`i',`=c'XB4])
					mata: `=m'[`i',`=c'PR3] = `=m'[`i',`=c'PR2] + `=m'[`i',`=c'XB3]/(`=m'[`i',`=c'XB1] + `=m'[`i',`=c'XB2] + `=m'[`i',`=c'XB3] + `=m'[`i',`=c'XB4])
					mata: `=m'[`i',`=c'PR4] = `=m'[`i',`=c'PR3] + `=m'[`i',`=c'XB4]/(`=m'[`i',`=c'XB1] + `=m'[`i',`=c'XB2] + `=m'[`i',`=c'XB3] + `=m'[`i',`=c'XB4])
						
				*Compare to RN
					mata {
						if		(mMOR[`i',`=OMC'-1] == 0 & `=m'[`i',`=c'PR1] != . & `=m'[`i',`=c'RN] < `=m'[`i',`=c'PR1]) `=m'[`i',`=c'OC] = `=o'[1,1]
						else if	(mMOR[`i',`=OMC'-1] == 0 & `=m'[`i',`=c'PR2] != . & `=m'[`i',`=c'RN] > `=m'[`i',`=c'PR1]  & `=m'[`i',`=c'RN] < `=m'[`i',`=c'PR2]) `=m'[`i',`=c'OC] = `=o'[1,2]
						else if	(mMOR[`i',`=OMC'-1] == 0 & `=m'[`i',`=c'PR3] != . & `=m'[`i',`=c'RN] > `=m'[`i',`=c'PR2]  & `=m'[`i',`=c'RN] < `=m'[`i',`=c'PR3]) `=m'[`i',`=c'OC] = `=o'[1,3]
						else if	(mMOR[`i',`=OMC'-1] == 0 & `=m'[`i',`=c'PR4] != . & `=m'[`i',`=c'RN] > `=m'[`i',`=c'PR3]  & `=m'[`i',`=c'RN] < `=m'[`i',`=c'PR4]) `=m'[`i',`=c'OC] = `=o'[1,4]
					}
			
				*Update outcome matrices (for those with mMOR == 0 only)
					mata {
						*mCore
							if 	(mMOR[`i',`=OMC'-1] == 0) mCore[`i',cCR] = `=m'[`i',`=c'OC]
						*mCR
							if 	(mMOR[`i',`=OMC'-1] == 0) mCR[`i',`=LX'+2] = `=m'[`i',`=c'OC]
					}
			}
