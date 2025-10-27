**********
	*SIM CR L1
**********

	*Patient matrix	
		mata: `=m' = mCore , mML
		mata: r`=m' = rmCore
		mata: c`=m' = cmCore \ cmML
		mata: `=m'[.,`=c'RN] = runiform(`=Obs',1)	
		mata: _matrix_list(`=m', r`=m', c`=m')
		
	*Create nL1
		mata: st_matrix("stL1_CR", oL1_CR)
		scalar nL1 = colsof(stL1_CR)	
	
	*Calculate outcome for nL1 == 4
		if nL1 == 4 {
			forvalues i=1/`=Obs'{
				*Set e(XB1) to 1
					mata: `=m'[`i',`=c'XB1] = 1	
				
				*Calculate e(XB2)
					*Age
						mata: `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + (`=m'[`i',cAge] * `=b'[1,10])
					*Male
						mata: `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + (`=m'[`i',cMale] * `=b'[1,11])
					*ECOG 
						mata {
							if 		(`=m'[`i',cECOGc] == 1) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,13] 
							else if (`=m'[`i',cECOGc] == 2) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,14] 
						}
					*ISS 
						mata {
							if 		(`=m'[`i',cISS] == 2) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,16] 
							else if (`=m'[`i',cISS] == 3) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,17]
						}
					*cons
						mata: `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,18]
					*Exponent
						mata: `=m'[`i',`=c'XB2] = exp(`=m'[`i',`=c'XB2])
				
				*Calculate e(XB3)
					*Age
						mata: `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + (`=m'[`i',cAge] * `=b'[1,19])
					*Male
						mata: `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + (`=m'[`i',cMale] * `=b'[1,20])
					*ECOG 
						mata {
							if 		(`=m'[`i',cECOGc] == 1) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,22] 
							else if (`=m'[`i',cECOGc] == 2) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,23] 
						}
					*ISS 
						mata {
							if 		(`=m'[`i',cISS] == 2) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,25] 
							else if (`=m'[`i',cISS] == 3) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,26]
						}
					*cons
						mata: `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,27]
					*Exponent
						mata: `=m'[`i',`=c'XB3] = exp(`=m'[`i',`=c'XB3])
				
				*Calculate e(XB4)
					*Age
						mata: `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + (`=m'[`i',cAge] * `=b'[1,28])
					*Male
						mata: `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + (`=m'[`i',cMale] * `=b'[1,29])
					*ECOG 
						mata {
							if 		(`=m'[`i',cECOGc] == 1) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,31] 
							else if (`=m'[`i',cECOGc] == 2) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,32] 
						}
					*ISS 
						mata {
							if 		(`=m'[`i',cISS] == 2) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,34] 
							else if (`=m'[`i',cISS] == 3) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,35]
						}
					*cons
						mata: `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,36]
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
							if 	(mMOR[`i',`=OMC'-1] != 0) mCore[`i',cCR] = .
						*mCR
							if 	(mMOR[`i',`=OMC'-1] == 0) mCR[`i',`=LX'+2] = `=m'[`i',`=c'OC]
						}
			}
		}
	
	*Calculate outcome for nL1 == 3
		if nL1 == 3 {
			forvalues i=1/`=Obs'{
				*Set e(XB1) to 1
					mata: `=m'[`i',`=c'XB1] = 1	
				
				*Calculate e(XB2)
					*Age
						mata: `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + (`=m'[`i',cAge] * `=b'[1,10])
					*Male
						mata: `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + (`=m'[`i',cMale] * `=b'[1,11])
					*ECOG 
						mata {
							if 		(`=m'[`i',cECOGc] == 1) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,13] 
							else if (`=m'[`i',cECOGc] == 2) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,14] 
						}
					*ISS 
						mata {
							if 		(`=m'[`i',cISS] == 2) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,16] 
							else if (`=m'[`i',cISS] == 3) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,17]
						}
					*cons
						mata: `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,18]
					*Exponent
						mata: `=m'[`i',`=c'XB2] = exp(`=m'[`i',`=c'XB2])
				
				*Calculate e(XB3)
					*Age
						mata: `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + (`=m'[`i',cAge] * `=b'[1,19])
					*Male
						mata: `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + (`=m'[`i',cMale] * `=b'[1,20])
					*ECOG 
						mata {
							if 		(`=m'[`i',cECOGc] == 1) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,22] 
							else if (`=m'[`i',cECOGc] == 2) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,23] 
						}
					*ISS 
						mata {
							if 		(`=m'[`i',cISS] == 2) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,25] 
							else if (`=m'[`i',cISS] == 3) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,26]
						}
					*cons
						mata: `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,27]
					*Exponent
						mata: `=m'[`i',`=c'XB3] = exp(`=m'[`i',`=c'XB3])
									
				*Calculate probabilities
					mata: `=m'[`i',`=c'PR1] = `=m'[`i',`=c'XB1]/(`=m'[`i',`=c'XB1] + `=m'[`i',`=c'XB2] + `=m'[`i',`=c'XB3] + `=m'[`i',`=c'XB4])
					mata: `=m'[`i',`=c'PR2] = `=m'[`i',`=c'PR1] + `=m'[`i',`=c'XB2]/(`=m'[`i',`=c'XB1] + `=m'[`i',`=c'XB2] + `=m'[`i',`=c'XB3] + `=m'[`i',`=c'XB4])
					mata: `=m'[`i',`=c'PR3] = `=m'[`i',`=c'PR2] + `=m'[`i',`=c'XB3]/(`=m'[`i',`=c'XB1] + `=m'[`i',`=c'XB2] + `=m'[`i',`=c'XB3] + `=m'[`i',`=c'XB4])
							
				*Compare to RN
					mata {
						if		(mMOR[`i',`=OMC'-1] == 0 & `=m'[`i',`=c'PR1] != . & `=m'[`i',`=c'RN] < `=m'[`i',`=c'PR1]) `=m'[`i',`=c'OC] = `=o'[1,1]
						else if	(mMOR[`i',`=OMC'-1] == 0 & `=m'[`i',`=c'PR2] != . & `=m'[`i',`=c'RN] > `=m'[`i',`=c'PR1]  & `=m'[`i',`=c'RN] < `=m'[`i',`=c'PR2]) `=m'[`i',`=c'OC] = `=o'[1,2]
						else if	(mMOR[`i',`=OMC'-1] == 0 & `=m'[`i',`=c'PR3] != . & `=m'[`i',`=c'RN] > `=m'[`i',`=c'PR2]  & `=m'[`i',`=c'RN] < `=m'[`i',`=c'PR3]) `=m'[`i',`=c'OC] = `=o'[1,3]
					}	
					
				*Update outcome matrices (for those with mMOR == 0 only)
					mata {
						*mCore
							if 	(mMOR[`i',`=OMC'-1] == 0) mCore[`i',cCR] = `=m'[`i',`=c'OC]
							if 	(mMOR[`i',`=OMC'-1] != 0) mCore[`i',cCR] = .
						*mCR
							if 	(mMOR[`i',`=OMC'-1] == 0) mCR[`i',`=LX'+2] = `=m'[`i',`=c'OC]
					}
			}
		}
	
	*Calculate outcome for nL1 == 2
		if nL1 == 2 {
			forvalues i=1/`=Obs'{
				*Set e(XB1) to 1
					mata: `=m'[`i',`=c'XB1] = 1	
				
								*Calculate e(XB2)
					*Age
						mata: `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + (`=m'[`i',cAge] * `=b'[1,10])
					*Male
						mata: `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + (`=m'[`i',cMale] * `=b'[1,11])
					*ECOG 
						mata {
							if 		(`=m'[`i',cECOGc] == 1) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,13] 
							else if (`=m'[`i',cECOGc] == 2) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,14] 
						}
					*ISS 
						mata {
							if 		(`=m'[`i',cISS] == 2) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,16] 
							else if (`=m'[`i',cISS] == 3) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,17]
						}
					*cons
						mata: `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,18]
					*Exponent
						mata: `=m'[`i',`=c'XB2] = exp(`=m'[`i',`=c'XB2])
									
				*Calculate probabilities
					mata: `=m'[`i',`=c'PR1] = `=m'[`i',`=c'XB1]/(`=m'[`i',`=c'XB1] + `=m'[`i',`=c'XB2] + `=m'[`i',`=c'XB3] + `=m'[`i',`=c'XB4])
					mata: `=m'[`i',`=c'PR2] = `=m'[`i',`=c'PR1] + `=m'[`i',`=c'XB2]/(`=m'[`i',`=c'XB1] + `=m'[`i',`=c'XB2] + `=m'[`i',`=c'XB3] + `=m'[`i',`=c'XB4])
							
				*Compare to RN
					mata {
						if		(mMOR[`i',`=OMC'-1] == 0 & `=m'[`i',`=c'PR1] != . & `=m'[`i',`=c'RN] < `=m'[`i',`=c'PR1]) `=m'[`i',`=c'OC] = `=o'[1,1]
						else if	(mMOR[`i',`=OMC'-1] == 0 & `=m'[`i',`=c'PR2] != . & `=m'[`i',`=c'RN] > `=m'[`i',`=c'PR1]  & `=m'[`i',`=c'RN] < `=m'[`i',`=c'PR2]) `=m'[`i',`=c'OC] = `=o'[1,2]
					}
					
				*Update outcome matrices (for those with mMOR == 0 only)
					mata {
						*mCore
							if 	(mMOR[`i',`=OMC'-1] == 0) mCore[`i',cCR] = `=m'[`i',`=c'OC]
							if 	(mMOR[`i',`=OMC'-1] != 0) mCore[`i',cCR] = .
						*mCR
							if 	(mMOR[`i',`=OMC'-1] == 0) mCR[`i',`=LX'+2] = `=m'[`i',`=c'OC]
						}
			}
		}
