**********
	*SIM CR L3
**********

	*Patient matrix	
		mata: `=m' = mCore , mML
		mata: r`=m' = rmCore
		mata: c`=m' = cmCore \ cmML
		mata: `=m'[.,`=c'RN] = runiform(`=Obs',1)	
		mata: _matrix_list(`=m', r`=m', c`=m')
		
	*Create nL3 and update nL
		mata: st_matrix("stL3_CR", oL3_CR)
		scalar nL3 = colsof(stL3_CR)
		scalar nL = nL + nL2

	*Calculate outcome for nL3 == 4
		if nL3 == 4 {
			forvalues i=1(1)`=Obs'{
				*Set e(XB1) to 1
					mata: `=m'[`i',`=c'XB1] = 1	
				
				*Calculate e(XB2)
					*Age
						mata: `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + (`=m'[`i',cAge] * `=b'[1,11+`=nL'])
					*Male
						mata: `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + (`=m'[`i',cMale] * `=b'[1,12+`=nL'])
					*ECOG 
						mata {
							if 		(`=m'[`i',cECOGc] == 1) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,14+`=nL'] 
							else if (`=m'[`i',cECOGc] == 2) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,15+`=nL'] 
							else if (`=m'[`i',cECOGc] == 3) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,16+`=nL'] 
						}
					*ISS 
						mata {
							if 		(`=m'[`i',cISS] == 2) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,18+`=nL1'] 
							else if (`=m'[`i',cISS] == 3) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,19+`=nL1']
						}
					*CR_L1
						forvalues j = 1/`=nL1'{
							mata {
								if	(mCR[`i',2] == oL1_CR[1,`j']) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,19+`=nL'+`j']
							}
						}
					*CR_L2
						forvalues j = 1/`=nL2'{
							mata {
								if	(mCR[`i',3] == oL2_CR[1,`j']) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,19+`=nL'+`j'+`=nL1']
							}
						}	
					*cons
						mata: `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,20+(`=nL')*2]
					*Exponent
						mata: `=m'[`i',`=c'XB2] = exp(`=m'[`i',`=c'XB2])
				
				*Calculate e(XB3)
					*Age
						mata: `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + (`=m'[`i',cAge] * `=b'[1,21+(`=nL')*2])
					*Male
						mata: `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + (`=m'[`i',cMale] * `=b'[1,22+(`=nL')*2])
					*ECOG 
						mata {
							if 		(`=m'[`i',cECOGc] == 1) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,24+(`=nL')*2]
							else if (`=m'[`i',cECOGc] == 2) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,25+(`=nL')*2]
							else if (`=m'[`i',cECOGc] == 3) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,26+(`=nL')*2]
						}
					*ISS 
						mata {
							if 		(`=m'[`i',cISS] == 2) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,28+(`=nL')*2]
							else if (`=m'[`i',cISS] == 3) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,29+(`=nL')*2]
						}
					*CR_L1
						forvalues j = 1/`=nL1'{
							mata {
								if	(mCR[`i',2] == oL1_CR[1,`j']) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,29+(`=nL')*2+`j']
							}
						}
					*CR_L2
						forvalues j = 1/`=nL2'{
							mata {
								if	(mCR[`i',3] == oL2_CR[1,`j']) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,29+(`=nL')*2+`j'+`=nL1']
							}
						}	
					*cons
						mata: `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,30+(`=nL')*3]
					*Exponent
						mata: `=m'[`i',`=c'XB3] = exp(`=m'[`i',`=c'XB3])
				
				*Calculate e(XB4)
					*Age
						mata: `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + (`=m'[`i',cAge] * `=b'[1,31+(`=nL')*3])
					*Male
						mata: `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + (`=m'[`i',cMale] * `=b'[1,32+(`=nL')*3])
					*ECOG 
						mata {
							if 		(`=m'[`i',cECOGc] == 1) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,34+(`=nL')*3] 
							else if (`=m'[`i',cECOGc] == 2) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,35+(`=nL')*3] 
							else if (`=m'[`i',cECOGc] == 3) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,36+(`=nL')*3] 
						}
					*ISS 
						mata {
							if 		(`=m'[`i',cISS] == 2) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,38+(`=nL')*3] 
							else if (`=m'[`i',cISS] == 3) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,39+(`=nL')*3]
						}
					*CR_L1
						forvalues j = 1/`=nL1'{
							mata {
								if	(mCR[`i',2] == oL1_CR[1,`j']) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,39+(`=nL')*3+`j']
							}
						}
					*CR_L2
						forvalues j = 1/`=nL2'{
							mata {
								if	(mCR[`i',3] == oL2_CR[1,`j']) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,39+(`=nL')*3+`j'+`=nL1']
							}
						}		
					*cons
						mata: `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,40+(`=nL')*4]
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
		
	*Calculate outcome for nL3 == 3
		if nL3 == 3 {
			forvalues i=1(1)`=Obs'{
				*Set e(XB1) to 1
					mata: `=m'[`i',`=c'XB1] = 1	
				
				*Calculate e(XB2)
					*Age
						mata: `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + (`=m'[`i',cAge] * `=b'[1,11+`=nL'])
					*Male
						mata: `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + (`=m'[`i',cMale] * `=b'[1,12+`=nL'])
					*ECOG 
						mata {
							if 		(`=m'[`i',cECOGc] == 1) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,14+`=nL'] 
							else if (`=m'[`i',cECOGc] == 2) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,15+`=nL'] 
							else if (`=m'[`i',cECOGc] == 3) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,16+`=nL'] 
						}
					*ISS 
						mata {
							if 		(`=m'[`i',cISS] == 2) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,18+`=nL1'] 
							else if (`=m'[`i',cISS] == 3) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,19+`=nL1']
						}
					*CR_L1
						forvalues j = 1/`=nL1'{
							mata {
								if	(mCR[`i',2] == oL1_CR[1,`j']) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,19+`=nL'+`j']
							}
						}
					*CR_L2
						forvalues j = 1/`=nL2'{
							mata {
								if	(mCR[`i',3] == oL2_CR[1,`j']) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,19+`=nL'+`j'+`=nL1']
							}
						}	
					*cons
						mata: `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,20+(`=nL')*2]
					*Exponent
						mata: `=m'[`i',`=c'XB2] = exp(`=m'[`i',`=c'XB2])
				
				*Calculate e(XB3)
					*Age
						mata: `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + (`=m'[`i',cAge] * `=b'[1,21+(`=nL')*2])
					*Male
						mata: `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + (`=m'[`i',cMale] * `=b'[1,22+(`=nL')*2])
					*ECOG 
						mata {
							if 		(`=m'[`i',cECOGc] == 1) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,24+(`=nL')*2]
							else if (`=m'[`i',cECOGc] == 2) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,25+(`=nL')*2]
							else if (`=m'[`i',cECOGc] == 3) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,26+(`=nL')*2]
						}
					*ISS 
						mata {
							if 		(`=m'[`i',cISS] == 2) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,28+(`=nL')*2]
							else if (`=m'[`i',cISS] == 3) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,29+(`=nL')*2]
						}
					*CR_L1
						forvalues j = 1/`=nL1'{
							mata {
								if	(mCR[`i',2] == oL1_CR[1,`j']) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,29+(`=nL')*2+`j']
							}
						}
					*CR_L2
						forvalues j = 1/`=nL2'{
							mata {
								if	(mCR[`i',3] == oL2_CR[1,`j']) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,29+(`=nL')*2+`j'+`=nL1']
							}
						}	
					*cons
						mata: `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,30+(`=nL')*3]
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
		
	*Calculate outcome for nL3 == 2	
		if nL3 == 2 {
			forvalues i=1(1)`=Obs'{
				*Set e(XB1) to 1
					mata: `=m'[`i',`=c'XB1] = 1	
				
				*Calculate e(XB2)
					*Age
						mata: `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + (`=m'[`i',cAge] * `=b'[1,11+`=nL'])
					*Male
						mata: `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + (`=m'[`i',cMale] * `=b'[1,12+`=nL'])
					*ECOG 
						mata {
							if 		(`=m'[`i',cECOGc] == 1) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,14+`=nL'] 
							else if (`=m'[`i',cECOGc] == 2) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,15+`=nL'] 
							else if (`=m'[`i',cECOGc] == 3) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,16+`=nL'] 
						}
					*ISS 
						mata {
							if 		(`=m'[`i',cISS] == 2) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,18+`=nL1'] 
							else if (`=m'[`i',cISS] == 3) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,19+`=nL1']
						}
					*CR_L1
						forvalues j = 1/`=nL1'{
							mata {
								if	(mCR[`i',2] == oL1_CR[1,`j']) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,19+`=nL'+`j']
							}
						}
					*CR_L2
						forvalues j = 1/`=nL2'{
							mata {
								if	(mCR[`i',3] == oL2_CR[1,`j']) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,19+`=nL'+`j'+`=nL1']
							}
						}	
					*cons
						mata: `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,20+(`=nL')*2]
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
