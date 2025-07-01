**********
	*SIM BCR
**********

	*Patient Matrix
		mata: `=m' = mCore , mOL
		mata: r`=m' = rmCore
		mata: c`=m' = cmCore \ cmOL
		mata: `=m'[.,`=c'RN] = runiform(`=Obs',1)			
		mata: _matrix_list(`=m', r`=m', c`=m')
			
	*Determine outcome
		forvalues i = 1/`=Obs'{
			*Calculate xb
				*Age
					mata: `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cAge] * `=b'[1,1])
				*Male 
					mata: `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cMale] * `=b'[1,2])
				*ECOGcc
					mata {
						if 		(`=m'[`i',cECOGc] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,4] 
						else if (`=m'[`i',cECOGc] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,5]  
					}
				*ISS 
					mata {
						if 		(`=m'[`i',cISS] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,7] 
						else if (`=m'[`i',cISS] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,8]
					}
				*CR
					mata {
						if (cols(`=o') >= 2) {
							if	(`=m'[`i',cCR] == `=o'[1,2]) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,10]
						}
						if (cols(`=o') >= 3) {
							if	(`=m'[`i',cCR] == `=o'[1,3]) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,11]
						}
						if (cols(`=o') >= 4) {
							if	(`=m'[`i',cCR] == `=o'[1,4]) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,12]
						}
					} 
				*CD
					*mata: `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cCD] * `=b'[1,9 + cols(`=o')])
				*pBCR
					mata {
						if		(`=m'[`i',cBCR] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,10 + cols(`=o')] 
						else if	(`=m'[`i',cBCR] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,11 + cols(`=o')] 
						else if	(`=m'[`i',cBCR] == 4) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,12 + cols(`=o')]
						else if	(`=m'[`i',cBCR] == 5) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,13 + cols(`=o')] 
						else if	(`=m'[`i',cBCR] == 6) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,14 + cols(`=o')] 
					}
			
			*Probabilities
				mata: `=m'[`i',`=c'PR1] = 1/(1+exp(`=m'[`i',`=c'XB] - `=b'[1,15 + cols(`=o')]))
				mata: `=m'[`i',`=c'PR2] = 1/(1+exp(`=m'[`i',`=c'XB] - `=b'[1,16 + cols(`=o')]))
				mata: `=m'[`i',`=c'PR3] = 1/(1+exp(`=m'[`i',`=c'XB] - `=b'[1,17 + cols(`=o')]))
				mata: `=m'[`i',`=c'PR4] = 1/(1+exp(`=m'[`i',`=c'XB] - `=b'[1,18 + cols(`=o')]))
				mata: `=m'[`i',`=c'PR5] = 1/(1+exp(`=m'[`i',`=c'XB] - `=b'[1,19 + cols(`=o')]))
				mata {
					if (`=m'[`i',`=c'PR5] != .) `=m'[`i',`=c'PR6] = 1
				}
			
			*BCR outcome
				mata {
					if		((`=m'[`i',`=c'RN] < `=m'[`i',`=c'PR1]) & (`=m'[`i',`=c'PR6] == 1)) `=m'[`i',`=c'OC] = 1
					else if ((`=m'[`i',`=c'RN] > `=m'[`i',`=c'PR1]) & (`=m'[`i',`=c'RN] < `=m'[`i',`=c'PR2])) `=m'[`i',`=c'OC] = 2
					else if ((`=m'[`i',`=c'RN] > `=m'[`i',`=c'PR2]) & (`=m'[`i',`=c'RN] < `=m'[`i',`=c'PR3])) `=m'[`i',`=c'OC] = 3
					else if ((`=m'[`i',`=c'RN] > `=m'[`i',`=c'PR3]) & (`=m'[`i',`=c'RN] < `=m'[`i',`=c'PR4])) `=m'[`i',`=c'OC] = 4
					else if ((`=m'[`i',`=c'RN] > `=m'[`i',`=c'PR4]) & (`=m'[`i',`=c'RN] < `=m'[`i',`=c'PR5])) `=m'[`i',`=c'OC] = 5
					else if (`=m'[`i',`=c'RN] > `=m'[`i',`=c'PR5]) `=m'[`i',`=c'OC] = 6
				}
					
			*Update outcome matrices
				mata {
					if 	(mMOR[`i',`=OMC'-1] == 0) mCore[`i',cBCR] = `=m'[`i',`=c'OC]					
					if 	(mMOR[`i',`=OMC'-1] == 0) mBCR[`i',`=LX'+1] = `=m'[`i',`=c'OC]
				}			
		}
