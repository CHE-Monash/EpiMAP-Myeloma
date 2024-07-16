**********
	*SIM BCR SCT
**********
		
	*Patient Matrix
		mata: `=m' = mCore , mOL
		mata: r`=m' = rmCore
		mata: c`=m' = cmCore \ cmOL
		mata: `=m'[.,`=c'RN] = runiform(`=Obs',1)			
		mata: _matrix_list(`=m', r`=m', c`=m')
		
	*Determine outcome
		forvalues i = 1/`=Obs' {
			*Calculate xb
				*Age
					mata {
						if 	(`=m'[`i',cSCT] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cAge] * `=b'[1,1])
					}
				*Male 
					mata {
						if 	(`=m'[`i',cSCT] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cMale] * `=b'[1,2])
					} 
				*ECOGcc 
					mata {
						if 		(`=m'[`i',cSCT] == 1 & `=m'[`i',cECOGc] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,4] 
						else if (`=m'[`i',cSCT] == 1 & `=m'[`i',cECOGc] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,5]  
					}
				*ISS 
					mata {
						if 		(`=m'[`i',cSCT] == 1 & `=m'[`i',cISS] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,7] 
						else if (`=m'[`i',cSCT] == 1 & `=m'[`i',cISS] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,8]
					}
				*pBCR
					mata {
						if		(`=m'[`i',cSCT] == 1 & `=m'[`i',cBCR] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,10] 
						else if	(`=m'[`i',cSCT] == 1 & `=m'[`i',cBCR] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,11] 
						else if	(`=m'[`i',cSCT] == 1 & `=m'[`i',cBCR] == 4) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,12]
						else if	(`=m'[`i',cSCT] == 1 & `=m'[`i',cBCR] == 5) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,13] 
					}
				
			*Probabilities	
				mata { 
					if (`=m'[`i',cSCT] == 1) `=m'[`i',`=c'PR1] = 1/(1+exp(`=m'[`i',`=c'XB] - `=b'[1,cols(`=b')-2]))
					if (`=m'[`i',cSCT] == 1) `=m'[`i',`=c'PR2] = 1/(1+exp(`=m'[`i',`=c'XB] - `=b'[1,cols(`=b')-1]))	
					if (`=m'[`i',cSCT] == 1) `=m'[`i',`=c'PR3] = 1/(1+exp(`=m'[`i',`=c'XB] - `=b'[1,cols(`=b')]))
					if (`=m'[`i',cSCT] == 1) `=m'[`i',`=c'PR4] = 1
				} 
			
			*BCR outcome for SCT patients
				mata {
					if		(`=m'[`i',cSCT] == 1 & (`=m'[`i',`=c'RN] < `=m'[`i',`=c'PR1])) `=m'[`i',`=c'OC] = 1
					else if (`=m'[`i',cSCT] == 1 & (`=m'[`i',`=c'RN] > `=m'[`i',`=c'PR1]) & (`=m'[`i',`=c'RN] < `=m'[`i',`=c'PR2])) `=m'[`i',`=c'OC] = 2
					else if (`=m'[`i',cSCT] == 1 & (`=m'[`i',`=c'RN] > `=m'[`i',`=c'PR2]) & (`=m'[`i',`=c'RN] < `=m'[`i',`=c'PR3])) `=m'[`i',`=c'OC] = 3
					else if (`=m'[`i',cSCT] == 1 & (`=m'[`i',`=c'RN] > `=m'[`i',`=c'PR3])) `=m'[`i',`=c'OC] = 4

				}		
		
			*Update outcome matrices
				mata {
					if 	(`=m'[`i',cSCT] == 1 & mMOR[`i',`=OMC'-1] == 0) mCore[`i',cBCR] = `=m'[`i',`=c'OC]
				*Set outcome to missing if patients die (otherwise they get carried through to Results)	- Turned off on 23/01/24 don't think needed anymore - tab BCR_SCT MOR_L1S shows no SCT == 1
					*if	(`=m'[`i',cSCT] == 1 & mMOR[`i',`=OMC'-1] == 1) `=m'[`i',`=c'OC] = .
				}
		}
