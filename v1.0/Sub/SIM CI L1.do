**********
	*SIM CI L1
**********

	*Patient Matrix
		mata: `=m' = mCore , mSU
		mata: r`=m' = rmCore
		mata: c`=m' = cmCore \ cmSU	
		mata: _matrix_list(`=m', r`=m', c`=m')
			
	*Determine outcome - SCT
		scalar b = "bL1_CI_S1"
		forvalues i = 1/`=Obs'{
			*Calculate xb
				*Age
					mata {
						if (`=m'[`i',cSCT] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cAge] * `=b'[1,1])
					}
				*Male 
					mata {
						if (`=m'[`i',cSCT] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cMale] * `=b'[1,2])
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
				*MNT 
					mata {
						if (`=m'[`i',cSCT] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cMNT] * `=b'[1,9])
					}	
				*BCR
					mata {
						if		(`=m'[`i',cSCT] == 1 & `=m'[`i',cBCR] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,11]
						else if	(`=m'[`i',cSCT] == 1 & `=m'[`i',cBCR] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,12] 
						else if	(`=m'[`i',cSCT] == 1 & `=m'[`i',cBCR] == 4) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,13] 
						* else if	(`=m'[`i',cSCT] == 1 & `=m'[`i',cBCR] == 5) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,15] Turned off 13/02/24 as SCT patients cannot have SD
						* else if	(`=m'[`i',cSCT] == 1 & `=m'[`i',cBCR] == 6) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,16] Turned off 13/02/24 as SCT patients cannot have PD
					}
				*cons
					mata {
						if (`=m'[`i',cSCT] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,14]
					}	
			
			*Calculate survival time
				mata: `=m'[`i',`=c'RN] = runiform(1, 1, 0, 1)
				mata {
					if 		(`=m'[`i',cSCT] == 1 & f`=b' == "exponential") 	`=m'[`i',`=c'OC] = (ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])
					else if (`=m'[`i',cSCT] == 1 & f`=b' == "weibull")  	`=m'[`i',`=c'OC] = ((ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])):^(1:/exp(`=b'[1,cols(`=b')]))
					else if (`=m'[`i',cSCT] == 1 & f`=b' == "gompertz")	 	`=m'[`i',`=c'OC] = (ln(1:-((`=b'[1,cols(`=b')]:*(ln(`=m'[`i',`=c'RN]))):/exp(`=m'[`i',`=c'XB])))):/`=b'[1,cols(`=b')]
				}
				
			*Curtail if beyond last observed in the data
				mata {
					if	(`=m'[`i',cSCT] == 1 & `=m'[`i',`=c'OC] > maxL1_CI_S1)	`=m'[`i',`=c'OC] = maxL1_CI_S1
				}
		}

	*Determine outcome - No SCT
		scalar b = "bL1_CI_S0"
		forvalues i = 1/`=Obs'{
			*Calculate xb
				*Age
					mata {
						if (`=m'[`i',cSCT] == 0) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cAge] * `=b'[1,1])
					}
				*Male 
					mata {
						if (`=m'[`i',cSCT] == 0) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cMale] * `=b'[1,2])
					}	
				*ECOGcc
					mata {
						if 		(`=m'[`i',cSCT] == 0 & `=m'[`i',cECOGc] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,4] 
						else if (`=m'[`i',cSCT] == 0 & `=m'[`i',cECOGc] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,5]  
					}
				*ISS 
					mata {
						if 		(`=m'[`i',cSCT] == 0 & `=m'[`i',cISS] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,7] 
						else if (`=m'[`i',cSCT] == 0 & `=m'[`i',cISS] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,8]
					}
				*MNT 
					mata {
						if (`=m'[`i',cSCT] == 0) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cMNT] * `=b'[1,9])
					}	
				*BCR
					mata {
						if		(`=m'[`i',cSCT] == 0 & `=m'[`i',cBCR] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,11] 
						else if	(`=m'[`i',cSCT] == 0 & `=m'[`i',cBCR] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,12]
						else if	(`=m'[`i',cSCT] == 0 & `=m'[`i',cBCR] == 4) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,13] 
						else if	(`=m'[`i',cSCT] == 0 & `=m'[`i',cBCR] == 5) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,14] 
						else if	(`=m'[`i',cSCT] == 0 & `=m'[`i',cBCR] == 6) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,15] 
					}
				*cons
					mata {
						if (`=m'[`i',cSCT] == 0) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,16]
					}
				
			*Calculate survival time
				mata: `=m'[`i',`=c'RN] = runiform(1, 1, 0, 1)
				mata {
					if 		(`=m'[`i',cSCT] == 0 & f`=b' == "exponential") 	`=m'[`i',`=c'OC] = (ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])
					else if (`=m'[`i',cSCT] == 0 & f`=b' == "weibull")  	`=m'[`i',`=c'OC] = ((ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])):^(1:/exp(`=b'[1,cols(`=b')]))
					else if (`=m'[`i',cSCT] == 0 & f`=b' == "gompertz")		`=m'[`i',`=c'OC] = (ln(1:-((`=b'[1,cols(`=b')]:*(ln(`=m'[`i',`=c'RN]))):/exp(`=m'[`i',`=c'XB])))):/`=b'[1,cols(`=b')]
				}
			
			*Curtail if beyond last observed in the data
				mata {
					if	(`=m'[`i',cSCT] == 0 & `=m'[`i',`=c'OC] > maxL1_CI_S0)	`=m'[`i',`=c'OC] = maxL1_CI_S0
				}
		}
				
	*Update outcome matrices (for those with mMOR == 0 only)
		forvalues i = 1/`=Obs'{
			mata {
				if 	(mMOR[`i',`=OMC'-1] == 0) mTNE[`i',`=OMC'] = `=m'[`i',`=c'OC]/365.25
				if 	(mMOR[`i',`=OMC'-1] == 0) mTSD[`i',`=OMC'+1] = mTSD[`i',`=OMC'] + mTNE[`i',`=OMC']
				if 	(mMOR[`i',`=OMC'-1] == 0) mCI[`i',`=LX'+1] = `=m'[`i',`=c'OC]
			}
		}
