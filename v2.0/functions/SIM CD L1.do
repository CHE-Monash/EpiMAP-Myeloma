**********
	*SIM CD L1
**********

	*Patient matrix	
		mata: `=m' = mCore , mSU
		mata: r`=m' = rmCore
		mata: c`=m' = cmCore \ cmSU
		mata: `=m'[.,`=c'RN] = runiform(`=Obs',1)	
		mata: _matrix_list(`=m', r`=m', c`=m')
	
	*Determine outcome - Fixed + ASCT - Spline 1
		scalar b = "bL1F1_CD_S1"
		forvalues i = 1/`=Obs' {
			mata {
				if (mMOR[`i',`=OMC'-1] == 0 & mState[`i',2] <= `=OMC') { // Alive & State filters
					if (`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1) { // Fixed + SCT filter
					
					*Calculate XB
						*Age
							`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cAge] * `=b'[1,1])
						*Age2
							`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cAge]^2 * `=b'[1,2])							
						*Male 
							`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cMale] * `=b'[1,3])
						*ECOG 
							if (`=m'[`i',cECOG] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,5] 
							if (`=m'[`i',cECOG] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,6] 
						*RISS 
							if (`=m'[`i',cRISS] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,8] 
							if (`=m'[`i',cRISS] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,9]
						*CR 
							if (cols(`=o') >= 2) {
								if	(`=m'[`i',cCR] == `=o'[1,2]) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,11]
							}
							if (cols(`=o') >= 3) {
								if	(`=m'[`i',cCR] == `=o'[1,3]) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,12]
							}
						*cons
							`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,cols(`=b')-1]
				
					*Calculate survival time
						if (f`=b' == "ereg") 		`=m'[`i',`=c'OC] = (ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])
						if (f`=b' == "weibull")  	`=m'[`i',`=c'OC] = ((ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])):^(1:/exp(`=b'[1,cols(`=b')]))
						if (f`=b' == "gompertz")	`=m'[`i',`=c'OC] = (ln(1-((`=b'[1,cols(`=b')]:*(ln(`=m'[`i',`=c'RN]))):/exp(`=m'[`i',`=c'XB])))):/`=b'[1,cols(`=b')]
				
					*Reset XB
						mata: `=m'[`i',`=c'XB] = 0	
					}
				}
			}
		}
		
	*Determine outcome - Fixed + ASCT - Spline 2
		scalar b = "bL1F1_CD_S2"
		forvalues i = 1/`=Obs' {
			mata {
				if (mMOR[`i',`=OMC'-1] == 0 & mState[`i',2] <= `=OMC') { // Alive & State filters
					if (`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1) { // Fixed + SCT filter
					
					*Calculate XB
						*Age
							`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cAge] * `=b'[1,1])
						*Age2
							`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cAge]^2 * `=b'[1,2])	
						*Male 
							`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cMale] * `=b'[1,3])
						*ECOG 
							if (`=m'[`i',cECOG] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,5] 
							if (`=m'[`i',cECOG] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,6]
						*RISS
							if (`=m'[`i',cRISS] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,8] 
							if (`=m'[`i',cRISS] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,9]
						*CR
							if (cols(`=o') >= 2) {
								if	(`=m'[`i',cCR] == `=o'[1,2]) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,11]
							}
							if (cols(`=o') >= 3) {
								if	(`=m'[`i',cCR] == `=o'[1,3]) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,12]
							}
						*cons
							`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,cols(`=b')-1]
							
					*Draw RN conditional on survival to Spline 1 cut off
						if (`=m'[`i',`=c'XB] != 0) 	`=m'[`i',`=c'RN] = runiform(1, 1, 0, exp(-(exp(`=m'[`i',`=c'XB]))*(L1F1_CD_C1^exp(`=b'[1,cols(`=b')]))))
						
					*Recalculate survival time for those with OC beyond splice cut off 1 only
						if (f`=b' == "ereg" & `=m'[`i',`=c'OC] > L1F1_CD_C1) 		`=m'[`i',`=c'OC] = (ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])
						if (f`=b' == "weibull" & `=m'[`i',`=c'OC] > L1F1_CD_C1)  	`=m'[`i',`=c'OC] = ((ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])):^(1:/exp(`=b'[1,cols(`=b')]))
						if (f`=b' == "gompertz" & `=m'[`i',`=c'OC] > L1F1_CD_C1)	`=m'[`i',`=c'OC] = (ln(1-((`=b'[1,cols(`=b')]:*(ln(`=m'[`i',`=c'RN]))):/exp(`=m'[`i',`=c'XB])))):/`=b'[1,cols(`=b')]
					
					*Reset XB
						mata: `=m'[`i',`=c'XB] = 0		
					}
				}
			}
		}
		
	*Determine outocme - Fixed + ASCT - Spline 3
		scalar b = "bL1F1_CD_S3"
		forvalues i = 1/`=Obs' {
			mata {
				if (mMOR[`i',`=OMC'-1] == 0 & mState[`i',2] <= `=OMC') { // Alive & State filters
					if (`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1) { // Fixed + SCT filter
					
					*Calculate XB
						*Age
							`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cAge] * `=b'[1,1])
						*Age2
							`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cAge]^2 * `=b'[1,2])	
						*Male 
							`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cMale] * `=b'[1,3])
						*ECOG
							if (`=m'[`i',cECOG] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,5] 
							if (`=m'[`i',cECOG] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,6]
						*RISS
							if (`=m'[`i',cRISS] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,8] 
							if (`=m'[`i',cRISS] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,9]
						*CR
							if (cols(`=o') >= 2) {
								if	(`=m'[`i',cCR] == `=o'[1,2]) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,11]
							}
							if (cols(`=o') >= 3) {
								if	(`=m'[`i',cCR] == `=o'[1,3]) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,12]
							}
						*cons
							`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,cols(`=b')-1]
								
					*Draw RN conditional on survival to Spline 2 cut off
						if (`=m'[`i',`=c'XB] != 0) 	`=m'[`i',`=c'RN] = runiform(1, 1, 0, exp(-(exp(`=m'[`i',`=c'XB]))*(L1F1_CD_C2^exp(`=b'[1,cols(`=b')]))))

						*Recalculate survival time for those with OC beyond splice cut off 2 only
						if (f`=b' == "ereg" & `=m'[`i',`=c'OC] > L1F1_CD_C2) 		`=m'[`i',`=c'OC] = (ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])
						if (f`=b' == "weibull" & `=m'[`i',`=c'OC] > L1F1_CD_C2)  	`=m'[`i',`=c'OC] = ((ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])):^(1:/exp(`=b'[1,cols(`=b')]))
						if (f`=b' == "gompertz" & `=m'[`i',`=c'OC] > L1F1_CD_C2)	`=m'[`i',`=c'OC] = (ln(1-((`=b'[1,cols(`=b')]:*(ln(`=m'[`i',`=c'RN]))):/exp(`=m'[`i',`=c'XB])))):/`=b'[1,cols(`=b')]
					}
				}
			}
		}
		
	*Determine outcome - Fixed - ASCT
		scalar b = "bL1F0_CD"
		forvalues i = 1/`=Obs' {
			mata {
				if (mMOR[`i',`=OMC'-1] == 0 & mState[`i',2] <= `=OMC') { // Alive & State filters
					if (`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 0) { // Fixed - SCT filter
					
					*Calculate XB
						*Age
							`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cAge] * `=b'[1,1])
						*Age2
							`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cAge]^2 * `=b'[1,2])								
						*Male
							`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cMale] * `=b'[1,3])
						*ECOG 
							if (`=m'[`i',cECOG] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,5] 
							if (`=m'[`i',cECOG] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,6] 
						*RISS
							if (`=m'[`i',cRISS] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,8] 
							if (`=m'[`i',cRISS] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,9]
						*CR
							if (cols(`=o') >= 2) {
								if	(`=m'[`i',cCR] == `=o'[1,2]) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,11]
							}
							if (cols(`=o') >= 3) {
								if	(`=m'[`i',cCR] == `=o'[1,3]) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,12]
							}
						*cons
							`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,cols(`=b')-1]
			
					*Calculate survival time
						if (f`=b' == "ereg") 		`=m'[`i',`=c'OC] = (ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])
						if (f`=b' == "weibull")  	`=m'[`i',`=c'OC] = ((ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])):^(1:/exp(`=b'[1,cols(`=b')]))
						if (f`=b' == "gompertz")	`=m'[`i',`=c'OC] = (ln(1-((`=b'[1,cols(`=b')]:*(ln(`=m'[`i',`=c'RN]))):/exp(`=m'[`i',`=c'XB])))):/`=b'[1,cols(`=b')]
					}
				}
			}
		}
			
	*Continuous Therapy
		mata:
			CR_L1_7 = 0  // Check if oL1_CR includes 7
			for (j=1; j<=cols(oL1_CR); j++) {
				if (oL1_CR[1, j] == 7) {
					CR_L1_7 = 1
				}
			}
			st_numscalar("CR_L1_7", CR_L1_7)
		end
		if(`=CR_L1_7' == 1){
			scalar b = "bL1C_CD"
			forvalues i = 1/`=Obs' {
				mata {
					if (mMOR[`i',`=OMC'-1] == 0 & mState[`i',2] <= `=OMC') { // Alive & State filters
						if (`=m'[`i',cCR] == 7) { // Continuous Therapy filter	
						
						*Calculate XB
							*Age
								`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cAge] * `=b'[1,1])
							*Age2
								`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cAge]^2 * `=b'[1,2])								
							*Male 
								`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cMale] * `=b'[1,3])
							*ECOG 
								if (`=m'[`i',cECOG] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,5] 
								if (`=m'[`i',cECOG] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,6] 
							*RISS 
								if (`=m'[`i',cRISS] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,8] 
								if (`=m'[`i',cRISS] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,9]
							*cons
								`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,cols(`=b')-1]
				
						*Calculate survival time
							if (f`=b' == "ereg") 		`=m'[`i',`=c'OC] = (ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])
							if (f`=b' == "weibull")  	`=m'[`i',`=c'OC] = ((ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])):^(1:/exp(`=b'[1,cols(`=b')]))
							if (f`=b' == "gompertz")	`=m'[`i',`=c'OC] = (ln(1-((`=b'[1,cols(`=b')]:*(ln(`=m'[`i',`=c'RN]))):/exp(`=m'[`i',`=c'XB])))):/`=b'[1,cols(`=b')]
						
						*Curtail if outcome beyond last observed in the data
							if (`=m'[`i',`=c'OC] != . & `=m'[`i',`=c'OC] > maxL1C_CD)	`=m'[`i',`=c'OC] = maxL1C_CD
						}
					}
				}
			}	
		}

	*Update outcome matrices
		forvalues i = 1/`=Obs' {
			mata {
				if (mState[`i',2] > `=OMC') `=m'[`i',`=c'OC] = mTNE[`i',`=OMC']*365.25 // Grab prevalent patient data
				if (mMOR[`i',`=OMC'-1] == 0) {
					mTNE[`i',`=OMC'] = `=m'[`i',`=c'OC]/365.25
					mTSD[`i',`=OMC'+1] = mTSD[`i',`=OMC'] + mTNE[`i',`=OMC']
					mCore[`i',cCD] = `=m'[`i',`=c'OC]
					mCD[`i',`=LX'+2] = `=m'[`i',`=c'OC]	
				}
			}
		}
		


