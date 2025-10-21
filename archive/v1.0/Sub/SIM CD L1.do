**********
	*SIM CD L1
**********

	*Patient matrix	
		mata: `=m' = mCore , mSU
		mata: r`=m' = rmCore
		mata: c`=m' = cmCore \ cmSU
		mata: `=m'[.,`=c'RN] = runiform(`=Obs',1)	
		mata: _matrix_list(`=m', r`=m', c`=m')
	
	*Determine outcome - Fixed + SCT - Spline 1
		scalar b = "bL1F1_CD_S1"
		forvalues i = 1/`=Obs'{
			*Calculate xb
				*Age
					mata {
						if		(`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cAge] * `=b'[1,1])
					}	
				*Male 
					mata {
						if 		(`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cMale] * `=b'[1,2])
					} 
				*ECOG 
					mata {
						if 		(`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1 & `=m'[`i',cECOGc] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,4] 
						else if (`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1 & `=m'[`i',cECOGc] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,5] 
					}
				*ISS 
					mata {
						if 		(`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1 & `=m'[`i',cISS] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,7] 
						else if (`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1 & `=m'[`i',cISS] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,8]
					}
				*CR
					mata {
						if		(`=m'[`i',cCR] == 4 & `=m'[`i',cSCT] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,10] 
						else if	(`=m'[`i',cCR] == 31 & `=m'[`i',cSCT] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,11] 
					} 
				*cons
					mata {
						if		(`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,12]
					} 	
		
			*Calculate survival time
				mata {
					if 		(`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1 & f`=b' == "ereg") 		`=m'[`i',`=c'OC] = (ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])
					else if (`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1 & f`=b' == "weibull")  	`=m'[`i',`=c'OC] = ((ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])):^(1:/exp(`=b'[1,cols(`=b')]))
					else if (`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1 & f`=b' == "gompertz")	`=m'[`i',`=c'OC] = (ln(1-((`=b'[1,cols(`=b')]:*(ln(`=m'[`i',`=c'RN]))):/exp(`=m'[`i',`=c'XB])))):/`=b'[1,cols(`=b')]
				}
		}
			
		*Reset XB
			forvalues i = 1/`=Obs'{
				mata: `=m'[`i',`=c'XB] = 0
			}
		
	*Determine outcome - Fixed + SCT - Spline 2	
		scalar b = "bL1F1_CD_S2"
		forvalues i = 1/`=Obs'{
			*Calculate xb
				*Age
					mata {
						if		(`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cAge] * `=b'[1,1])
					}	
				*Male 
					mata {
						if 		(`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cMale] * `=b'[1,2])
					} 
				*ECOG 
					mata {
						if 		(`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1 & `=m'[`i',cECOGc] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,4] 
						else if (`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1 & `=m'[`i',cECOGc] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,5]  
					}
				*ISS 
					mata {
						if 		(`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1 & `=m'[`i',cISS] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,7] 
						else if (`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1 & `=m'[`i',cISS] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,8]
					}
				*CR
					mata {
						if		(`=m'[`i',cCR] == 4 & `=m'[`i',cSCT] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,10] 
						else if	(`=m'[`i',cCR] == 31 & `=m'[`i',cSCT] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,11] 
					} 
				*cons
					mata {
						if		(`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,12]
					} 
					
			*Draw RN conditional on survival to splice cut off 1
				mata {
					if (`=m'[`i',`=c'XB] != 0) 	`=m'[`i',`=c'RN] = runiform(1, 1, 0, exp(-(exp(`=m'[`i',`=c'XB]))*(L1F1_CD_C1^exp(`=b'[1,13]))))
				} 

			*Recalculate survival time for those with OC beyond splice cut off 1 only
				mata {
					if 		(`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1 & f`=b' == "ereg" & `=m'[`i',`=c'OC] > L1F1_CD_C1) 		`=m'[`i',`=c'OC] = (ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])
					else if (`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1 & f`=b' == "weibull" & `=m'[`i',`=c'OC] > L1F1_CD_C1)  	`=m'[`i',`=c'OC] = ((ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])):^(1:/exp(`=b'[1,cols(`=b')]))
					else if (`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1 & f`=b' == "gompertz" & `=m'[`i',`=c'OC] > L1F1_CD_C1)	`=m'[`i',`=c'OC] = (ln(1-((`=b'[1,cols(`=b')]:*(ln(`=m'[`i',`=c'RN]))):/exp(`=m'[`i',`=c'XB])))):/`=b'[1,cols(`=b')]
				}
		}
			
		*Reset XB
			forvalues i = 1/`=Obs'{
				mata: `=m'[`i',`=c'XB] = 0
			}
		
	*Determine outcome - Fixed + SCT - Spline 3
		scalar b = "bL1F1_CD_S3"
		forvalues i = 1/`=Obs'{
			*Calculate xb
				*Age
					mata {
						if		(`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cAge] * `=b'[1,1])
					}	
				*Male 
					mata {
						if 		(`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cMale] * `=b'[1,2])
					} 
				*ECOG 
					mata {
						if 		(`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1 & `=m'[`i',cECOGc] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,4] 
						else if (`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1 & `=m'[`i',cECOGc] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,5]
					}
				*ISS 
					mata {
						if 		(`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1 & `=m'[`i',cISS] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,7] 
						else if (`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1 & `=m'[`i',cISS] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,8]
					}
				*CR
					mata {
						if		(`=m'[`i',cCR] == 4 & `=m'[`i',cSCT] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,10] 
						else if	(`=m'[`i',cCR] == 31 & `=m'[`i',cSCT] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,11] 
					} 
				*cons
					mata {
						if		(`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,12]
					} 
					
			*Draw RN conditional on survival to splice cut off 2
				mata {
					if (`=m'[`i',`=c'XB] != 0) 	`=m'[`i',`=c'RN] = runiform(1, 1, 0, exp(-(exp(`=m'[`i',`=c'XB]))*(L1F1_CD_C2^exp(`=b'[1,13]))))
				} 

			*Recalculate survival time for those with OC beyond splice cut off 2 only
				mata {
					if 		(`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1 & f`=b' == "ereg" & `=m'[`i',`=c'OC] > L1F1_CD_C2) 		`=m'[`i',`=c'OC] = (ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])
					else if (`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1 & f`=b' == "weibull" & `=m'[`i',`=c'OC] > L1F1_CD_C2)  	`=m'[`i',`=c'OC] = ((ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])):^(1:/exp(`=b'[1,cols(`=b')]))
					else if (`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 1 & f`=b' == "gompertz" & `=m'[`i',`=c'OC] > L1F1_CD_C2)	`=m'[`i',`=c'OC] = (ln(1-((`=b'[1,cols(`=b')]:*(ln(`=m'[`i',`=c'RN]))):/exp(`=m'[`i',`=c'XB])))):/`=b'[1,cols(`=b')]
				}
		}
		
	*Determine outcome - Fixed - SCT
		scalar b = "bL1F0_CD"
		forvalues i = 1/`=Obs'{
			*Calculate xb
				*Age
					mata {
						if		(`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 0) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cAge] * `=b'[1,1])
					}	
				*Male 
					mata {
						if 		(`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 0) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cMale] * `=b'[1,2])
					} 
				*ECOG 
					mata {
						if 		(`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 0 & `=m'[`i',cECOGc] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,4] 
						else if (`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 0 & `=m'[`i',cECOGc] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,5] 
					}
				*ISS 
					mata {
						if 		(`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 0 & `=m'[`i',cISS] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,7] 
						else if (`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 0 & `=m'[`i',cISS] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,8]
					}
				*CR
					mata {
						if		(`=m'[`i',cCR] == 4 & `=m'[`i',cSCT] == 0) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,10] 
						else if	(`=m'[`i',cCR] == 31 & `=m'[`i',cSCT] == 0) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,11] 
					} 
				*cons
					mata {
						if		(`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 0) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,12]
					} 	
		
			*Calculate survival time
				mata {
					if 		(`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 0 & f`=b' == "ereg") 		`=m'[`i',`=c'OC] = (ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])
					else if (`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 0 & f`=b' == "weibull")  	`=m'[`i',`=c'OC] = ((ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])):^(1:/exp(`=b'[1,cols(`=b')]))
					else if (`=m'[`i',cCR] != 7 & `=m'[`i',cSCT] == 0 & f`=b' == "gompertz")	`=m'[`i',`=c'OC] = (ln(1-((`=b'[1,cols(`=b')]:*(ln(`=m'[`i',`=c'RN]))):/exp(`=m'[`i',`=c'XB])))):/`=b'[1,cols(`=b')]
				}
		}	
			
/*	*Determine outcome - Continuous Therapy
		scalar b = "bL1C_CD"
		forvalues i = 1/`=Obs'{
			*Calculate xb
				*Age
					mata {
						if (`=m'[`i',cCR] == 7) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cAge] * `=b'[1,1])
					} 
				*Male 
					mata {
						if (`=m'[`i',cCR] == 7) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cMale] * `=b'[1,2])
					} 
				*ECOG 
					mata {
						if 		(`=m'[`i',cCR] == 7 & `=m'[`i',cECOGc] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,4] 
						else if (`=m'[`i',cCR] == 7 & `=m'[`i',cECOGc] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,5] 
						else if (`=m'[`i',cCR] == 7 & `=m'[`i',cECOGc] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,6] 
					}
				*ISS 
					mata {
						if 		(`=m'[`i',cCR] == 7 & `=m'[`i',cISS] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,8] 
						else if (`=m'[`i',cCR] == 7 & `=m'[`i',cISS] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,9]
					}
				*cons
					mata {
						if (`=m'[`i',cCR] == 7) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,10]
					} 	
		
			*Calculate survival time
				mata {
					if 		(`=m'[`i',cCR] == 7 & f`=b' == "ereg") 		`=m'[`i',`=c'OC] = (ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])
					else if (`=m'[`i',cCR] == 7 & f`=b' == "weibull")  	`=m'[`i',`=c'OC] = ((ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])):^(1:/exp(`=b'[1,cols(`=b')]))
					else if (`=m'[`i',cCR] == 7 & f`=b' == "gompertz")	`=m'[`i',`=c'OC] = (ln(1-((`=b'[1,11]:*(ln(`=m'[`i',`=c'RN]))):/exp(`=m'[`i',`=c'XB])))):/`=b'[1,11]
				}
		}
	
*/	
	
		*Update outcome matrices	
			forvalues i = 1/`=Obs'{
				mata {
					if 	(mMOR[`i',`=OMC'-1] == 0) mTNE[`i',`=OMC'] = `=m'[`i',`=c'OC]/365.25
					if 	(mMOR[`i',`=OMC'-1] == 0) mTSD[`i',`=OMC'+1] = mTSD[`i',`=OMC'] + mTNE[`i',`=OMC']
					if 	(mMOR[`i',`=OMC'-1] == 0) mCore[`i',cCD] = `=m'[`i',`=c'OC]
					if 	(mMOR[`i',`=OMC'-1] != 0) mCore[`i',cCD] = .
					if 	(mMOR[`i',`=OMC'-1] == 0) mCD[`i',`=LX'+2] = `=m'[`i',`=c'OC]
				}
			}
