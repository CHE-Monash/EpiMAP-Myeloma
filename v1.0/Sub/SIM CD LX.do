**********
	*SIM CD LX
**********

	*Patient Matrix
		mata: `=m' = mCore , mSU
		mata: r`=m' = rmCore
		mata: c`=m' = cmCore \ cmSU
		mata: `=m'[.,`=c'RN] = runiform(`=Obs', 1, 0, 1)
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
					*mata{
						*if	(`=m'[`i',cCR] == 4) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,9] 
					*} 
				*pBCRc
					mata {
						if		(`=m'[`i',cBCR] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,11] 
						else if	(`=m'[`i',cBCR] == 5) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,12] 
					} 
				*cons
					mata: `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,13] 
		
			*Calculate survival time
				mata {
					if 		(f`=b' == "ereg") 		`=m'[`i',`=c'OC] = (ln(`=m'[`i',`=c'RN]))/-exp(`=m'[`i',`=c'XB])
					else if (f`=b' == "weibull")  	`=m'[`i',`=c'OC] = ((ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])):^(1:/exp(`=b'[1,cols(`=b')]))
					else if (f`=b' == "gompertz")	`=m'[`i',`=c'OC] = (ln(1-((`=b'[1,cols(`=b')]:*(ln(`=m'[`i',`=c'RN])))/exp(`=m'[`i',`=c'XB]))))/`=b'[1,cols(`=b')]
				}
		
			*Curtail if beyond last observed in the data
				mata {
					if	(`=Line' == 3 & `=m'[`i',`=c'OC] > maxL4_CD)	`=m'[`i',`=c'OC] = maxL4_CD // Simluation `=Line' is 3 when calculating CD_L4
					if	(`=Line' == 4 & `=m'[`i',`=c'OC] > maxL5_CD)	`=m'[`i',`=c'OC] = maxL5_CD // Simluation `=Line' is 4 when calculating CD_L5
				}
				
			*Update outcome matrices
				mata {
					if 	(mMOR[`i',`=OMC'-1] == 0) mTNE[`i',`=OMC'] = `=m'[`i',`=c'OC]/365.25
					if 	(mMOR[`i',`=OMC'-1] == 0) mTSD[`i',`=OMC'+1] = mTSD[`i',`=OMC'] + mTNE[`i',`=OMC']
					if 	(mMOR[`i',`=OMC'-1] == 0) mCore[`i',cCD] = `=m'[`i',`=c'OC]
					if 	(mMOR[`i',`=OMC'-1] != 0) mCore[`i',cCD] = .
					if 	(mMOR[`i',`=OMC'-1] == 0) mCD[`i',`=LX'+2] = `=m'[`i',`=c'OC]
				}
		}		
