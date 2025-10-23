**********
	*EpiMAP Myeloma - Process Data
**********

capture program drop process_data
program process_data

	di as text "Processing Simulated Data"

	*Create mSum in Mata 
		mata: mSum = mCore , mAge , mOS , mTNE , mTSD , mMOR , mOC , mTXR , mTXD , mBCR , mTFI , mState, mSCT
	
	*Convert mSum to stSum
		mata: st_matrix("stSum", mSum)
		drop _all
		
	*Convert stSum to variables
		svmat double stSum
	
	*Name variables
		local varnames ID Age Male ECOGcc RISS SCT MNT CR CD BCR ///
			Age_DN Age_L1S Age_L1E Age_L2S Age_L2E Age_L3S Age_L3E Age_L4S Age_L4E Age_L5S Age_L5E Age_L6S Age_L6E Age_L7S Age_L7E Age_L8S Age_L8E Age_L9S Age_L9E ///
			OS_DN OS_L1S OS_L1E OS_L2S OS_L2E OS_L3S OS_L3E OS_L4S OS_L4E OS_L5S OS_L5E OS_L6S OS_L6E OS_L7S OS_L7E OS_L8S OS_L8E OS_L9S OS_L9E ///
			TNE_DN TNE_L1S TNE_L1E TNE_L2S TNE_L2E TNE_L3S TNE_L3E TNE_L4S TNE_L4E TNE_L5S TNE_L5E TNE_L6S TNE_L6E TNE_L7S TNE_L7E TNE_L8S TNE_L8E TNE_L9S TNE_L9E ///
			TSD_DN TSD_L1S TSD_L1E TSD_L2S TSD_L2E TSD_L3S TSD_L3E TSD_L4S TSD_L4E TSD_L5S TSD_L5E TSD_L6S TSD_L6E TSD_L7S TSD_L7E TSD_L8S TSD_L8E TSD_L9S TSD_L9E ///
			MOR_DN MOR_L1S MOR_L1E MOR_L2S MOR_L2E MOR_L3S MOR_L3E MOR_L4S MOR_L4E MOR_L5S MOR_L5E MOR_L6S MOR_L6E MOR_L7S MOR_L7E MOR_L8S MOR_L8E MOR_L9S MOR_L9E ///
			OC_TIME OC_MORT ///
			TXR_L1 TXR_L2 TXR_L3 TXR_L4 TXR_L5 TXR_L6 TXR_L7 TXR_L8 TXR_L9 ///
			TXD_L1 TXD_L2 TXD_L3 TXD_L4 TXD_L5 TXD_L6 TXD_L7 TXD_L8 TXD_L9 ///
			BCR_L1 BCR_L2 BCR_L3 BCR_L4 BCR_L5 BCR_L6 BCR_L7 BCR_L8 BCR_L9 BCR_SCT ///
			TFI_L1 TFI_L2 TFI_L3 TFI_L4 TFI_L5 TFI_L6 TFI_L7 TFI_L8 TFI_L9 /// 
			State DateDN ///
			SCT_DN SCT_L1
		
		local varlength : word count `varnames'
		
		forvalues i = 1/`varlength'{
			local currentvar : word `i' of `varnames'
			rename stSum`i' `currentvar'
		}		
	
		format DateDN %td
	
	*Drop unnecessary variables
		drop CR CD BCR
		order ID Age Male ECOGcc RISS SCT MNT
		
	*Label
		label values State State_lbl
	
	*Generate Dates
		qui {
			gen DateL1S = DateDN + (TNE_DN*30.4375)
			gen DateL1E = DateL1S + (TNE_L1S*30.4375)
			gen DateL2S = DateL1E + (TNE_L1E*30.4375)
			gen DateL2E = DateL2S + (TNE_L2S*30.4375)
			gen DateL3S = DateL2E + (TNE_L2E*30.4375)
			gen DateL3E = DateL3S + (TNE_L3S*30.4375)
			gen DateL4S = DateL3E + (TNE_L3E*30.4375)
			gen DateL4E = DateL4S + (TNE_L4S*30.4375)
			gen DateL5S = DateL4E + (TNE_L4E*30.4375)
			gen DateL5E = DateL5S + (TNE_L5S*30.4375)
			gen DateL6S = DateL5E + (TNE_L5E*30.4375)
			gen DateL6E = DateL6S + (TNE_L6S*30.4375)
			gen DateL7S = DateL6E + (TNE_L6E*30.4375)
			gen DateL7E = DateL7S + (TNE_L7S*30.4375)
			gen DateL8S = DateL7E + (TNE_L7E*30.4375)
			gen DateL8E = DateL8S + (TNE_L8S*30.4375)
			gen DateL9S = DateL8E + (TNE_L8E*30.4375)
			gen DateL9E = DateL9S + (TNE_L9S*30.4375)
			gen DateSCT = DateL1E + 1 if(SCT == 1) // Fix DateSCT 1 day after DateL1E
			gen DateMOR = DateDN + (OC_TIME*30.4375)
		}
		format Date* %td
	
	*Generate Years
		qui {
		gen YearDN = yofd(DateDN)
		gen YearL1 = yofd(DateL1S)
		gen YearL2 = yofd(DateL2S)
		gen YearL3 = yofd(DateL3S)
		gen YearL4 = yofd(DateL4S)
		gen YearL5 = yofd(DateL5S)
		gen YearL6 = yofd(DateL6S)
		gen YearL7 = yofd(DateL7S)
		gen YearL8 = yofd(DateL8S)
		gen YearL9 = yofd(DateL9S)
		gen YearSCT = yofd(DateSCT)
		gen YearMOR = yofd(DateMOR)
		}


end
