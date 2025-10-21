**********
	*SIM CL
**********

	*L1 Chemo Regimens
		gen CR_L1 = 0 if(Event0 == 10)
		local CR_L1 4 31 // 4 - Bort/Cycl/Dexa, 31 - Bort/Lena/Dexa, 7 - Lena/Dexa
		foreach l of local CR_L1 {		
			replace CR_L1 = `l' if(Event0 == 10 & Regimen == `l')
		}
		bysort ID (CR_L1): replace CR_L1 = CR_L1[_n-1] if(CR_L1 == .)
		
	*L2 Chemo Regimens
		gen CR_L2 = 0 if(Event0 == 20)
		local CR_L2 7 80 // 7 - Lena/Dexa, 80 - Dara/Bort/Dexa
		foreach l of local CR_L2 {
			replace CR_L2 = `l' if(Event0 == 20 & Regimen == `l')
		}
		bysort ID (CR_L2): replace CR_L2 = CR_L2[_n-1] if(CR_L2 == .)
		
	*L3 Chemo Regimens
		gen CR_L3 = 0 if(Event0 == 30)
		local CR_L3 // 7 - Lena/Dexa, 49 - Carf/Dexa
		foreach l of local CR_L3 {
			replace CR_L3 = `l' if(Event0 == 30 & Regimen == `l')
		}
		bysort ID (CR_L3): replace CR_L3 = CR_L3[_n-1] if(CR_L3 == .)
		
	*L4 Chemo Regimens
		gen CR_L4 = 0 if(Event0 == 40)
		local CR_L4 // 7 - Lena/Dexa, 56 - Poma/Dexa
		foreach l of local CR_L4 {
			replace CR_L4 = `l' if(Event0 == 40 & Regimen == `l')
		}
		bysort ID (CR_L4): replace CR_L4 = CR_L4[_n-1] if(CR_L4 == .)
		
	*L5 Chemo Regimens
		gen CR_L5 = 0 if(Event0 == 50)
		local CR_L5 // 49 - Carf/Dexa, 56 - Poma/Dexa
		foreach l of local CR_L5 {
			replace CR_L5 = `l' if(Event0 == 50 & Regimen == `l')
		}
		bysort ID (CR_L5): replace CR_L5 = CR_L5[_n-1] if(CR_L5 == .)
		
	*L6 Chemo Regimens
		gen CR_L6 = 0 if(Event0 == 60)
		local CR_L6 // 49 - Carf/Dexa, 56 - Poma/Dexa
		foreach l of local CR_L6 {
			replace CR_L6 = `l' if(Event0 == 60 & Regimen == `l')
		}
		bysort ID (CR_L6): replace CR_L6 = CR_L6[_n-1] if(CR_L6 == .)
