cd "D:\User\Files\Rproject\Pandas_threshold"
forvalues n = 1/1{
	forvalues m = 1/4 {
		forvalues i = 1/10 {
			forvalues j= 1/10 {
			    local ThresholdFile Threshold\Pandas_Sample`i'_run`j'_presence`m'.xlsx
				putexcel set `ThresholdFile', modify sheet(Sheet 1)
				local PETable_file PE_data\PE_Table\Pandas_Sample`i'_run`j'_presence`m'.csv
				import delimited `PETable_file', clear
				local k = 1
				foreach var of varlist glm gam gbm mars maxnet rf {
					local k = `k' + 1
					capture threshold `var'ratio if `var'ratio>0.005, threshvar(`var') trim(10) regionvars(`var') noconstant
					if _rc == 0 {
						putexcel (J`k') = (e(thresholds)[1, 2])
					}
				}
			}
		}
	}	
}
