
program define olleyp, rclass

    * Author: Martim Leitão, Católica Lisbon and Center of Economics for Prosperity
	* Date: 11/06/2022
	* Modified: added rclass and return scalars so r() values are accessible after the call

	version 17.0
	syntax varlist(min=4 max=4) [ ,BASE(real 1) ///
								   END(real 1)]

	tokenize `varlist'

	di as blue "Time Variable Chosen: `1' "
	di as blue "Unit of Analysis Chosen: `2' "
	di as blue "Productivity Variable Chosen: `3' "
	di as blue "Weight Variable Chosen: `4'"

	di as green "Olley Pakes Productivity Growth Decomposition (in %):"

	// Total Productivity Growth

preserve

	keep `4' `3' `2' `1'

	qui: keep if `1'==`base' | `1'== `end'

	qui: bys `1': egen sumVA=total(`4')

	qui: bys `1': gen share_jt = `4'/sumVA

	// Compute yearly aggregate productivity and change in A_t

	qui: gen lnz_jt=share_jt*`3'

	qui: bys `1': egen A_t=total(lnz_jt)
	gcollapse (mean) A_t, by(`1')
	sort `1'
	qui: gen growth_overall = (A_t - A_t[_n-1])*100
	qui: replace growth_overall = growth_overall[2] if _n==1

	tempfile all
	qui: save `all'

restore

	// For surviving

preserve

	keep `4' `3' `2' `1'

	qui: keep if `1'== `base' | `1'== `end'

	qui: bys `1': egen sumVA=total(`4')

	qui: bys `1': gen share_jt = `4'/sumVA

	qui: bys `2': egen nvalues=nvals(`1')
	qui: drop if nvalues!=2

	// Compute average prod in each subperiod for surving

	qui: bys `1': egen unweigthed_mean = mean(`3')

	// Compute share for surviving

	qui: bys `1': egen share_surv=total(share_jt)

	// Compute prod of surviving

	qui: bys `1': gen prod_js=(1/share_surv)*share_jt*`3'
	qui: bys `1': egen A_st=total(prod_js)

	collapse (mean) A_st unweigthed_mean, by(`1')
	sort `1'
	qui: gen growth_surviving = (A_st - A_st[_n-1])*100
	qui: gen change_unweigthed = (unweigthed_mean-unweigthed_mean[_n-1])*100
	qui: gen reallocation = growth_surviving - change_unweigthed

	qui: replace growth_surviving = growth_surviving[2] if _n==1
	qui: replace change_unweigthed = change_unweigthed[2] if _n==1
	qui: replace reallocation = reallocation[2] if _n==1

	tempfile surviving
	qui: save `surviving'

restore

	// For Newborn


// Share of Newborns

preserve
	keep `4' `3' `2' `1'

	qui: bys `1': egen sumVA=total(`4')

	qui: bys `1': gen share_jt = `4'/sumVA

	qui: keep if `1'==`end'| `1'==`base'

	qui: bys `2': egen nvalues=nvals(`1')
	qui: drop if nvalues==2
	qui: drop if `1'==`base'

	// Compute weights for newborn in t

	qui: bys `1': egen share_newborn_t=total(share_jt)

	// Compute prod of newborn

	qui: bys `1': gen prod_jnb=(1/share_newborn_t)*share_jt*`3'
	qui: bys `1': egen prod_nb=total(prod_jnb)

	collapse (mean) share_newborn_t prod_nb, by(`1')

	tempfile newborn
	qui: save `newborn'

restore

					// Final Table Results

preserve
	qui: use `all', clear
	qui: merge 1:1 `1' using `surviving'
	qui: drop _merge
	qui: merge 1:1 `1' using `newborn'
	qui: drop _merge
	qui: egen tot_share_new = mean(share_newborn_t)
	drop share_newborn_t
	ren tot_share_new share_newborn_t

	qui: gen growth_newborn = share_newborn_t*(prod_nb-A_st)*100

	qui: keep if `1'==`end'
	keep growth_overall growth_surviving growth_newborn reallocation change_unweigthed
	qui: gen growth_exiting = growth_overall - (growth_surviving + growth_newborn)
	ren (growth_overall growth_surviving growth_newborn growth_exiting reallocation change_unweigthed) (Overall Surving Newborn Exiting Reallocation Distributional)
	list
	return scalar overall        = Overall[1]
	return scalar incumbents     = Surving[1]
	return scalar distributional = Distributional[1]
	return scalar reallocation   = Reallocation[1]
	return scalar newborns       = Newborn[1]
	return scalar exiters        = Exiting[1]
	erase `newborn'
	erase `all'
	erase `surviving'

restore

end
