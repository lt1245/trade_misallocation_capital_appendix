* ============================================================
* SETUP: packages, paths, directory creation
* Replication package location: //share_krtk/Kozos/ref_tmcmi
* ============================================================

* ===== REQUIRED PACKAGES =====
capture ssc install matsort, replace
capture ssc install gtools, replace
capture ssc install egenmore, replace
capture ssc install require, replace
capture ssc install reghdfe, replace
capture ssc install eventdd, replace

capture ssc install ftools,replace

/* *===== PASTE OLLEY-PAKES ADO TO PERSONAL ADO PATH =====
capture confirm file "olleyp.ado"

if _rc {
    display as error "olleyp.ado not found in `c(pwd)'"
    exit 601
}


copy "olleyp.ado" "`c(sysdir_personal)'olleyp.ado", replace
*/


* ===== HELPER: create a directory and all its parents =====
capture program drop mkdirs
program define mkdirs
    args path
    local path = subinstr(`"`path'"', "\", "/", .)

    * keep a leading "/" (root) or "//" (UNC share) intact
    local prefix ""
    if substr("`path'", 1, 2) == "//" {
        local prefix "//"
        local path   = substr("`path'", 3, .)
    }
    else if substr("`path'", 1, 1) == "/" {
        local prefix "/"
        local path   = substr("`path'", 2, .)
    }

    local sofar "`prefix'"
    local rest  "`path'"
    while `"`rest'"' != "" {
        gettoken piece rest : rest, parse("/")
        if `"`piece'"' == "/" {
            local sofar "`sofar'/"
        }
        else {
            local sofar "`sofar'`piece'"
            capture mkdir "`sofar'"
        }
    }
end
/*
* ===== PATHS =====
* Create the package root, then anchor everything to an ABSOLUTE path so
* the folder macros keep working after we cd into Empirics.
mkdirs "../replication_package"
cd "../replication_package"
*/

global root     = "/share_krtk/Kozos/ref_tmcmi"
global empirics "$root/Empirics"
global figures  "$root/Figures"
global tables   "$root/Tables"

mkdirs "$empirics"
mkdirs "$figures"
mkdirs "$tables"

cd "$empirics"



* ===== APPENDING THE FIRM PANEL (run once) =====
global nav "//share_krtk/Adattar/NAV_merleg"
*
use "$nav/2000/apeh_2000.dta", clear
forvalues y = 2001/2017 {
    append using "$nav/`y'/apeh_`y'.dta"
}
save "$empirics/firm_data_appended.dta", replace

*Assuming the append above already ran successfully:
use "$empirics/firm_data_appended.dta", clear


* ===== VARIABLE CONSTRUCTION =====
gen value_added = ereduzem + kecs + rszem
replace immat = 0 if immat ==.
gen capital = targyie + immat
gen rovkot1 = rovkot
replace rovkot1 = 0 if (hoskot!=. & rovkot ==.)
gen hoskot1 = hoskot
replace hoskot1 = 0 if (rovkot!=. & hoskot ==.)
gen equity= eszk - rovkot1 - hoskot1
gen log_equity = log(equity)
gen no_access_for_credit = 0 if (forsh>0 & forsh!=.)
replace no_access_for_credit =1 if no_access_for_credit==.
gen access_for_credit = 1 if (forsh>0 & forsh!=.)
replace access_for_credit =0 if no_access_for_credit==.
scalar alpha1 = 0.35
gen log_ARPL = log(value_added) - log(letszam)
gen log_ARPK = log(value_added) - log(capital)
gen year = ev

gen log_TFPR = alpha1 * log_ARPK  + (1 - alpha1) * log_ARPL
* Industry classification
gen industry = teaor1998_2
replace industry= teaor2008_2 if industry==.
replace industry= teaor2003_2 if industry==.
bysort sorszam: egen industry_initial = max(teaor1998_2)

gen ind10 = .
label define ind10 ///
 1 "Mezőgazdaság" ///
 2 "Bányászat" ///
 3 "Feldolgozóipar" ///
 4 "Energia, közmű" ///
 5 "Építőipar" ///
 6 "Kereskedelem" ///
 7 "Szállítás, raktár" ///
 8 "Pénzügy" ///
 9 "Üzleti szolgáltatások" ///
10 "Közszféra, egyéb"
label values ind10 ind10
replace ind10 = 1 if inrange(industry,1,5)
replace ind10 = 2 if inrange(industry,10,14) & ev <= 2007
replace ind10 = 2 if inrange(industry,5,9)   & ev >= 2008
replace ind10 = 3 if inrange(industry,15,37) & ev <= 2007
replace ind10 = 3 if inrange(industry,10,33) & ev >= 2008
replace ind10 = 4 if inrange(industry,40,41) & ev <= 2007
replace ind10 = 4 if inrange(industry,35,39) & ev >= 2008
replace ind10 = 5 if industry == 45 & ev <= 2007
replace ind10 = 5 if inrange(industry,41,43) & ev >= 2008
replace ind10 = 6 if inrange(industry,50,52) & ev <= 2007
replace ind10 = 6 if inrange(industry,45,47) & ev >= 2008
replace ind10 = 7 if inrange(industry,60,64) & ev <= 2007
replace ind10 = 7 if inrange(industry,49,53) & ev >= 2008
replace ind10 = 8 if inrange(industry,65,67) & ev <= 2007
replace ind10 = 8 if inrange(industry,64,66) & ev >= 2008
replace ind10 = 9 if inrange(industry,70,74) & ev <= 2007
replace ind10 = 9 if inrange(industry,69,82) & ev >= 2008
replace ind10 = 10 if ind10 == .

gen exportshare = arbevexp/arbevert
replace exportshare = 0 if (exportshare<0)
replace exportshare = 1 if (exportshare>1 & exportshare!=.)
replace exportshare = 0 if (exportshare==. & value_added>0& value_added!=.)
gen exportdummy = 1 if exportshare>0
replace exportdummy =0 if exportshare==0
gen export_value = exportshare*value_added
gen debt = eszk - equity
gen leverage = 100*debt/eszk
gen log_rev = log(arbevert)
gen log_value_added = log(value_added)
gen log_capital = log(capital)
gen log_letszam = log(letszam)
drop if sorszam ==.
drop if year ==.

* ===== SAMPLE RESTRICTIONS =====
drop if year>2008
bysort sorszam:egen flag_value_added = max(value_added<=0 | value_added==.)
bysort sorszam:egen flag_capital = max(capital<=0 | capital==.)
bysort sorszam:egen flag_equity = max(equity<=0 | equity==.)
bysort sorszam:egen flag_letszam = max(letszam<=1 | letszam==.)
drop if flag_value_added==1
drop if flag_capital==1
drop if flag_equity==1
drop if flag_letszam==1

duplicates drop sorszam year, force
bysort sorszam:egen flag_2000 = max(year ==2000)
bysort sorszam:egen flag_2001 = max(year ==2001)
bysort sorszam:egen flag_2002 = max(year ==2002)
bysort sorszam:egen flag_2003 = max(year ==2003)
bysort sorszam:egen flag_2004 = max(year ==2004)
bysort sorszam:egen flag_2005 = max(year ==2005)
bysort sorszam:egen flag_2006 = max(year ==2006)
bysort sorszam:egen flag_2007 = max(year ==2007)
bysort sorszam:egen flag_2008 = max(year ==2008)
gen flag_incumbent = 1 if(flag_2000==1 &flag_2001==1 &flag_2002==1 &flag_2003==1 &flag_2004==1 &flag_2005==1 &flag_2006==1 &flag_2007==1 &flag_2008==1)
replace flag_incumbent = 0 if flag_incumbent==.
drop if(flag_2000==1 & flag_2008==1 & (flag_2001 + flag_2002 + flag_2003+ flag_2004+ flag_2005+ flag_2006+ flag_2007<7))
drop if(flag_2000==1 & flag_2007==1 & (flag_2001 + flag_2002 + flag_2003+ flag_2004+ flag_2005+ flag_2006<6))
drop if(flag_2001==1 & flag_2008==1 & (flag_2007 + flag_2002 + flag_2003+ flag_2004+ flag_2005+ flag_2006<6))
drop if(flag_2000==1 & flag_2006==1 & (flag_2001 + flag_2002 + flag_2003+ flag_2004+ flag_2005<5))
drop if(flag_2001==1 & flag_2007==1 & (flag_2002 + flag_2003+ flag_2004+ flag_2005+ flag_2006<5))
drop if(flag_2002==1 & flag_2008==1 & (flag_2007 + flag_2003+ flag_2004+ flag_2005+ flag_2006<5))
drop if(flag_2000==1 & flag_2005==1 & (flag_2001 + flag_2002 + flag_2003+ flag_2004<4))
drop if(flag_2001==1 & flag_2006==1 & (flag_2002 + flag_2003+ flag_2004+ flag_2005<4))
drop if(flag_2002==1 & flag_2007==1 & (flag_2003+ flag_2004+ flag_2005+ flag_2006<4))
drop if(flag_2003==1 & flag_2008==1 & (flag_2007+ flag_2004+ flag_2005+ flag_2006<4))
drop if(flag_2000==1 & flag_2004==1 & (flag_2001 + flag_2002 + flag_2003<3))
drop if(flag_2001==1 & flag_2005==1 & (flag_2002 + flag_2003+ flag_2004<3))
drop if(flag_2002==1 & flag_2006==1 & (flag_2003+ flag_2004+ flag_2005<3))
drop if(flag_2003==1 & flag_2007==1 & (flag_2006+ flag_2004+ flag_2005<3))
drop if(flag_2004==1 & flag_2008==1 & (flag_2006+ flag_2007+ flag_2005<3))
drop if(flag_2000==1 & flag_2003==1 & (flag_2001 + flag_2002<2))
drop if(flag_2001==1 & flag_2004==1 & (flag_2002 + flag_2003<2))
drop if(flag_2002==1 & flag_2005==1 & (flag_2003+ flag_2004<2))
drop if(flag_2003==1 & flag_2006==1 & (flag_2004+ flag_2005<2))
drop if(flag_2004==1 & flag_2007==1 & (flag_2006+ flag_2005<2))
drop if(flag_2005==1 & flag_2008==1 & (flag_2006+ flag_2007<2))
drop if(flag_2000==1 & flag_2002==1 & (flag_2001<1))
drop if(flag_2001==1 & flag_2003==1 & (flag_2002<1))
drop if(flag_2002==1 & flag_2004==1 & (flag_2003<1))
drop if(flag_2003==1 & flag_2005==1 & (flag_2004<1))
drop if(flag_2004==1 & flag_2006==1 & (flag_2005<1))
drop if(flag_2005==1 & flag_2007==1 & (flag_2006<1))
drop if(flag_2006==1 & flag_2008==1 & (flag_2007<1))


* ===================================================================
* TABLE data_transition + TABLE data_regression_simple
* (main.tex Appendix lines ~990-1008 and ~1051-1077)
* Values printed below are copy-pasted into the hardcoded tables in main.tex
* ===================================================================
preserve
scalar final_year = 2008
gen period_considered =1 if year ==2000
replace period_considered =2 if year ==final_year
keep if (year==2000 | year==final_year)
gen flag_final_year = flag_2008
xtset sorszam period_considered

gen exportentrants= D.exportdummy
gen no_access_for_credit_change= D.no_access_for_credit
gen capital_growth = D.log_capital
gen value_added_growth = D.log_value_added
gen log_TFPR_growth = D.log_TFPR
gen leverage_change = D.leverage
bysort sorszam:egen flag_exporter_ever = max(exportdummy ==1)
bysort sorszam:egen flag_exporterentrant_ever = max(exportentrants ==1)
bysort sorszam:egen flag_exportexit_ever = max(exportentrants ==-1)
gen exportincumbent = 1 if (flag_exporter_ever ==1 & flag_exporterentrant_ever==0 & flag_exportexit_ever==0 & flag_2000 == 1 & flag_final_year == 1)
replace exportincumbent =0 if exportincumbent==.
gen flag_exporterentrant = 1 if flag_exporterentrant_ever ==1
replace flag_exporterentrant = 1 if (flag_exporterentrant ==. & flag_exporter_ever ==1 & flag_2000 ==0)
replace flag_exporterentrant= 0 if flag_exporterentrant==.
gen flag_exportexit = 1 if flag_exportexit_ever==1
replace flag_exportexit = 1 if (flag_exportexit ==. & flag_exporter_ever ==1  & flag_final_year ==0)
replace flag_exportexit = 0 if flag_exportexit==.
bysort sorszam:egen flag_foreign_owned_ever = max(no_access_for_credit ==0)
bysort sorszam:egen flag_foreign_credit_acquired = max(no_access_for_credit_change ==-1)
bysort sorszam:egen flag_foreign_credit_lost = max(no_access_for_credit_change ==1)
gen no_access_credit_ever = 1 if (flag_foreign_owned_ever ==0 )
replace no_access_credit_ever =0 if no_access_credit_ever==.
gen access_credit_always = 1 if (flag_foreign_owned_ever ==1 & flag_foreign_credit_acquired==0 & flag_foreign_credit_lost==0)
replace access_credit_always =0 if access_credit_always==.

gen exportdummy_foreign = exportdummy if no_access_for_credit==0
replace exportdummy_foreign=0 if exportdummy_foreign==.
gen exportdummy_domestic = exportdummy if no_access_for_credit==1
replace exportdummy_domestic=0 if exportdummy_domestic==.
gen domesticdummy_foreign = 1 if (no_access_for_credit==0 & exportdummy == 0)
replace domesticdummy_foreign=0 if domesticdummy_foreign==.
gen domesticdummy_domestic = 1 if (no_access_for_credit==1 & exportdummy == 0)
replace domesticdummy_domestic=0 if domesticdummy_domestic==.

gen tabulatable_dummy = 1 if (exportdummy == 1 & no_access_for_credit==1 & period_considered ==1 & flag_final_year==1 & exportincumbent ==1 &no_access_credit_ever==1)
replace tabulatable_dummy =2 if (exportdummy == 0 & no_access_for_credit==1& period_considered ==1 & flag_final_year==1 & flag_exporterentrant_ever ==0 &no_access_credit_ever==1)
replace tabulatable_dummy =3 if (exportdummy == 1 & no_access_for_credit==0 & period_considered ==1 & flag_final_year==1 & exportincumbent ==1 &access_credit_always==1)
replace tabulatable_dummy =4 if (exportdummy == 0 & no_access_for_credit==0 & period_considered ==1 & flag_final_year==1 & flag_exporterentrant_ever ==0 &access_credit_always==1)
replace tabulatable_dummy =5 if (exportdummy == 1 & no_access_for_credit==1 & period_considered ==1 & flag_final_year==1 & flag_exportexit ==1 &no_access_credit_ever==1)
replace tabulatable_dummy =6 if (exportdummy == 0 & no_access_for_credit==1 & period_considered ==1 & flag_final_year==1 & flag_exporterentrant ==1 &no_access_credit_ever==1)
replace tabulatable_dummy =7 if (exportdummy == 1 & no_access_for_credit==0 & period_considered ==1 & flag_final_year==1 & flag_exportexit ==1 &access_credit_always==1)
replace tabulatable_dummy =8 if (exportdummy == 0 & no_access_for_credit==0 & period_considered ==1 & flag_final_year==1 & flag_exporterentrant ==1 &access_credit_always==1)
replace tabulatable_dummy =9 if (exportdummy == 1 & no_access_for_credit==1 & period_considered ==1 & flag_final_year==1 & exportincumbent ==1 &flag_foreign_credit_acquired==1)
replace tabulatable_dummy =10 if (exportdummy == 0 & no_access_for_credit==1 & period_considered ==1 & flag_final_year==1 & flag_exporterentrant_ever ==0 &flag_foreign_credit_acquired==1)
replace tabulatable_dummy =11 if (exportdummy == 1 & no_access_for_credit==0 & period_considered ==1 & flag_final_year==1 & exportincumbent ==1 &flag_foreign_credit_lost==1)
replace tabulatable_dummy =12 if (exportdummy == 0 & no_access_for_credit==0 & period_considered ==1 & flag_final_year==1 & flag_exporterentrant_ever ==0 &flag_foreign_credit_lost==1)
replace tabulatable_dummy =13 if (exportdummy == 1 & no_access_for_credit==1 & period_considered ==1 & flag_final_year==1 & flag_exportexit ==1 &flag_foreign_credit_acquired==1)
replace tabulatable_dummy =14 if (exportdummy == 0 & no_access_for_credit==1 & period_considered ==1 & flag_final_year==1 & flag_exporterentrant ==1 &flag_foreign_credit_acquired==1)
replace tabulatable_dummy =15 if (exportdummy == 1 & no_access_for_credit==0 & period_considered ==1 & flag_final_year==1 & flag_exportexit ==1 &flag_foreign_credit_lost==1)
replace tabulatable_dummy =16 if (exportdummy == 0 & no_access_for_credit==0 & period_considered ==1 & flag_final_year==1 & flag_exporterentrant ==1 &flag_foreign_credit_lost==1)
replace tabulatable_dummy =17 if (exportdummy == 1 & no_access_for_credit==1 & period_considered ==1 & flag_final_year==0)
replace tabulatable_dummy =18 if (exportdummy == 0 & no_access_for_credit==1 & period_considered ==1 & flag_final_year==0)
replace tabulatable_dummy =19 if (exportdummy == 1 & no_access_for_credit==0 & period_considered ==1 & flag_final_year==0)
replace tabulatable_dummy =20 if (exportdummy == 0 & no_access_for_credit==0 & period_considered ==1 & flag_final_year==0)
replace tabulatable_dummy =21 if (exportdummy == 1 & no_access_for_credit==1 & period_considered ==2 & flag_2000==0)
replace tabulatable_dummy =22 if (exportdummy == 0 & no_access_for_credit==1 & period_considered ==2 & flag_2000==0)
replace tabulatable_dummy =23 if (exportdummy == 1 & no_access_for_credit==0 & period_considered ==2 & flag_2000==0)
replace tabulatable_dummy =24 if (exportdummy == 0 & no_access_for_credit==0 & period_considered ==2 & flag_2000==0)

tabulate tabulatable_dummy, matcell(distribution_values)

gen tabulatable_dummy_sum = 1 if (exportdummy == 1 & no_access_for_credit==0 & period_considered ==1)
replace tabulatable_dummy_sum =2 if (exportdummy == 0 & no_access_for_credit==0& period_considered ==1)
replace tabulatable_dummy_sum =3 if (exportdummy == 1 & no_access_for_credit==1 & period_considered ==1)
replace tabulatable_dummy_sum =4 if (exportdummy == 0 & no_access_for_credit==1 & period_considered ==1)
replace tabulatable_dummy_sum =5 if (exportdummy == 1 & no_access_for_credit==0 & period_considered ==2)
replace tabulatable_dummy_sum =6 if (exportdummy == 0 & no_access_for_credit==0 & period_considered ==2)
replace tabulatable_dummy_sum =7 if (exportdummy == 1 & no_access_for_credit==1 & period_considered ==2)
replace tabulatable_dummy_sum =8 if (exportdummy == 0 & no_access_for_credit==1 & period_considered ==2)
gen tabulatable_dummy_exit =0 if (flag_2000 ==0 &flag_final_year ==1)
replace tabulatable_dummy_exit =1 if (flag_2000 ==1 &flag_final_year ==0)
gen tabulatable_dummy_total =0 if (period_considered == 1)
replace tabulatable_dummy_total =1 if ( period_considered == 2)
tabulate tabulatable_dummy_sum, matcell(distribution_values_sum)
tabulate tabulatable_dummy_exit, matcell(distribution_values_exit)
tabulate tabulatable_dummy_total, matcell(distribution_values_total)

* --- Table data_transition (main.tex ~line 990) ---
di "=== TABLE data_transition rows ==="
di distribution_values[3,1] "& " distribution_values[7,1] "& " distribution_values[11,1] "& " distribution_values[15,1] "& " distribution_values[19,1] "& " distribution_values_sum[1,1] "\\ "
di distribution_values[8,1] "& " distribution_values[4,1] "& " distribution_values[16,1] "& " distribution_values[12,1] "& " distribution_values[20,1] "& " distribution_values_sum[2,1] "\\ "
di distribution_values[9,1] "& " distribution_values[13,1] "& " distribution_values[1,1] "& " distribution_values[5,1] "& " distribution_values[17,1] "& "  distribution_values_sum[3,1] "\\ "
di distribution_values[14,1] "& " distribution_values[10,1] "& " distribution_values[6,1] "& " distribution_values[2,1] "& " distribution_values[18,1] "& " distribution_values_sum[4,1] "\\ "
di distribution_values[23,1] "& " distribution_values[24,1] "& " distribution_values[21,1] "& " distribution_values[22,1] "& " distribution_values_exit[1,1] "\\ "
di distribution_values_sum[5,1]  "& " distribution_values_sum[6,1] "& " distribution_values_sum[7,1] "& " distribution_values_sum[8,1] "& " distribution_values_exit[2,1] "& " (distribution_values_total[2,1]+distribution_values_exit[2,1]) "\\ "
scalar tot2008 = distribution_values_total[2,1]+distribution_values_exit[2,1]
file open dtfile using "$tables/Table_data_transition.tex", write replace
file write dtfile "\begin{table}" _n
file write dtfile "\caption{ Transition matrix}" _n
file write dtfile "    \centering" _n
file write dtfile "     \scalebox{0.45}{\begin{tabular}{cccccc|c}" _n
file write dtfile "2000 vs. 2008 & Exporter w foreign owner & Non-exporters w foreign owner & Exporter w/o foreign owner & Non-exporters w/o foreign owner & Exits & Total in 2000\\" _n
file write dtfile "\hline" _n
file write dtfile "Exporter w foreign owner & " %9.0g (distribution_values[3,1]) " & " %9.0g (distribution_values[7,1]) " & " %9.0g (distribution_values[11,1]) " & " %9.0g (distribution_values[15,1]) " & " %9.0g (distribution_values[19,1]) " & " %9.0g (distribution_values_sum[1,1]) " \\" _n
file write dtfile "Non-exporters w foreign owner & " %9.0g (distribution_values[8,1]) " & " %9.0g (distribution_values[4,1]) " & " %9.0g (distribution_values[16,1]) " & " %9.0g (distribution_values[12,1]) " & " %9.0g (distribution_values[20,1]) " & " %9.0g (distribution_values_sum[2,1]) " \\" _n
file write dtfile "Exporter w/o foreign owner & " %9.0g (distribution_values[9,1]) " & " %9.0g (distribution_values[13,1]) " & " %9.0g (distribution_values[1,1]) " & " %9.0g (distribution_values[5,1]) " & " %9.0g (distribution_values[17,1]) " & " %9.0g (distribution_values_sum[3,1]) " \\" _n
file write dtfile "Non-exporters w/o foreign owner & " %9.0g (distribution_values[14,1]) " & " %9.0g (distribution_values[10,1]) " & " %9.0g (distribution_values[6,1]) " & " %9.0g (distribution_values[2,1]) " & " %9.0g (distribution_values[18,1]) " & " %9.0g (distribution_values_sum[4,1]) " \\" _n
file write dtfile "Entrants & " %9.0g (distribution_values[23,1]) " & " %9.0g (distribution_values[24,1]) " & " %9.0g (distribution_values[21,1]) " & " %9.0g (distribution_values[22,1]) " & - & " %9.0g (distribution_values_exit[1,1]) " \\" _n
file write dtfile "\hline" _n
file write dtfile "Total in 2008 & " %9.0g (distribution_values_sum[5,1]) " & " %9.0g (distribution_values_sum[6,1]) " & " %9.0g (distribution_values_sum[7,1]) " & " %9.0g (distribution_values_sum[8,1]) " & " %9.0g (distribution_values_exit[2,1]) " & " %9.0g (tot2008) " \\" _n
file write dtfile "        \end{tabular}} \caption*{ \scriptsize \textit{Note:} Each row corresponds to the status in 2000, and each column correspond to the status in 2008. Note that the total number of firms are lower than in Table \ref{table:data_micro}, as some firms enter after 2000 and exit before 2008.}" _n
file write dtfile "        \label{table:data_transition}" _n
file write dtfile "    \end{table}" _n
file close dtfile

* --- Table data_transition_probability (main.tex ~line 995) ---
di "=== TABLE data_transition_probability rows ==="
matrix prob_trans = (distribution_values[3,1]/distribution_values_sum[1,1],distribution_values[7,1]/distribution_values_sum[1,1],distribution_values[11,1]/distribution_values_sum[1,1],distribution_values[15,1]/distribution_values_sum[1,1],distribution_values[19,1]/distribution_values_sum[1,1], distribution_values_sum[1,1]/(distribution_values_total[2,1]+distribution_values_exit[2,1]) \ distribution_values[8,1]/distribution_values_sum[2,1],distribution_values[4,1]/distribution_values_sum[2,1],distribution_values[16,1]/distribution_values_sum[2,1],distribution_values[12,1]/distribution_values_sum[2,1],distribution_values[20,1]/distribution_values_sum[2,1],distribution_values_sum[2,1]/(distribution_values_total[2,1]+distribution_values_exit[2,1]) \ distribution_values[9,1]/distribution_values_sum[3,1],distribution_values[13,1]/distribution_values_sum[3,1],distribution_values[1,1]/distribution_values_sum[3,1],distribution_values[5,1]/distribution_values_sum[3,1],distribution_values[17,1]/distribution_values_sum[3,1],distribution_values_sum[3,1]/(distribution_values_total[2,1]+distribution_values_exit[2,1]) \ distribution_values[14,1]/distribution_values_sum[4,1],distribution_values[10,1]/distribution_values_sum[4,1],distribution_values[6,1]/distribution_values_sum[4,1],distribution_values[2,1]/distribution_values_sum[4,1],distribution_values[18,1]/distribution_values_sum[4,1],distribution_values_sum[4,1]/(distribution_values_total[2,1]+distribution_values_exit[2,1]) \ distribution_values[23,1]/distribution_values_exit[1,1] ,distribution_values[24,1]/distribution_values_exit[1,1] ,distribution_values[21,1]/distribution_values_exit[1,1] ,distribution_values[22,1]/distribution_values_exit[1,1] ,0,distribution_values_exit[1,1] /(distribution_values_total[2,1]+distribution_values_exit[2,1]) \  distribution_values_sum[5,1]/ (distribution_values_total[2,1]+distribution_values_exit[2,1]) ,distribution_values_sum[6,1]/ (distribution_values_total[2,1]+distribution_values_exit[2,1]) ,distribution_values_sum[7,1]/ (distribution_values_total[2,1]+distribution_values_exit[2,1]) ,distribution_values_sum[8,1]/ (distribution_values_total[2,1]+distribution_values_exit[2,1]) ,distribution_values_exit[2,1]/ (distribution_values_total[2,1]+distribution_values_exit[2,1]),1 )
matrix prob_trans_100 = 100 * prob_trans
mata : Mrounded=round(st_matrix("prob_trans_100"),1.0)
mata : st_matrix("prob_trans_100",Mrounded)
di "Exporter w foreign owner & " prob_trans_100[1,1] "& "prob_trans_100[1,2] "& " prob_trans_100[1,3] "& " prob_trans_100[1,4] "& " prob_trans_100[1,5] "& " prob_trans_100[1,6] "\\ "
di "Non-exporters w foreign owner & " prob_trans_100[2,1] "& "prob_trans_100[2,2] "& " prob_trans_100[2,3] "& " prob_trans_100[2,4] "& " prob_trans_100[2,5] "& " prob_trans_100[2,6] "\\ "
di "Exporter w/o foreign owner & " prob_trans_100[3,1] "& "prob_trans_100[3,2] "& " prob_trans_100[3,3] "& " prob_trans_100[3,4] "& " prob_trans_100[3,5] "& " prob_trans_100[3,6] "\\ "
di "Non-exporters w/o foreign owner & " prob_trans_100[4,1] "& "prob_trans_100[4,2] "& " prob_trans_100[4,3] "& " prob_trans_100[4,4] "& " prob_trans_100[4,5] "& " prob_trans_100[4,6] "\\ "
di "Entrants & " prob_trans_100[5,1] "& "prob_trans_100[5,2] "& " prob_trans_100[5,3] "& " prob_trans_100[5,4] "& " prob_trans_100[5,5] "& " prob_trans_100[5,6] "\\ "
di "Total in 2008 & " prob_trans_100[6,1] "& "prob_trans_100[6,2] "& " prob_trans_100[6,3] "& " prob_trans_100[6,4] "& " prob_trans_100[6,5] "& " prob_trans_100[6,6] "\\ "
file open dpfile using "$tables/Table_data_transition_probability.tex", write replace
file write dpfile "\begin{table}" _n
file write dpfile "\caption{ Transition probability matrix}" _n
file write dpfile "\label{table:data_transition_probability}" _n
file write dpfile "    \centering" _n
file write dpfile "     \scalebox{0.45}{\begin{tabular}{cccccc|c}" _n
file write dpfile "2000 vs. 2008 & Exporter w foreign owner & Non-exporters w foreign owner & Exporter w/o foreign owner & Non-exporters w/o foreign owner & Exits & Total in 2000\\" _n
file write dpfile "\hline" _n
file write dpfile "Exporter w foreign owner & " %9.0g (prob_trans_100[1,1]) " & " %9.0g (prob_trans_100[1,2]) " & " %9.0g (prob_trans_100[1,3]) " & " %9.0g (prob_trans_100[1,4]) " & " %9.0g (prob_trans_100[1,5]) " & " %9.0g (prob_trans_100[1,6]) " \\" _n
file write dpfile "Non-exporters w foreign owner & " %9.0g (prob_trans_100[2,1]) " & " %9.0g (prob_trans_100[2,2]) " & " %9.0g (prob_trans_100[2,3]) " & " %9.0g (prob_trans_100[2,4]) " & " %9.0g (prob_trans_100[2,5]) " & " %9.0g (prob_trans_100[2,6]) " \\" _n
file write dpfile "Exporter w/o foreign owner & " %9.0g (prob_trans_100[3,1]) " & " %9.0g (prob_trans_100[3,2]) " & " %9.0g (prob_trans_100[3,3]) " & " %9.0g (prob_trans_100[3,4]) " & " %9.0g (prob_trans_100[3,5]) " & " %9.0g (prob_trans_100[3,6]) " \\" _n
file write dpfile "Non-exporters w/o foreign owner & " %9.0g (prob_trans_100[4,1]) " & " %9.0g (prob_trans_100[4,2]) " & " %9.0g (prob_trans_100[4,3]) " & " %9.0g (prob_trans_100[4,4]) " & " %9.0g (prob_trans_100[4,5]) " & " %9.0g (prob_trans_100[4,6]) " \\" _n
file write dpfile "Entrants & " %9.0g (prob_trans_100[5,1]) " & " %9.0g (prob_trans_100[5,2]) " & " %9.0g (prob_trans_100[5,3]) " & " %9.0g (prob_trans_100[5,4]) " & " %9.0g (prob_trans_100[5,5]) " & " %9.0g (prob_trans_100[5,6]) " \\" _n
file write dpfile "\hline" _n
file write dpfile "Total in 2008 & " %9.0g (prob_trans_100[6,1]) " & " %9.0g (prob_trans_100[6,2]) " & " %9.0g (prob_trans_100[6,3]) " & " %9.0g (prob_trans_100[6,4]) " & " %9.0g (prob_trans_100[6,5]) " & " %9.0g (prob_trans_100[6,6]) " \\" _n
file write dpfile "        \end{tabular}}\caption*{\scriptsize \textit{Note:} This table presents the empirical transition probability of Hungarian firms, with each row corresponding to the status in 2000, and each column correspond to the status in 2008. Totals add up in 2000 and in 2008 as the difference is accounted by entrants and exits. Firms that do not exist in either 2000 or in 2008 are omitted. Transitions within this 9 year period are not represented in this matrix.}" _n
file write dpfile "    \end{table}" _n
file close dpfile

* --- Table data_regression_simple (main.tex ~line 1051) ---
di "=== TABLE data_regression_simple rows ==="
reg capital_growth no_access_for_credit
mata : Mrounded=round(st_matrix("e(b)"),0.0001)
mata : st_matrix("coeff_simple_reg",Mrounded)
matrix _V_dlK_simple = e(V)
scalar dlK_dom     = 100*(coeff_simple_reg[1,1]+coeff_simple_reg[1,2])
scalar dlK_for     = 100*coeff_simple_reg[1,2]
scalar dlK_dom_se  = round(10*sqrt(e(V)[1,1]+e(V)[2,2]),0.01)
scalar dlK_for_se  = round(10*sqrt(e(V)[2,2]),0.01)
di "$\Delta log K $&  \multicolumn{2}{c}{" dlK_dom "}&\multicolumn{2}{c}{" dlK_for "}\\"
di "s.e.&  \multicolumn{2}{c}{(" dlK_dom_se ")}&\multicolumn{2}{c}{(" dlK_for_se ")}\\"

reg capital_growth exportdummy_domestic domesticdummy_domestic exportdummy_foreign domesticdummy_foreign, noconstant
mata : Mrounded=round(st_matrix("e(b)"),0.0001)
mata : st_matrix("coeff_simple_reg_exporting",Mrounded)
scalar dlK_expdom  = 100*coeff_simple_reg_exporting[1,1]
scalar dlK_nexpdom = 100*coeff_simple_reg_exporting[1,2]
scalar dlK_expfor  = 100*coeff_simple_reg_exporting[1,3]
scalar dlK_nexpfor = 100*coeff_simple_reg_exporting[1,4]
scalar dlK_expdom_se  = 10*round(sqrt(e(V)[1,1]),0.01)
scalar dlK_nexpdom_se = 10*round(sqrt(e(V)[2,2]),0.01)
scalar dlK_expfor_se  = 10*round(sqrt(e(V)[3,3]),0.01)
scalar dlK_nexpfor_se = 10*round(sqrt(e(V)[4,4]),0.01)
di "$\Delta log K $&  " dlK_expdom "&" dlK_nexpdom "&" dlK_expfor "&" dlK_nexpfor "\\"
di "s.e.&  (" dlK_expdom_se ")&(" dlK_nexpdom_se ")&(" dlK_expfor_se ")&(" dlK_nexpfor_se ")\\"

reg leverage_change no_access_for_credit
mata : Mrounded=round(st_matrix("e(b)"),0.01)
mata : st_matrix("coeff_simple_reg_lev",Mrounded)
scalar lev_dom    = coeff_simple_reg_lev[1,1]+coeff_simple_reg_lev[1,2]
scalar lev_for    = coeff_simple_reg_lev[1,2]
scalar lev_dom_se = round(sqrt(e(V)[1,1]+e(V)[2,2]),0.01)
scalar lev_for_se = round(sqrt(e(V)[2,2]),0.01)
di "$\Delta \frac{Debt}{Asset} $&  \multicolumn{2}{c}{" lev_dom "}&\multicolumn{2}{c}{" lev_for "}\\"
di "s.e.&  \multicolumn{2}{c}{(" lev_dom_se ")}&\multicolumn{2}{c}{(" lev_for_se ")}\\"

reg leverage_change exportdummy_domestic domesticdummy_domestic exportdummy_foreign domesticdummy_foreign, noconstant
mata : Mrounded=round(st_matrix("e(b)"),0.01)
mata : st_matrix("coeff_simple_reg_exporting_lev",Mrounded)
scalar lev_expdom  = coeff_simple_reg_exporting_lev[1,1]
scalar lev_nexpdom = coeff_simple_reg_exporting_lev[1,2]
scalar lev_expfor  = coeff_simple_reg_exporting_lev[1,3]
scalar lev_nexpfor = coeff_simple_reg_exporting_lev[1,4]
scalar lev_expdom_se  = round(sqrt(e(V)[1,1]),0.01)
scalar lev_nexpdom_se = round(sqrt(e(V)[2,2]),0.01)
scalar lev_expfor_se  = round(sqrt(e(V)[3,3]),0.01)
scalar lev_nexpfor_se = round(sqrt(e(V)[4,4]),0.01)
di "$\Delta \frac{Debt}{Asset} $&  " lev_expdom "&" lev_nexpdom "&" lev_expfor "&" lev_nexpfor "\\"
di "s.e.&  (" lev_expdom_se ")&(" lev_nexpdom_se ")&(" lev_expfor_se ")&(" lev_nexpfor_se ")\\"
restore

* --- Table data_transition_smaller (main.tex ~line 1023) ---
di "=== TABLE data_transition_smaller ==="
scalar fo_ww   = distribution_values[1,1]+distribution_values[2,1]+distribution_values[5,1]+distribution_values[6,1]
scalar fo_wfw  = distribution_values[11,1]+distribution_values[12,1]+distribution_values[15,1]+distribution_values[16,1]
scalar fo_went = distribution_values[21,1]+distribution_values[22,1]
scalar fo_fww  = distribution_values[9,1]+distribution_values[10,1]+distribution_values[13,1]+distribution_values[14,1]
scalar fo_fwf  = distribution_values[3,1]+distribution_values[4,1]+distribution_values[7,1]+distribution_values[8,1]
scalar fo_fent = distribution_values[23,1]+distribution_values[24,1]
scalar fo_exitw = distribution_values[17,1]+distribution_values[18,1]
scalar fo_exitf = distribution_values[19,1]+distribution_values[20,1]
scalar ex_nn   = distribution_values[2,1]+distribution_values[4,1]+distribution_values[10,1]+distribution_values[12,1]
scalar ex_yn   = distribution_values[5,1]+distribution_values[7,1]+distribution_values[13,1]+distribution_values[15,1]
scalar ex_nent = distribution_values[22,1]+distribution_values[24,1]
scalar ex_ny   = distribution_values[6,1]+distribution_values[8,1]+distribution_values[14,1]+distribution_values[16,1]
scalar ex_yy   = distribution_values[1,1]+distribution_values[3,1]+distribution_values[9,1]+distribution_values[11,1]
scalar ex_yent = distribution_values[21,1]+distribution_values[23,1]
scalar ex_exitn = distribution_values[18,1]+distribution_values[20,1]
scalar ex_exity = distribution_values[17,1]+distribution_values[19,1]
file open dsfile using "$tables/Table_data_transition_smaller.tex", write replace
file write dsfile "\begin{table}" _n
file write dsfile "\caption{ Transition matrix}" _n
file write dsfile "    \centering" _n
file write dsfile "     \scalebox{0.70}{\begin{tabular}{ccccc |cc ccccc}" _n
file write dsfile "   \multicolumn{2}{c}{Foreign owner} & \multicolumn{2}{c}{2000} & &&   \multicolumn{2}{c}{Export status} & \multicolumn{2}{c}{2000} &\\" _n
file write dsfile " Year  & & without & with & entrants  &&  Year&  & No & Yes & entrants \\" _n
file write dsfile "    \hline" _n
file write dsfile "2008 & without & " %9.0g (fo_ww) " & " %9.0g (fo_wfw) " & " %9.0g (fo_went) " && 2008 & No  & " %9.0g (ex_nn) " & " %9.0g (ex_yn) " & " %9.0g (ex_nent) " \\" _n
file write dsfile " & with & " %9.0g (fo_fww) " & " %9.0g (fo_fwf) " & " %9.0g (fo_fent) "  &&  & Yes & " %9.0g (ex_ny) " & " %9.0g (ex_yy) " & " %9.0g (ex_yent) " \\" _n
file write dsfile "& exit & " %9.0g (fo_exitw) " & " %9.0g (fo_exitf) " & - &&  &exit & " %9.0g (ex_exitn) " & " %9.0g (ex_exity) " & -\\" _n
file write dsfile "\hline" _n
file write dsfile "    \end{tabular}}" _n
file write dsfile " \parbox{0.7\textwidth}{   \caption*{ \scriptsize \textit{Note:} The transition of the number of firms for foreign ownership and export status. For a more detailed transition matrix across each group, see Tables \ref{table:data_transition} and \ref{table:data_transition_probability} }. }" _n
file write dsfile "    \label{table:data_transition_smaller}" _n
file write dsfile "\end{table}" _n
file close dsfile

* --- Table data_regression_simple file export ---
file open dsrfile using "$tables/Table_data_regression_simple.tex", write replace
file write dsrfile "\begin{table}" _n
file write dsrfile "\caption{Capital and leverage growth between 2000 and 2008}" _n
file write dsrfile "\label{table:data_regression_simple}" _n
file write dsrfile "    \centering" _n
file write dsrfile "    \scalebox{0.9}{" _n
file write dsrfile "\def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi}" _n
file write dsrfile "\begin{tabular}{l*{5}{c}}" _n
file write dsrfile "\hline\hline" _n
file write dsrfile "& \multicolumn{2}{c}{Domestic owned} & \multicolumn{2}{c}{Foreign owned}\\" _n
file write dsrfile "\hline" _n
file write dsrfile "\$\Delta \log K \$ & \multicolumn{2}{c}{" %9.0g (dlK_dom) "} & \multicolumn{2}{c}{" %9.0g (dlK_for) "}\\" _n
file write dsrfile "s.e. & \multicolumn{2}{c}{(" %9.0g (dlK_dom_se) ")} & \multicolumn{2}{c}{(" %9.0g (dlK_for_se) ")}\\" _n
file write dsrfile "\$\Delta \frac{Debt}{Asset}\$ & \multicolumn{2}{c}{" %9.0g (lev_dom) "} & \multicolumn{2}{c}{" %9.0g (lev_for) "}\\" _n
file write dsrfile "s.e. & \multicolumn{2}{c}{(" %9.0g (lev_dom_se) ")} & \multicolumn{2}{c}{(" %9.0g (lev_for_se) ")}\\" _n
file write dsrfile " & Exporter & Non-exporter & Exporter & Non-exporter  \\" _n
file write dsrfile "\$\Delta \log K \$ & " %9.0g (dlK_expdom) " & " %9.0g (dlK_nexpdom) " & " %9.0g (dlK_expfor) " & " %9.0g (dlK_nexpfor) " \\" _n
file write dsrfile "s.e. & (" %9.0g (dlK_expdom_se) ") & (" %9.0g (dlK_nexpdom_se) ") & (" %9.0g (dlK_expfor_se) ") & (" %9.0g (dlK_nexpfor_se) ") \\" _n
file write dsrfile "\$\Delta \frac{Debt}{Asset}\$ & " %9.0g (lev_expdom) " & " %9.0g (lev_nexpdom) " & " %9.0g (lev_expfor) " & " %9.0g (lev_nexpfor) " \\" _n
file write dsrfile "s.e. & (" %9.0g (lev_expdom_se) ") & (" %9.0g (lev_nexpdom_se) ") & (" %9.0g (lev_expfor_se) ") & (" %9.0g (lev_nexpfor_se) ") \\" _n
file write dsrfile "\hline" _n
file write dsrfile "\hline" _n
file write dsrfile "\end{tabular}" _n
file write dsrfile "}" _n
file write dsrfile "\parbox{0.65\textwidth}{\caption*{\scriptsize \textit{Note:} Capital and leverage growth for incumbent firms, with or without foreign owners. The averages and standard errors are obtained from a simple OLS on dummy variables. Firms are also split based on export status. Average capital growth was 68 \%, leverage declined by 5.58 p.p.} }" _n
file write dsrfile "\end{table}" _n
file close dsrfile


* ===== PANEL SETUP (full-year panel) =====
xtset sorszam year
gen exportentrants= D.exportdummy
gen no_access_for_credit_change= D.no_access_for_credit
gen capital_growth = D.log_capital
gen value_added_growth = D.log_value_added
gen log_TFPR_growth = D.log_TFPR
gen leverage_change = D.leverage
bysort sorszam:egen flag_exporter_ever = max(exportdummy ==1)
bysort sorszam:egen flag_exporterentrant_ever = max(exportentrants ==1)
bysort sorszam:egen flag_exportexit_ever = max(exportentrants ==-1)
gen exportincumbent = 1 if (flag_exporter_ever ==1 & flag_exporterentrant_ever==0 & flag_exportexit_ever==0 & flag_2000 == 1 & flag_2008 == 1)
replace exportincumbent =0 if exportincumbent==.
gen flag_exporterentrant = 1 if flag_exporterentrant_ever ==1
replace flag_exporterentrant = 1 if (flag_exporterentrant ==. & flag_exporter_ever ==1 & flag_2000 ==0)
replace flag_exporterentrant= 0 if flag_exporterentrant==.
gen flag_exportexit = 1 if flag_exportexit_ever==1
replace flag_exportexit = 1 if (flag_exportexit ==. & flag_exporter_ever ==1  & flag_2008 ==0)
replace flag_exportexit = 0 if flag_exportexit==.
bysort sorszam:egen flag_foreign_owned_ever = max(no_access_for_credit ==0)
bysort sorszam:egen flag_foreign_credit_acquired = max(no_access_for_credit_change ==-1)
bysort sorszam:egen flag_foreign_credit_lost = max(no_access_for_credit_change ==1)
gen no_access_credit_ever = 1 if (flag_foreign_owned_ever ==0 )
replace no_access_credit_ever =0 if no_access_credit_ever==.
gen access_credit_always = 1 if (flag_foreign_owned_ever ==1 & flag_foreign_credit_acquired==0 & flag_foreign_credit_lost==0)
replace access_credit_always =0 if access_credit_always==.

gen exportdummy_foreign = exportdummy if no_access_for_credit==0
replace exportdummy_foreign=0 if exportdummy_foreign==.
gen exportdummy_domestic = exportdummy if no_access_for_credit==1
replace exportdummy_domestic=0 if exportdummy_domestic==.
gen domesticdummy_foreign = 1 if (no_access_for_credit==0 & exportdummy == 0)
replace domesticdummy_foreign=0 if domesticdummy_foreign==.
gen domesticdummy_domestic = 1 if (no_access_for_credit==1 & exportdummy == 0)
replace domesticdummy_domestic=0 if domesticdummy_domestic==.
egen unique_firms = group(sorszam)

gen value_added_millions = value_added/1000
gen capital_millions = capital/1000
gen equity_millions = equity/1000
gen debt_millions = debt/1000

* ===== DiD SETUP =====
gen treated = 1 if (no_access_credit_ever==1 & year>2001)
replace treated = 0 if (no_access_credit_ever==1 & year<2002)
replace treated = 0 if (flag_foreign_credit_lost==1 )
bysort sorszam: egen treated1 = max(treated)
gen reform = 2002 if treated1==1
gen TimeToTreat=year-reform

bys ev: egen med_equity = median(equity)
gen median_equity = equity > med_equity
bys ev: egen median_lev= median(leverage)
gen median_leverage = leverage > median_lev
gen median_equity_2000 = median_equity if year ==2000
bysort sorszam: egen median_equity_2000a = max(median_equity_2000)

bys ev: egen med_equity_export1 = median(equity) if exportincumbent==1
gen median_equity_export1 = equity > med_equity_export1
gen median_equity_2000_export1 = median_equity_export1 if (year ==2000 & exportincumbent==1)
bysort sorszam: egen median_equity_2000a_export1 = max(median_equity_2000_export1)


* ===================================================================
* TABLE data_micro (main.tex ~line 182)
* Descriptive statistics of Hungarian firms 2000-2008
* Exports to Tables/Table_data_micro.tex
* ===================================================================
di "=== TABLE data_micro: Descriptive statistics ==="
preserve
    gen value_added_m = value_added/1000
    gen capital_m     = capital/1000
    gen equity_m      = equity/1000

    * Firm counts per group (collapse to one obs per firm using period_considered==1 or flag_2000)
    * We count each unique firm once using the full multi-year panel
    * Group indicators at firm level (same for all years for a given firm)
    gen g_exit    = (flag_2008==0)
    gen g_entrant = (flag_2000==0)

    * Counts of unique firms and obs per group
    * Use tabstat to get means and se's for each group
    * --- All ---
    tabstat value_added_m capital_m equity_m letszam, statistics(mean semean n) save
    matrix r_all = r(StatTotal)
    * --- Incumbents ---
    tabstat value_added_m capital_m equity_m letszam if flag_incumbent==1, statistics(mean semean n) save
    matrix r_incumb = r(StatTotal)
    * --- Exits ---
    tabstat value_added_m capital_m equity_m letszam if g_exit==1, statistics(mean semean n) save
    matrix r_exit = r(StatTotal)
    * --- Entrants ---
    tabstat value_added_m capital_m equity_m letszam if g_entrant==1, statistics(mean semean n) save
    matrix r_entrant = r(StatTotal)
    * --- Exporters: least once ---
    tabstat value_added_m capital_m equity_m letszam if flag_exporter_ever==1, statistics(mean semean n) save
    matrix r_expleast = r(StatTotal)
    * --- Exporters: always ---
    tabstat value_added_m capital_m equity_m letszam if exportincumbent==1, statistics(mean semean n) save
    matrix r_expalways = r(StatTotal)
    * --- Exporters: entrants ---
    tabstat value_added_m capital_m equity_m letszam if flag_exporterentrant==1, statistics(mean semean n) save
    matrix r_expent = r(StatTotal)
    * --- Exporters: exits ---
    tabstat value_added_m capital_m equity_m letszam if flag_exportexit==1, statistics(mean semean n) save
    matrix r_expexit = r(StatTotal)
    * --- Foreign: least once ---
    tabstat value_added_m capital_m equity_m letszam if flag_foreign_owned_ever==1, statistics(mean semean n) save
    matrix r_forleast = r(StatTotal)
    * --- Foreign: always ---
    tabstat value_added_m capital_m equity_m letszam if access_credit_always==1, statistics(mean semean n) save
    matrix r_foralways = r(StatTotal)
    * --- Foreign: never ---
    tabstat value_added_m capital_m equity_m letszam if flag_foreign_owned_ever==0, statistics(mean semean n) save
    matrix r_fornever = r(StatTotal)

    * Unique firm counts per group (tag first obs per firm, then count)
    bysort sorszam: gen byte _firmtag = (_n==1)
    count if _firmtag==1
    scalar n_all = r(N)
    count if _firmtag==1 & flag_incumbent==1
    scalar n_incumb = r(N)
    count if _firmtag==1 & g_exit==1
    scalar n_exit = r(N)
    count if _firmtag==1 & g_entrant==1
    scalar n_entrant = r(N)
    count if _firmtag==1 & flag_exporter_ever==1
    scalar n_expleast = r(N)
    count if _firmtag==1 & exportincumbent==1
    scalar n_expalways = r(N)
    count if _firmtag==1 & flag_exporterentrant==1
    scalar n_expent = r(N)
    count if _firmtag==1 & flag_exportexit==1
    scalar n_expexit = r(N)
    count if _firmtag==1 & flag_foreign_owned_ever==1
    scalar n_forleast = r(N)
    count if _firmtag==1 & access_credit_always==1
    scalar n_foralways = r(N)
    count if _firmtag==1 & flag_foreign_owned_ever==0
    scalar n_fornever = r(N)
    drop _firmtag

    * Write table
    file open dmfile using "$tables/Table_data_micro.tex", write replace
    file write dmfile "\begin{table}[ht!]" _n
    file write dmfile "\caption{Descriptive statistics of Hungarian firms between 2000-2008}" _n
    file write dmfile " \scalebox{0.6}{" _n
    file write dmfile "\begin{tabular}{lcccc| cccc | ccc}" _n
    file write dmfile "            & \multicolumn{4}{c|}{All firms} & \multicolumn{4}{c|}{Exporters} & \multicolumn{3}{c}{Foreign owned} \\" _n
    file write dmfile "             & All & Incumbents  & Exits  & Entrants & Least once & Always & Entrants  &Exits & Least once & Always & Never\\ \cline{1-12}" _n
    file write dmfile "Firms & " %9.0g (n_all) " & " %9.0g (n_incumb) " & " %9.0g (n_exit) " & " %9.0g (n_entrant) " & " %9.0g (n_expleast) " & " %9.0g (n_expalways) " & " %9.0g (n_expent) " & " %9.0g (n_expexit) " & " %9.0g (n_forleast) " & " %9.0g (n_foralways) " & " %9.0g (n_fornever) "\\" _n
    file write dmfile "\% of total: & 100 & " %9.0g (round(100*n_incumb/n_all)) " & " %9.0g (round(100*n_exit/n_all)) " & " %9.0g (round(100*n_entrant/n_all)) " & " %9.0g (round(100*n_expleast/n_all)) " & " %9.0g (round(100*n_expalways/n_all)) " & " %9.0g (round(100*n_expent/n_all)) " & " %9.0g (round(100*n_expexit/n_all)) " & " %9.0g (round(100*n_forleast/n_all)) " & " %9.0g (round(100*n_foralways/n_all)) " & " %9.0g (round(100*n_fornever/n_all)) "\\" _n
    file write dmfile "\% of exporters: & - & - & - & - & 100 & " %9.0g (round(100*n_expalways/n_expleast)) " & " %9.0g (round(100*n_expent/n_expleast)) " & " %9.0g (round(100*n_expexit/n_expleast)) " & - & - & \\" _n
    file write dmfile "Observations & " %9.0g (r_all[3,1]) " & " %9.0g (r_incumb[3,1]) " & " %9.0g (r_exit[3,1]) " & " %9.0g (r_entrant[3,1]) " & " %9.0g (r_expleast[3,1]) " & " %9.0g (r_expalways[3,1]) " & " %9.0g (r_expent[3,1]) " & " %9.0g (r_expexit[3,1]) " & " %9.0g (r_forleast[3,1]) " & " %9.0g (r_foralways[3,1]) " & " %9.0g (r_fornever[3,1]) "\\" _n
    file write dmfile "Value added & & & & & & & & & & &\\" _n
    file write dmfile "\hspace{3mm} Mean & " %9.0g (round(r_all[1,1])) " & " %9.0g (round(r_incumb[1,1])) " & " %9.0g (round(r_exit[1,1])) " & " %9.0g (round(r_entrant[1,1])) " & " %9.0g (round(r_expleast[1,1])) " & " %9.0g (round(r_expalways[1,1])) " & " %9.0g (round(r_expent[1,1])) " & " %9.0g (round(r_expexit[1,1])) " & " %9.0g (round(r_forleast[1,1])) " & " %9.0g (round(r_foralways[1,1])) " & " %9.0g (round(r_fornever[1,1])) "\\" _n
    file write dmfile "\hspace{3mm} s.d. & " %9.0g (round(r_all[2,1])) " & " %9.0g (round(r_incumb[2,1])) " & " %9.0g (round(r_exit[2,1])) " & " %9.0g (round(r_entrant[2,1])) " & " %9.0g (round(r_expleast[2,1])) " & " %9.0g (round(r_expalways[2,1])) " & " %9.0g (round(r_expent[2,1])) " & " %9.0g (round(r_expexit[2,1])) " & " %9.0g (round(r_forleast[2,1])) " & " %9.0g (round(r_foralways[2,1])) " & " %9.0g (round(r_fornever[2,1])) "\\" _n
    file write dmfile "Capital & & & & & & & & & & &\\" _n
    file write dmfile "\hspace{3mm} Mean & " %9.0g (round(r_all[1,2])) " & " %9.0g (round(r_incumb[1,2])) " & " %9.0g (round(r_exit[1,2])) " & " %9.0g (round(r_entrant[1,2])) " & " %9.0g (round(r_expleast[1,2])) " & " %9.0g (round(r_expalways[1,2])) " & " %9.0g (round(r_expent[1,2])) " & " %9.0g (round(r_expexit[1,2])) " & " %9.0g (round(r_forleast[1,2])) " & " %9.0g (round(r_foralways[1,2])) " & " %9.0g (round(r_fornever[1,2])) "\\" _n
    file write dmfile "\hspace{3mm} s.d. & " %9.0g (round(r_all[2,2])) " & " %9.0g (round(r_incumb[2,2])) " & " %9.0g (round(r_exit[2,2])) " & " %9.0g (round(r_entrant[2,2])) " & " %9.0g (round(r_expleast[2,2])) " & " %9.0g (round(r_expalways[2,2])) " & " %9.0g (round(r_expent[2,2])) " & " %9.0g (round(r_expexit[2,2])) " & " %9.0g (round(r_forleast[2,2])) " & " %9.0g (round(r_foralways[2,2])) " & " %9.0g (round(r_fornever[2,2])) "\\" _n
    file write dmfile "Equity & & & & & & & & & & &\\" _n
    file write dmfile "\hspace{3mm} Mean & " %9.0g (round(r_all[1,3])) " & " %9.0g (round(r_incumb[1,3])) " & " %9.0g (round(r_exit[1,3])) " & " %9.0g (round(r_entrant[1,3])) " & " %9.0g (round(r_expleast[1,3])) " & " %9.0g (round(r_expalways[1,3])) " & " %9.0g (round(r_expent[1,3])) " & " %9.0g (round(r_expexit[1,3])) " & " %9.0g (round(r_forleast[1,3])) " & " %9.0g (round(r_foralways[1,3])) " & " %9.0g (round(r_fornever[1,3])) "\\" _n
    file write dmfile "\hspace{3mm} s.d. & " %9.0g (round(r_all[2,3])) " & " %9.0g (round(r_incumb[2,3])) " & " %9.0g (round(r_exit[2,3])) " & " %9.0g (round(r_entrant[2,3])) " & " %9.0g (round(r_expleast[2,3])) " & " %9.0g (round(r_expalways[2,3])) " & " %9.0g (round(r_expent[2,3])) " & " %9.0g (round(r_expexit[2,3])) " & " %9.0g (round(r_forleast[2,3])) " & " %9.0g (round(r_foralways[2,3])) " & " %9.0g (round(r_fornever[2,3])) "\\" _n
    file write dmfile "Employment & & & & & & & & & & &\\" _n
    file write dmfile "\hspace{3mm} Mean & " %9.0g (round(r_all[1,4])) " & " %9.0g (round(r_incumb[1,4])) " & " %9.0g (round(r_exit[1,4])) " & " %9.0g (round(r_entrant[1,4])) " & " %9.0g (round(r_expleast[1,4])) " & " %9.0g (round(r_expalways[1,4])) " & " %9.0g (round(r_expent[1,4])) " & " %9.0g (round(r_expexit[1,4])) " & " %9.0g (round(r_forleast[1,4])) " & " %9.0g (round(r_foralways[1,4])) " & " %9.0g (round(r_fornever[1,4])) "\\" _n
    file write dmfile "\hspace{3mm} s.d. & " %9.0g (round(r_all[2,4])) " & " %9.0g (round(r_incumb[2,4])) " & " %9.0g (round(r_exit[2,4])) " & " %9.0g (round(r_entrant[2,4])) " & " %9.0g (round(r_expleast[2,4])) " & " %9.0g (round(r_expalways[2,4])) " & " %9.0g (round(r_expent[2,4])) " & " %9.0g (round(r_expexit[2,4])) " & " %9.0g (round(r_forleast[2,4])) " & " %9.0g (round(r_foralways[2,4])) " & " %9.0g (round(r_fornever[2,4])) "\\" _n
    file write dmfile "\cline{1-12}" _n
    file write dmfile `"\end{tabular}} \caption*{ \scriptsize \textit{Note:} Export and foreign-owned status are based on export sales and on foreign ownership of equity, respectively. Quantities displayed are means in millions of Hungarian Forints (HUF). "Entrants" and "Exits" describe firms that entered or left the sample. "Always" refers to incumbent firms that always exported or had foreign owners between 2000 and 2008. Appendix \ref{appendix:firm_data} contains additional details, in particular, the transition matrix and the overlap across these groups.}"' _n
    file write dmfile "\label{table:data_micro}" _n
    file write dmfile "\vspace{-.6cm}" _n
    file write dmfile "\end{table}" _n
    file close dmfile
restore

* ===================================================================
* FIGURE overlap (main.tex Appendix ~line 1027)
* DiD_capital_overlap.pdf, DiD_value_added_overlap.pdf,
* DiD_equity_overlap.pdf, DiD_employment_overlap.pdf
* ===================================================================
su log_capital if year ==2000 & treated1 == 0
su log_capital if year ==2000 & treated1 == 1
local start = 1
local bin = 1
local opts start(`start') width(`bin')
local normal_begin = 1
local normal_end = 20
local range ra(`normal_begin' `normal_end')
local colors red blue black
local call
levelsof treated1
tokenize "`r(levels)'"
local nlevels : word count `r(levels)'
forval j = 1/`nlevels' {
    su log_capital if year ==2000 & treated1 == ``j''
    scalar mu`j' = r(mean)
    scalar sd`j' = r(sd)
    local color : word `j' of `colors'
    local call `call' || histogram log_capital if year ==2000 & treated1 == ``j'', bcolor(`color'%20) `opts'
    local call `call' || function normalden(x, mu`j', sd`j') , `range' lcolor(`color')
}
twoway `call'  xla(`normal_begin'(1)`normal_end')  legend(order(3 "Treated" 1 "Control") pos(12))  ytitle(Density) xtitle("log(capital)")
graph export "$figures/DiD_capital_overlap.pdf", as(pdf) name("Graph") replace

local start = 1
local bin = 1
local opts start(`start') width(`bin')
local normal_begin = 1
local normal_end = 20
local range ra(`normal_begin' `normal_end')
local colors red blue black
local call
levelsof treated1
tokenize "`r(levels)'"
local nlevels : word count `r(levels)'
forval j = 1/`nlevels' {
    su log_value_added if year ==2000 & treated1 == ``j''
    scalar mu`j' = r(mean)
    scalar sd`j' = r(sd)
    local color : word `j' of `colors'
    local call `call' || histogram log_value_added if year ==2000 & treated1 == ``j'', bcolor(`color'%20) `opts'
    local call `call' || function normalden(x, mu`j', sd`j') , `range' lcolor(`color')
}
twoway `call'  xla(`normal_begin'(1)`normal_end')  legend(order(3 "Treated" 1 "Control") pos(12))  ytitle(Density) xtitle("log(value_added)")
graph export "$figures/DiD_value_added_overlap.pdf", as(pdf) name("Graph") replace

su equity if year ==2000 & treated1 == 0
su equity if year ==2000 & treated1 == 1
local start =  0.0
local bin = 1
local opts start(`start') width(`bin')
local normal_begin = 1
local normal_end = 20
local range ra(`normal_begin' `normal_end')
local colors red blue black
local call
levelsof treated1
tokenize "`r(levels)'"
local nlevels : word count `r(levels)'
forval j = 1/`nlevels' {
    su log_equity if year ==2000 & treated1 == ``j''
    scalar mu`j' = r(mean)
    scalar sd`j' = r(sd)
    local color : word `j' of `colors'
    local call `call' || histogram log_equity if year ==2000 & treated1 == ``j'', bcolor(`color'%20) `opts'
    local call `call' || function normalden(x, mu`j', sd`j') , `range' lcolor(`color')
}
twoway `call'  xla(`normal_begin'(1)`normal_end')  legend(order(3 "Treated" 1 "Control") pos(12))  ytitle(Density) xtitle("log(equity)")
graph export "$figures/DiD_equity_overlap.pdf", as(pdf) name("Graph") replace

su letszam if year ==2000 & treated1 == 0
su letszam if year ==2000 & treated1 == 1
local start = 0.68
local bin = 1
local opts start(`start') width(`bin')
local normal_begin = 1
local normal_end = 11
local range ra(`normal_begin' `normal_end')
local colors red blue black
local call
levelsof treated1
tokenize "`r(levels)'"
local nlevels : word count `r(levels)'
forval j = 1/`nlevels' {
    su log_letszam if year ==2000 & treated1 == ``j''
    scalar mu`j' = r(mean)
    scalar sd`j' = r(sd)
    local color : word `j' of `colors'
    local call `call' || histogram log_letszam if year ==2000 & treated1 == ``j'', bcolor(`color'%20) `opts'
    local call `call' || function normalden(x, mu`j', sd`j') , `range' lcolor(`color')
}
twoway `call'  xla(`normal_begin'(1)`normal_end')  legend(order(3 "Treated" 1 "Control") pos(12))  ytitle(Density) xtitle("log(employment)")
graph export "$figures/DiD_employment_overlap.pdf", as(pdf) name("Graph") replace


* ===================================================================
* DiD REGRESSIONS — produces:
*  Figure DiD_combined (main.tex Section 2, ~line 219):
*    DiD_capital_incumbents_ARPL.eps, DiD_capital_exporters_ARPL.eps
*  Figure DiD_combined_noARPL (main.tex Appendix ~line 1097):
*    DiD_capital_incumbents.eps, DiD_capital_exporters.eps
*  Table data_regression_did (main.tex Appendix ~line 1105): collect layout
* ===================================================================
collect create ex6, replace
preserve
keep if(flag_incumbent==1)
collect _r_b _r_se, tag(model[(1)]): eventdd log_capital log_ARPL, timevar(TimeToTreat) ci(rcap) method(hdfe, cluster(sorszam) absorb( no_access_credit_ever alapitas regiokod ind10 year )) graph_op( title("") scheme(s1mono) xlabel(-2 "2000" -1 "2001" 0 "2002" 1 "2003" 2 "2004" 3 "2005" 4 "2006" 5 "2007" 6 "2008") ylabel(-0.1 "-10" 0 "0" 0.1 "10" 0.2 "20"))
matrix results_incumbent_ARPL = e(leads) \ e(lags)
scalar arpl_b_1 = _b[log_ARPL]
scalar arpl_se_1 = _se[log_ARPL]
scalar eN_1 = e(N)
collect p_d=r(p), tag(model[(1)]): estat leads
scalar p_1 = r(p)
collect _r_b _r_se, tag(model[(2)]): eventdd log_capital, timevar(TimeToTreat) ci(rcap) method(hdfe, cluster(sorszam) absorb( no_access_credit_ever alapitas regiokod ind10 year )) graph_op( title("") scheme(s1mono) xlabel(-2 "2000" -1 "2001" 0 "2002" 1 "2003" 2 "2004" 3 "2005" 4 "2006" 5 "2007" 6 "2008") ylabel(-0.1 "-10" 0 "0" 0.1 "10" 0.2 "20"))
scalar eN_2 = e(N)
collect p_d=r(p), tag(model[(2)]): estat leads
scalar p_2 = r(p)
matrix results_incumbent = e(leads) \ e(lags)
restore

svmat results_incumbent_ARPL, names(event_incumbent_ARPL)
svmat results_incumbent, names(event_incumbent)

preserve
keep if(flag_incumbent==1 & flag_exporter_ever==0)
collect _r_b _r_se, tag(model[(3)]):  eventdd log_capital log_ARPL, timevar(TimeToTreat) ci(rcap) method(hdfe, cluster(sorszam) absorb( no_access_credit_ever alapitas regiokod ind10 year )) graph_op( title("") scheme(s1mono) xlabel(-2 "2000" -1 "2001" 0 "2002" 1 "2003" 2 "2004" 3 "2005" 4 "2006" 5 "2007" 6 "2008") ylabel(-0.1 "-10" 0 "0" 0.1 "10" 0.2 "20"))
matrix results_nexp_ARPL = e(leads) \ e(lags)
scalar arpl_b_3 = _b[log_ARPL]
scalar arpl_se_3 = _se[log_ARPL]
scalar eN_3 = e(N)
collect p_d=r(p), tag(model[(3)]): estat leads
scalar p_3 = r(p)
collect _r_b _r_se, tag(model[(4)]):  eventdd log_capital, timevar(TimeToTreat) ci(rcap) method(hdfe, cluster(sorszam) absorb( no_access_credit_ever alapitas regiokod ind10 year )) graph_op( title("") scheme(s1mono) xlabel(-2 "2000" -1 "2001" 0 "2002" 1 "2003" 2 "2004" 3 "2005" 4 "2006" 5 "2007" 6 "2008") ylabel(-0.1 "-10" 0 "0" 0.1 "10" 0.2 "20"))
scalar eN_4 = e(N)
collect p_d=r(p), tag(model[(4)]): estat leads
scalar p_4 = r(p)
matrix results_nexp = e(leads) \ e(lags)
restore

svmat results_nexp_ARPL, names(event_nexp_ARPL)
svmat results_nexp, names(event_nexp)

preserve
keep if(exportincumbent==1)
collect _r_b _r_se, tag(model[(5)]):  eventdd log_capital log_ARPL, timevar(TimeToTreat) ci(rcap) method(hdfe, cluster(sorszam) absorb( no_access_credit_ever alapitas regiokod ind10 year )) graph_op( title("") scheme(s1mono) xlabel(-2 "2000" -1 "2001" 0 "2002" 1 "2003" 2 "2004" 3 "2005" 4 "2006" 5 "2007" 6 "2008") ylabel(-0.1 "-10" 0 "0" 0.1 "10" 0.2 "20"))
matrix results_exp_ARPL = e(leads) \ e(lags)
scalar arpl_b_5 = _b[log_ARPL]
scalar arpl_se_5 = _se[log_ARPL]
scalar eN_5 = e(N)
collect p_d=r(p), tag(model[(5)]): estat leads
scalar p_5 = r(p)
collect _r_b _r_se, tag(model[(6)]): eventdd log_capital, timevar(TimeToTreat) ci(rcap) method(hdfe, cluster(sorszam) absorb( no_access_credit_ever alapitas regiokod ind10 year )) graph_op( title("") scheme(s1mono) xlabel(-2 "2000" -1 "2001" 0 "2002" 1 "2003" 2 "2004" 3 "2005" 4 "2006" 5 "2007" 6 "2008") ylabel(-0.1 "-10" 0 "0" 0.1 "10" 0.2 "20"))
scalar eN_6 = e(N)
collect p_d=r(p), tag(model[(6)]): estat leads
scalar p_6 = r(p)
matrix results_exp = e(leads) \ e(lags)
restore
svmat results_exp, names(event_exp)
svmat results_exp_ARPL, names(event_exp_ARPL)

preserve
keep if(exportincumbent==1 & median_equity_2000a==0)
collect _r_b _r_se, tag(model[(7)]):  eventdd log_capital log_ARPL, timevar(TimeToTreat) ci(rcap) method(hdfe, cluster(sorszam) absorb( no_access_credit_ever alapitas regiokod ind10 year )) graph_op( title("") scheme(s1mono) xlabel(-2 "2000" -1 "2001" 0 "2002" 1 "2003" 2 "2004" 3 "2005" 4 "2006" 5 "2007" 6 "2008") ylabel(-0.1 "-10" 0 "0" 0.1 "10" 0.2 "20"))
matrix results_exp_belowmed_ARPL = e(leads) \ e(lags)
scalar arpl_b_7 = _b[log_ARPL]
scalar arpl_se_7 = _se[log_ARPL]
scalar eN_7 = e(N)
collect p_d=r(p), tag(model[(7)]): estat leads
scalar p_7 = r(p)
collect _r_b _r_se, tag(model[(8)]): eventdd log_capital, timevar(TimeToTreat) ci(rcap) method(hdfe, cluster(sorszam) absorb( no_access_credit_ever alapitas regiokod ind10 year )) graph_op( title("") scheme(s1mono) xlabel(-2 "2000" -1 "2001" 0 "2002" 1 "2003" 2 "2004" 3 "2005" 4 "2006" 5 "2007" 6 "2008") ylabel(-0.1 "-10" 0 "0" 0.1 "10" 0.2 "20"))
scalar eN_8 = e(N)
collect p_d=r(p), tag(model[(8)]): estat leads
scalar p_8 = r(p)
matrix results_exp_belowmed = e(leads) \ e(lags)
restore
svmat results_exp_belowmed_ARPL, names(event_exp_belowmed_ARPL)
svmat results_exp_belowmed, names(event_exp_belowmed)

preserve
keep if(exportincumbent==1 & median_equity_2000a==1)
collect _r_b _r_se, tag(model[(9)]):  eventdd log_capital log_ARPL, timevar(TimeToTreat) ci(rcap) method(hdfe, cluster(sorszam) absorb( no_access_credit_ever alapitas regiokod ind10 year )) graph_op( title("") scheme(s1mono) xlabel(-2 "2000" -1 "2001" 0 "2002" 1 "2003" 2 "2004" 3 "2005" 4 "2006" 5 "2007" 6 "2008") ylabel(-0.1 "-10" 0 "0" 0.1 "10" 0.2 "20"))
matrix results_exp_abovemed_ARPL = e(leads) \ e(lags)
scalar arpl_b_9 = _b[log_ARPL]
scalar arpl_se_9 = _se[log_ARPL]
scalar eN_9 = e(N)
collect p_d=r(p), tag(model[(9)]): estat leads
scalar p_9 = r(p)
collect _r_b _r_se, tag(model[(10)]): eventdd log_capital, timevar(TimeToTreat) ci(rcap) method(hdfe, cluster(sorszam) absorb( no_access_credit_ever alapitas regiokod ind10 year )) graph_op( title("") scheme(s1mono) xlabel(-2 "2000" -1 "2001" 0 "2002" 1 "2003" 2 "2004" 3 "2005" 4 "2006" 5 "2007" 6 "2008") ylabel(-0.1 "-10" 0 "0" 0.1 "10" 0.2 "20"))
scalar eN_10 = e(N)
collect p_d=r(p), tag(model[(10)]): estat leads
scalar p_10 = r(p)
matrix results_exp_abovemed = e(leads) \ e(lags)
restore
svmat results_exp_abovemed_ARPL, names(event_exp_abovemed_ARPL)
svmat results_exp_abovemed, names(event_exp_abovemed)

* Set event year axis variables
replace event_exp_abovemed1 = 2001 in 1
replace event_exp_abovemed1 = 2000 in 2
replace event_exp_abovemed1 = 2002 in 3
replace event_exp_abovemed1 = 2003 in 4
replace event_exp_abovemed1 = 2004 in 5
replace event_exp_abovemed1 = 2005 in 6
replace event_exp_abovemed1 = 2006 in 7
replace event_exp_abovemed1 = 2007 in 8
replace event_exp_abovemed1 = 2008 in 9

replace event_exp_belowmed1 = 2001.0 in 1
replace event_exp_belowmed1 = 2000.0 in 2
replace event_exp_belowmed1 = 2002.0 in 3
replace event_exp_belowmed1 = 2003.0 in 4
replace event_exp_belowmed1 = 2004.0 in 5
replace event_exp_belowmed1 = 2005.0 in 6
replace event_exp_belowmed1 = 2006.0 in 7
replace event_exp_belowmed1 = 2007.0 in 8
replace event_exp_belowmed1 = 2008.0 in 9

replace event_incumbent1 = 2001.0 in 1
replace event_incumbent1 = 2000.0 in 2
replace event_incumbent1 = 2002.0 in 3
replace event_incumbent1 = 2003.0 in 4
replace event_incumbent1 = 2004.0 in 5
replace event_incumbent1 = 2005.0 in 6
replace event_incumbent1 = 2006.0 in 7
replace event_incumbent1 = 2007.0 in 8
replace event_incumbent1 = 2008.0 in 9

replace event_exp_abovemed_ARPL1 = 2001 in 1
replace event_exp_abovemed_ARPL1 = 2000 in 2
replace event_exp_abovemed_ARPL1 = 2002 in 3
replace event_exp_abovemed_ARPL1 = 2003 in 4
replace event_exp_abovemed_ARPL1 = 2004 in 5
replace event_exp_abovemed_ARPL1 = 2005 in 6
replace event_exp_abovemed_ARPL1 = 2006 in 7
replace event_exp_abovemed_ARPL1 = 2007 in 8
replace event_exp_abovemed_ARPL1 = 2008 in 9

replace event_exp_belowmed_ARPL1 = 2001.0 in 1
replace event_exp_belowmed_ARPL1 = 2000.0 in 2
replace event_exp_belowmed_ARPL1 = 2002.0 in 3
replace event_exp_belowmed_ARPL1 = 2003.0 in 4
replace event_exp_belowmed_ARPL1 = 2004.0 in 5
replace event_exp_belowmed_ARPL1 = 2005.0 in 6
replace event_exp_belowmed_ARPL1 = 2006.0 in 7
replace event_exp_belowmed_ARPL1 = 2007.0 in 8
replace event_exp_belowmed_ARPL1 = 2008.0 in 9

replace event_incumbent_ARPL1 = 2001.0 in 1
replace event_incumbent_ARPL1 = 2000.0 in 2
replace event_incumbent_ARPL1 = 2002.0 in 3
replace event_incumbent_ARPL1 = 2003.0 in 4
replace event_incumbent_ARPL1 = 2004.0 in 5
replace event_incumbent_ARPL1 = 2005.0 in 6
replace event_incumbent_ARPL1 = 2006.0 in 7
replace event_incumbent_ARPL1 = 2007.0 in 8
replace event_incumbent_ARPL1 = 2008.0 in 9

sort event_exp_abovemed1
graph set window fontface "Times New Roman"

* --- Figure DiD_combined_noARPL (main.tex Appendix ~line 1097) ---
twoway rarea event_incumbent4  event_incumbent2 event_exp_abovemed1, color(gs15) || line event_incumbent3 event_exp_abovemed1, mc(black) ms(S)  lwidth(1.5)  ///
|| line event_exp3 event_exp_belowmed1, mc(red) ms(S) lpattern(dash) lwidth(1.5)    ///
|| line event_nexp3 event_exp_abovemed1,  mc(black) lpattern(dot) lwidth(2.5) ///
legend(order(2 "All" 3 "Exporters" 4 "Non-exporters ") row(1)) ylabel(-0.2 "-20%" -0.1 "-10%" 0 "0%" 0.1 "10%" 0.2 "20%" 0.3 "30%" 0.4 "40%" 0.5 "50%") xtitle("") xla(2000/2008) scheme(s1mono) yscale(range (-0.2,0.5)) ///
xli(2001, lc(gs8) lw(thin)) xtitle("") ytitle("log(K)") yla(, ang(h))
graph export "$figures/DiD_capital_incumbents.eps", as(eps) name("Graph") preview(off) replace

twoway rarea  event_exp4  event_exp2 event_exp_abovemed1, color(gs15) ||line event_exp3 event_exp_belowmed1, mc(black) lpattern(dash) lwidth(1.5)  || line event_exp_abovemed3 event_exp_abovemed1, mc(red) ms(S) lwidth(1.5)    ///
 || line event_exp_belowmed3 event_exp_abovemed1, mc(black)  lpattern(dot) lwidth(2.5) ///
legend(order(2 "All" 3 "Above median equity" 4 "Below median equity" ) row(1))  ylabel(-0.2 "-20%" -0.1 "-10%" 0 "0%" 0.1 "10%" 0.2 "20%" 0.3 "30%" 0.4 "40%" 0.5 "50%") xtitle("") xla(2000/2008) scheme(s1mono) yscale(range (-0.2,0.5)) ///
xli(2001, lc(gs8) lw(thin)) xtitle("") ytitle("log(K)") yla(, ang(h))
graph export "$figures/DiD_capital_exporters.eps", as(eps) name("Graph") preview(off) replace

* --- Figure DiD_combined (main.tex Section 2, ~line 219) ---
twoway rarea event_incumbent_ARPL4  event_incumbent_ARPL2 event_incumbent_ARPL1, color(gs15) || line event_incumbent_ARPL3 event_incumbent_ARPL1, mc(black) ms(S) lpattern(dash) lwidth(1.5)  ///
|| line event_exp_ARPL3 event_incumbent_ARPL1, mc(red) ms(S) lwidth(1.5)    ///
|| line event_nexp_ARPL3 event_incumbent_ARPL1,  mc(black) lpattern(dot) lwidth(2.5) ///
legend(order(2 "All" 3 "Exporters" 4 "Non-exporters ") row(1)) ylabel(-0.2 "-20%" -0.1 "-10%" 0 "0%" 0.1 "10%" 0.2 "20%" 0.3 "30%" 0.4 "40%" 0.5 "50%") xtitle("") xla(2000/2008) scheme(s1mono) yscale(range (-0.2,0.5)) ///
xli(2001, lc(gs8) lw(thin)) xtitle("") ytitle("log(K)") yla(, ang(h))
graph export "$figures/DiD_capital_incumbents_ARPL.eps", as(eps) name("Graph") preview(off) replace

twoway rarea  event_exp_ARPL4  event_exp_ARPL2 event_exp_abovemed_ARPL1, color(gs15) ||line event_exp_ARPL3 event_exp_belowmed_ARPL1, mc(black)  lwidth(1.5)  || line event_exp_abovemed_ARPL3 event_exp_abovemed_ARPL1, mc(red) lpattern(dash) ms(S) lwidth(1.5)    ///
 || line event_exp_belowmed_ARPL3 event_exp_abovemed_ARPL1, mc(black)  lpattern(dot) lwidth(2.5) ///
legend(order(2 "All" 3 "Above median equity" 4 "Below median equity" ) row(1))  ylabel(-0.2 "-20%" -0.1 "-10%" 0 "0%" 0.1 "10%" 0.2 "20%" 0.3 "30%" 0.4 "40%" 0.5 "50%") xtitle("") xla(2000/2008) scheme(s1mono) yscale(range (-0.2,0.5)) ///
xli(2001, lc(gs8) lw(thin)) xtitle("") ytitle("log(K)") yla(, ang(h))
graph export "$figures/DiD_capital_exporters_ARPL.eps", as(eps) name("Graph") preview(off) replace

* --- Table data_regression_did (main.tex Appendix ~line 1105) ---
collect style cell, nformat(%5.2f)
collect style cell result[_r_se], sformat("(%s)")
collect layout (colname#result p_d) (model)
* Write table using stored scalars and matrices (col3=coef, col2=lo CI, col4=hi CI)
* SE derived from 95% CI: se = (upper - lower) / 3.92
file open didfile using "$tables/Table_data_regression_did.tex", write replace
file write didfile "    \begin{table}" _n
file write didfile "    \caption{Capital growth difference in difference}" _n
file write didfile "    \label{table:data_regression_did}" _n
file write didfile "        \centering" _n
file write didfile "\scalebox{0.5}{" _n
file write didfile "\def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi}" _n
file write didfile "\begin{tabular}{lcccc | cccccc}" _n
file write didfile "\hline\hline" _n
file write didfile "  &\multicolumn{4}{c}{Incumbents}&\multicolumn{6}{c}{Export incumbents} \\" _n
file write didfile "    &\multicolumn{2}{c}{All}&\multicolumn{2}{c}{Never-exporters} &\multicolumn{2}{c}{All}&\multicolumn{2}{c}{Below median equity} &\multicolumn{2}{c}{Above median wealth}\\" _n
file write didfile "\hline" _n
file write didfile "\$\beta_{4,ARPL} \$ & " %5.2f (arpl_b_1) " & - & " %5.2f (arpl_b_3) " & - & " %5.2f (arpl_b_5) " & - & " %5.2f (arpl_b_7) " & - & " %5.2f (arpl_b_9) " & - \\" _n
file write didfile "& (" %5.2f (arpl_se_1) ") & & (" %5.2f (arpl_se_3) ") & & (" %5.2f (arpl_se_5) ") & & (" %5.2f (arpl_se_7) ") & & (" %5.2f (arpl_se_9) ") & \\" _n
file write didfile "\$\beta_{3, t= 2000} \$& " %5.2f (results_incumbent_ARPL[2,3]) " & " %5.2f (results_incumbent[2,3]) " & " %5.2f (results_nexp_ARPL[2,3]) " & " %5.2f (results_nexp[2,3]) " & " %5.2f (results_exp_ARPL[2,3]) " & " %5.2f (results_exp[2,3]) " & " %5.2f (results_exp_belowmed_ARPL[2,3]) " & " %5.2f (results_exp_belowmed[2,3]) " & " %5.2f (results_exp_abovemed_ARPL[2,3]) " & " %5.2f (results_exp_abovemed[2,3]) " \\" _n
file write didfile "& (" %5.2f ((results_incumbent_ARPL[2,4]-results_incumbent_ARPL[2,2])/3.92) ") & (" %5.2f ((results_incumbent[2,4]-results_incumbent[2,2])/3.92) ") & (" %5.2f ((results_nexp_ARPL[2,4]-results_nexp_ARPL[2,2])/3.92) ") & (" %5.2f ((results_nexp[2,4]-results_nexp[2,2])/3.92) ") & (" %5.2f ((results_exp_ARPL[2,4]-results_exp_ARPL[2,2])/3.92) ") & (" %5.2f ((results_exp[2,4]-results_exp[2,2])/3.92) ") & (" %5.2f ((results_exp_belowmed_ARPL[2,4]-results_exp_belowmed_ARPL[2,2])/3.92) ") & (" %5.2f ((results_exp_belowmed[2,4]-results_exp_belowmed[2,2])/3.92) ") & (" %5.2f ((results_exp_abovemed_ARPL[2,4]-results_exp_abovemed_ARPL[2,2])/3.92) ") & (" %5.2f ((results_exp_abovemed[2,4]-results_exp_abovemed[2,2])/3.92) ") \\" _n
file write didfile "[1em]" _n
file write didfile "\$\beta_{3, t= 2002} \$& " %5.2f (results_incumbent_ARPL[3,3]) " & " %5.2f (results_incumbent[3,3]) " & " %5.2f (results_nexp_ARPL[3,3]) " & " %5.2f (results_nexp[3,3]) " & " %5.2f (results_exp_ARPL[3,3]) " & " %5.2f (results_exp[3,3]) " & " %5.2f (results_exp_belowmed_ARPL[3,3]) " & " %5.2f (results_exp_belowmed[3,3]) " & " %5.2f (results_exp_abovemed_ARPL[3,3]) " & " %5.2f (results_exp_abovemed[3,3]) " \\" _n
file write didfile "& (" %5.2f ((results_incumbent_ARPL[3,4]-results_incumbent_ARPL[3,2])/3.92) ") & (" %5.2f ((results_incumbent[3,4]-results_incumbent[3,2])/3.92) ") & (" %5.2f ((results_nexp_ARPL[3,4]-results_nexp_ARPL[3,2])/3.92) ") & (" %5.2f ((results_nexp[3,4]-results_nexp[3,2])/3.92) ") & (" %5.2f ((results_exp_ARPL[3,4]-results_exp_ARPL[3,2])/3.92) ") & (" %5.2f ((results_exp[3,4]-results_exp[3,2])/3.92) ") & (" %5.2f ((results_exp_belowmed_ARPL[3,4]-results_exp_belowmed_ARPL[3,2])/3.92) ") & (" %5.2f ((results_exp_belowmed[3,4]-results_exp_belowmed[3,2])/3.92) ") & (" %5.2f ((results_exp_abovemed_ARPL[3,4]-results_exp_abovemed_ARPL[3,2])/3.92) ") & (" %5.2f ((results_exp_abovemed[3,4]-results_exp_abovemed[3,2])/3.92) ") \\" _n
file write didfile "[1em]" _n
file write didfile "\$\beta_{3, t= 2003} \$& " %5.2f (results_incumbent_ARPL[4,3]) " & " %5.2f (results_incumbent[4,3]) " & " %5.2f (results_nexp_ARPL[4,3]) " & " %5.2f (results_nexp[4,3]) " & " %5.2f (results_exp_ARPL[4,3]) " & " %5.2f (results_exp[4,3]) " & " %5.2f (results_exp_belowmed_ARPL[4,3]) " & " %5.2f (results_exp_belowmed[4,3]) " & " %5.2f (results_exp_abovemed_ARPL[4,3]) " & " %5.2f (results_exp_abovemed[4,3]) " \\" _n
file write didfile "& (" %5.2f ((results_incumbent_ARPL[4,4]-results_incumbent_ARPL[4,2])/3.92) ") & (" %5.2f ((results_incumbent[4,4]-results_incumbent[4,2])/3.92) ") & (" %5.2f ((results_nexp_ARPL[4,4]-results_nexp_ARPL[4,2])/3.92) ") & (" %5.2f ((results_nexp[4,4]-results_nexp[4,2])/3.92) ") & (" %5.2f ((results_exp_ARPL[4,4]-results_exp_ARPL[4,2])/3.92) ") & (" %5.2f ((results_exp[4,4]-results_exp[4,2])/3.92) ") & (" %5.2f ((results_exp_belowmed_ARPL[4,4]-results_exp_belowmed_ARPL[4,2])/3.92) ") & (" %5.2f ((results_exp_belowmed[4,4]-results_exp_belowmed[4,2])/3.92) ") & (" %5.2f ((results_exp_abovemed_ARPL[4,4]-results_exp_abovemed_ARPL[4,2])/3.92) ") & (" %5.2f ((results_exp_abovemed[4,4]-results_exp_abovemed[4,2])/3.92) ") \\" _n
file write didfile "[1em]" _n
file write didfile "\$\beta_{3, t= 2004} \$& " %5.2f (results_incumbent_ARPL[5,3]) " & " %5.2f (results_incumbent[5,3]) " & " %5.2f (results_nexp_ARPL[5,3]) " & " %5.2f (results_nexp[5,3]) " & " %5.2f (results_exp_ARPL[5,3]) " & " %5.2f (results_exp[5,3]) " & " %5.2f (results_exp_belowmed_ARPL[5,3]) " & " %5.2f (results_exp_belowmed[5,3]) " & " %5.2f (results_exp_abovemed_ARPL[5,3]) " & " %5.2f (results_exp_abovemed[5,3]) " \\" _n
file write didfile "& (" %5.2f ((results_incumbent_ARPL[5,4]-results_incumbent_ARPL[5,2])/3.92) ") & (" %5.2f ((results_incumbent[5,4]-results_incumbent[5,2])/3.92) ") & (" %5.2f ((results_nexp_ARPL[5,4]-results_nexp_ARPL[5,2])/3.92) ") & (" %5.2f ((results_nexp[5,4]-results_nexp[5,2])/3.92) ") & (" %5.2f ((results_exp_ARPL[5,4]-results_exp_ARPL[5,2])/3.92) ") & (" %5.2f ((results_exp[5,4]-results_exp[5,2])/3.92) ") & (" %5.2f ((results_exp_belowmed_ARPL[5,4]-results_exp_belowmed_ARPL[5,2])/3.92) ") & (" %5.2f ((results_exp_belowmed[5,4]-results_exp_belowmed[5,2])/3.92) ") & (" %5.2f ((results_exp_abovemed_ARPL[5,4]-results_exp_abovemed_ARPL[5,2])/3.92) ") & (" %5.2f ((results_exp_abovemed[5,4]-results_exp_abovemed[5,2])/3.92) ") \\" _n
file write didfile "[1em]" _n
file write didfile "\$\beta_{3, t= 2005} \$& " %5.2f (results_incumbent_ARPL[6,3]) " & " %5.2f (results_incumbent[6,3]) " & " %5.2f (results_nexp_ARPL[6,3]) " & " %5.2f (results_nexp[6,3]) " & " %5.2f (results_exp_ARPL[6,3]) " & " %5.2f (results_exp[6,3]) " & " %5.2f (results_exp_belowmed_ARPL[6,3]) " & " %5.2f (results_exp_belowmed[6,3]) " & " %5.2f (results_exp_abovemed_ARPL[6,3]) " & " %5.2f (results_exp_abovemed[6,3]) " \\" _n
file write didfile "& (" %5.2f ((results_incumbent_ARPL[6,4]-results_incumbent_ARPL[6,2])/3.92) ") & (" %5.2f ((results_incumbent[6,4]-results_incumbent[6,2])/3.92) ") & (" %5.2f ((results_nexp_ARPL[6,4]-results_nexp_ARPL[6,2])/3.92) ") & (" %5.2f ((results_nexp[6,4]-results_nexp[6,2])/3.92) ") & (" %5.2f ((results_exp_ARPL[6,4]-results_exp_ARPL[6,2])/3.92) ") & (" %5.2f ((results_exp[6,4]-results_exp[6,2])/3.92) ") & (" %5.2f ((results_exp_belowmed_ARPL[6,4]-results_exp_belowmed_ARPL[6,2])/3.92) ") & (" %5.2f ((results_exp_belowmed[6,4]-results_exp_belowmed[6,2])/3.92) ") & (" %5.2f ((results_exp_abovemed_ARPL[6,4]-results_exp_abovemed_ARPL[6,2])/3.92) ") & (" %5.2f ((results_exp_abovemed[6,4]-results_exp_abovemed[6,2])/3.92) ") \\" _n
file write didfile "[1em]" _n
file write didfile "\$\beta_{3, t= 2006} \$& " %5.2f (results_incumbent_ARPL[7,3]) " & " %5.2f (results_incumbent[7,3]) " & " %5.2f (results_nexp_ARPL[7,3]) " & " %5.2f (results_nexp[7,3]) " & " %5.2f (results_exp_ARPL[7,3]) " & " %5.2f (results_exp[7,3]) " & " %5.2f (results_exp_belowmed_ARPL[7,3]) " & " %5.2f (results_exp_belowmed[7,3]) " & " %5.2f (results_exp_abovemed_ARPL[7,3]) " & " %5.2f (results_exp_abovemed[7,3]) " \\" _n
file write didfile "& (" %5.2f ((results_incumbent_ARPL[7,4]-results_incumbent_ARPL[7,2])/3.92) ") & (" %5.2f ((results_incumbent[7,4]-results_incumbent[7,2])/3.92) ") & (" %5.2f ((results_nexp_ARPL[7,4]-results_nexp_ARPL[7,2])/3.92) ") & (" %5.2f ((results_nexp[7,4]-results_nexp[7,2])/3.92) ") & (" %5.2f ((results_exp_ARPL[7,4]-results_exp_ARPL[7,2])/3.92) ") & (" %5.2f ((results_exp[7,4]-results_exp[7,2])/3.92) ") & (" %5.2f ((results_exp_belowmed_ARPL[7,4]-results_exp_belowmed_ARPL[7,2])/3.92) ") & (" %5.2f ((results_exp_belowmed[7,4]-results_exp_belowmed[7,2])/3.92) ") & (" %5.2f ((results_exp_abovemed_ARPL[7,4]-results_exp_abovemed_ARPL[7,2])/3.92) ") & (" %5.2f ((results_exp_abovemed[7,4]-results_exp_abovemed[7,2])/3.92) ") \\" _n
file write didfile "[1em]" _n
file write didfile "\$\beta_{3, t= 2007} \$& " %5.2f (results_incumbent_ARPL[8,3]) " & " %5.2f (results_incumbent[8,3]) " & " %5.2f (results_nexp_ARPL[8,3]) " & " %5.2f (results_nexp[8,3]) " & " %5.2f (results_exp_ARPL[8,3]) " & " %5.2f (results_exp[8,3]) " & " %5.2f (results_exp_belowmed_ARPL[8,3]) " & " %5.2f (results_exp_belowmed[8,3]) " & " %5.2f (results_exp_abovemed_ARPL[8,3]) " & " %5.2f (results_exp_abovemed[8,3]) " \\" _n
file write didfile "& (" %5.2f ((results_incumbent_ARPL[8,4]-results_incumbent_ARPL[8,2])/3.92) ") & (" %5.2f ((results_incumbent[8,4]-results_incumbent[8,2])/3.92) ") & (" %5.2f ((results_nexp_ARPL[8,4]-results_nexp_ARPL[8,2])/3.92) ") & (" %5.2f ((results_nexp[8,4]-results_nexp[8,2])/3.92) ") & (" %5.2f ((results_exp_ARPL[8,4]-results_exp_ARPL[8,2])/3.92) ") & (" %5.2f ((results_exp[8,4]-results_exp[8,2])/3.92) ") & (" %5.2f ((results_exp_belowmed_ARPL[8,4]-results_exp_belowmed_ARPL[8,2])/3.92) ") & (" %5.2f ((results_exp_belowmed[8,4]-results_exp_belowmed[8,2])/3.92) ") & (" %5.2f ((results_exp_abovemed_ARPL[8,4]-results_exp_abovemed_ARPL[8,2])/3.92) ") & (" %5.2f ((results_exp_abovemed[8,4]-results_exp_abovemed[8,2])/3.92) ") \\" _n
file write didfile "[1em]" _n
file write didfile "\$\beta_{3, t= 2008} \$& " %5.2f (results_incumbent_ARPL[9,3]) " & " %5.2f (results_incumbent[9,3]) " & " %5.2f (results_nexp_ARPL[9,3]) " & " %5.2f (results_nexp[9,3]) " & " %5.2f (results_exp_ARPL[9,3]) " & " %5.2f (results_exp[9,3]) " & " %5.2f (results_exp_belowmed_ARPL[9,3]) " & " %5.2f (results_exp_belowmed[9,3]) " & " %5.2f (results_exp_abovemed_ARPL[9,3]) " & " %5.2f (results_exp_abovemed[9,3]) " \\" _n
file write didfile "& (" %5.2f ((results_incumbent_ARPL[9,4]-results_incumbent_ARPL[9,2])/3.92) ") & (" %5.2f ((results_incumbent[9,4]-results_incumbent[9,2])/3.92) ") & (" %5.2f ((results_nexp_ARPL[9,4]-results_nexp_ARPL[9,2])/3.92) ") & (" %5.2f ((results_nexp[9,4]-results_nexp[9,2])/3.92) ") & (" %5.2f ((results_exp_ARPL[9,4]-results_exp_ARPL[9,2])/3.92) ") & (" %5.2f ((results_exp[9,4]-results_exp[9,2])/3.92) ") & (" %5.2f ((results_exp_belowmed_ARPL[9,4]-results_exp_belowmed_ARPL[9,2])/3.92) ") & (" %5.2f ((results_exp_belowmed[9,4]-results_exp_belowmed[9,2])/3.92) ") & (" %5.2f ((results_exp_abovemed_ARPL[9,4]-results_exp_abovemed_ARPL[9,2])/3.92) ") & (" %5.2f ((results_exp_abovemed[9,4]-results_exp_abovemed[9,2])/3.92) ") \\" _n
file write didfile "\hline" _n
file write didfile "pre-trend p-values & " %5.2f (p_1) " & " %5.2f (p_2) " & " %5.2f (p_3) " & " %5.2f (p_4) " & " %5.2f (p_5) " & " %5.2f (p_6) " & " %5.2f (p_7) " & " %5.2f (p_8) " & " %5.2f (p_9) " & " %5.2f (p_10) " \\" _n
file write didfile "\(N\) & " %15.0fc (eN_1) " & " %15.0fc (eN_2) " & " %15.0fc (eN_3) " & " %15.0fc (eN_4) " & " %15.0fc (eN_5) " & " %15.0fc (eN_6) " & " %15.0fc (eN_7) " & " %15.0fc (eN_8) " & " %15.0fc (eN_9) " & " %15.0fc (eN_10) " \\" _n
file write didfile "\hline\hline" _n
file write didfile "\end{tabular}" _n
file write didfile "}" _n
file write didfile "\parbox{0.75\textwidth}{\caption*{\scriptsize \textit{Note:} Coefficients \$\beta_3\$ and \$\beta_{4,ARPL}\$ for the regression in equation \ref{eq:DiD_reg}. Controls include each firm's industry code in 2000 (10 categories), regional dummies (8 categories), and age (24 distinct start dates). Standard errors are in parentheses and are clustered at the firm level.  } }" _n
file write didfile "\end{table}" _n
file close didfile


* ===================================================================
* FIGURE DiD_combined_placebo (main.tex Appendix ~line 1085)
* DiD_capital_incumbents_ARPL_placebo.eps, DiD_capital_exporters_ARPL_placebo.eps
* Source: Tetenyi_fake_did.do
* ===================================================================
gen treated2 = 1 if (no_access_credit_ever==1 & year>2004)
replace treated2 = 0 if (no_access_credit_ever==1 & year<2005)
replace treated2 = 0 if (flag_foreign_credit_lost==1 )
bysort sorszam: egen treated3 = max(treated2)
gen reform_fake = 2005 if treated3==1
gen TimeToTreat_fake=year-reform_fake

collect create exh6, replace
preserve
keep if(flag_incumbent==1)
collect _r_b _r_se, tag(model[(1)]): eventdd log_capital log_ARPL, timevar(TimeToTreat_fake) ci(rcap) method(hdfe, cluster(sorszam) absorb( no_access_credit_ever alapitas regiokod ind10 year )) graph_op( title("") scheme(s1mono) xlabel(-5 "2000" -4 "2001" -3 "2002" -2 "2003" -1 "2004" 0 "2005" 1 "2006" 2 "2007" 3 "2008") ylabel(-0.1 "-10" 0 "0" 0.1 "10" 0.2 "20"))
matrix results_incumbent_ARPL_fake = e(leads) \ e(lags)
collect p_d=r(p), tag(model[(1)]): estat leads
collect _r_b _r_se, tag(model[(2)]): eventdd log_capital, timevar(TimeToTreat_fake) ci(rcap) method(hdfe, cluster(sorszam) absorb( no_access_credit_ever alapitas regiokod ind10 year )) graph_op( title("") scheme(s1mono) xlabel(-5 "2000" -4 "2001" -3 "2002" -2 "2003" -1 "2004" 0 "2005" 1 "2006" 2 "2007" 3 "2008")  ylabel(-0.1 "-10" 0 "0" 0.1 "10" 0.2 "20"))
collect p_d=r(p), tag(model[(2)]): estat leads
matrix results_incumbent_fake = e(leads) \ e(lags)
estat leads
estat lags
restore

svmat results_incumbent_ARPL_fake, names(event_incumbent_ARPL_fake)
svmat results_incumbent_fake, names(event_incumbent_fake)

preserve
keep if(flag_incumbent==1 & flag_exporter_ever==0)
collect _r_b _r_se, tag(model[(3)]):  eventdd log_capital log_ARPL, timevar(TimeToTreat_fake) ci(rcap) method(hdfe, cluster(sorszam) absorb( no_access_credit_ever alapitas regiokod ind10 year )) graph_op( title("") scheme(s1mono) xlabel(-5 "2000" -4 "2001" -3 "2002" -2 "2003" -1 "2004" 0 "2005" 1 "2006" 2 "2007" 3 "2008")  ylabel(-0.1 "-10" 0 "0" 0.1 "10" 0.2 "20"))
matrix results_nexp_ARPL_fake = e(leads) \ e(lags)
collect p_d=r(p), tag(model[(3)]): estat leads
collect _r_b _r_se, tag(model[(4)]):  eventdd log_capital, timevar(TimeToTreat_fake) ci(rcap) method(hdfe, cluster(sorszam) absorb( no_access_credit_ever alapitas regiokod ind10 year )) graph_op( title("") scheme(s1mono) xlabel(-5 "2000" -4 "2001" -3 "2002" -2 "2003" -1 "2004" 0 "2005" 1 "2006" 2 "2007" 3 "2008")  ylabel(-0.1 "-10" 0 "0" 0.1 "10" 0.2 "20"))
collect p_d=r(p), tag(model[(4)]): estat leads
matrix results_nexp_fake = e(leads) \ e(lags)
estat leads
estat lags
restore

svmat results_nexp_ARPL_fake, names(event_nexp_ARPL_fake)
svmat results_nexp_fake, names(event_nexp_fake)

preserve
keep if(exportincumbent==1)
collect _r_b _r_se, tag(model[(5)]):  eventdd log_capital log_ARPL, timevar(TimeToTreat_fake) ci(rcap) method(hdfe, cluster(sorszam) absorb( no_access_credit_ever alapitas regiokod ind10 year )) graph_op( title("") scheme(s1mono) xlabel(-5 "2000" -4 "2001" -3 "2002" -2 "2003" -1 "2004" 0 "2005" 1 "2006" 2 "2007" 3 "2008")  ylabel(-0.1 "-10" 0 "0" 0.1 "10" 0.2 "20"))
matrix results_exp_ARPL_fake = e(leads) \ e(lags)
collect p_d=r(p), tag(model[(5)]): estat leads
collect _r_b _r_se, tag(model[(6)]): eventdd log_capital, timevar(TimeToTreat_fake) ci(rcap) method(hdfe, cluster(sorszam) absorb( no_access_credit_ever alapitas regiokod ind10 year )) graph_op( title("") scheme(s1mono) xlabel(-5 "2000" -4 "2001" -3 "2002" -2 "2003" -1 "2004" 0 "2005" 1 "2006" 2 "2007" 3 "2008")  ylabel(-0.1 "-10" 0 "0" 0.1 "10" 0.2 "20"))
collect p_d=r(p), tag(model[(6)]): estat leads
matrix results_exp_fake = e(leads) \ e(lags)
estat leads
estat lags
restore
svmat results_exp_fake, names(event_exp_fake)
svmat results_exp_ARPL_fake, names(event_exp_ARPL_fake)

preserve
keep if(exportincumbent==1 & median_equity_2000a==0)
collect _r_b _r_se, tag(model[(7)]):  eventdd log_capital log_ARPL, timevar(TimeToTreat_fake) ci(rcap) method(hdfe, cluster(sorszam) absorb( no_access_credit_ever alapitas regiokod ind10 year )) graph_op( title("") scheme(s1mono) xlabel(-5 "2000" -4 "2001" -3 "2002" -2 "2003" -1 "2004" 0 "2005" 1 "2006" 2 "2007" 3 "2008")  ylabel(-0.1 "-10" 0 "0" 0.1 "10" 0.2 "20"))
matrix results_exp_belowmed_ARPL_fake = e(leads) \ e(lags)
collect p_d=r(p), tag(model[(7)]): estat leads
collect _r_b _r_se, tag(model[(8)]): eventdd log_capital, timevar(TimeToTreat_fake) ci(rcap) method(hdfe, cluster(sorszam) absorb( no_access_credit_ever alapitas regiokod ind10 year )) graph_op( title("") scheme(s1mono) xlabel(-5 "2000" -4 "2001" -3 "2002" -2 "2003" -1 "2004" 0 "2005" 1 "2006" 2 "2007" 3 "2008")  ylabel(-0.1 "-10" 0 "0" 0.1 "10" 0.2 "20"))
collect p_d=r(p), tag(model[(8)]): estat leads
matrix results_exp_belowmed_fake = e(leads) \ e(lags)
estat leads
estat lags
restore
svmat results_exp_belowmed_ARPL_fake, names(event_exp_belowmed_ARPL_fake)
svmat results_exp_belowmed_fake, names(event_exp_belowmed_fake)

preserve
keep if(exportincumbent==1 & median_equity_2000a==1)
collect _r_b _r_se, tag(model[(9)]):  eventdd log_capital log_ARPL, timevar(TimeToTreat_fake) ci(rcap) method(hdfe, cluster(sorszam) absorb( no_access_credit_ever alapitas regiokod ind10 year )) graph_op( title("") scheme(s1mono) xlabel(-5 "2000" -4 "2001" -3 "2002" -2 "2003" -1 "2004" 0 "2005" 1 "2006" 2 "2007" 3 "2008")  ylabel(-0.1 "-10" 0 "0" 0.1 "10" 0.2 "20"))
matrix results_exp_abovemed_ARPL_fake = e(leads) \ e(lags)
collect p_d=r(p), tag(model[(9)]): estat leads
collect _r_b _r_se, tag(model[(10)]): eventdd log_capital, timevar(TimeToTreat_fake) ci(rcap) method(hdfe, cluster(sorszam) absorb( no_access_credit_ever alapitas regiokod ind10 year )) graph_op( title("") scheme(s1mono) xlabel(-5 "2000" -4 "2001" -3 "2002" -2 "2003" -1 "2004" 0 "2005" 1 "2006" 2 "2007" 3 "2008")  ylabel(-0.1 "-10" 0 "0" 0.1 "10" 0.2 "20"))
collect p_d=r(p), tag(model[(10)]): estat leads
matrix results_exp_abovemed_fake = e(leads) \ e(lags)
estat leads
estat lags
restore
svmat results_exp_abovemed_ARPL_fake, names(event_exp_abovemed_ARPL_fake)
svmat results_exp_abovemed_fake, names(event_exp_abovemed_fake)

replace event_exp_abovemed_fake1 = 2004 in 1
replace event_exp_abovemed_fake1 = 2003 in 2
replace event_exp_abovemed_fake1 = 2002 in 3
replace event_exp_abovemed_fake1 = 2001 in 4
replace event_exp_abovemed_fake1 = 2000 in 5
replace event_exp_abovemed_fake1 = 2005 in 6
replace event_exp_abovemed_fake1 = 2006 in 7
replace event_exp_abovemed_fake1 = 2007 in 8
replace event_exp_abovemed_fake1 = 2008 in 9

replace event_exp_belowmed_fake1 = 2004 in 1
replace event_exp_belowmed_fake1 = 2003 in 2
replace event_exp_belowmed_fake1 = 2002 in 3
replace event_exp_belowmed_fake1 = 2001 in 4
replace event_exp_belowmed_fake1 = 2000 in 5
replace event_exp_belowmed_fake1 = 2005 in 6
replace event_exp_belowmed_fake1 = 2006 in 7
replace event_exp_belowmed_fake1 = 2007 in 8
replace event_exp_belowmed_fake1 = 2008 in 9

replace event_incumbent_fake1 = 2004 in 1
replace event_incumbent_fake1 = 2003 in 2
replace event_incumbent_fake1 = 2002 in 3
replace event_incumbent_fake1 = 2001 in 4
replace event_incumbent_fake1 = 2000 in 5
replace event_incumbent_fake1 = 2005 in 6
replace event_incumbent_fake1 = 2006 in 7
replace event_incumbent_fake1 = 2007 in 8
replace event_incumbent_fake1 = 2008 in 9

replace event_exp_abovemed_ARPL_fake1 = 2004 in 1
replace event_exp_abovemed_ARPL_fake1 = 2003 in 2
replace event_exp_abovemed_ARPL_fake1 = 2002 in 3
replace event_exp_abovemed_ARPL_fake1 = 2001 in 4
replace event_exp_abovemed_ARPL_fake1 = 2000 in 5
replace event_exp_abovemed_ARPL_fake1 = 2005 in 6
replace event_exp_abovemed_ARPL_fake1 = 2006 in 7
replace event_exp_abovemed_ARPL_fake1 = 2007 in 8
replace event_exp_abovemed_ARPL_fake1 = 2008 in 9

replace event_exp_belowmed_ARPL_fake1 = 2004 in 1
replace event_exp_belowmed_ARPL_fake1 = 2003 in 2
replace event_exp_belowmed_ARPL_fake1 = 2002 in 3
replace event_exp_belowmed_ARPL_fake1 = 2001 in 4
replace event_exp_belowmed_ARPL_fake1 = 2000 in 5
replace event_exp_belowmed_ARPL_fake1 = 2005 in 6
replace event_exp_belowmed_ARPL_fake1 = 2006 in 7
replace event_exp_belowmed_ARPL_fake1 = 2007 in 8
replace event_exp_belowmed_ARPL_fake1 = 2008 in 9

replace event_incumbent_ARPL_fake1 = 2004 in 1
replace event_incumbent_ARPL_fake1 = 2003 in 2
replace event_incumbent_ARPL_fake1 = 2002 in 3
replace event_incumbent_ARPL_fake1 = 2001 in 4
replace event_incumbent_ARPL_fake1 = 2000 in 5
replace event_incumbent_ARPL_fake1 = 2005 in 6
replace event_incumbent_ARPL_fake1 = 2006 in 7
replace event_incumbent_ARPL_fake1 = 2007 in 8
replace event_incumbent_ARPL_fake1 = 2008 in 9

sort event_exp_abovemed_fake1
graph set window fontface "Times New Roman"

twoway rarea event_incumbent_ARPL_fake4  event_incumbent_ARPL_fake2 event_incumbent_ARPL_fake1, color(gs15) || line event_incumbent_ARPL_fake3 event_incumbent_ARPL_fake1, mc(black) ms(S) lpattern(dash) lwidth(1.5)  ///
|| line event_exp_ARPL_fake3 event_incumbent_ARPL_fake1, mc(red) ms(S) lwidth(1.5)    ///
|| line event_nexp_ARPL_fake3 event_incumbent_ARPL_fake1,  mc(black) lpattern(dot) lwidth(2.5) ///
legend(order(2 "All" 3 "Exporters" 4 "Non-exporters ") row(1)) ylabel(-0.2 "-20%" -0.1 "-10%" 0 "0%" 0.1 "10%" 0.2 "20%" 0.3 "30%" 0.4 "40%" 0.5 "50%") xtitle("") xla(2000/2008) scheme(s1mono) yscale(range (-0.2,0.5)) ///
xli(2004, lc(gs8) lw(thin)) xtitle("") ytitle("log(K)") yla(, ang(h))
graph export "$figures/DiD_capital_incumbents_ARPL_placebo.eps", as(eps) name("Graph") preview(off) replace

twoway rarea  event_exp_ARPL_fake4  event_exp_ARPL_fake2 event_exp_abovemed_ARPL_fake1, color(gs15) ||line event_exp_ARPL_fake3 event_exp_belowmed_ARPL_fake1, mc(black)  lwidth(1.5)  || line event_exp_abovemed_ARPL_fake3 event_exp_abovemed_ARPL_fake1, mc(red) lpattern(dash) ms(S) lwidth(1.5)    ///
 || line event_exp_belowmed_ARPL_fake3 event_exp_abovemed_ARPL_fake1, mc(black)  lpattern(dot) lwidth(2.5) ///
legend(order(2 "All" 3 "Above median equity" 4 "Below median equity" ) row(1))  ylabel(-0.2 "-20%" -0.1 "-10%" 0 "0%" 0.1 "10%" 0.2 "20%" 0.3 "30%" 0.4 "40%" 0.5 "50%") xtitle("") xla(2000/2008) scheme(s1mono) yscale(range (-0.2,0.5)) ///
xli(2004, lc(gs8) lw(thin)) xtitle("") ytitle("log(K)") yla(, ang(h))
graph export "$figures/DiD_capital_exporters_ARPL_placebo.eps", as(eps) name("Graph") preview(off) replace


* ===================================================================
* FIGURE DiD_combined_manufacturing (main.tex Appendix ~line 1179)
* DiD_capital_incumbents_manuf.eps, DiD_capital_incumbents_nomanuf.eps
* ===================================================================
collect create exa6, replace
preserve
keep if(flag_incumbent==1 & ind10 ==3  )
collect _r_b _r_se, tag(model[(1)]): eventdd log_capital log_ARPL, timevar(TimeToTreat) ci(rcap) method(hdfe, cluster(sorszam) absorb( no_access_credit_ever alapitas regiokod ind10 year )) graph_op( title("") scheme(s1mono) xlabel(-2 "2000" -1 "2001" 0 "2002" 1 "2003" 2 "2004" 3 "2005" 4 "2006" 5 "2007" 6 "2008") ylabel(-0.1 "-10" 0 "0" 0.1 "10" 0.2 "20"))
matrix results_incumbent_manuf = e(leads) \ e(lags)
collect p_d=r(p), tag(model[(1)]): estat leads
restore

svmat results_incumbent_manuf, names(event_incumbent_manuf)

preserve
keep if(flag_incumbent==1 & ind10 !=3  )
collect _r_b _r_se, tag(model[(2)]): eventdd log_capital log_ARPL, timevar(TimeToTreat) ci(rcap) method(hdfe, cluster(sorszam) absorb( no_access_credit_ever alapitas regiokod ind10 year )) graph_op( title("") scheme(s1mono) xlabel(-2 "2000" -1 "2001" 0 "2002" 1 "2003" 2 "2004" 3 "2005" 4 "2006" 5 "2007" 6 "2008") ylabel(-0.1 "-10" 0 "0" 0.1 "10" 0.2 "20"))
matrix results_incumbent_nomanuf = e(leads) \ e(lags)
collect p_d=r(p), tag(model[(2)]): estat leads
restore
svmat results_incumbent_nomanuf, names(event_incumbent_nomanuf)

preserve
keep if(flag_incumbent==1 & ind10 ==3  & flag_exporter_ever==0)
collect _r_b _r_se, tag(model[(3)]): eventdd log_capital log_ARPL, timevar(TimeToTreat) ci(rcap) method(hdfe, cluster(sorszam) absorb( no_access_credit_ever alapitas regiokod ind10 year )) graph_op( title("") scheme(s1mono) xlabel(-2 "2000" -1 "2001" 0 "2002" 1 "2003" 2 "2004" 3 "2005" 4 "2006" 5 "2007" 6 "2008") ylabel(-0.1 "-10" 0 "0" 0.1 "10" 0.2 "20"))
matrix results_nexp_manuf = e(leads) \ e(lags)
collect p_d=r(p), tag(model[(3)]): estat leads
restore

svmat results_nexp_manuf, names(event_nexp_manuf)

preserve
keep if(flag_incumbent==1 & ind10 !=3  & flag_exporter_ever==0)
collect _r_b _r_se, tag(model[(4)]): eventdd log_capital log_ARPL, timevar(TimeToTreat) ci(rcap) method(hdfe, cluster(sorszam) absorb( no_access_credit_ever alapitas regiokod ind10 year )) graph_op( title("") scheme(s1mono) xlabel(-2 "2000" -1 "2001" 0 "2002" 1 "2003" 2 "2004" 3 "2005" 4 "2006" 5 "2007" 6 "2008") ylabel(-0.1 "-10" 0 "0" 0.1 "10" 0.2 "20"))
matrix results_nexp_nomanuf = e(leads) \ e(lags)
collect p_d=r(p), tag(model[(4)]): estat leads
restore
svmat results_nexp_nomanuf, names(event_nexp_nomanuf)

preserve
keep if(ind10 ==3  & exportincumbent==1)
collect _r_b _r_se, tag(model[(5)]): eventdd log_capital log_ARPL, timevar(TimeToTreat) ci(rcap) method(hdfe, cluster(sorszam) absorb( no_access_credit_ever alapitas regiokod ind10 year )) graph_op( title("") scheme(s1mono) xlabel(-2 "2000" -1 "2001" 0 "2002" 1 "2003" 2 "2004" 3 "2005" 4 "2006" 5 "2007" 6 "2008") ylabel(-0.1 "-10" 0 "0" 0.1 "10" 0.2 "20"))
matrix results_exp_manuf = e(leads) \ e(lags)
collect p_d=r(p), tag(model[(5)]): estat leads
restore

svmat results_exp_manuf, names(event_exp_manuf)

preserve
keep if(ind10 !=3  & exportincumbent==1)
collect _r_b _r_se, tag(model[(6)]): eventdd log_capital log_ARPL, timevar(TimeToTreat) ci(rcap) method(hdfe, cluster(sorszam) absorb( no_access_credit_ever alapitas regiokod ind10 year )) graph_op( title("") scheme(s1mono) xlabel(-2 "2000" -1 "2001" 0 "2002" 1 "2003" 2 "2004" 3 "2005" 4 "2006" 5 "2007" 6 "2008") ylabel(-0.1 "-10" 0 "0" 0.1 "10" 0.2 "20"))
matrix results_exp_nomanuf = e(leads) \ e(lags)
collect p_d=r(p), tag(model[(6)]): estat leads
restore
svmat results_exp_nomanuf, names(event_exp_nomanuf)

replace event_incumbent_manuf1 = 2001.0 in 1
replace event_incumbent_manuf1 = 2000.0 in 2
replace event_incumbent_manuf1 = 2002.0 in 3
replace event_incumbent_manuf1 = 2003.0 in 4
replace event_incumbent_manuf1 = 2004.0 in 5
replace event_incumbent_manuf1 = 2005.0 in 6
replace event_incumbent_manuf1 = 2006.0 in 7
replace event_incumbent_manuf1 = 2007.0 in 8
replace event_incumbent_manuf1 = 2008.0 in 9

replace event_incumbent_nomanuf1 = 2001.0 in 1
replace event_incumbent_nomanuf1 = 2000.0 in 2
replace event_incumbent_nomanuf1 = 2002.0 in 3
replace event_incumbent_nomanuf1 = 2003.0 in 4
replace event_incumbent_nomanuf1 = 2004.0 in 5
replace event_incumbent_nomanuf1 = 2005.0 in 6
replace event_incumbent_nomanuf1 = 2006.0 in 7
replace event_incumbent_nomanuf1 = 2007.0 in 8
replace event_incumbent_nomanuf1 = 2008.0 in 9

sort event_incumbent_manuf1
graph set window fontface "Times New Roman"

twoway rarea event_incumbent_manuf4  event_incumbent_manuf2 event_incumbent_manuf1, color(gs15) || line event_incumbent_manuf3 event_incumbent_manuf1, mc(black) ms(S)  lwidth(1.5)  ///
|| line event_exp_manuf3 event_incumbent_manuf1, mc(red) ms(S) lpattern(dash) lwidth(1.5)    ///
|| line event_nexp_manuf3 event_incumbent_manuf1,  mc(black) lpattern(dot) lwidth(2.5) ///
legend(order(2 "All" 3 "Exporters" 4 "Non-exporters ") row(1)) ylabel(-0.2 "-20%" -0.1 "-10%" 0 "0%" 0.1 "10%" 0.2 "20%" 0.3 "30%" 0.4 "40%" 0.5 "50%") xtitle("") xla(2000/2008) scheme(s1mono) yscale(range (-0.2,0.5)) ///
xli(2001, lc(gs8) lw(thin)) xtitle("") ytitle("log(K)") yla(, ang(h))
graph export "$figures/DiD_capital_incumbents_manuf.eps", as(eps) name("Graph") preview(off) replace

twoway rarea event_incumbent_nomanuf4  event_incumbent_nomanuf2 event_incumbent_nomanuf1, color(gs15) || line event_incumbent_nomanuf3 event_incumbent_nomanuf1, mc(black) ms(S)  lwidth(1.5)  ///
|| line event_exp_nomanuf3 event_incumbent_nomanuf1, mc(red) ms(S) lpattern(dash) lwidth(1.5)    ///
|| line event_nexp_nomanuf3 event_incumbent_nomanuf1,  mc(black) lpattern(dot) lwidth(2.5) ///
legend(order(2 "All" 3 "Exporters" 4 "Non-exporters ") row(1)) ylabel(-0.2 "-20%" -0.1 "-10%" 0 "0%" 0.1 "10%" 0.2 "20%" 0.3 "30%" 0.4 "40%" 0.5 "50%") xtitle("") xla(2000/2008) scheme(s1mono) yscale(range (-0.2,0.5)) ///
xli(2001, lc(gs8) lw(thin)) xtitle("") ytitle("log(K)") yla(, ang(h))
graph export "$figures/DiD_capital_incumbents_nomanuf.eps", as(eps) name("Graph") preview(off) replace


* ===================================================================
* TABLE industry_shares_did (main.tex ~line 1153)
* Industry shares in 2000 and 2008 for selected sectors
* Source: Tetenyi_main_v5.do lines 1501-1568
* Exports to Tables/Table_industry_shares.tex; values printed to log
* ===================================================================
di "=== TABLE industry_shares_did: Industry shares by year ==="
preserve
    collapse (sum) value_added capital export_value exportdummy, by(year ind10)
    bys year: egen total_value_added = total(value_added)
    bys year: egen total_capital = total(capital)
    bys year: egen total_export_value = total(export_value)
    bys year: egen total_exportdummy = total(exportdummy)
    gen share_value_added = value_added/total_value_added
    gen share_capital = capital/total_capital
    gen share_export_value = export_value/total_export_value
    gen share_exportdummy = exportdummy/total_exportdummy
    keep if inlist(ind10, 3, 6, 7, 9, 10)
    keep if inlist(year, 2000, 2008)
    collapse (mean) share_value_added share_capital share_export_value share_exportdummy, by(year ind10)
    reshape wide share_value_added share_capital share_export_value share_exportdummy, i(ind10) j(year)
    gen order = .
    replace order = 1 if ind10 == 3
    replace order = 2 if ind10 == 6
    replace order = 3 if ind10 == 7
    replace order = 4 if ind10 == 9
    replace order = 5 if ind10 == 10
    sort order
    format share* %9.3f
    list ind10 share_value_added2000 share_capital2000 share_export_value2000 share_exportdummy2000 ///
               share_value_added2008 share_capital2008 share_export_value2008 share_exportdummy2008
    local name1 "Manufacturing"
    local name2 "Retail"
    local name3 "Transportation"
    local name4 "Business services"
    local name5 "Public services"
    file open myfile using "$tables/Table_industry_shares.tex", write replace
    file write myfile "\begin{table}[htbp]" _n
    file write myfile "\centering" _n
    file write myfile "\caption{Industry Shares in 2000 and 2008}" _n
    file write myfile "\label{tab:industry_shares_did}" _n
    file write myfile "\scalebox{0.8}{" _n
    file write myfile "\begin{tabular}{lcccccccc}" _n
    file write myfile "\toprule" _n
    file write myfile " & \multicolumn{4}{c}{2000} & \multicolumn{4}{c}{2008} \\" _n
    file write myfile "\cmidrule(lr){2-5} \cmidrule(lr){6-9}" _n
    file write myfile "Industry & VA & Capital & Exports & Exporters & VA & Capital & Exports & Exporters \\" _n
    file write myfile "\midrule" _n
    forvalues i = 1/5 {
        local va00  = string(share_value_added2000[`i'],  "%9.3f")
        local cap00 = string(share_capital2000[`i'],      "%9.3f")
        local exp00 = string(share_export_value2000[`i'], "%9.3f")
        local exd00 = string(share_exportdummy2000[`i'],  "%9.3f")
        local va08  = string(share_value_added2008[`i'],  "%9.3f")
        local cap08 = string(share_capital2008[`i'],      "%9.3f")
        local exp08 = string(share_export_value2008[`i'], "%9.3f")
        local exd08 = string(share_exportdummy2008[`i'],  "%9.3f")
        local nm `name`i''
        file write myfile "`nm' & `va00' & `cap00' & `exp00' & `exd00' & `va08' & `cap08' & `exp08' & `exd08' \\" _n
    }
    file write myfile "\bottomrule" _n
    file write myfile "\end{tabular}}" _n
    file write myfile "\parbox{0.90\textwidth}{\caption*{\scriptsize \textit{Note: Only the 5 sectors that are responsible for more than 5\% of the export value or exporter firms in 2000 or in 2008 are represented out of the 10 industries.}   } }" _n
    file write myfile "\end{table}" _n
    file close myfile
restore

* ===================================================================
* TABLE Table23 (main.tex ~line 1188)
* Value-added fixed-effect regression for calibrating productivity process
* Source: Tetenyi_main_v5.do line 976
* Results printed to log; copy values into hardcoded table in main.tex
* ===================================================================
di "=== TABLE Table23: Value-added fixed-effect regression ==="
xtset sorszam year
xtreg log_value_added L.log_value_added if flag_incumbent==1, fe
di "Lag coef:  " _b[L.log_value_added] "  s.e.: " _se[L.log_value_added]
di "Constant:  " _b[_cons]              "  s.e.: " _se[_cons]
di "N obs:     " e(N)
di "N groups:  " e(N_g)
di "sigma_u:   " e(sigma_u)
di "sigma_e:   " e(sigma_e)
di "R2 within: " e(r2_w)
di "R2 between:" e(r2_b)
di "R2 overall:" e(r2_o)
scalar t23_b_lag   = _b[L.log_value_added]
scalar t23_se_lag  = _se[L.log_value_added]
scalar t23_b_cons  = _b[_cons]
scalar t23_se_cons = _se[_cons]
scalar t23_n       = e(N)
scalar t23_ng      = e(N_g)
scalar t23_sigu    = e(sigma_u)
scalar t23_sige    = e(sigma_e)
scalar t23_r2w     = e(r2_w)
scalar t23_r2b     = e(r2_b)
scalar t23_r2o     = e(r2_o)
file open t23file using "$tables/Table_Table23.tex", write replace
file write t23file "\begin{table}[h!]" _n
file write t23file "\caption{Value added regression}" _n
file write t23file "\centering" _n
file write t23file "\scalebox{0.75}{" _n
file write t23file "\def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi}" _n
file write t23file "\begin{tabular}{l*{3}{c}}" _n
file write t23file "\hline\hline" _n
file write t23file "            &  Log(Value added) & s.e. \\" _n
file write t23file "\hline" _n
file write t23file "Lag of Log(Value added) & " %9.3f (t23_b_lag) "\sym{***} & (" %9.4f (t23_se_lag) ")\\" _n
file write t23file "Constant & " %9.3f (t23_b_cons) "\sym{***} & (" %9.4f (t23_se_cons) ")\\" _n
file write t23file "\hline" _n
file write t23file "Number of observations & " %9.0g (t23_n) " & \\" _n
file write t23file "Number of firms & " %9.0g (t23_ng) " & \\" _n
file write t23file "Standard deviation of individual effects & " %9.2f (t23_sigu) " & \\" _n
file write t23file "Standard deviation of error term & " %9.2f (t23_sige) " & \\" _n
file write t23file " Within  \$R^2\$ & " %9.4f (t23_r2w) " & \\" _n
file write t23file " Between \$R^2\$ & " %9.4f (t23_r2b) " & \\" _n
file write t23file " Overall \$R^2\$ & " %9.4f (t23_r2o) " & \\" _n
file write t23file "\hline" _n
file write t23file "\hline" _n
file write t23file "\end{tabular}" _n
file write t23file "}" _n
file write t23file "\parbox{0.62\textwidth}{\caption*{\scriptsize \textit{Note:} The fixed effect lagged dependent variable regression on firms that always existed between 2000 and 2008. As there are no industries in the model, 0.82 is the relevant value for the calibration.}}" _n
file write t23file "\label{table:Table23}" _n
file write t23file "\end{table}" _n
file close t23file

* ===================================================================
* TABLE data_op_decomp (main.tex Appendix ~line 1226)
* Dynamic Olley-Pakes decomposition for ARPK, 2001-2008
* Source: Tetenyi_alternatives.do lines 109-136
* Results stored in scalars, exported to Tables/Table_data_op_decomp.tex
* ===================================================================
di "=== TABLE data_op_decomp: Dynamic Olley-Pakes decomposition for ARPK ==="

di "Full sample:"
olleyp year sorszam log_ARPK capital, base(2001) end(2008)
scalar op_overall_full      = r(overall)
scalar op_incumbents_full   = r(incumbents)
scalar op_distrib_full      = r(distributional)
scalar op_realloc_full      = r(reallocation)
scalar op_newborn_full      = r(newborns)
scalar op_exit_full         = r(exiters)

di "Firms with foreign owners:"
preserve
keep if (no_access_credit_ever ==0)
olleyp year sorszam log_ARPK capital, base(2001) end(2008)
scalar op_overall_for       = r(overall)
scalar op_incumbents_for    = r(incumbents)
scalar op_distrib_for       = r(distributional)
scalar op_realloc_for       = r(reallocation)
scalar op_newborn_for       = r(newborns)
scalar op_exit_for          = r(exiters)
restore

di "Firms without foreign owners:"
preserve
keep if (no_access_credit_ever ==1)
olleyp year sorszam log_ARPK capital, base(2001) end(2008)
scalar op_overall_dom       = r(overall)
scalar op_incumbents_dom    = r(incumbents)
scalar op_distrib_dom       = r(distributional)
scalar op_realloc_dom       = r(reallocation)
scalar op_newborn_dom       = r(newborns)
scalar op_exit_dom          = r(exiters)
restore

di "Exporters with foreign owners:"
preserve
keep if (no_access_credit_ever ==0 & exportdummy==1)
olleyp year sorszam log_ARPK capital, base(2001) end(2008)
scalar op_overall_expfor    = r(overall)
scalar op_incumbents_expfor = r(incumbents)
scalar op_distrib_expfor    = r(distributional)
scalar op_realloc_expfor    = r(reallocation)
scalar op_newborn_expfor    = r(newborns)
scalar op_exit_expfor       = r(exiters)
restore

di "Exporters without foreign owners:"
preserve
keep if (no_access_credit_ever ==1 & exportdummy==1)
olleyp year sorszam log_ARPK capital, base(2001) end(2008)
scalar op_overall_expdom    = r(overall)
scalar op_incumbents_expdom = r(incumbents)
scalar op_distrib_expdom    = r(distributional)
scalar op_realloc_expdom    = r(reallocation)
scalar op_newborn_expdom    = r(newborns)
scalar op_exit_expdom       = r(exiters)
restore

di "Non-exporters with foreign owners:"
preserve
keep if (no_access_credit_ever ==0 & exportdummy==0)
olleyp year sorszam log_ARPK capital, base(2001) end(2008)
scalar op_overall_nexpfor    = r(overall)
scalar op_incumbents_nexpfor = r(incumbents)
scalar op_distrib_nexpfor    = r(distributional)
scalar op_realloc_nexpfor    = r(reallocation)
scalar op_newborn_nexpfor    = r(newborns)
scalar op_exit_nexpfor       = r(exiters)
restore

di "Non-exporters without foreign owners:"
preserve
keep if (no_access_credit_ever ==1 & exportdummy==0)
olleyp year sorszam log_ARPK capital, base(2001) end(2008)
scalar op_overall_nexpdom    = r(overall)
scalar op_incumbents_nexpdom = r(incumbents)
scalar op_distrib_nexpdom    = r(distributional)
scalar op_realloc_nexpdom    = r(reallocation)
scalar op_newborn_nexpdom    = r(newborns)
scalar op_exit_nexpdom       = r(exiters)
restore

file open opfile using "$tables/Table_data_op_decomp.tex", write replace
file write opfile "\begin{table}[h!]" _n
file write opfile "\caption{Dynamic Olley-Pakes decomposition for ARPK growth, \% , 2001-2008}" _n
file write opfile "\label{table:data_op_decomp}" _n
file write opfile "    \centering" _n
file write opfile "\scalebox{0.7}{\begin{tabular}{l*{7}{c}}" _n
file write opfile "\hline\hline" _n
file write opfile "&  Overall &   Incumbent &  Distributional &  Reallocation  &  Newborn  &   Exiting  \\" _n
file write opfile "\hline" _n
file write opfile "Full sample & " %9.0g (round(op_overall_full)) " & " %9.0g (round(op_incumbents_full)) " & " %9.0g (round(op_distrib_full)) " & " %9.0g (round(op_realloc_full)) " & " %9.0g (round(op_newborn_full)) " & " %9.0g (round(op_exit_full)) " \\" _n
file write opfile "Firms with foreign owners & " %9.0g (round(op_overall_for)) " & " %9.0g (round(op_incumbents_for)) " & " %9.0g (round(op_distrib_for)) " & " %9.0g (round(op_realloc_for)) " & " %9.0g (round(op_newborn_for)) " & " %9.0g (round(op_exit_for)) " \\" _n
file write opfile "Firms without foreign owners & " %9.0g (round(op_overall_dom)) " & " %9.0g (round(op_incumbents_dom)) " & " %9.0g (round(op_distrib_dom)) " & " %9.0g (round(op_realloc_dom)) " & " %9.0g (round(op_newborn_dom)) " & " %9.0g (round(op_exit_dom)) " \\" _n
file write opfile "Exporters with foreign owners & " %9.0g (round(op_overall_expfor)) " & " %9.0g (round(op_incumbents_expfor)) " & " %9.0g (round(op_distrib_expfor)) " & " %9.0g (round(op_realloc_expfor)) " & " %9.0g (round(op_newborn_expfor)) " & " %9.0g (round(op_exit_expfor)) " \\" _n
file write opfile "Exporters without foreign owners & " %9.0g (round(op_overall_expdom)) " & " %9.0g (round(op_incumbents_expdom)) " & " %9.0g (round(op_distrib_expdom)) " & " %9.0g (round(op_realloc_expdom)) " & " %9.0g (round(op_newborn_expdom)) " & " %9.0g (round(op_exit_expdom)) " \\" _n
file write opfile "Non-exporters with foreign owners & " %9.0g (round(op_overall_nexpfor)) " & " %9.0g (round(op_incumbents_nexpfor)) " & " %9.0g (round(op_distrib_nexpfor)) " & " %9.0g (round(op_realloc_nexpfor)) " & " %9.0g (round(op_newborn_nexpfor)) " & " %9.0g (round(op_exit_nexpfor)) " \\" _n
file write opfile "Non-exporters without foreign owners & " %9.0g (round(op_overall_nexpdom)) " & " %9.0g (round(op_incumbents_nexpdom)) " & " %9.0g (round(op_distrib_nexpdom)) " & " %9.0g (round(op_realloc_nexpdom)) " & " %9.0g (round(op_newborn_nexpdom)) " & " %9.0g (round(op_exit_nexpdom)) " \\" _n
file write opfile "\hline" _n
file write opfile "\hline" _n
file write opfile "\end{tabular}}" _n
file write opfile "\caption*{  \parbox{0.95\textwidth}{\scriptsize \textit{Note:}  Capital used as weights \$s_j = \frac{k_j}{K}\$ used to decompose the overall growth of average revenue productivity of capital, \$ARPK = \sum_{j \in J} ARPK_j s_j\$ with \$ARPK_j = \frac{VA_j}{K_j}\$. Overall refers to all firms while Incumbent, Newborn and Exiting correspond to separating the set of indices \$J\$ into firms that are either incumbents in, enter to, or exit from the sub-sample year on year. A firm can be incumbent in one year, exiting in another one. Incumbent ARPK gains can be due to firm-level distributional shifts or as a result of capital reallocation. I provide further details about the interpretation of each column in Appendix \ref{appendix:firm_data}}}" _n
file write opfile "\end{table}" _n
file close opfile

* ===================================================================
* MOMENTS FOR MIXED TABLES (Table1, Table_regs, Table5)
* These three tables combine Hungarian firm-level moments (computed
* below) with moments from the Julia quantitative model. The output
* file Hungary_firm_level_data.txt documents which row/column of each
* table comes from which source, and reports all firm-level values.
* ===================================================================
di "=== MOMENTS FOR MIXED TABLES (Table1, Table_regs, Table5) ==="

preserve

sort sorszam year

* Lagged export status and share — only valid if consecutive years (no panel gaps)
bysort sorszam (year): gen lag_exportdummy = exportdummy[_n-1] if year - year[_n-1] == 1
bysort sorszam (year): gen lag_exportshare = exportshare[_n-1] if year - year[_n-1] == 1

* Consecutive export duration: how many years in a row has firm been exporting up to year t
bysort sorszam (year): gen _cstart = (exportdummy==1) & (exportdummy[_n-1]==0 | _n==1)
replace _cstart = 0 if exportdummy==0
bysort sorszam (year): gen _spell = sum(_cstart) if exportdummy==1
bysort sorszam _spell (year): gen cum_dur = _n if exportdummy==1 & !missing(_spell)


* ===== TABLE 1: Calibrated parameters and moments =====
* Firm dynamics rows: f_ex, sigma_z, rho_z
* sigma_z (0.82) and rho_z (0.43) already come from Table23 regression above.
* f_ex target: fraction of all firms that ever enter exports (n_expent/n_all)
* n_expent and n_all are defined in the data_micro section above.

scalar t1_entry_rate = round(100 * n_expent / n_all, 0.1)

* ===== TABLE_regs: Exporting dynamics, Data columns (1)-(4) =====
* LPM of current export status on log value added, lagged status/share.
* Industry and year FEs (no industry for model 4). SE clustered by industry.

qui: reghdfe exportdummy log_value_added, absorb(year industry) vce(cluster industry)
scalar tr_va1    = round(_b[log_value_added],  0.0001)
scalar tr_se1    = round(_se[log_value_added], 0.0001)
scalar tr_N1     = e(N)
scalar tr_r2_1   = round(e(r2_a), 0.01)

qui: reghdfe exportdummy log_value_added lag_exportdummy, absorb(year industry) vce(cluster industry)
scalar tr_va2    = round(_b[log_value_added],  0.0001)
scalar tr_se2    = round(_se[log_value_added], 0.0001)
scalar tr_ld2    = round(_b[lag_exportdummy],  0.0001)
scalar tr_ld_se2 = round(_se[lag_exportdummy], 0.0001)
scalar tr_N2     = e(N)
scalar tr_r2_2   = round(e(r2_a), 0.01)

qui: reghdfe exportdummy log_value_added lag_exportdummy lag_exportshare, absorb(year industry) vce(cluster industry)
scalar tr_va3    = round(_b[log_value_added],  0.0001)
scalar tr_se3    = round(_se[log_value_added], 0.0001)
scalar tr_ld3    = round(_b[lag_exportdummy],  0.0001)
scalar tr_ld_se3 = round(_se[lag_exportdummy], 0.0001)
scalar tr_ls3    = round(_b[lag_exportshare],  0.0001)
scalar tr_ls_se3 = round(_se[lag_exportshare], 0.0001)
scalar tr_N3     = e(N)
scalar tr_r2_3   = round(e(r2_a), 0.01)

qui: reghdfe exportdummy log_value_added lag_exportdummy lag_exportshare, absorb(year) vce(cluster industry)
scalar tr_va4    = round(_b[log_value_added],  0.0001)
scalar tr_se4    = round(_se[log_value_added], 0.0001)
scalar tr_ld4    = round(_b[lag_exportdummy],  0.0001)
scalar tr_ld_se4 = round(_se[lag_exportdummy], 0.0001)
scalar tr_ls4    = round(_b[lag_exportshare],  0.0001)
scalar tr_ls_se4 = round(_se[lag_exportshare], 0.0001)
scalar tr_N4     = e(N)
scalar tr_r2_4   = round(e(r2_a), 0.01)

* ===== TABLE 5: Impact of trade liberalization, Data columns (2001 and 2008) =====

foreach y in 2001 2008 {

    qui: sum capital if year==`y' & exportdummy==1
    local capex = r(sum)
    qui: sum capital if year==`y'
    scalar t5_pct_cap_`y' = round(100*`capex'/r(sum), 0.1)

    qui: sum letszam if year==`y' & exportdummy==1
    local labex = r(sum)
    qui: sum letszam if year==`y'
    scalar t5_pct_lab_`y' = round(100*`labex'/r(sum), 0.1)

    qui: sum log_ARPK if year==`y' & exportdummy==0
    scalar t5_sd_arpk_d_`y' = round(r(sd), 0.01)
    qui: sum log_ARPK if year==`y' & exportdummy==1
    scalar t5_sd_arpk_x_`y' = round(r(sd), 0.01)

    qui: sum exportshare if year==`y'
    scalar t5_expsales_`y' = round(100*r(mean), 0.1)

    qui: sum exportdummy if year==`y'
    scalar t5_ext_`y' = round(100*r(mean), 0.1)

    qui: sum exportshare if year==`y' & exportdummy==1
    scalar t5_int_`y' = round(100*r(mean), 0.1)

    qui: sum arbevert if year==`y' & exportdummy==1
    local mean_ex = r(mean)
    qui: sum arbevert if year==`y'
    scalar t5_premium_`y' = round(100*(`mean_ex'/r(mean) - 1), 0.1)

    qui: count if lag_exportdummy==0 & exportdummy==1 & year==`y'
    local nstart = r(N)
    qui: count if lag_exportdummy==0 & year==`y' & !missing(lag_exportdummy)
    scalar t5_starter_`y' = round(100*`nstart'/r(N), 0.1)

    qui: sum lag_exportdummy if year==`y' & exportdummy==0 & !missing(lag_exportdummy)
    scalar t5_stopper_`y' = round(100*r(mean), 0.1)

    qui: sum capital if year==`y' & exportdummy==0
    scalar t5_cap_d_`y'   = round(r(mean), 1)
    qui: sum capital if year==`y' & exportdummy==1
    scalar t5_cap_x_`y'   = round(r(mean), 1)
    qui: sum capital if year==`y'
    scalar t5_cap_all_`y' = round(r(mean), 1)
}

* Override 2008 capital averages using DiD event-study growth rates.
* results_X_ARPL = e(leads) \ e(lags) from eventdd log_capital log_ARPL (9 rows:
*   base 2001 in row 1, lead 2000 in row 2, lags 2002-2008 in rows 3-9).
*   Column layout: [year | ll | b | ul]. Row 9, col 3 = 2008 point estimate.
*   ARPL-controlled spec matches what was eyeballed from Figure 2 for main.tex.
scalar t5_cap_d_2008   = round(t5_cap_d_2001   * exp(results_nexp_ARPL[9,3]),       1)
scalar t5_cap_x_2008   = round(t5_cap_x_2001   * exp(results_exp_ARPL[9,3]),        1)
scalar t5_cap_all_2008 = round(t5_cap_all_2001  * exp(results_incumbent_ARPL[9,3]), 1)

qui: sum cum_dur if year==2001 & exportdummy==1
scalar t5_dur_2001 = round(r(mean), 0.01)
qui: sum cum_dur if year==2008 & exportdummy==1
scalar t5_dur_2008 = round(r(mean), 0.01)

restore

* ===== Write Hungary_firm_level_data.txt =====
file open hfl using "$tables/Hungary_firm_level_data.txt", write replace

file write hfl "======================================================" _n
file write hfl " Hungary Firm-Level Moments for Mixed Tables" _n
file write hfl " Source: firm_data_appended.dta via empirics_replication.do" _n
file write hfl " Julia model moments are indicated separately per table." _n
file write hfl "======================================================" _n _n

file write hfl "TABLE 1: Calibrated parameters and moments" _n
file write hfl "------------------------------------------------------" _n
file write hfl "Firm dynamics rows — Data column from firm_appended.dta:" _n
file write hfl "  f_ex    | Entry rate to exports  | Data = " %5.1f (t1_entry_rate) " % | Model target = 20 %" _n
file write hfl "           (fraction of all firms that ever enter exports: n_expent/n_all)" _n
file write hfl "  sigma_z | SD of value added      | Data = 0.82  (from Table23 FE regression)" _n
file write hfl "  rho_z   | AR(1) of value added   | Data = 0.431 (from Table23 FE regression)" _n _n
file write hfl "Financial Development and Trade rows: macro data only (BIS, ECB, WB, TiVA)." _n
file write hfl "Model column (Data and Model): Julia output." _n _n

file write hfl "TABLE_regs: Exporting dynamics in model vs. data" _n
file write hfl "------------------------------------------------------" _n
file write hfl "Entry cost / Per-period cost columns: Julia model output." _n
file write hfl "Data columns (1)-(4): LPM regressions on firm_appended.dta." _n
file write hfl "  Dep. var: exportdummy. Regressors: log(VA), lag(export_d), lag(export_share)." _n
file write hfl "  FEs: year always; industry for specs (1)-(3), not (4). SE clustered by industry." _n _n
file write hfl "  Data (1)  N = " %9.0fc (tr_N1) "  Adj.R2 = " %4.2f (tr_r2_1) _n
file write hfl "    Log value added:  " %7.4f (tr_va1) "  (" %7.4f (tr_se1) ")" _n _n
file write hfl "  Data (2)  N = " %9.0fc (tr_N2) "  Adj.R2 = " %4.2f (tr_r2_2) _n
file write hfl "    Log value added:  " %7.4f (tr_va2) "  (" %7.4f (tr_se2) ")" _n
file write hfl "    Prev export status: " %7.4f (tr_ld2) "  (" %7.4f (tr_ld_se2) ")" _n _n
file write hfl "  Data (3)  N = " %9.0fc (tr_N3) "  Adj.R2 = " %4.2f (tr_r2_3) _n
file write hfl "    Log value added:  " %7.4f (tr_va3) "  (" %7.4f (tr_se3) ")" _n
file write hfl "    Prev export status: " %7.4f (tr_ld3) "  (" %7.4f (tr_ld_se3) ")" _n
file write hfl "    Prev export share: " %7.4f (tr_ls3) "  (" %7.4f (tr_ls_se3) ")" _n _n
file write hfl "  Data (4)  N = " %9.0fc (tr_N4) "  Adj.R2 = " %4.2f (tr_r2_4) _n
file write hfl "    Log value added:  " %7.4f (tr_va4) "  (" %7.4f (tr_se4) ")" _n
file write hfl "    Prev export status: " %7.4f (tr_ld4) "  (" %7.4f (tr_ld_se4) ")" _n
file write hfl "    Prev export share: " %7.4f (tr_ls4) "  (" %7.4f (tr_ls_se4) ")" _n _n

file write hfl "TABLE 5: Impact of trade liberalization on firms" _n
file write hfl "------------------------------------------------------" _n
file write hfl "Model columns (None / Trade / Trade & capital): Julia output." _n
file write hfl "Data columns (2001, 2008): firm_appended.dta moments below." _n
file write hfl "Note: firm counts and capital sizes in the table are normalized so that" _n
file write hfl "  2001 data = model Trade column (= 100 for capital, 134/75 for firms)." _n
file write hfl "  Raw capital averages below need that scaling before entering the table." _n
file write hfl "  Productivity loss rows: model only, no data counterpart." _n
file write hfl "  Distribution of exporters: needs model wealth/productivity thresholds." _n _n
file write hfl "  % capital used by exporters       | 2001 = " %5.1f (t5_pct_cap_2001) " | 2008 = " %5.1f (t5_pct_cap_2008) _n
file write hfl "  % labor used by exporters         | 2001 = " %5.1f (t5_pct_lab_2001) " | 2008 = " %5.1f (t5_pct_lab_2008) _n
file write hfl "  Avg export duration (consec yrs)  | 2001 = " %5.2f (t5_dur_2001) " | 2008 = " %5.2f (t5_dur_2008) _n
file write hfl "  SD of ARPK — non-exporters        | 2001 = " %5.2f (t5_sd_arpk_d_2001) " | 2008 = " %5.2f (t5_sd_arpk_d_2008) _n
file write hfl "  SD of ARPK — exporters            | 2001 = " %5.2f (t5_sd_arpk_x_2001) " | 2008 = " %5.2f (t5_sd_arpk_x_2008) _n
file write hfl "  Export-sales ratio (%)            | 2001 = " %5.1f (t5_expsales_2001) " | 2008 = " %5.1f (t5_expsales_2008) _n
file write hfl "  Extensive margin (% exporters)    | 2001 = " %5.1f (t5_ext_2001) " | 2008 = " %5.1f (t5_ext_2008) _n
file write hfl "  Intensive margin (exp/tot sales%) | 2001 = " %5.1f (t5_int_2001) " | 2008 = " %5.1f (t5_int_2008) _n
file write hfl "  Exporter premium (%)              | 2001 = " %7.1f (t5_premium_2001) " | 2008 = " %7.1f (t5_premium_2008) _n
file write hfl "  Starter rate (%)                  | 2001 = " %5.1f (t5_starter_2001) " | 2008 = " %5.1f (t5_starter_2008) _n
file write hfl "  Stopper rate (%)                  | 2001 = " %5.1f (t5_stopper_2001) " | 2008 = " %5.1f (t5_stopper_2008) _n
file write hfl "  Avg capital non-exporters (raw)   | 2001 = " %9.1f (t5_cap_d_2001) " | 2008 = " %9.1f (t5_cap_d_2008) _n
file write hfl "  Avg capital exporters (raw)       | 2001 = " %9.1f (t5_cap_x_2001) " | 2008 = " %9.1f (t5_cap_x_2008) _n
file write hfl "  Avg capital all firms (raw)       | 2001 = " %9.1f (t5_cap_all_2001) " | 2008 = " %9.1f (t5_cap_all_2008) _n _n
file write hfl "======================================================" _n

file close hfl
di "Hungary_firm_level_data.txt written to $tables"
