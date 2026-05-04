* ============================================================
* figures_15_18_table13.do
*
* Reproduces Figures 15-18 and Table 13 from main.tex
*
* Figures 15-17: Choropleth maps of European countries showing
*   - Figure 15: Relative change in Import/GDP, 1992-2008
*   - Figure 16: Change in TFP, 1992-2008
*   - Figure 17: Credit to private sector / GDP in 1992
* Figure 18: Chinn-Ito capital openness index by country group (time series)
* Table 13:  Panel FE regression of log(TFP) on trade and financial variables
*            (hardcoded in main.tex lines 893-914; do-file reproduces numbers)
*
* ============================================================
* REQUIRED STATA PACKAGES
*   spmap    -- choropleth maps  (ssc install spmap)
*   estout   -- regression tables (ssc install estout)
* ============================================================
* DATA FILES TO COPY INTO Data/ FOLDER
*   1. complete_country_panel_map.dta
*         Pre-merged European panel (EU + EFTA, ~25 countries, 1990-2008)
*         containing: id (spmap polygon ID), _merge (map merge flag),
*         country (ISO3), year, import_share (=100*import/GDP from PWT),
*         GFDDDI14 (credit/GDP from IMF), ka_open (Chinn-Ito index),
*         rtfpna (TFP from PWT 9.0), diff_tfp_growth_past17 (17-yr TFP change x100),
*         diff_import_share_past17 (17-yr change in import share, percentage points).
*
*   2. Europe_coord.dta
*         European shapefile polygon vertex coordinates for spmap.
*         Generated from Europe.shp via: shp2dta using Europe, ...
*
*   3. importshare_WB.dta
*         World Bank import/GDP shares by country-year (string variable
*         import_share_WB in percentage units, e.g. "42.3").
*         Merged on: country (ISO3), year.
*
*   4. complete_country_panel.dta
*         Full world panel (100+ countries, 1970-2014) merging:
*         Penn World Table 9.0, IMF financial development (GFDDDI14),
*         Chinn-Ito index (ka_open). Used for Figure 18 and Table 13.
*         Key variables: country (ISO3), year, rtfpna, csh_m (imports/GDP,
*         negative convention), GFDDDI14, ka_open.
* ============================================================

set more off
capture ssc install spmap,   replace
capture ssc install estout,  replace

* ===== PATHS =====
local data    "C:\Users\bpu063010\Downloads\replication_package\Data"
local figures "C:\Users\bpu063010\Downloads\replication_package\Figures"
local tables  "C:\Users\bpu063010\Downloads\replication_package\Tables"

* ============================================================
* FIGURES 15-17: European choropleth maps
* Source data: complete_country_panel_map.dta + importshare_WB.dta
*              + Europe_coord.dta (coordinates for spmap)
* ============================================================
set scheme s1mono
use "`data'/complete_country_panel_map.dta", clear

* --- Merge World Bank import shares ---
* import_share_WB in importshare_WB.dta is a string percentage (e.g. "42.3")
* Use m:1 since the map dataset has multiple polygon rows per country-year
merge m:1 country year using "`data'/importshare_WB.dta", ///
    keepusing(import_share_WB) keep(master match) nogen

* Convert string to numeric fraction (matching pwt_full_do.do lines 163-167)
gen double import_share_WB_num = real(import_share_WB) / 100   // fraction [0-1]
drop import_share_WB
rename import_share_WB_num import_share_WB

* Re-declare panel after merge (country_n encoding may have changed)
drop if country == "ALB"     // Albania lacks shapefile match; drop to avoid gaps
capture drop country_n
encode country, gen(country_n)
xtset country_n year

* --- Construct relative import share change (Figure 15) ---
* Ratio: import share in 2008 / import share in 1992 (17-year lag).
* Primary source: World Bank fraction. Fall back to PWT-based import_share/100 if missing.
gen double import_share_PWT = import_share / 100   // PWT: import_share is %-points
gen double diff_import_share_WBpast17 = import_share_WB   / L17.import_share_WB
replace    diff_import_share_WBpast17 = import_share_PWT  / L17.import_share_PWT ///
    if diff_import_share_WBpast17 == .

* --- Norway diagnostic ---
di "Norway (NOR) debug:"
list country year _merge import_share_WB import_share_PWT ///
    diff_import_share_WBpast17 if country == "NOR" & inlist(year, 1991, 1992, 2008), noobs

* ---------------------------------------------------------------
* DIAGNOSTICS: countries in the shapefile with no data (shown as
* empty polygons on the map)
* ---------------------------------------------------------------
di ""
di "=== Countries on map with NO DATA ==="
di ""
di "Figure 15 (Import/GDP ratio 2008 -- missing):"
list country if _merge == 3 & year == 2008 & diff_import_share_WBpast17 == ., ///
    noobs clean
di ""
di "Figure 16 (TFP change 2008 -- missing):"
list country if _merge == 3 & year == 2008 & diff_tfp_growth_past17 == ., ///
    noobs clean
di ""
di "Figure 17 (Credit/GDP 1992 -- missing):"
list country if _merge == 3 & year == 1992 & GFDDDI14 == ., ///
    noobs clean
di "---------------------------------------------------------------"
di ""

format diff_import_share_WBpast17 %5.2f

* --- Figure 15: Relative import/GDP change 1992-2008 ---
* Filter ratio >= 1 inline so 1992 rows are not dropped (needed for Figure 17)
spmap diff_import_share_WBpast17 using "`data'/Europe_coord.dta" ///
    if _merge == 3 & year == 2008 & diff_import_share_WBpast17 != ., ///
    id(id) fcolor(Blues) clnumber(5) ///
    legend(position(2) ///
        label(2 "+20%") label(3 "+37%") label(4 "+62%") ///
        label(5 "+81%") label(6 "+90%") size(huge)) ///
    ocolor(black ..)
graph export "`figures'/Figure15.eps", replace as(eps)

* --- Figure 16: Change in TFP 1992-2008 ---
* diff_tfp_growth_past17 = 100*(rtfpna_t - rtfpna_{t-17}) is pre-computed in .dta
* Compute integer breakpoints from data range
quietly summarize diff_tfp_growth_past17 if _merge == 3 & year == 2008
local vmin  = floor(`r(min)')
local vmax  = ceil(`r(max)')
local step  = ceil((`vmax' - `vmin') / 5)
local b0 = `vmin'
local b1 = `b0' + `step'
local b2 = `b1' + `step'
local b3 = `b2' + `step'
local b4 = `b3' + `step'
local b5 = `vmax'

spmap diff_tfp_growth_past17 using "`data'/Europe_coord.dta" ///
    if _merge == 3 & year == 2008, ///
    id(id) fcolor(Blues) ///
    clmethod(custom) clbreaks(`b0' `b1' `b2' `b3' `b4' `b5') ///
    legend(position(7) size(medium) ///
        label(2 "`b0'% to `b1'%") label(3 "`b1'% to `b2'%") ///
        label(4 "`b2'% to `b3'%") label(5 "`b3'% to `b4'%") ///
        label(6 "`b4'% to `b5'%")) ///
    ocolor(white ..)
graph export "`figures'/Figure16.eps", replace as(eps)

* --- Figure 17: Credit to private sector / GDP in 1992 ---
* GFDDDI14 = domestic credit to private sector, % GDP (IMF Financial Development Index)
spmap GFDDDI14 using "`data'/Europe_coord.dta" ///
    if _merge == 3 & year == 1992, ///
    id(id) fcolor(Blues) clmethod(custom) clbreaks(0 40 60 148) ///
    legtitle("") ///
    legend(position(7) ///
        label(2 "Below 40%") label(3 "Below 60%") label(4 "Above 60%") ///
        size(huge)) ///
    ocolor(black ..)
graph export "`figures'/Figure17.eps", replace as(eps)

* ============================================================
* FIGURE 18: Chinn-Ito index by country group (time series)
* Source data: complete_country_panel.dta (full world panel)
* ============================================================
use "`data'/complete_country_panel.dta", clear

* ----- Define country groups -----
* Definitions from paper text (see main.tex Section 2 and Appendix A):
*   South: Spain, Italy, Portugal, Greece
*   NMS  : CEE EU members (joined 2004+): CZE, EST, HUN, LVA, LTU,
*           POL, SVK, SVN, ROU, BGR, HRV
*   Core : all remaining EU + EFTA members
gen country_group = ""
replace country_group = "South" if inlist(country, "ESP","ITA","PRT","GRC")
replace country_group = "NMS"   if inlist(country, "CZE","EST","HUN","LVA","LTU") | ///
                                    inlist(country, "POL","SVK","SVN","ROU","BGR","HRV")
replace country_group = "Core"  if inlist(country, "AUT","BEL","DEU","DNK","FIN") | ///
                                    inlist(country, "FRA","GBR","IRL","ISL","LUX") | ///
                                    inlist(country, "NLD","NOR","SWE","CHE")
keep if country_group != ""
keep if year >= 1985 & year <= 2010

* Unweighted average of Chinn-Ito index within each group-year
collapse (mean) ka_open, by(country_group year)

* Reshape to wide for twoway plot
reshape wide ka_open, i(year) j(country_group) string

label variable ka_openCore  "Core"
label variable ka_openNMS   "NMS"
label variable ka_openSouth "South"

twoway ///
    (line ka_openCore  year, lpattern(solid)    lcolor(navy)        lwidth(medthick)) ///
    (line ka_openSouth year, lpattern(longdash) lcolor(cranberry)   lwidth(medthick)) ///
    (line ka_openNMS   year, lpattern(shortdash) lcolor(dkorange)   lwidth(medthick)), ///
    legend(order(1 "Core" 2 "South" 3 "NMS") position(6) rows(1) size(medsmall)) ///
    xtitle("") ytitle("Chinn-Ito index (unweighted avg.)") ///
    xline(2001, lpattern(shortdash) lcolor(gray) lwidth(thin)) ///
    note("Vertical line: 2001 Hungarian capital market reform.") ///
    scheme(s1mono)
graph export "`figures'/Figure18.eps", replace as(eps)

* ============================================================
* TABLE 13: Panel FE regression --- TFP and trade
* Source: class_panel.do (original specification)
* Dataset: complete_country_panel_map.dta (European panel, NOT the world panel)
* Sample: 27 EU/EFTA countries, 1990+, MLT dropped
*   log(TFP_it) = b0 + b1*log(Import/GDP)_it + b2*log(Credit/GDP)_it
*               + b3*(log(Import/GDP) x log(Credit/GDP))_it
*               + b4*CMI_it + b5*(log(Import/GDP) x CMI)_it
*               + country FE + year FE + e_it
* Expected (hardcoded in main.tex): b1=0.184, b2=0.185, b3=0.1061,
*   b4=-0.0343, b5=-0.0889; N=3983; Country and time FE.
* ============================================================
use "`data'/complete_country_panel.dta", clear

* --- Variable construction (class_panel.do uses csh_m_neg = -csh_m) ---
capture drop log_tfp
capture drop log_import
capture drop findev
gen double log_tfp    = log(rtfpna)
gen double log_import = log(-csh_m)
gen double findev     = log(GFDDDI14 / 100)

capture drop country_n
encode country, gen(country_n)
xtset country_n year

* --- Drop observations with missing key variables ---
keep if !missing(log_tfp, log_import, findev, ka_open)

* --- Regression (matches Table 13 in main.tex) ---
eststo clear
eststo table13: xtreg log_tfp log_import findev c.log_import#c.findev ///
    ka_open c.log_import#c.ka_open i.year, fe

* --- Extract coefficients and standard errors ---
local N_obs = e(N)

local b1 = _b[log_import]
local b2 = _b[findev]
local b3 = _b[c.log_import#c.findev]
local b4 = _b[ka_open]
local b5 = _b[c.log_import#c.ka_open]

local s1 = _se[log_import]
local s2 = _se[findev]
local s3 = _se[c.log_import#c.findev]
local s4 = _se[ka_open]
local s5 = _se[c.log_import#c.ka_open]

* --- Significance stars (two-tailed, using normal approximation for large N) ---
foreach v in 1 2 3 4 5 {
    local t`v' = `b`v'' / `s`v''
    if abs(`t`v'') >= 2.576 local star`v' "\sym{***}"
    else if abs(`t`v'') >= 1.960 local star`v' "\sym{**}"
    else if abs(`t`v'') >= 1.645 local star`v' "\sym{*}"
    else local star`v' ""
}

* --- Format numbers (4 significant figures, trailing zeros removed) ---
foreach v in 1 2 3 4 5 {
    local fb`v' = string(`b`v'', "%9.4g")
    local fs`v' = string(`s`v'', "%9.4g")
}

* --- Write LaTeX table matching main.tex structure exactly ---
* Use _char(36) for dollar signs to prevent Stata global macro expansion
* (e.g. $CMI would expand as a global macro if written as a string)

file open ftex using "`tables'/Table13.tex", write replace
file write ftex "\begin{table}[h!]" _n
file write ftex "\caption{TFP and trade}" _n
file write ftex "            \centering" _n
file write ftex "\scalebox{0.9}{" _n
file write ftex "\def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi}" _n
file write ftex "\begin{tabular}{l*{7}{c}}" _n
file write ftex "\hline\hline" _n
file write ftex "            & " _char(36) "\log(\frac{Import}{GDP})  " _char(36) "& " _n
file write ftex _char(36) "\log(\frac{Credit}{GDP})" _char(36) " " _n
file write ftex "&" _char(36) "\log(\frac{Import}{GDP}) \times \log(\frac{Credit}{GDP})" _char(36) "& " _n
file write ftex _char(36) "CMI " _char(36) " & " _char(36) "\log(\frac{Import}{GDP}) \times CMI" _char(36) " \\" _n
file write ftex "\hline" _n
file write ftex "Log(TFP) &`fb1'`star1' & `fb2'`star2'&     `fb3'`star3'&      `fb4'`star4'&      `fb5'`star5'\\" _n
file write ftex "s.e. &(`fs1')&(`fs2')&(`fs3')&    (`fs4')&    (`fs5')\\" _n
file write ftex "\hline" _n
file write ftex "\hline" _n
file write ftex "\multicolumn{6}{c}{\footnotesize Standard errors in parentheses. N = `N_obs', Country and time FE}" _n
file write ftex "\end{tabular}" _n
file write ftex "}" _n
file write ftex "" _n
file write ftex "" _n
file write ftex "\label{table:Table13}" _n
file write ftex "        \end{table}" _n
file close ftex

di ""
di "========================================================"
di "DONE. Files written:"
di "  `figures'/Figure15.eps   -- Relative import share change map"
di "  `figures'/Figure16.eps   -- TFP change map"
di "  `figures'/Figure17.eps   -- Credit/GDP 1992 map"
di "  `figures'/Figure18.eps   -- Chinn-Ito index by country group"
di "  `tables'/Table13.tex     -- Regression table (matches main.tex format)"
di "========================================================"
di ""
di "Table 13 coefficients (compare to main.tex hardcoded values):"
di "  log(Import/GDP):                   `fb1' `star1'  (se=`fs1')"
di "  log(Credit/GDP):                   `fb2' `star2'  (se=`fs2')"
di "  log(Import/GDP) x log(Credit/GDP): `fb3' `star3'  (se=`fs3')"
di "  CMI:                               `fb4' `star4'  (se=`fs4')"
di "  log(Import/GDP) x CMI:             `fb5' `star5'  (se=`fs5')"
di "  N = `N_obs'"
