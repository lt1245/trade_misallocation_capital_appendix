# =============================================================================
# figures_tables.jl
#
# Reproduces every figure and table that appears in main.tex and is generated
# by the Julia model code.  Sections are named after their main.tex label so
# that each output can be located immediately.
#
# OUTPUTS NOT PRODUCED (present in original plot files but absent from main.tex):
#   Figure5d.pdf   Figure8a.pdf   Figure8b.pdf   Figure11b.pdf   Figure12j.pdf
#   Figure19a/b    Figure20a/b    Figure21a/b    Figure22a/b/c
#
# HOW TO RESUME
# =============
#   1. Re-run prerequisites (or load saved state) so all SS variables are in scope.
#   2. Set each PROGRESS FLAG below to 1 for every step already completed.
#   3. Run this file; it continues from the first flag that equals 0.
# =============================================================================

using Plots
using StatsBase
using StatsPlots
using GeometryTypes
using GLM
using Printf
using XLSX

# Use GR backend and point it at the Windows system font directory so that
# serif fonts (Times New Roman etc.) are available.
gr()
if Sys.iswindows()
    ENV["GKS_FONT_DIRS"] = "C:/Windows/Fonts"
end
# Apply margins so descenders are not clipped by GR.
default(fontfamily = "Computer Modern",
        bottom_margin = 12Plots.mm,
        left_margin   = 12Plots.mm,
        right_margin  = 6Plots.mm,
        top_margin    = 6Plots.mm)

# Utility functions (copied from plot.jl)
fpct(x) = string(round(Int, round(100*x, digits=2, RoundNearestTiesAway), RoundNearestTiesAway))
function movmean(array::Array{Float64,1}, window::Int64)
    array_size = size(array)[1]
    array_smooth = zeros(array_size)
    for i = 1:array_size
        if i < window
            array_smooth[i] = sum(array[1:i]) / i
        else
            array_smooth[i] = sum(array[(i-window+1):i]) / window
        end
    end
    return array_smooth
end

function mat_creator(y_value)
    y_mat = zeros(36, 200)
    for x_index = 1:36
        for y_index = 1:200
            y_mat[x_index, y_index] = y_value[200*(x_index-1) + y_index]
        end
    end
    return y_mat
end

# =============================================================================
# PROGRESS FLAGS  (set each to 1 once the step has finished successfully)
# Steps run in the order listed here.
# =============================================================================
step_baseline_ss      = 0   # Step 1  — main_script.jl baseline SS
step_dev_ss           = 0   # Step 2  — developed_steady_states.jl (needed before transitions)
step_transitions      = 0   # Step 3  — transition cases 0-4, 6-7 (baseline + dev)
step_save_reg_vars    = 0   # Step 3b — baseline regression inputs (*_BL aliases)
step_asym_trans       = 0   # Step 3c — transition cases 2-4 with asymmetric_liberalization=1
step_fc_ss            = 0   # Step 4  — fixed_cost_steady_states.jl FC SS
step_fc_trans         = 0   # Step 5  — transition cases 18-20 (FC model)

# =============================================================================
# STEP 1 — Baseline steady states
# Provides: all _initial, _closed_CM_open_trade, _open_CM_open_trade,
#           _open_CM_closed_trade variable families; Baseline_parameter;
#           open_CM_open_trade_parameter; prices_*; aggregates_*.
# Required by: Table1, Table2, Figure4, Table4, Table5, Table_regs,
#              Figure5, Figure6, Figure7, Figure9, Figure10, Figure11,
#              Figure12, Table8, Table9, Table10, Table11, Table4b.
# =============================================================================
if step_baseline_ss == 0
    include("main_script.jl")
    step_baseline_ss = 1
    println("Step 1 done: baseline steady states loaded.")
end

# =============================================================================
# STEP 2 — Developed-economy steady states
# Must run before Step 3 because transition cases 6 and 7 require
# prices_initial_dev and other *_dev variables defined in this file.
# Provides: all _dev variable families, Developed_parameter.
# Required by: transition cases 6-7, tab:combined/Table10, tab:combined/Table11,
#              tab:Table4b.
# =============================================================================
if step_dev_ss == 0
    include("developed_steady_states.jl")
    step_dev_ss = 1
    println("Step 2 done: developed-economy steady states loaded.")
end

# =============================================================================
# STEP 3 — Baseline transition dynamics
# Loads: transition.jl cases 0-4 (NMS baseline) and 6-7 (dev economy)
# Case 0: initial → open_CM_open_trade
# Case 1: initial → closed_CM_open_trade
# Case 2: initial → open_CM_open_trade_initial  (from trade-only SS)
# Case 3: initial → open_CMdelayed_open_trade    (delayed CM, 10-yr lag)
# Case 4: initial → open_CM_closed_trade         (CM only, no trade)
# Case 6: dev initial → dev closed_CM_open_trade
# Case 7: dev initial → dev open_CM_open_trade
# Provides: V_saved_store_*, capital_supply_future_*, K_x_ratio_store_*,
#           TFP_store_*, total_consumption_store_*, coeff_store_*_dev, etc.
# Required by: Table8, Figure12, Figure12b, Figure12c, Table9, Table4b,
#              Table10, Table11, Table_regs (welfare calcs).
# =============================================================================
if step_transitions == 0
    load_solution = 1; save_solution = 0; asymmetric_liberalization = 0
    case_final = 0; include("transition.jl")
    case_final = 1; include("transition.jl")
    case_final = 2; include("transition.jl")
    case_final = 3; include("transition.jl")
    case_final = 4; include("transition.jl")
    case_final = 6; include("transition.jl")
    case_final = 7; include("transition.jl")
    step_transitions = 1
    println("Step 3 done: baseline transition dynamics loaded.")
end

# =============================================================================
# STEP 3b — Compute baseline regression inputs from main_script.jl outputs.
# The variables worker_past_dist, Q_trans, etc. are LOCAL to the solver and
# never set at global scope.  We reconstruct them here from the named outputs
# of Residual_stst_detailed (open_CM_open_trade, country 1 = NMS).
# Q_trans_BL is rebuilt by calling Q_transition with the converged coeff.
# Required by: [tab:Table_regs] columns 1-2.
# =============================================================================
if step_save_reg_vars == 0
    # ------------------------------------------------------------------
    # Grid size and stationary distribution slices for country 1
    # ------------------------------------------------------------------
    ns_fine_BL    = open_CM_open_trade_parameter.Country_spec_p[1][4]
    distr_c1      = Vector(distr_current_open_CM_open_trade[:, 1])
    worker_past_dist_BL   = distr_c1[1:ns_fine_BL]
    domestic_past_dist_BL = distr_c1[(ns_fine_BL+1):(2*ns_fine_BL)]
    exporter_past_dist_BL = distr_c1[(2*ns_fine_BL+1):(3*ns_fine_BL)]

    # ------------------------------------------------------------------
    # future_occupation_fine_avgs[j,jj]: probability-weighted probability
    # that a type-j agent transitions to type jj (each element is a vector
    # of length ns_fine, indexed by asset-productivity state)
    # ------------------------------------------------------------------
    _iid_prob_BL  = open_CM_open_trade_parameter.iid_cost_prob
    _niid_BL      = length(_iid_prob_BL)
    _focc_BL      = Array(future_occupation_fine_open_CM_open_trade[:, :, 1, :])
    future_occupation_fine_avgs_BL = Array{Vector{Float64}}(undef, 3, 3)
    for _j = 1:3, _jj = 1:3
        future_occupation_fine_avgs_BL[_j, _jj] = zeros(ns_fine_BL)
        for _jc = 1:_niid_BL
            future_occupation_fine_avgs_BL[_j, _jj] .+=
                _iid_prob_BL[_jc] .* (_focc_BL[:, _j, _jc] .== _jj)
        end
    end

    # ------------------------------------------------------------------
    # Derived firm-type flow variables (same formulas as in
    # country_residual_detailed, lines 1679-1687 of code_new.jl)
    # ------------------------------------------------------------------
    incumbents_exporter_BL             = exporter_past_dist_BL .* future_occupation_fine_avgs_BL[3, 3]
    entrants_domestic_from_workers_BL  = worker_past_dist_BL   .* future_occupation_fine_avgs_BL[1, 2]
    incumbents_domestic_BL             = domestic_past_dist_BL .* future_occupation_fine_avgs_BL[2, 2]
    exit_exporting_to_domestic_BL      = exporter_past_dist_BL .* future_occupation_fine_avgs_BL[3, 2]
    worker_pop_tmp_BL                  = worker_pop_open_CM_open_trade[1]

    # ------------------------------------------------------------------
    # Revenue variables on fine grid (country 1)
    # ------------------------------------------------------------------
    rev_xx_fine_BL = Vector(rev_xx_fine_store_open_CM_open_trade[:, 1])
    rev_dx_fine_BL = Vector(rev_dx_fine_store_open_CM_open_trade[:, 1])
    rev_d_fine_BL  = Vector(rev_d_fine_store_open_CM_open_trade[:, 1])

    # ------------------------------------------------------------------
    # Q_trans_BL: steady-state Markov transition matrix (row-stochastic).
    # Rebuilt by calling Q_transition with the converged coefficients from
    # main_script.jl.  All intermediate variables are scoped to this let
    # block to avoid polluting the global namespace.
    # ------------------------------------------------------------------
    Q_trans_BL = let
        _country = 1
        (_β, _α, _δ, _θ, _α₁, _α₂, _σ, _α₁eff, _α₂eff, _ω, _L,
         _FM, _FMₓ, _F, _Fₓ, _Exit, _Exitₓ, _iid_val, _iid_prob,
         _niid, _country_no, _τ, _a_min, _a_max, _fspace_a, _fspace_a_fine,
         _agrid, _bank_cost, _bounds, _CSp, _s_cell, _ns_cell, _sf_cell,
         _nsf_cell, _Phiz_cell, _Phizf_cell, _Phi_cell, _Phiaug_cell,
         _Pkron_cell, _Pkron1_cell, _Pkronf_cell, _expeg_cell,
         _ns_tmp, _nsf_tmp, _openness) = local_parameters(open_CM_open_trade_parameter)

        (_W, _r, _out_final, _pf, _, _) = price_reshaper(
            prices_open_CM_open_trade, _openness, _country_no, _bounds)

        (_, _s, _ns, _sf, _nsf, _Phiz, _Phizf, _Phi, _Phiaug,
         _Pk, _Pk1, _Pkf, _, _, _imat, _imatf,
         _apf, _foccf, _conf, _coeff, _, _pfprev, _pfcurr,
         _, _Wloc, _θloc, _, _τloc, _rloc,
         _bcloc, _, _R, _cd, _cx, _cz, _czbar, _czf, _czfbar,
         _ones, _, _, _βloc, _, _, _, _, _, _Qt,
         _onesf) = country_local_initialize(
            _country, _s_cell, _ns_cell, _sf_cell, _nsf_cell,
            _Phiz_cell, _Phizf_cell, _Phi_cell, _Phiaug_cell,
            _Pkron_cell, _Pkron1_cell, _Pkronf_cell,
            _σ, _niid, _pf, _W, _r, _out_final,
            _θ, _L, _τ, _country_no, _bank_cost, _δ, _ω, _β)

        # Plug in the converged value-function coefficients
        _coeff_conv = Array(coeff_final_open_CM_open_trade[:, :, _country])

        # Firm profits on fine grid (needed by income_creator)
        (_Pd, _Px, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _) =
            true_profit(_pfprev, _R, _Wloc, _sf,
                _cx, _cd, _czf, _czfbar,
                _θloc, _τloc, _σ, _α₁eff, _α₂eff, _α₁, _α₂, _onesf)

        # Income matrix and maximum savings constraint on fine grid
        (_imatf, _max_xf) = income_creator(_Wloc, _onesf, _Pd, _Px,
            _FM, _FMₓ, _country, _rloc, _sf, _iid_val, _F, _Fₓ,
            _niid, _imatf, _nsf, _a_min, _a_max)

        # Build transition matrix — argument order for positions 16-18
        # matches the call in country_residual_detailed (cons, a_prime, focc)
        (_Qtp, _, _, _, _) = Q_transition(_coeff_conv,
            _nsf, _onesf, _niid, _iid_prob, _imatf,
            _Pkf, _a_min, _Phizf, _βloc,
            _fspace_a, _fspace_a_fine, _pfcurr, _max_xf, _Pk1,
            _conf, _apf, _foccf, _Qt)

        # Q_transition returns column-stochastic Q'; transpose → row-stochastic
        sparse(_Qtp')
    end

    step_save_reg_vars = 1
    println("Step 3b done: baseline regression variables saved as *_BL aliases.")
end

# =============================================================================
# STEP 3c — Asymmetric trade-liberalization transition dynamics
# Re-runs transition cases 2, 3, 4 with asymmetric_liberalization = 1.
# NMS reduces tariffs over 11 years; Core reduces over 4 years.
# Outputs: all *_asym variable families
#   (TFP_store_*_asym, K_x_ratio_store_*_asym, V_saved_store_*_asym, etc.)
# Required by: [fig:Figure12_cons], [fig:Figure12_TFP], [fig:Figure12_ARPK],
#              [fig:Figure12_asym], [table:Table_asym].
# =============================================================================
if step_asym_trans == 0
    load_solution = 1; save_solution = 0; asymmetric_liberalization = 1
    case_final = 2; include("transition.jl")   # → *_closed_CM_open_trade_asym
    case_final = 3; include("transition.jl")   # → *_open_CM_open_trade_asym
    case_final = 4; include("transition.jl")   # → *_open_CMdelayed_open_trade_asym (sets T_cap=10)
    asymmetric_liberalization = 0              # restore global flag
    step_asym_trans = 1
    println("Step 3c done: asymmetric transition dynamics loaded.")
end

# =============================================================================
# STEP 4 — Fixed-cost model steady states
# Provides: all _FC variable families, Baseline_parameter_FC,
#           open_CM_open_trade_parameter_FC.
# Required by: tab:calibration_fixed_cost, tab:nontargeted_fixed_cost,
#              tab:results_fixed_cost, tab:results_fixed_cost_micro,
#              tab:Table_regs (columns 3-4).
# =============================================================================
if step_fc_ss == 0
    include("fixed_cost_steady_states.jl")
    step_fc_ss = 1
    println("Step 4 done: FC steady states loaded.")
end

# =============================================================================
# STEP 5 — FC model transition dynamics
# Loads: transition.jl cases 18 (trade-only FC), 19 (full FC), 20 (delayed FC)
# Provides: V_saved_store_*_FC, capital_supply_future_*_FC, etc.
# Required by: tab:calibration_fixed_cost, tab:results_fixed_cost,
#              tab:results_fixed_cost_micro, tab:Table_regs (columns 3-4).
# =============================================================================
if step_fc_trans == 0
    load_solution = 1; save_solution = 0; asymmetric_liberalization = 0
    case_final = 18; include("transition.jl")
    case_final = 19; include("transition.jl")
    case_final = 20; include("transition.jl")
    step_fc_trans = 1
    println("Step 5 done: FC transition dynamics loaded.")
end

# =============================================================================
# SHARED WELFARE CALCULATIONS
# Must run after Steps 1-2. Produces scalar welfare CEVs and V_trans_* vectors
# used by Table8, Table9, and several figures.
# Source: plot.jl lines 39-160
# =============================================================================

V_saved_initial = reshape(V_saved_fine_initial[:,:,1], size(V_saved_fine_initial)[1]*3)
V_saved_closed_CM_open_trade = reshape(V_saved_fine_closed_CM_open_trade[:,:,1], size(V_saved_fine_initial)[1]*3)
V_saved_open_CM_open_trade = reshape(V_saved_fine_open_CM_open_trade[:,:,1], size(V_saved_fine_initial)[1]*3)

welfare_init_closed_CM_open_trade_stst = sum(current_distr_store_initial[:,1] .* (exp.((V_saved_closed_CM_open_trade - V_saved_initial) * (1.0 - Baseline_parameter.β[1])) .- 1))
welfare_init_open_CM_open_trade_stst   = sum(current_distr_store_initial[:,1] .* (exp.((V_saved_open_CM_open_trade  - V_saved_initial) * (1.0 - Baseline_parameter.β[1])) .- 1))

V_trans_closed_CM_open_trade       = reshape(V_saved_store_closed_CM_open_trade[:,:,1,2],       size(V_saved_fine_initial)[1]*3)
V_trans_open_CM_open_trade         = reshape(V_saved_store_open_CM_open_trade[:,:,1,2],         size(V_saved_fine_initial)[1]*3)
V_trans_open_CM_open_trade_initial = reshape(V_saved_store_open_CM_open_trade_initial[:,:,1,2], size(V_saved_fine_initial)[1]*3)
V_trans_open_CM_closed_trade       = reshape(V_saved_store_open_CM_closed_trade[:,:,1,2],       size(V_saved_fine_initial)[1]*3)
V_trans_open_CMdelayed_open_trade  = reshape(V_saved_store_open_CMdelayed_open_trade[:,:,1,2],  size(V_saved_fine_initial)[1]*3)
V_trans_10yr_closed_CM_open_trade  = reshape(V_saved_store_closed_CM_open_trade[:,:,1,11],      size(V_saved_fine_initial)[1]*3)
V_trans_10yr_open_CM_open_trade    = reshape(V_saved_store_open_CM_open_trade[:,:,1,11],        size(V_saved_fine_initial)[1]*3)

welfare_change_clCM_otrade_trans =
    sum(current_distr_store_initial[:,1] .* (exp.((V_trans_closed_CM_open_trade - V_saved_initial) * (1.0 - Baseline_parameter.β[1])) .- 1))
welfare_change_oCM_otrade_trans =
    sum(current_distr_store_initial[:,1] .* (exp.((V_trans_open_CM_open_trade  - V_saved_initial) * (1.0 - Baseline_parameter.β[1])) .- 1))
welfare_change_oCM_cltrade_trans =
    sum(current_distr_store_initial[:,1] .* (exp.((V_trans_open_CM_closed_trade - V_saved_initial) * (1.0 - Baseline_parameter.β[1])) .- 1))
welfare_change_oCM_otrade_delayed10_trans =
    sum(current_distr_store_initial[:,1] .* (exp.((V_trans_open_CMdelayed_open_trade - V_saved_initial) * (1.0 - Baseline_parameter.β[1])) .- 1))
welfare_change_oCM_from_otrade_trans =
    sum(current_distr_store_closed_CM_open_trade[:,1] .* (exp.((V_trans_open_CM_open_trade_initial - V_saved_closed_CM_open_trade) * (1.0 - Baseline_parameter.β[1])) .- 1))

# Foreign-country welfare
V_saved_initial_F               = reshape(V_saved_fine_initial[:,:,2],               size(V_saved_fine_initial)[1]*3)
V_saved_closed_CM_open_trade_F  = reshape(V_saved_fine_closed_CM_open_trade[:,:,2],  size(V_saved_fine_initial)[1]*3)
V_saved_open_CM_open_trade_F    = reshape(V_saved_fine_open_CM_open_trade[:,:,2],    size(V_saved_fine_initial)[1]*3)

V_trans_closed_CM_open_trade_F            = reshape(V_saved_store_closed_CM_open_trade[:,:,2,2],       size(V_saved_fine_initial)[1]*3)
V_trans_open_CM_open_trade_F              = reshape(V_saved_store_open_CM_open_trade[:,:,2,2],         size(V_saved_fine_initial)[1]*3)
V_trans_open_CM_open_trade_initial_F      = reshape(V_saved_store_open_CM_open_trade_initial[:,:,2,2], size(V_saved_fine_initial)[1]*3)
V_trans_open_CM_closed_trade_F            = reshape(V_saved_store_open_CM_closed_trade[:,:,2,2],       size(V_saved_fine_initial)[1]*3)
V_trans_open_CMdelayed_open_trade_F       = reshape(V_saved_store_open_CMdelayed_open_trade[:,:,2,2],  size(V_saved_fine_initial)[1]*3)

welfare_init_closed_CM_open_trade_stst_F =
    sum(current_distr_store_initial[:,2] ./ Baseline_parameter.L[2] .* (exp.((V_saved_closed_CM_open_trade_F - V_saved_initial_F) * (1.0 - Baseline_parameter.β[2])) .- 1))
welfare_init_open_CM_open_trade_stst_F =
    sum(current_distr_store_initial[:,2] ./ Baseline_parameter.L[2] .* (exp.((V_saved_open_CM_open_trade_F  - V_saved_initial_F) * (1.0 - Baseline_parameter.β[2])) .- 1))
welfare_change_clCM_otrade_trans_F =
    sum(current_distr_store_initial[:,2] ./ Baseline_parameter.L[2] .* (exp.((V_trans_closed_CM_open_trade_F - V_saved_initial_F) * (1.0 - Baseline_parameter.β[2])) .- 1))
welfare_change_oCM_otrade_trans_F =
    sum(current_distr_store_initial[:,2] ./ Baseline_parameter.L[2] .* (exp.((V_trans_open_CM_open_trade_F  - V_saved_initial_F) * (1.0 - Baseline_parameter.β[2])) .- 1))
welfare_change_oCM_cltrade_trans_F =
    sum(current_distr_store_initial[:,2] ./ Baseline_parameter.L[2] .* (exp.((V_trans_open_CM_closed_trade_F - V_saved_initial_F) * (1.0 - Baseline_parameter.β[2])) .- 1))
welfare_change_oCM_otrade_delayed10_trans_F =
    sum(current_distr_store_initial[:,2] ./ Baseline_parameter.L[2] .* (exp.((V_trans_open_CMdelayed_open_trade_F - V_saved_initial_F) * (1.0 - Baseline_parameter.β[2])) .- 1))

welfare_clCM_otrade_totalEU =
    (welfare_change_clCM_otrade_trans + Baseline_parameter.L[2] * welfare_change_clCM_otrade_trans_F) / (Baseline_parameter.L[2] + 1)
welfare_oCM_otrade_totalEU =
    (welfare_change_oCM_otrade_trans + Baseline_parameter.L[2] * welfare_change_oCM_otrade_trans_F) / (Baseline_parameter.L[2] + 1)
welfare_delayed_totalEU =
    (welfare_change_oCM_otrade_delayed10_trans + Baseline_parameter.L[2] * welfare_change_oCM_otrade_delayed10_trans_F) / (Baseline_parameter.L[2] + 1)

# Distribution slices (reused across multiple sections)
current_exporter_initial  = current_distr_store_initial[convert(Int64, 2/3*size(current_distr_store_closed_CM_open_trade)[1]+1):end, 1]
current_domestic_initial  = current_distr_store_initial[7201:14400, 1]
current_exporter_closed_CM_open_trade = current_distr_store_closed_CM_open_trade[convert(Int64, 2/3*size(current_distr_store_closed_CM_open_trade)[1]+1):end, 1]
current_domestic_closed_CM_open_trade = current_distr_store_closed_CM_open_trade[7201:14400, 1]
current_exporter_open_CM_open_trade   = current_distr_store_open_CM_open_trade[convert(Int64, 2/3*size(current_distr_store_closed_CM_open_trade)[1]+1):end, 1]
current_domestic_open_CM_open_trade   = current_distr_store_open_CM_open_trade[7201:14400, 1]

current_exporter_closed_CM_open_trade_mat = mat_creator(current_exporter_closed_CM_open_trade)
current_domestic_closed_CM_open_trade_mat = mat_creator(current_domestic_closed_CM_open_trade)
current_producer_closed_CM_open_trade_mat = current_exporter_closed_CM_open_trade_mat + current_domestic_closed_CM_open_trade_mat
current_exporter_open_CM_open_trade_mat   = mat_creator(current_exporter_open_CM_open_trade)
current_domestic_open_CM_open_trade_mat   = mat_creator(current_domestic_open_CM_open_trade)
current_producer_open_CM_open_trade_mat   = current_exporter_open_CM_open_trade_mat + current_domestic_open_CM_open_trade_mat

println("Shared welfare calculations done.")

# =============================================================================
# [fig:Figure1]  Evolution of Hungary's foreign credit, imports, and tariffs
# main.tex: \label{fig:Figure1}
# Source: ../Data/Hungary_macrodata.xlsx
# Prerequisites: none (reads directly from Excel)
# =============================================================================

let
    # Series A: Foreign Credit to NFCs (quarterly 1987–2019)
    # Col A = date "YYYY-MM-DD", Col F = "Relevant NFC Foreign credit"
    nfc_raw = XLSX.readdata("../Data/Hungary_macrodata.xlsx", "Foregin credit to NFC", "A1:F222")
    fig1_credit = String[]
    for i in 3:size(nfc_raw, 1)
        d = nfc_raw[i, 1]
        v = nfc_raw[i, 6]
        (d isa String && !ismissing(v)) || continue
        yr = parse(Int, d[1:4])
        yr < 1987 && continue
        mo = parse(Int, d[6:7])
        q  = mo == 3 ? "$(yr).00" :
             mo == 6 ? "$(yr).25" :
             mo == 9 ? "$(yr).5"  :
                       "$(yr).75"
        push!(fig1_credit, "$(q)& $(Float64(v))\t")
    end

    # Series B: Import Share (annual, midpoints 1991.5–2017.5)
    # Column A = year (integer), Column B = import share % (may use comma decimal)
    imp_raw = XLSX.readdata("../Data/Hungary_macrodata.xlsx", "Importshares", "A1:B28")
    fig1_import = String[]
    for i in 2:size(imp_raw, 1)
        yr  = imp_raw[i, 1]
        val = imp_raw[i, 2]
        if !ismissing(yr) && !ismissing(val)
            y = isa(yr, Number) ? Int(yr) : parse(Int, string(yr))
            v = isa(val, Number) ? Float64(val) : parse(Float64, replace(string(val), "," => "."))
            push!(fig1_import, @sprintf("%.2f &%.5f", y + 0.5, v))
        end
    end

    # Series C: Avg. Tariff Rate on Imports
    # Implied rate = import tariff revenue / total import value, from "Hungary Imports" sheet
    # col 1 = year, col 11 = rate as fraction (1991–2003 only; 0 from 2004 after EU accession)
    tar_raw = XLSX.readdata("../Data/Hungary_macrodata.xlsx", "Hungary Imports", "A1:K16")
    tariff_dict = Dict{Int,Float64}()
    for i in axes(tar_raw, 1)[2:end]
        yr  = tar_raw[i, 1]
        rat = tar_raw[i, 11]
        !ismissing(yr) && !ismissing(rat) || continue
        tariff_dict[Int(yr)] = Float64(rat) * 100
    end
    fig1_tariff_lines = String[]
    for yy in 1991.5:1.0:2017.5
        t = get(tariff_dict, Int(round(yy - 0.5)), 0.0)
        push!(fig1_tariff_lines, t == 0.0 ? @sprintf("%.2f &0", yy) : @sprintf("%.2f &%.7f", yy, t))
    end

    open("../Figures/Figure1.tikz", "w") do io
        write(io, "\\begin{tikzpicture}\n\n")
        write(io, "% Quarterly series of Foreign Credit to firms\n")
        write(io, "\\pgfplotstableread[col sep=&, header=true]{\n")
        write(io, "description& A\n")
        foreach(l -> write(io, l * "\n"), fig1_credit)
        write(io, "}\\datatableentry\n\n")
        write(io, "% Import Share\n")
        write(io, "\\pgfplotstableread[col sep=&, header=true]{\n")
        write(io, "description& B\n")
        foreach(l -> write(io, l * "\n"), fig1_import)
        write(io, "}\\datatariff\n\n")
        write(io, "% Avg. Tariff Rate on Imports\n")
        write(io, "\\pgfplotstableread[col sep=&, header=true]{\n")
        write(io, "description& C\n")
        foreach(l -> write(io, l * "\n"), fig1_tariff_lines)
        write(io, "}\\dataeu\n\n")
        write(io, "\\begin{axis}[\n")
        write(io, "  name=mainaxis,\n")
        write(io, "  clip=false,\n")
        write(io, "  width=\\textwidth,\n")
        write(io, "  height=9cm,\n")
        write(io, "  xlabel={},\n")
        write(io, "  ylabel={Foreign Credit \\& Import Share (\\%)},\n")
        write(io, "  ylabel style={xshift=0pt, yshift=-30pt},\n")
        write(io, "  xmin=1990, xmax=2018,\n")
        write(io, "  ymin=0, ymax=90,\n")
        write(io, "  xtick={1990,1993,1996,1999,2002,2005,2008,2011,2014,2017},\n")
        write(io, "  xticklabel={\\pgfmathprintnumber[assume math mode=true,1000 sep={}]{\\tick}},\n")
        write(io, "  xticklabel style={rotate=0, font=\\small},\n")
        write(io, "  yticklabel={\\pgfmathprintnumber{\\tick}\\%},\n")
        write(io, "  yticklabel style={font=\\small},\n")
        write(io, "  label style={font=\\bfseries\\small},\n")
        write(io, "  tick align=outside,\n")
        write(io, "  tick style={black},\n")
        write(io, "  grid=both,\n")
        write(io, "  grid style={dashed, gray!30},\n")
        write(io, "  axis y line*=left,\n")
        write(io, "  legend columns=1,\n")
        write(io, "  legend style={\n")
        write(io, "    font=\\small,\n")
        write(io, "    at={(1.00,0.02)},\n")
        write(io, "    anchor=south east,\n")
        write(io, "    draw=black,\n")
        write(io, "    fill=white,\n")
        write(io, "    cells={align=left}\n")
        write(io, "  },\n")
        write(io, "]\n\n")
        write(io, "% Plot 1: Foreign Credit to Firms (Quarterly)\n")
        write(io, "\\addplot [blue, thick, restrict x to domain=1990:2018]\n")
        write(io, "    table [x=description, y=A] {\\datatableentry};\n")
        write(io, "\\addlegendentry{Foreign Credit}\n\n")
        write(io, "% Plot 2: Import Share\n")
        write(io, "\\addplot [red!70!black, thick, mark=*, mark options={scale=1, fill=white}] \n")
        write(io, "    table [x=description, y=B] {\\datatariff};\n")
        write(io, "\\addlegendentry{Import Share}\n\n")
        write(io, "% Plot 3: Avg. Tariff Rate on EU Imports (scaled to left axis)\n")
        write(io, "\\addplot [green!50!black, thick, dashed, mark=square*, mark options={scale=1}] \n")
        write(io, "    table [x=description, y expr=\\thisrow{C}*90/12] {\\dataeu};\n")
        write(io, "\\addlegendentry{Import Tariff}\n\n")
        write(io, "\\end{axis}\n\n")
        write(io, "% Second y-axis for Tariff only\n")
        write(io, "\\begin{axis}[\n")
        write(io, "  width=\\textwidth,\n")
        write(io, "  height=9cm,\n")
        write(io, "  ylabel={Import Tariff (\\%)},\n")
        write(io, "  ylabel style={green!50!black},\n")
        write(io, "  xmin=1990, xmax=2018,\n")
        write(io, "  ymin=0, ymax=12,\n")
        write(io, "  yticklabel={\\pgfmathprintnumber{\\tick}\\%},\n")
        write(io, "  yticklabel style={font=\\small, green!50!black},\n")
        write(io, "  label style={font=\\bfseries\\small},\n")
        write(io, "  tick align=outside,\n")
        write(io, "  tick style={black},\n")
        write(io, "  axis y line*=right,\n")
        write(io, "  axis x line=none,\n")
        write(io, "  every y tick/.append style={xshift=-0pt},\n")
        write(io, "  every y tick label/.append style={xshift=-30pt},\n")
        write(io, "]\n\n")
        write(io, "\\end{axis}\n\n")
        write(io, "\\end{tikzpicture}")
    end
    println("Figure1.tikz written.")
end

# =============================================================================
# Load Hungarian firm-level moments from Stata-generated text file.
# Run empirics_replication.do first to produce ../Tables/Hungary_firm_level_data.txt.
# =============================================================================
const hfd = let
    d = Dict{String,Float64}()
    spec = 0
    for ln in readlines("../Tables/Hungary_firm_level_data.txt")
        # Dual-year rows: "... | 2001 = X | 2008 = Y"
        m2 = match(r"\|\s*2001 =\s*([\d.]+)\s*\|\s*2008 =\s*([\d.]+)", ln)
        if m2 !== nothing
            v1, v2 = parse(Float64, m2[1]), parse(Float64, m2[2])
            if     occursin("% capital used", ln)
                d["pct_cap_2001"]   = v1/100; d["pct_cap_2008"]   = v2/100
            elseif occursin("% labor used", ln)
                d["pct_lab_2001"]   = v1/100; d["pct_lab_2008"]   = v2/100
            elseif occursin("Avg export duration", ln)
                d["dur_2001"]       = v1;     d["dur_2008"]       = v2
            elseif occursin("ARPK", ln) && occursin("non-exporters", ln)
                d["sd_arpk_d_2001"] = v1;     d["sd_arpk_d_2008"] = v2
            elseif occursin("ARPK", ln) && occursin("exporters", ln)
                d["sd_arpk_x_2001"] = v1;     d["sd_arpk_x_2008"] = v2
            elseif occursin("Export-sales ratio", ln)
                d["expsales_2001"]  = v1/100; d["expsales_2008"]  = v2/100
            elseif occursin("Extensive margin", ln)
                d["ext_2001"]       = v1/100; d["ext_2008"]       = v2/100
            elseif occursin("Intensive margin", ln)
                d["int_2001"]       = v1/100; d["int_2008"]       = v2/100
            elseif occursin("Exporter premium", ln)
                d["premium_2001"]   = v1/100; d["premium_2008"]   = v2/100
            elseif occursin("Starter rate", ln)
                d["starter_2001"]   = v1/100; d["starter_2008"]   = v2/100
            elseif occursin("Stopper rate", ln)
                d["stopper_2001"]   = v1/100; d["stopper_2008"]   = v2/100
            elseif occursin("capital non-exporters", ln)
                d["cap_d_2001"]     = v1;     d["cap_d_2008"]     = v2
            elseif occursin("capital exporters", ln)
                d["cap_x_2001"]     = v1;     d["cap_x_2008"]     = v2
            elseif occursin("capital all firms", ln)
                d["cap_all_2001"]   = v1;     d["cap_all_2008"]   = v2
            end
            continue
        end
        # Single-value Table 1 rows
        m = match(r"f_ex\s*\|.*Data =\s*([\d.]+)", ln)
        if m !== nothing; d["f_ex_pct"] = parse(Float64, m[1]); continue; end
        m = match(r"sigma_z\s*\|.*Data =\s*([\d.]+)", ln)
        if m !== nothing; d["sigma_z"] = parse(Float64, m[1]); continue; end
        m = match(r"rho_z\s*\|.*Data =\s*([\d.]+)", ln)
        if m !== nothing; d["rho_z"] = parse(Float64, m[1]); continue; end
        # Table_regs Data(k) header: "  Data (k)  N = X  Adj.R2 = Y"
        m = match(r"Data \((\d)\)\s+N =\s*([\d,]+)\s+Adj\.R2 =\s*([\d.]+)", ln)
        if m !== nothing
            spec = parse(Int, m[1])
            d["tr_N$spec"]   = parse(Float64, replace(m[2], "," => ""))
            d["tr_r2_$spec"] = parse(Float64, m[3])
            continue
        end
        # Table_regs regression coefficients and SEs
        if spec > 0
            m = match(r"Log value added:\s*([-\d.]+)\s+\(\s*([\d.]+)\)", ln)
            if m !== nothing
                d["tr_va$spec"] = parse(Float64, m[1])
                d["tr_se$spec"] = parse(Float64, m[2]); continue
            end
            m = match(r"Prev export status:\s*([-\d.]+)\s+\(\s*([\d.]+)\)", ln)
            if m !== nothing
                d["tr_ld$spec"]    = parse(Float64, m[1])
                d["tr_ld_se$spec"] = parse(Float64, m[2]); continue
            end
            m = match(r"Prev export share:\s*([-\d.]+)\s+\(\s*([\d.]+)\)", ln)
            if m !== nothing
                d["tr_ls$spec"]    = parse(Float64, m[1])
                d["tr_ls_se$spec"] = parse(Float64, m[2]); continue
            end
        end
    end
    d
end
println("  Loaded Hungary firm-level data ($(length(hfd)) entries) from Hungary_firm_level_data.txt")

# =============================================================================
# [table:Table1]  Calibrated parameters and moments
# main.tex: \label{table:Table1} — "Calibrated parameters and moments"
# Source: plot.jl lines 162–175
# Prerequisites: Step 1 (main_script.jl) only — no transition paths needed.
#
# Columns of table_1_results (8×3):
#   col 1 = calibrated parameter value
#   col 2 = data target (data moment)
#   col 3 = model moment (what the model produces for the same statistic)
#
# Rows:
#   1  θ       borrowing tightness          → total credit to NFCs / GDP
#   2  β*      core discount factor         → bank lending rate r*
#   3  β       NMS discount factor          → foreign credit / total credit
#   4  τ₀      initial import trade cost    → initial import / GDP
#   5  τ∞      final import trade cost      → final import / GDP
#   6  f_ex    avg export entry cost        → entry rate to exports / exporter pop
#   7  σ_z     s.d. productivity innovation → s.d. of log value added growth
#   8  ρ_z     AR(1) productivity           → autocorrelation of log value added
# =============================================================================

total_credit_open_CM_open_trade =
    -(domestic_firm_debt_open_CM_open_trade[1] + exporter_firm_debt_open_CM_open_trade[1]) /
    nomGDP_open_CM_open_trade[1]
domestic_credit_open_CM_open_trade =
    (worker_bond_holding_open_CM_open_trade[1] +
     domestic_bond_holding_open_CM_open_trade[1] +
     exporter_bond_holding_open_CM_open_trade[1]) ./
    nomGDP_open_CM_open_trade[1]
foreign_credit_open_CM_open_trade       = total_credit_open_CM_open_trade - domestic_credit_open_CM_open_trade
foreign_credit_share_open_CM_open_trade = foreign_credit_open_CM_open_trade / total_credit_open_CM_open_trade

let _nfc = XLSX.readdata("../Data/Hungary_macrodata.xlsx", "Foregin credit to NFC", "A1:G222")
    global _data_total_credit_pct = round(Int, Float64(_nfc[153, 7]))  # col G: corrected Credit/GDP
    global _data_fgn_share_pct    = round(Int, 100 * Float64(_nfc[153, 5]))  # col E: Foreign ratio
end

table_1_results = zeros(8, 3)
table_1_results[1,:] = [open_CM_open_trade_parameter.θ[1],  _data_total_credit_pct/100, total_credit_open_CM_open_trade]
table_1_results[2,:] = [open_CM_open_trade_parameter.β[2],  0.05, prices_open_CM_open_trade[6]]
table_1_results[3,:] = [open_CM_open_trade_parameter.β[1],  _data_fgn_share_pct/100, foreign_credit_share_open_CM_open_trade]
table_1_results[4,:] = [Baseline_parameter.τ[1],            0.21, import_share_initial[1]]
table_1_results[5,:] = [open_CM_open_trade_parameter.τ[1],  0.42, import_share_open_CM_open_trade[1]]
table_1_results[6,:] = [open_CM_open_trade_parameter.Fₓ[1], 0.27,
    entry_share_to_exporter_open_CM_open_trade[1] / exporter_pop_open_CM_open_trade[1]]
table_1_results[7,:] = [open_CM_open_trade_parameter.σₛ[1], 0.82, sd_growth_rev_open_CM_open_trade[1]]
table_1_results[8,:] = [open_CM_open_trade_parameter.ρ[1],  0.43, autocorr_rev_open_CM_open_trade[1]]

println("[table:Table1] table_1_results computed.")
println(round.(table_1_results, digits=3))

# --- Save Tables/Table1.tex ---
# Model moments (col 3): rows 1-6 shown as integers (×100 = %); rows 7-8 as decimals.
# Parameter values (col 1): displayed as decimals, except f_ex shown as integer × 100 + "%".
mkpath("../Tables")
let
    p = table_1_results[:,1]   # calibrated parameter values
    m = table_1_results[:,3]   # model moments

    # Value-column formatting
    val_theta  = @sprintf("%.2f", round(p[1], digits=2, RoundNearestTiesAway))
    val_betaSt = @sprintf("%.2f", round(p[2], digits=2, RoundNearestTiesAway))
    val_beta   = @sprintf("%.2f", round(p[3], digits=2, RoundNearestTiesAway))
    val_tau0   = @sprintf("%.2f", round(1 + p[4], digits=2, RoundNearestTiesAway))
    val_tauInf = @sprintf("%.2f", round(1 + p[5], digits=2, RoundNearestTiesAway))
    val_fex    = string(round(Int, round(100 * p[6], digits=1, RoundNearestTiesAway), RoundNearestTiesAway)) * " \\%"   # e.g. "450 \%"
    val_sigma  = @sprintf("%.2f", round(p[7], digits=2, RoundNearestTiesAway))
    val_rho    = @sprintf("%.2f", round(p[8], digits=2, RoundNearestTiesAway))

    # Model-column formatting (percentages for rows 1-6, decimals for 7-8)
    mod1 = string(round(Int, 100 * m[1]))
    mod2 = string(round(Int, 100 * m[2]))
    mod3 = string(round(Int, 100 * m[3]))
    mod4 = string(round(Int, 100 * m[4]))
    mod5 = string(round(Int, 100 * m[5]))
    mod6 = string(round(Int, 100 * m[6]))
    mod7 = @sprintf("%.2f", round(m[7], digits=2, RoundNearestTiesAway))
    mod8 = @sprintf("%.2f", round(m[8], digits=2, RoundNearestTiesAway))

    tex = """\\begin{table}[!ht]
\\centering
\\caption{Calibrated parameters and moments}
\\label{table:Table1}
\\scalebox{0.55}{
\\begin{tabular}{lccccc}
\\\\[-1.8ex]\\hline
\\hline \\\\[-1.8ex]
Parameter & Value & Target & Source \\& Year & Data & Model \\\\
\\hline \\\\[-1.8ex]
\\textbf{Financial Development} & & & & & \\\\[1ex]
 Borrowing tightness, \$\\theta\$ & $(val_theta) & \$\\text{Total Credit to nonfinancials}\$, \$\\%\\text{GDP}\$ & BIS 2008 & $(_data_total_credit_pct) & $(mod1) \\\\
 Core discount factor, \$\\beta^*\$ & $(val_betaSt) & Bank lending rate in Germany \$r^*\$ & ECB 2008 January & 5 & $(mod2) \\\\
 NMS discount factor, \$\\beta\$ & $(val_beta) & \$\\text{Foreign Credit to nonfinancials}\$, \$\\%\\text{Total credit}\$ & BIS 2008 & $(_data_fgn_share_pct) & $(mod3) \\\\[1ex]
\\textbf{Trade} & & & & & \\\\[1ex]
 Initial import trade cost, \$\\tau_0\$ & $(val_tau0) & Initial \$\\frac{\\text{Import}}{\\text{GDP}}\$ & WB 1991/TiVA 1995 & 21 & $(mod4) \\\\
 Final import trade cost, \$\\tau_{\\infty}\$ & $(val_tauInf) & Final \$\\frac{\\text{Import}}{\\text{GDP}}\$ & WB 2008/TiVA 2008 & 42 & $(mod5) \\\\[1ex]
\\textbf{Firm dynamics} & & & & & \\\\[1ex]
 Avg.\\ export entry cost, \$f_{ex}\$ & \$$(val_fex)\$ & Entry rate to exports & Firm level, Hungary & $(round(Int, hfd["f_ex_pct"])) & $(mod6) \\\\
 s.d.\\ of LN productivity innovation, \$\\sigma_{z}\$ & $(val_sigma) & s.d.\\ value added & Firm level, Hungary & $(hfd["sigma_z"]) & $(mod7) \\\\
 AR(1) of LN productivity innovation, \$\\rho_{z}\$ & $(val_rho) & Auto-correlation of value added & Firm level, Hungary & $(@sprintf("%.2f", hfd["rho_z"])) & $(mod8) \\\\
\\hline
\\hline
\\end{tabular}
}
\\parbox{0.9\\textwidth}{\\caption*{\\scriptsize \\textit{Note:} Simulated moments of the model matched with moments in the data. Construction of the moments is described in detail in Appendix \\ref{appendix:data_sources}. All parameters are calibrated under the assumption that both trade and capital market integration are completed by 2008, except for \$\\tau_0\$, which is calibrated to target the pre-reform economy as closely as possible.}}
\\end{table}
"""
    open("../Tables/Table1.tex", "w") do io
        write(io, tex)
    end
    println("  Saved → ../Tables/Table1.tex")
end

# =============================================================================
# [table:Table2]  Non-targeted / externally calibrated parameters
# main.tex: \label{table:Table2} — appears just after Table 1
# Source: plot.jl lines 177–183
# Prerequisites: Step 1 only.
#
# Columns of table_2_results (6×2):
#   col 1 = parameter value
#   col 2 = data or auxiliary moment (where applicable; else 0)
#
# Rows:
#   1  L[1]  NMS labour endowment
#   2  L[2]  Core labour endowment
#   3  σ     CES elasticity of substitution
#   4  θ[2]  Core borrowing tightness
#   5  δ     capital depreciation rate
#   6  F_x[2] foreign export entry cost — ratio compared with actual wage cost
# =============================================================================

table_2_results = zeros(6, 2)
table_2_results[1,:] = [open_CM_open_trade_parameter.L[1], 0]
table_2_results[2,:] = [open_CM_open_trade_parameter.L[2], 0]
table_2_results[3,:] = [open_CM_open_trade_parameter.σ[1], 0]
table_2_results[4,:] = [open_CM_open_trade_parameter.θ[2], 0]
table_2_results[5,:] = [open_CM_open_trade_parameter.δ[1], 0]
table_2_results[6,:] = [open_CM_open_trade_parameter.Fₓ[2],
    prices_open_CM_open_trade[4]^(1/4) / prices_open_CM_open_trade[3]^(1/4) /
    prices_open_CM_open_trade[1] * prices_open_CM_open_trade[2] *
    prices_open_CM_open_trade[5]]

println("[table:Table2] table_2_results computed.")
println(round.(table_2_results, digits=4))

# --- Save Tables/Table2.tex ---
# Only row 6 col 2 is model-computed (the market-size ratio ≈ 6).
# All other cells are parameter values or fixed text.
let
    p = table_2_results[:,1]   # parameter values
    ratio = table_2_results[6,2]   # [Y*/Y]^(1/σ) * w*/w * P  ≈ 6

    # Value-column formatting
    val_L1    = string(round(Int, p[1]))             # NMS pop  (= 1)
    val_L2    = string(round(Int, p[2]))             # Core pop (= 4)
    val_sigma = string(round(Int, p[3]))             # CES σ    (= 4)
    val_theta = @sprintf("%.2f", round(p[4], digits=2, RoundNearestTiesAway))               # Core θ*  (= 0.86)
    val_delta = @sprintf("%.2f", round(p[5], digits=2, RoundNearestTiesAway))               # δ        (= 0.06)
    val_fex   = @sprintf("%.0f \\%%", round(100 * p[6], digits=0, RoundNearestTiesAway))   # f_ex* as % (= 75 \%)
    val_ratio = @sprintf("%.0f", round(ratio, digits=0, RoundNearestTiesAway))              # market-size ratio (= 6)

    tex = """\\begin{table}[ht]
\\centering
\\caption{Preassigned parameters}
\\label{table:Table2}
\\scalebox{0.6}{
\\begin{tabular}{lccc}
\\\\[-1.8ex]\\hline
\\hline \\\\[-1.8ex]
Parameter & Value & Source/Target & Comments \\\\
\\hline \\\\[-1.8ex]
NMS population, \$L\$ & $(val_L1) & - & Normalization \\\\
Core population, \$L^*\$ & $(val_L2) & UN 1989 & Population ratio, Core vs.\\ NMS \\\\
Elasticity of substitution, \$\\sigma\$ & $(val_sigma) & \\cite{trade_elasticity} & Trade, not substitution \\\\
Core borrowing tightness, \$\\theta^*\$ & $(val_theta) & \\cite{10.1257/aer.104.2.422} & Developed countries (Korean) firm data \\\\
Depreciation, \$\\delta\$ & $(val_delta) & \\cite{10.1257/aer.104.2.422} & - \\\\
Avg.\\ export entry cost, \$f_{ex}^*\$ & \$$(val_fex)\$ & \$f_{ex} \\times\$ Market size proportional to cost & \$\\Big[\\frac{Y^*}{Y}\\Big]^{\\frac{1}{\\sigma}} \\frac{w^*}{w} P = $(val_ratio)\$ \\\\
\\hline
\\hline
\\end{tabular}
}
\\parbox{0.88\\textwidth}{\\caption*{\\scriptsize \\textit{Note:} The Core economy is assumed to be more developed, productive, and larger than the NMS economy, but is still comparable. Parameters controlling Core variable export costs or firm productivity are set the same as in NMS, except for export entry cost, which corrects for the size differences across markets.}}
\\end{table}
"""
    open("../Tables/Table2.tex", "w") do io
        write(io, tex)
    end
    println("  Saved → ../Tables/Table2.tex")
end

# =============================================================================
# [fig:Figure4]  Decision rules — capital utilisation and export decision map
# main.tex: \label{fig:Figure4}
# Source: plot.jl lines 187–210
# Prerequisites: Step 1 (main_script.jl).
# Outputs: ../Figures/Figure4a.pdf  (capital utilisation rate map)
#          ../Figures/Figure4b.pdf  (export status decision map)
# =============================================================================

mkpath("../Figures")

# Figure4a — capital utilisation (k_used / k_opt) for exporters
y_value = k_x_fine_open_CM_open_trade[:,1] ./ k_opt_x_fine_open_CM_open_trade[:,1]
y_mat = mat_creator(y_value)
y_mat_applied = 100.0 * y_mat[10:36, 1:100]
contour(y_mat_applied', fill=true, levels=5, c=cgrad(:Blues_6, rev=true),
    xlabel="Productivity", ylabel="Wealth", axis=nothing,
    tickfontsize=14, xguidefontsize=18, yguidefontsize=18,
    legendfontsize=18)
savefig("../Figures/Figure4a.pdf")

# Figure4b — export status decision map
y_value = future_occupation_fine_open_CM_open_trade[:,3,1,1]
entry_dec = (future_occupation_fine_open_CM_open_trade[:,2,1,1] .== 3.0)
y_value[entry_dec] .= 4.0
making_losses_indexes = (making_losses_exporter_open_CM_open_trade[:,1]) .* (y_value .== 3.0)
y_value[making_losses_indexes .== 1.0] .= 2.0
y_mat = mat_creator(y_value)
y_mat = y_mat .- 1.0
y_mat = y_mat * 10
y_mat_applied = y_mat[1:36, 1:30]
contourf(y_mat_applied', levels=5, c=cgrad(:Blues_6, categorical=true, rev=false),
    xlabel="Productivity", ylabel="Wealth", axis=nothing, colorbar=nothing, lw=0,
    tickfontsize=14, xguidefontsize=18, yguidefontsize=18,
    legendfontsize=18)
savefig("../Figures/Figure4b.pdf")

println("[fig:Figure4] Figure4a.pdf and Figure4b.pdf saved → ../Figures/")

# =============================================================================
# [tab:Table4]  Trade liberalization under closed and integrated capital markets
# main.tex: \label{tab:Table4}
# Source: plot.jl lines 395–430
# Prerequisites: Steps 1-2 (transitions for welfare row).
# Columns (6): None / Trade Only / Trade+Capital | 1991 / 2001 / 2008
# Rows (selected from 20-row array):
#   1  TFP (relative)           17  Credit/GDP
#   2  sd MRPK                   8  Welfare
#   3  Output (relative)         9  Top 10% income share
#   5  Consumption (relative)   10  Real wage index
#   6  Capital (relative)       11  Interest rate premium
#   4  Relative output per cap  12  Import/GDP
#  15  Share of exporters        13  NFA/GDP
#  14  Entrepreneurship rate    16  Real exchange rate
# Output: ../Tables/Table4.tex
# =============================================================================

price_union_wide_initial =
    prices_initial[5] * prices_initial[5] * GDP_initial[1] /
    (prices_initial[5] * GDP_initial[1] + GDP_initial[2]) +
    1.0 * GDP_initial[2] / (prices_initial[5] * GDP_initial[1] + GDP_initial[2])
price_union_wide_closed_CM_open_trade =
    prices_closed_CM_open_trade[5] * prices_closed_CM_open_trade[5] * GDP_closed_CM_open_trade[1] /
    (prices_closed_CM_open_trade[5] * GDP_closed_CM_open_trade[1] + GDP_closed_CM_open_trade[2]) +
    1.0 * GDP_closed_CM_open_trade[2] /
    (prices_closed_CM_open_trade[5] * GDP_closed_CM_open_trade[1] + GDP_closed_CM_open_trade[2])
price_union_wide_open_CM_open_trade =
    prices_open_CM_open_trade[5] * prices_open_CM_open_trade[5] * GDP_open_CM_open_trade[1] /
    (prices_open_CM_open_trade[5] * GDP_open_CM_open_trade[1] + GDP_open_CM_open_trade[2]) +
    1.0 * GDP_open_CM_open_trade[2] /
    (prices_open_CM_open_trade[5] * GDP_open_CM_open_trade[1] + GDP_open_CM_open_trade[2])

table_4_results = zeros(20, 6)
table_4_results[1,:] = [TFP_initial[1]/TFP_initial[1], TFP_closed_CM_open_trade[1]/TFP_initial[1],
    TFP_open_CM_open_trade[1]/TFP_initial[1], 100, 116, 132]
table_4_results[2,:] = [sd_MRPK_initial[1], sd_MRPK_closed_CM_open_trade[1],
    sd_MRPK_open_CM_open_trade[1], 0, 1.27, 1.51]
table_4_results[3,:] = [GDP_initial[1]/GDP_initial[1], GDP_closed_CM_open_trade[1]/GDP_initial[1],
    GDP_open_CM_open_trade[1]/GDP_initial[1], 100, 158, 249]
table_4_results[4,:] = [
    prices_initial[5]*GDP_initial[1]/(prices_initial[5]*GDP_initial[1]+GDP_initial[2]) *
        price_union_wide_initial / open_CM_open_trade_parameter.L[1] *
        (open_CM_open_trade_parameter.L[1] + open_CM_open_trade_parameter.L[2]),
    prices_closed_CM_open_trade[5]*GDP_closed_CM_open_trade[1]/
        (prices_closed_CM_open_trade[5]*GDP_closed_CM_open_trade[1]+GDP_closed_CM_open_trade[2]) *
        price_union_wide_closed_CM_open_trade / closed_CM_open_trade_parameter.L[1] *
        (closed_CM_open_trade_parameter.L[1] + closed_CM_open_trade_parameter.L[2]),
    prices_open_CM_open_trade[5]*GDP_open_CM_open_trade[1]/
        (prices_open_CM_open_trade[5]*GDP_open_CM_open_trade[1]+GDP_open_CM_open_trade[2]) *
        price_union_wide_open_CM_open_trade / open_CM_open_trade_parameter.L[1] *
        (open_CM_open_trade_parameter.L[1] + open_CM_open_trade_parameter.L[2]),
    0.54, 0.57, 0.64]
table_4_results[5,:] = [total_consumption_initial[1]/total_consumption_initial[1],
    total_consumption_closed_CM_open_trade[1]/total_consumption_initial[1],
    total_consumption_open_CM_open_trade[1]/total_consumption_initial[1], 100, 116, 133]
table_4_results[6,:] = [capital_demand_initial[1]/capital_demand_initial[1],
    capital_demand_closed_CM_open_trade[1]/capital_demand_initial[1],
    capital_demand_open_CM_open_trade[1]/capital_demand_initial[1], 100, 124, 154]
table_4_results[8,:] = [0, welfare_change_clCM_otrade_trans, welfare_change_oCM_otrade_trans, 0, 0, 0]
table_4_results[9,:] = [p90_income_initial[1], p90_income_closed_CM_open_trade[1],
    p90_income_open_CM_open_trade[1], 0.24, 0.3, 0.34]
table_4_results[10,:] = [prices_initial[1]/prices_initial[5]/(prices_initial[1]/prices_initial[5]),
    prices_closed_CM_open_trade[1]/prices_closed_CM_open_trade[5]/(prices_initial[1]/prices_initial[5]),
    prices_open_CM_open_trade[1]/prices_open_CM_open_trade[5]/(prices_initial[1]/prices_initial[5]),
    0, 100, 144]
table_4_results[11,:] = [prices_initial[6]-prices_initial[7],
    prices_closed_CM_open_trade[6]-prices_closed_CM_open_trade[7],
    prices_open_CM_open_trade[6]-prices_open_CM_open_trade[6], 0.21, 0.05, 0.03]
table_4_results[12,:] = [import_share_initial[1], import_share_closed_CM_open_trade[1],
    import_share_open_CM_open_trade[1], 0.21, 0.42, 0.42]
table_4_results[13,:] = [0, 0, NFA_open_CM_open_trade[1]/nomGDP_open_CM_open_trade[1],
    0.01, -0.05, -0.06]
table_4_results[14,:] = [(domestic_pop_initial[1]+exporter_pop_initial[1]),
    (domestic_pop_closed_CM_open_trade[1]+exporter_pop_closed_CM_open_trade[1]),
    (domestic_pop_open_CM_open_trade[1]+exporter_pop_open_CM_open_trade[1]), 0, 0, 0.07]
table_4_results[15,:] = [exporter_pop_initial[1]/(domestic_pop_initial[1]+exporter_pop_initial[1]),
    exporter_pop_closed_CM_open_trade[1]/(domestic_pop_closed_CM_open_trade[1]+exporter_pop_closed_CM_open_trade[1]),
    exporter_pop_open_CM_open_trade[1]/(domestic_pop_open_CM_open_trade[1]+exporter_pop_open_CM_open_trade[1]),
    0, 0.22, 0.21]
table_4_results[16,:] = [prices_initial[5], prices_closed_CM_open_trade[5],
    prices_open_CM_open_trade[5], 0, prices_closed_CM_open_trade[5],
    prices_closed_CM_open_trade[5]*1.41]

total_credit_initial =
    -(domestic_firm_debt_initial[1] + exporter_firm_debt_initial[1]) / nomGDP_initial[1]
domestic_credit_initial =
    (worker_bond_holding_initial[1] + domestic_bond_holding_initial[1] +
     exporter_bond_holding_initial[1]) / nomGDP_initial[1]
foreign_credit_initial = total_credit_initial - domestic_credit_initial
foreign_credit_share_initial = foreign_credit_initial / total_credit_initial

total_credit_closed_CM_open_trade =
    -(domestic_firm_debt_closed_CM_open_trade[1] + exporter_firm_debt_closed_CM_open_trade[1]) /
    nomGDP_closed_CM_open_trade[1]
domestic_credit_closed_CM_open_trade =
    (worker_bond_holding_closed_CM_open_trade[1] + domestic_bond_holding_closed_CM_open_trade[1] +
     exporter_bond_holding_closed_CM_open_trade[1]) / nomGDP_closed_CM_open_trade[1]
foreign_credit_closed_CM_open_trade =
    total_credit_closed_CM_open_trade - domestic_credit_closed_CM_open_trade
foreign_credit_share_closed_CM_open_trade =
    foreign_credit_closed_CM_open_trade / total_credit_closed_CM_open_trade
table_4_results[17,:] = [total_credit_initial, total_credit_closed_CM_open_trade,
    total_credit_open_CM_open_trade, 0.43, 0.54, 0.62]

println("[tab:Table4] table_4_results (20×6) computed.")

# --- Save Tables/Table4.tex ---
let
    r = table_4_results
    # Columns: None | Trade Only | Trade+Capital || 1991 | 2001 | 2008
    fmt_2f(x)   = @sprintf("%.2f", round(x, digits=2, RoundNearestTiesAway))
    fmt_3f(x)   = @sprintf("%.3f", round(x, digits=3, RoundNearestTiesAway))
    fmt_0f(x)   = @sprintf("%.0f", round(x, digits=0, RoundNearestTiesAway))

    tex = """\\begin{table}[!ht]
\\centering
\\caption{Trade liberalization under closed and integrated capital markets}
\\label{table:Table4}
\\scalebox{0.90}{
\\begin{tabular}{l ccc | ccc c}
\\\\[-1.8ex]\\hline
\\hline
 & \\multicolumn{3}{c|}{Model} & \\multicolumn{3}{c}{Hungarian Data} \\\\
Integration & None & Trade & Trade and capital & 1991 & 2001 & 2008 \\\\
\\hline
\\textbf{Productivity} & & & & & & \\\\
TFP & $(fpct(r[1,1])) & $(fpct(r[1,2])) & $(fpct(r[1,3])) & $(round(Int,r[1,4])) & $(round(Int,r[1,5])) & $(round(Int,r[1,6])) \\\\
s.d. of ARPK & $(fmt_2f(r[2,1])) & $(fmt_2f(r[2,2])) & $(fmt_2f(r[2,3])) & - & $(fmt_2f(r[2,5])) & $(fmt_2f(r[2,6])) \\\\[1ex]
\\textbf{Aggregates} & & & & & & \\\\
Output & $(fpct(r[3,1])) & $(fpct(r[3,2])) & $(fpct(r[3,3])) & $(round(Int,r[3,4])) & $(round(Int,r[3,5])) & $(round(Int,r[3,6])) \\\\
Relative Output pc. & $(fpct(r[4,1])) & $(fpct(r[4,2])) & $(fpct(r[4,3])) & $(fpct(r[4,4])) & $(fpct(r[4,5])) & $(fpct(r[4,6])) \\\\
Consumption & $(@sprintf("%.1f", round(100*r[5,1], digits=1, RoundNearestTiesAway))) & $(@sprintf("%.1f", round(100*r[5,2], digits=1, RoundNearestTiesAway))) & $(@sprintf("%.1f", round(100*r[5,3], digits=1, RoundNearestTiesAway))) & $(round(Int,r[5,4],RoundNearestTiesAway)) & $(round(Int,r[5,5],RoundNearestTiesAway)) & $(round(Int,r[5,6],RoundNearestTiesAway)) \\\\
Capital & $(fpct(r[6,1])) & $(fpct(r[6,2])) & $(fpct(r[6,3])) & $(round(Int,r[6,4])) & $(round(Int,r[6,5])) & $(round(Int,r[6,6])) \\\\[1ex]
\\textbf{Welfare \\& Inequality} & & & & & & \\\\
CE Welfare change & 0 & $(@sprintf("%.1f", round(round(100*r[8,2], digits=2, RoundNearestTiesAway), digits=1, RoundNearestTiesAway))) & $(@sprintf("%.1f", round(round(100*r[8,3], digits=2, RoundNearestTiesAway), digits=1, RoundNearestTiesAway)))* \\\\
Top \$10\\%\$ income share & $(fpct(r[9,1])) & $(fpct(r[9,2])) & $(fpct(r[9,3])) & $(fpct(r[9,4])) & $(fpct(r[9,5])) & $(fpct(r[9,6])) \\\\[1ex]
\\textbf{Factor prices} & & & & & & \\\\
Real wage & $(fpct(r[10,1])) & $(fpct(r[10,2])) & $(fpct(r[10,3])) & - & $(round(Int,r[10,5])) & $(round(Int,r[10,6])) \\\\
Interest rate premium & $(fpct(r[11,1])) & $(fpct(r[11,2])) & $(fpct(r[11,3])) & $(fpct(r[11,4])) & $(fpct(r[11,5])) & $(fpct(r[11,6])) \\\\[1ex]
\\textbf{Trade} & & & & & & \\\\
\$\\frac{\\text{Import}}{\\text{GDP}}\$ & $(fpct(r[12,1])) & $(fpct(r[12,2])) & $(fpct(r[12,3])) & $(fpct(r[12,4])) & $(fpct(r[12,5])) & $(fpct(r[12,6])) \\\\
\$\\frac{\\text{NFA}}{\\text{GDP}}\$ & $(fpct(r[13,1])) & $(fpct(r[13,2])) & $(fpct(r[13,3])) & $(fpct(r[13,4])) & $(fpct(r[13,5])) & $(fpct(r[13,6])) \\\\[0.2ex]
Entrepreneurship rate & $(fpct(r[14,1])) & $(fpct(r[14,2])) & $(fpct(r[14,3])) & - & - & $(fpct(r[14,6])) \\\\
Share of exporters & $(fpct(r[15,1])) & $(fpct(r[15,2])) & $(fpct(r[15,3])) & - & $(fpct(r[15,5])) & $(fpct(r[15,6])) \\\\
Real exchange rate & $(fpct(r[16,1])) & $(fpct(r[16,2])) & $(fpct(r[16,3])) & - & $(fpct(r[16,5])) & $(fpct(r[16,6])) \\\\
\$\\frac{\\text{Credit}}{\\text{GDP}}\$ & $(fpct(r[17,1])) & $(fpct(r[17,2])) & $(fpct(r[17,3])) & $(fpct(r[17,4])) & $(fpct(r[17,5])) & $(fpct(r[17,6])) \\\\
\\hline
\\hline
\\end{tabular}
}
\\parbox{0.90\\textwidth}{\\caption*{\\scriptsize \\textit{Note:} The baseline economy without integration refers to the hypothetical state of the economy before reforms. The middle column shows the post-trade-liberalization, pre-capital-market-integration state of the economy, and the left column indicates the state of the economy after both trade and capital market integration, largely corresponding to 1991, 2001, and 2008 states of the NMS economy. Consumption-equivalent (CE) welfare change accounts for transition dynamics, assuming immediate integration of capital markets. Postponing by 10 years reduces gains to 12.2\\%. Starting from liberalized trade yields welfare gains of 6.2\\%. s.d. of ARPK stands for the standard deviation of the log of the average revenue product of capital. Relative output per capita is measured as NMS's GDP per capita relative to Core GDP per capita. Data counterparts are for Hungary only, with variables converted to the same initial values as in the model (when the value is 100 or for the real exchange rate).}}
\\end{table}
"""
    open("../Tables/Table4.tex", "w") do io
        write(io, tex)
    end
    println("  Saved → ../Tables/Table4.tex")
end

# =============================================================================
# [tab:Table5]  Impact of trade liberalization on firms
# main.tex: \label{tab:Table5}
# Source: plot.jl lines 433–522 (table_5_results) and 634–660 (table_6_results)
# Prerequisites: Step 1.
# Columns (5): None / Trade Only / Trade+Capital | 2001 / 2008
# Output: ../Tables/Table5.tex
# =============================================================================

# table_5_results (22×5)
table_5_results = zeros(22, 5)
table_5_results[1,:] = [1,
    domestic_pop_closed_CM_open_trade[1]/domestic_pop_initial[1],
    domestic_pop_open_CM_open_trade[1]/domestic_pop_initial[1], 0.75,
    42781/29958 * domestic_pop_closed_CM_open_trade[1]/domestic_pop_initial[1]]
table_5_results[2,:] = [1,
    exporter_pop_closed_CM_open_trade[1]/exporter_pop_initial[1],
    exporter_pop_open_CM_open_trade[1]/exporter_pop_initial[1], 1.34,
    11228/8414 * exporter_pop_closed_CM_open_trade[1]/exporter_pop_initial[1]]
table_5_results[3,1] = K_x_initial[1]/(K_x_initial[1]+K_d_initial[1])
table_5_results[3,2] = K_x_closed_CM_open_trade[1]/(K_x_closed_CM_open_trade[1]+K_d_closed_CM_open_trade[1])
table_5_results[3,3] = K_x_open_CM_open_trade[1]/(K_x_open_CM_open_trade[1]+K_d_open_CM_open_trade[1])
table_5_results[3,4] = hfd["pct_cap_2001"]; table_5_results[3,5] = hfd["pct_cap_2008"]
table_5_results[4,1] = L_x_initial[1]/(L_x_initial[1]+L_d_initial[1])
table_5_results[4,2] = L_x_closed_CM_open_trade[1]/(L_x_closed_CM_open_trade[1]+L_d_closed_CM_open_trade[1])
table_5_results[4,3] = L_x_open_CM_open_trade[1]/(L_x_open_CM_open_trade[1]+L_d_open_CM_open_trade[1])
table_5_results[4,4] = hfd["pct_lab_2001"]; table_5_results[4,5] = hfd["pct_lab_2008"]
table_5_results[5,:] = [exporter_pop_initial[1]/entry_share_to_exporter_initial[1],
    exporter_pop_closed_CM_open_trade[1]/entry_share_to_exporter_closed_CM_open_trade[1],
    exporter_pop_open_CM_open_trade[1]/entry_share_to_exporter_open_CM_open_trade[1],
    hfd["dur_2001"], hfd["dur_2008"]]
table_5_results[6,1:3] = [1,
    sum(k_d_fine_closed_CM_open_trade[:,1].*current_domestic_closed_CM_open_trade)/
        domestic_pop_closed_CM_open_trade[1] * domestic_pop_initial[1] /
        sum(k_d_fine_initial[:,1].*current_domestic_initial),
    sum(k_d_fine_open_CM_open_trade[:,1].*current_domestic_open_CM_open_trade)/
        domestic_pop_open_CM_open_trade[1] * domestic_pop_initial[1] /
        sum(k_d_fine_initial[:,1].*current_domestic_initial)]
table_5_results[6,4] = table_5_results[6,2]
table_5_results[6,5] = (hfd["cap_d_2008"] / hfd["cap_d_2001"]) * table_5_results[6,4]
table_5_results[7,1:3] = [1,
    sum(k_x_fine_closed_CM_open_trade[:,1].*current_exporter_closed_CM_open_trade)/
        exporter_pop_closed_CM_open_trade[1] * exporter_pop_initial[1] /
        sum(k_x_fine_initial[:,1].*current_exporter_initial),
    sum(k_x_fine_open_CM_open_trade[:,1].*current_exporter_open_CM_open_trade)/
        exporter_pop_open_CM_open_trade[1] * exporter_pop_initial[1] /
        sum(k_x_fine_initial[:,1].*current_exporter_initial)]
table_5_results[7,4] = table_5_results[7,2]
table_5_results[7,5] = (hfd["cap_x_2008"] / hfd["cap_x_2001"]) * table_5_results[7,4]
table_5_results[8,1:3] = [1,
    sum(k_x_fine_closed_CM_open_trade[:,1].*current_exporter_closed_CM_open_trade +
        k_d_fine_closed_CM_open_trade[:,1].*current_domestic_closed_CM_open_trade) /
        (domestic_pop_closed_CM_open_trade[1]+exporter_pop_closed_CM_open_trade[1]) *
        (domestic_pop_initial[1]+exporter_pop_initial[1]) /
        sum(k_x_fine_initial[:,1].*current_exporter_initial + k_d_fine_initial[:,1].*current_domestic_initial),
    sum(k_x_fine_open_CM_open_trade[:,1].*current_exporter_open_CM_open_trade +
        k_d_fine_open_CM_open_trade[:,1].*current_domestic_open_CM_open_trade) /
        (domestic_pop_open_CM_open_trade[1]+exporter_pop_open_CM_open_trade[1]) *
        (domestic_pop_initial[1]+exporter_pop_initial[1]) /
        sum(k_x_fine_initial[:,1].*current_exporter_initial + k_d_fine_initial[:,1].*current_domestic_initial)]
table_5_results[8,4] = table_5_results[8,2]
table_5_results[8,5] = (hfd["cap_all_2008"] / hfd["cap_all_2001"]) * table_5_results[8,4]
table_5_results[9,:]  = [sd_MRPK_d_initial[1], sd_MRPK_d_closed_CM_open_trade[1],
    sd_MRPK_d_open_CM_open_trade[1], hfd["sd_arpk_d_2001"], hfd["sd_arpk_d_2008"]]
table_5_results[10,:] = [sd_MRPK_x_initial[1], sd_MRPK_x_closed_CM_open_trade[1],
    sd_MRPK_x_open_CM_open_trade[1], hfd["sd_arpk_x_2001"], hfd["sd_arpk_x_2008"]]

export_revenue_total_initial =
    sum(rev_xx_fine_store_initial[:,1] .* current_exporter_initial[:,1])
domestic_revenue_total_initial =
    sum(rev_dx_fine_store_initial[:,1] .* current_exporter_initial[:,1] +
        rev_d_fine_store_initial[:,1] .* current_domestic_initial[:,1])
export_sales_ratio_initial = export_revenue_total_initial /
    (export_revenue_total_initial + domestic_revenue_total_initial)

export_revenue_total_closed_CM_open_trade =
    sum(rev_xx_fine_store_closed_CM_open_trade[:,1] .* current_exporter_closed_CM_open_trade[:,1])
domestic_revenue_total_closed_CM_open_trade =
    sum(rev_dx_fine_store_closed_CM_open_trade[:,1] .* current_exporter_closed_CM_open_trade[:,1] +
        rev_d_fine_store_closed_CM_open_trade[:,1] .* current_domestic_closed_CM_open_trade[:,1])
export_sales_ratio_closed_CM_open_trade = export_revenue_total_closed_CM_open_trade /
    (export_revenue_total_closed_CM_open_trade + domestic_revenue_total_closed_CM_open_trade)

export_revenue_total_open_CM_open_trade =
    sum(rev_xx_fine_store_open_CM_open_trade[:,1] .* current_exporter_open_CM_open_trade[:,1])
domestic_revenue_total_open_CM_open_trade =
    sum(rev_dx_fine_store_open_CM_open_trade[:,1] .* current_exporter_open_CM_open_trade[:,1] +
        rev_d_fine_store_open_CM_open_trade[:,1] .* current_domestic_open_CM_open_trade[:,1])
export_sales_ratio_open_CM_open_trade = export_revenue_total_open_CM_open_trade /
    (export_revenue_total_open_CM_open_trade + domestic_revenue_total_open_CM_open_trade)

table_5_results[11,:] = [export_sales_ratio_initial, export_sales_ratio_closed_CM_open_trade,
    export_sales_ratio_open_CM_open_trade, hfd["expsales_2001"], hfd["expsales_2008"]]
table_5_results[12,1:3] = table_4_results[15,1:3]
table_5_results[12,4:5] = table_4_results[15,5:6]

exs_initial = rev_xx_fine_store_initial[:,1] ./
    (rev_xx_fine_store_initial[:,1] + rev_dx_fine_store_initial[:,1])
total_revenue_exporters_initial = rev_xx_fine_store_initial[:,1] + rev_dx_fine_store_initial[:,1]
intensive_margin_initial = sum(exs_initial .* total_revenue_exporters_initial) /
    sum(total_revenue_exporters_initial)

exs_closed_CM_open_trade = rev_xx_fine_store_closed_CM_open_trade[:,1] ./
    (rev_xx_fine_store_closed_CM_open_trade[:,1] + rev_dx_fine_store_closed_CM_open_trade[:,1])
total_revenue_exporters_closed_CM_open_trade =
    rev_xx_fine_store_closed_CM_open_trade[:,1] + rev_dx_fine_store_closed_CM_open_trade[:,1]
intensive_margin_closed_CM_open_trade =
    sum(exs_closed_CM_open_trade .* total_revenue_exporters_closed_CM_open_trade) /
    sum(total_revenue_exporters_closed_CM_open_trade)

exs_open_CM_open_trade = rev_xx_fine_store_open_CM_open_trade[:,1] ./
    (rev_xx_fine_store_open_CM_open_trade[:,1] + rev_dx_fine_store_open_CM_open_trade[:,1])
total_revenue_exporters_open_CM_open_trade =
    rev_xx_fine_store_open_CM_open_trade[:,1] + rev_dx_fine_store_open_CM_open_trade[:,1]
intensive_margin_open_CM_open_trade =
    sum(exs_open_CM_open_trade .* total_revenue_exporters_open_CM_open_trade) /
    sum(total_revenue_exporters_open_CM_open_trade)

table_5_results[13,:] = [intensive_margin_initial, intensive_margin_closed_CM_open_trade,
    intensive_margin_open_CM_open_trade, hfd["int_2001"], hfd["int_2008"]]
table_5_results[14,:] = [
    sum(total_revenue_exporters_initial)/exporter_pop_initial[1] /
        sum(rev_d_fine_store_initial[:,1]) * domestic_pop_initial[1],
    sum(total_revenue_exporters_closed_CM_open_trade)/exporter_pop_closed_CM_open_trade[1] /
        sum(rev_d_fine_store_closed_CM_open_trade[:,1]) * domestic_pop_closed_CM_open_trade[1],
    sum(total_revenue_exporters_open_CM_open_trade)/exporter_pop_open_CM_open_trade[1] /
        sum(rev_d_fine_store_open_CM_open_trade[:,1]) * domestic_pop_open_CM_open_trade[1],
    hfd["premium_2001"], hfd["premium_2008"]]
table_5_results[15,:] = [
    entry_share_to_exporter_initial[1]/(entry_share_to_exporter_initial[1]+domestic_pop_initial[1]),
    entry_share_to_exporter_closed_CM_open_trade[1]/
        (exporter_pop_closed_CM_open_trade[1]+domestic_pop_closed_CM_open_trade[1]),
    entry_share_to_exporter_open_CM_open_trade[1]/
        (exporter_pop_open_CM_open_trade[1]+domestic_pop_open_CM_open_trade[1]),
    hfd["starter_2001"], hfd["starter_2008"]]
table_5_results[16,:] = [
    exit_exporting_to_work_sum_initial[1]/exporter_pop_initial[1],
    exit_exporting_to_work_sum_closed_CM_open_trade[1]/exporter_pop_closed_CM_open_trade[1],
    exit_exporting_to_work_sum_open_CM_open_trade[1]/exporter_pop_open_CM_open_trade[1],
    hfd["stopper_2001"], hfd["stopper_2008"]]
table_5_results[17,:] = [Misallocation_within_d_initial[1], Misallocation_within_d_closed_CM_open_trade[1],
    Misallocation_within_d_open_CM_open_trade[1], 0, 0]
table_5_results[18,:] = [Misallocation_within_x_initial[1], Misallocation_within_x_closed_CM_open_trade[1],
    Misallocation_within_x_open_CM_open_trade[1], 0, 0]

# table_6_results (4×3) — distribution of exporters by wealth and productivity
current_exporter_initial_mat   = mat_creator(current_exporter_initial)
current_domestic_initial_mat   = mat_creator(current_domestic_initial)
current_producer_initial_mat   = current_exporter_initial_mat + current_domestic_initial_mat
y_start = 5; x_start = 18
table_6_results = zeros(4, 3)
table_6_results[4,1] = sum(current_exporter_initial_mat[(x_start+1):end,(y_start+1):end])/sum(current_exporter_initial)
table_6_results[3,1] = sum(current_exporter_initial_mat[1:x_start,(y_start+1):end])/sum(current_exporter_initial)
table_6_results[2,1] = sum(current_exporter_initial_mat[(x_start+1):end,1:y_start])/sum(current_exporter_initial)
table_6_results[1,1] = sum(current_exporter_initial_mat[1:x_start,1:y_start])/sum(current_exporter_initial)
current_exporter_closed_CM_mat = mat_creator(current_exporter_closed_CM_open_trade)
table_6_results[4,2] = sum(current_exporter_closed_CM_mat[(x_start+1):end,(y_start+1):end])/sum(current_exporter_closed_CM_open_trade)
table_6_results[3,2] = sum(current_exporter_closed_CM_mat[1:x_start,(y_start+1):end])/sum(current_exporter_closed_CM_open_trade)
table_6_results[2,2] = sum(current_exporter_closed_CM_mat[(x_start+1):end,1:y_start])/sum(current_exporter_closed_CM_open_trade)
table_6_results[1,2] = sum(current_exporter_closed_CM_mat[1:x_start,1:y_start])/sum(current_exporter_closed_CM_open_trade)
current_exporter_open_CM_mat   = mat_creator(current_exporter_open_CM_open_trade)
table_6_results[4,3] = sum(current_exporter_open_CM_mat[(x_start+1):end,(y_start+1):end])/sum(current_exporter_open_CM_open_trade)
table_6_results[3,3] = sum(current_exporter_open_CM_mat[1:x_start,(y_start+1):end])/sum(current_exporter_open_CM_open_trade)
table_6_results[2,3] = sum(current_exporter_open_CM_mat[(x_start+1):end,1:y_start])/sum(current_exporter_open_CM_open_trade)
table_6_results[1,3] = sum(current_exporter_open_CM_mat[1:x_start,1:y_start])/sum(current_exporter_open_CM_open_trade)

println("[tab:Table5] table_5_results (22×5) and table_6_results (4×3) computed.")

# --- Save Tables/Table5.tex ---
let
    r5 = table_5_results
    r6 = table_6_results
    fmt2(x) = @sprintf("%.2f", round(x, digits=2, RoundNearestTiesAway))
    fmt3(x) = @sprintf("%.3f", round(x, digits=3, RoundNearestTiesAway))

    tex = """\\begin{table}[!ht]
\\centering
\\caption{Impact of trade liberalization on firms}
\\label{table:Table5}
\\scalebox{0.90}{
\\begin{tabular}{l ccc | ccc c}
\\\\[-1.8ex]\\hline
\\hline
 & \\multicolumn{3}{c|}{Model} & \\multicolumn{2}{c}{Data} \\\\
Integration & None & Trade & Trade and capital & 2001 & 2008 \\\\
\\hline
\\textbf{Extensive margin} & & & & & \\\\
Non-exporting firms & $(fpct(r5[1,1])) & $(fpct(r5[1,2])) & $(fpct(r5[1,3])) & $(fpct(r5[1,4])) & $(fpct(r5[1,5])) \\\\
Exporting firms & $(fpct(r5[2,1])) & $(fpct(r5[2,2])) & $(fpct(r5[2,3])) & $(fpct(r5[2,4])) & $(fpct(r5[2,5])) \\\\[1ex]
\\textbf{Intensive margin} & & & & & \\\\
\$ \\%\$ of capital used by exporters & $(fpct(r5[3,1])) & $(fpct(r5[3,2])) & $(fpct(r5[3,3])) & $(fpct(r5[3,4])) & $(fpct(r5[3,5])) \\\\
\$ \\%\$ of labor used by exporters & $(fpct(r5[4,1])) & $(fpct(r5[4,2])) & $(fpct(r5[4,3])) & $(fpct(r5[4,4])) & $(fpct(r5[4,5])) \\\\
Avg. duration (years) of export status & $(@sprintf("%.1f", round(r5[5,1], digits=1, RoundNearestTiesAway))) & $(@sprintf("%.1f", round(r5[5,2], digits=1, RoundNearestTiesAway))) & $(@sprintf("%.1f", round(r5[5,3], digits=1, RoundNearestTiesAway))) & $(@sprintf("%.1f", round(r5[5,4], digits=1, RoundNearestTiesAway))) & $(@sprintf("%.1f", round(r5[5,5], digits=1, RoundNearestTiesAway))) \\\\
Average capital size of non-exporters & $(fpct(r5[6,1])) & $(fpct(r5[6,2])) & $(fpct(r5[6,3])) & $(fpct(r5[6,4])) & $(fpct(r5[6,5])) \\\\
Average capital size of exporters & $(fpct(r5[7,1])) & $(fpct(r5[7,2])) & $(fpct(r5[7,3])) & $(fpct(r5[7,4])) & $(fpct(r5[7,5])) \\\\
Average capital size & $(fpct(r5[8,1])) & $(fpct(r5[8,2])) & $(fpct(r5[8,3])) & $(fpct(r5[8,4])) & $(fpct(r5[8,5])) \\\\[1ex]
\\textbf{Standard deviation of ARPK} & & & & & \\\\
 Non-exporting & $(fmt2(r5[9,1])) & $(fmt2(r5[9,2])) & $(fmt2(r5[9,3])) & $(@sprintf("%.1f", r5[9,4])) & $(fmt2(r5[9,5])) \\\\
 Exporter & $(fmt2(r5[10,1])) & $(fmt2(r5[10,2])) & $(fmt2(r5[10,3])) & $(fmt2(r5[10,4])) & $(fmt2(r5[10,5])) \\\\[1ex]
\\textbf{Export sales decomposition} & & & & & \\\\
Export-sales ratio & $(fpct(r5[11,1])) & $(fpct(r5[11,2])) & $(fpct(r5[11,3])) & $(fpct(r5[11,4])) & $(fpct(r5[11,5])) \\\\
Extensive margin & $(fpct(r5[12,1])) & $(fpct(r5[12,2])) & $(fpct(r5[12,3])) & $(fpct(r5[12,4])) & $(fpct(r5[12,5])) \\\\
Intensive margin & $(fpct(r5[13,1])) & $(fpct(r5[13,2])) & $(fpct(r5[13,3])) & $(fpct(r5[13,4])) & $(fpct(r5[13,5])) \\\\
Exporter premium & $(fpct(r5[14,1])) & $(fpct(r5[14,2])) & $(fpct(r5[14,3])) & $(fpct(r5[14,4])) & $(fpct(r5[14,5])) \\\\
Starter rate & $(fpct(r5[15,1])) & $(fpct(r5[15,2])) & $(fpct(r5[15,3])) & $(fpct(r5[15,4])) & $(fpct(r5[15,5])) \\\\
Stopper rate & $(fpct(r5[16,1])) & $(fpct(r5[16,2])) & $(fpct(r5[16,3])) & $(fpct(r5[16,4])) & $(fpct(r5[16,5])) \\\\
\\textbf{Productivity loss} & & & & & \\\\
 Non-exporting & $(@sprintf("%.1f", round(100*r5[17,1], digits=1, RoundNearestTiesAway))) & $(@sprintf("%.1f", round(100*r5[17,2], digits=1, RoundNearestTiesAway))) & $(@sprintf("%.1f", round(100*r5[17,3], digits=1, RoundNearestTiesAway))) \\\\
 Exporter & $(@sprintf("%.1f", round(100*r5[18,1], digits=1, RoundNearestTiesAway))) & $(@sprintf("%.1f", round(100*r5[18,2], digits=1, RoundNearestTiesAway))) & $(@sprintf("%.1f", round(100*r5[18,3], digits=1, RoundNearestTiesAway))) \\\\
\\textbf{Distribution of exporters} & & & & & \\\\
 Low wealth and low productivity & $(fpct(r6[1,1])) & $(fpct(r6[1,2])) & $(fpct(r6[1,3])) & \\multirow{2}{*}{13} & \\multirow{2}{*}{5} \\\\
Low wealth and high productivity & $(fpct(r6[2,1])) & $(fpct(r6[2,2])) & $(fpct(r6[2,3])) \\\\
High wealth and low productivity & $(fpct(r6[3,1])) & $(fpct(r6[3,2])) & $(fpct(r6[3,3])) & \\multirow{2}{*}{87} & \\multirow{2}{*}{95} \\\\
High wealth and high productivity & $(fpct(r6[4,1])) & $(fpct(r6[4,2])) & $(fpct(r6[4,3])) \\\\
\\hline
\\hline
\\end{tabular}
}
\\parbox{0.95\\textwidth}{\\caption*{\\scriptsize \\textit{Note:} Standard deviation of average capital productivity is calculated within export status \\(\\{d, ex\\}\\). Productivity loss within export status is defined as the percentage deviation of \\(TFP_{i}\\) from \\(TFP_i^{e}\\), defined as the productivity of firms with the same export status \\(\\{d, ex\\}\\) that would occur if all firms were able to obtain their unconstrained input choice. For the distribution of exporters, each column sums to 1. The productivity threshold is the mean entrepreneurial productivity across all households, including non-firm operators. The wealth threshold is the median wealth of all entrepreneurs in a closed capital-market, liberalized-trade steady-state economy, for comparison with Figure \\ref{fig:DiD_combined}. Export sales decomposition follows \\cite{annurev:/content/journals/10.1146/annurev-economics-090919-025159} table 1 for both the data and the model-simulated data. The number of firms and capital size are normalized to their model counterparts, with 2001 values scaled to steady-state values after trade liberalization without capital market integration.}}
\\end{table}
"""
    open("../Tables/Table5.tex", "w") do io
        write(io, tex)
    end
    println("  Saved → ../Tables/Table5.tex")
end

# =============================================================================
# [tab:Table_regs]  GLM regressions: entry cost vs per-period cost model
# main.tex: \label{tab:Table_regs}
# Source: plot.jl lines 524–616 (baseline/entry-cost) + plot_FC.jl lines 223–315 (FC/per-period)
# Prerequisites: Steps 1-4 (needs both BL and FC SS variables).
#
# Uses *_BL aliases for the entry-cost columns (saved in Step 2b before
# fixed_cost_steady_states.jl overwrote the originals).
# Uses direct variable names for FC columns (which after Step 3 hold FC values).
#
# Columns (8):
#   Model: Entry cost (1) | Entry cost (2) | Per-period cost (1) | Per-period cost (2)
#   Data:  (1) | (2) | (3) | (4)
# Rows: Log value added coef / s.e. / Prev. export status coef / s.e. / N / Adj. R²
# Output: ../Tables/Table_regs.tex
# =============================================================================

# --- Baseline (entry-cost) GLM regressions using *_BL variables ---
prices_BL          = prices_open_CM_open_trade   # prices not overwritten
parameters_BL      = open_CM_open_trade_parameter
country_BL         = 1

prob_distr_BL = [worker_past_dist_BL .* future_occupation_fine_avgs_BL[1,3]
    domestic_past_dist_BL .* future_occupation_fine_avgs_BL[2,3]
    incumbents_exporter_BL
    entrants_domestic_from_workers_BL
    incumbents_domestic_BL
    exit_exporting_to_domestic_BL] ./ (1 - worker_pop_tmp_BL)
prob_distr_BL = convert(Array, prob_distr_BL)[:,1]

weights_BL = Weights(prob_distr_BL)
N_sample_BL  = 408554
N_sample2_BL = 340429

exportstatus_BL = zeros(size(prob_distr_BL))
exportstatus_BL[1:21600] .= 1.0
export_share_BL = zeros(size(prob_distr_BL))
export_share_BL[1:ns_fine_BL] = rev_xx_fine_BL ./ (rev_xx_fine_BL + rev_dx_fine_BL)
export_share_BL[(ns_fine_BL+1):(2*ns_fine_BL)] = rev_xx_fine_BL ./ (rev_xx_fine_BL + rev_dx_fine_BL)
export_share_BL[(2*ns_fine_BL+1):(3*ns_fine_BL)] = rev_xx_fine_BL ./ (rev_xx_fine_BL + rev_dx_fine_BL)
value_added_BL = zeros(size(prob_distr_BL))
value_added_BL[1:ns_fine_BL] = rev_xx_fine_BL + rev_dx_fine_BL
value_added_BL[(ns_fine_BL+1):(2*ns_fine_BL)] = rev_xx_fine_BL + rev_dx_fine_BL
value_added_BL[(2*ns_fine_BL+1):(3*ns_fine_BL)] = rev_xx_fine_BL + rev_dx_fine_BL
value_added_BL[(3*ns_fine_BL+1):(4*ns_fine_BL)] = rev_d_fine_BL
value_added_BL[(4*ns_fine_BL+1):(5*ns_fine_BL)] = rev_d_fine_BL
value_added_BL[(5*ns_fine_BL+1):(6*ns_fine_BL)] = rev_d_fine_BL
log_value_added_BL = log.(value_added_BL)

Q_trans_forward_BL = Q_trans_BL * Q_trans_BL
sum_prob_nonex_to_ex_BL = sum(Q_trans_forward_BL[(ns_fine_BL+1):(2*ns_fine_BL),(2*ns_fine_BL+1):(3*ns_fine_BL)], dims=2)
sum_prob_ex_to_ex_BL    = sum(Q_trans_forward_BL[(2*ns_fine_BL+1):(3*ns_fine_BL),(2*ns_fine_BL+1):(3*ns_fine_BL)], dims=2)

future_export_status_BL = zeros(size(prob_distr_BL))
future_export_status_BL[1:ns_fine_BL] = sum_prob_ex_to_ex_BL
future_export_status_BL[(ns_fine_BL+1):(2*ns_fine_BL)] = sum_prob_ex_to_ex_BL
future_export_status_BL[(2*ns_fine_BL+1):(3*ns_fine_BL)] = sum_prob_ex_to_ex_BL
future_export_status_BL[(3*ns_fine_BL+1):(4*ns_fine_BL)] = sum_prob_nonex_to_ex_BL
future_export_status_BL[(4*ns_fine_BL+1):(5*ns_fine_BL)] = sum_prob_nonex_to_ex_BL
future_export_status_BL[(5*ns_fine_BL+1):(6*ns_fine_BL)] = sum_prob_nonex_to_ex_BL

nume_simu = 100
m1BL_coef      = zeros(nume_simu, 2); m1BL_se = zeros(nume_simu, 2); m1BL_r2 = zeros(nume_simu)
m2BL_coef      = zeros(nume_simu, 3); m2BL_se = zeros(nume_simu, 3); m2BL_r2 = zeros(nume_simu)

for ii = 1:nume_simu
    Atmp_BL = [exportstatus_BL export_share_BL log_value_added_BL future_export_status_BL]
    indx = sample(axes(Atmp_BL, 1), weights_BL, N_sample_BL)
    sd = Atmp_BL[indx, :]
    df = DataFrame(log_value_added=sd[:,3], exportstatus=sd[:,1],
                   export_share=sd[:,2], future_export_status=sd[:,4])
    mdl = lm(@formula(future_export_status ~ 1 + log_value_added), df)
    m1BL_coef[ii,:] = coef(mdl); m1BL_se[ii,:] = stderror(mdl)
    n = nobs(mdl); m1BL_r2[ii] = 1 - (1 - r2(mdl))*(n-1)/(n-3)
end
for ii = 1:nume_simu
    Atmp_BL = [exportstatus_BL export_share_BL log_value_added_BL future_export_status_BL]
    indx2 = sample(axes(Atmp_BL, 1), weights_BL, N_sample2_BL)
    sd2 = Atmp_BL[indx2, :]
    df2 = DataFrame(log_value_added=sd2[:,3], exportstatus=sd2[:,1],
                    export_share=sd2[:,2], future_export_status=sd2[:,4])
    mdl2 = lm(@formula(future_export_status ~ 1 + log_value_added + exportstatus), df2)
    m2BL_coef[ii,:] = coef(mdl2); m2BL_se[ii,:] = stderror(mdl2)
    n2 = nobs(mdl2); m2BL_r2[ii] = 1 - (1 - r2(mdl2))*(n2-1)/(n2-4)
end

# --- FC (per-period cost) GLM regressions ---
# Derive all needed variables from fixed_cost_steady_states.jl outputs
# (open_CM_open_trade_FC, country 1 = NMS).
ns_fine_FC    = open_CM_open_trade_parameter_FC.Country_spec_p[1][4]
distr_c1_FC   = Vector(distr_current_open_CM_open_trade_FC[:, 1])
worker_past_dist_FC   = distr_c1_FC[1:ns_fine_FC]
domestic_past_dist_FC = distr_c1_FC[(ns_fine_FC+1):(2*ns_fine_FC)]
exporter_past_dist_FC = distr_c1_FC[(2*ns_fine_FC+1):(3*ns_fine_FC)]

_iid_prob_FC  = open_CM_open_trade_parameter_FC.iid_cost_prob
_niid_FC      = length(_iid_prob_FC)
_focc_FC      = Array(future_occupation_fine_open_CM_open_trade_FC[:, :, 1, :])
_focc_avgs_FC = Array{Vector{Float64}}(undef, 3, 3)
for _j = 1:3, _jj = 1:3
    _focc_avgs_FC[_j, _jj] = zeros(ns_fine_FC)
    for _jc = 1:_niid_FC
        _focc_avgs_FC[_j, _jj] .+= _iid_prob_FC[_jc] .* (_focc_FC[:, _j, _jc] .== _jj)
    end
end

incumbents_exporter_FC            = exporter_past_dist_FC .* _focc_avgs_FC[3, 3]
entrants_domestic_from_workers_FC = worker_past_dist_FC   .* _focc_avgs_FC[1, 2]
incumbents_domestic_FC            = domestic_past_dist_FC .* _focc_avgs_FC[2, 2]
exit_exporting_to_domestic_FC     = exporter_past_dist_FC .* _focc_avgs_FC[3, 2]
worker_pop_tmp_FC                 = worker_pop_open_CM_open_trade_FC[1]

rev_xx_fine_FC = Vector(rev_xx_fine_store_open_CM_open_trade_FC[:, 1])
rev_dx_fine_FC = Vector(rev_dx_fine_store_open_CM_open_trade_FC[:, 1])
rev_d_fine_FC  = Vector(rev_d_fine_store_open_CM_open_trade_FC[:, 1])

Q_trans_FC = let
    _country = 1
    (_β, _α, _δ, _θ, _α₁, _α₂, _σ, _α₁eff, _α₂eff, _ω, _L,
     _FM, _FMₓ, _F, _Fₓ, _Exit, _Exitₓ, _iid_val, _iid_prob,
     _niid, _country_no, _τ, _a_min, _a_max, _fspace_a, _fspace_a_fine,
     _agrid, _bank_cost, _bounds, _CSp, _s_cell, _ns_cell, _sf_cell,
     _nsf_cell, _Phiz_cell, _Phizf_cell, _Phi_cell, _Phiaug_cell,
     _Pkron_cell, _Pkron1_cell, _Pkronf_cell, _expeg_cell,
     _ns_tmp, _nsf_tmp, _openness) = local_parameters(open_CM_open_trade_parameter_FC)

    (_W, _r, _out_final, _pf, _, _) = price_reshaper(
        prices_open_CM_open_trade_FC, _openness, _country_no, _bounds)

    (_, _s, _ns, _sf, _nsf, _Phiz, _Phizf, _Phi, _Phiaug,
     _Pk, _Pk1, _Pkf, _, _, _imat, _imatf,
     _apf, _foccf, _conf, _coeff, _, _pfprev, _pfcurr,
     _, _Wloc, _θloc, _, _τloc, _rloc,
     _bcloc, _, _R, _cd, _cx, _cz, _czbar, _czf, _czfbar,
     _ones, _, _, _βloc, _, _, _, _, _, _Qt,
     _onesf) = country_local_initialize(
        _country, _s_cell, _ns_cell, _sf_cell, _nsf_cell,
        _Phiz_cell, _Phizf_cell, _Phi_cell, _Phiaug_cell,
        _Pkron_cell, _Pkron1_cell, _Pkronf_cell,
        _σ, _niid, _pf, _W, _r, _out_final,
        _θ, _L, _τ, _country_no, _bank_cost, _δ, _ω, _β)

    _coeff_conv = Array(coeff_final_open_CM_open_trade_FC[:, :, _country])

    (_Pd, _Px, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _) =
        true_profit(_pfprev, _R, _Wloc, _sf,
            _cx, _cd, _czf, _czfbar,
            _θloc, _τloc, _σ, _α₁eff, _α₂eff, _α₁, _α₂, _onesf)

    (_imatf, _max_xf) = income_creator(_Wloc, _onesf, _Pd, _Px,
        _FM, _FMₓ, _country, _rloc, _sf, _iid_val, _F, _Fₓ,
        _niid, _imatf, _nsf, _a_min, _a_max)

    (_Qtp, _, _, _, _) = Q_transition(_coeff_conv,
        _nsf, _onesf, _niid, _iid_prob, _imatf,
        _Pkf, _a_min, _Phizf, _βloc,
        _fspace_a, _fspace_a_fine, _pfcurr, _max_xf, _Pk1,
        _conf, _apf, _foccf, _Qt)

    collect(sparse(_Qtp'))
end

prob_distr_FC = [worker_past_dist_FC .* _focc_avgs_FC[1, 3]
    domestic_past_dist_FC .* _focc_avgs_FC[2, 3]
    incumbents_exporter_FC
    entrants_domestic_from_workers_FC
    incumbents_domestic_FC
    exit_exporting_to_domestic_FC] ./ (1 - worker_pop_tmp_FC)
prob_distr_FC = convert(Array, prob_distr_FC)[:, 1]

weights_FC = Weights(prob_distr_FC)
N_sample_FC  = 408554
N_sample2_FC = 340429

exportstatus_FC = zeros(size(prob_distr_FC))
exportstatus_FC[1:21600] .= 1.0
export_share_FC = zeros(size(prob_distr_FC))
export_share_FC[1:ns_fine_FC] = rev_xx_fine_FC ./ (rev_xx_fine_FC + rev_dx_fine_FC)
export_share_FC[(ns_fine_FC+1):(2*ns_fine_FC)] = rev_xx_fine_FC ./ (rev_xx_fine_FC + rev_dx_fine_FC)
export_share_FC[(2*ns_fine_FC+1):(3*ns_fine_FC)] = rev_xx_fine_FC ./ (rev_xx_fine_FC + rev_dx_fine_FC)
value_added_FC = zeros(size(prob_distr_FC))
value_added_FC[1:ns_fine_FC] = rev_xx_fine_FC + rev_dx_fine_FC
value_added_FC[(ns_fine_FC+1):(2*ns_fine_FC)] = rev_xx_fine_FC + rev_dx_fine_FC
value_added_FC[(2*ns_fine_FC+1):(3*ns_fine_FC)] = rev_xx_fine_FC + rev_dx_fine_FC
value_added_FC[(3*ns_fine_FC+1):(4*ns_fine_FC)] = rev_d_fine_FC
value_added_FC[(4*ns_fine_FC+1):(5*ns_fine_FC)] = rev_d_fine_FC
value_added_FC[(5*ns_fine_FC+1):(6*ns_fine_FC)] = rev_d_fine_FC
log_value_added_FC = log.(value_added_FC)

Q_trans_forward_FC = Q_trans_FC * Q_trans_FC
sum_prob_nonex_to_ex_FC = sum(Q_trans_forward_FC[(ns_fine_FC+1):(2*ns_fine_FC),(2*ns_fine_FC+1):(3*ns_fine_FC)], dims=2)
sum_prob_ex_to_ex_FC    = sum(Q_trans_forward_FC[(2*ns_fine_FC+1):(3*ns_fine_FC),(2*ns_fine_FC+1):(3*ns_fine_FC)], dims=2)

future_export_status_FC = zeros(size(prob_distr_FC))
future_export_status_FC[1:ns_fine_FC] = sum_prob_ex_to_ex_FC
future_export_status_FC[(ns_fine_FC+1):(2*ns_fine_FC)] = sum_prob_ex_to_ex_FC
future_export_status_FC[(2*ns_fine_FC+1):(3*ns_fine_FC)] = sum_prob_ex_to_ex_FC
future_export_status_FC[(3*ns_fine_FC+1):(4*ns_fine_FC)] = sum_prob_nonex_to_ex_FC
future_export_status_FC[(4*ns_fine_FC+1):(5*ns_fine_FC)] = sum_prob_nonex_to_ex_FC
future_export_status_FC[(5*ns_fine_FC+1):(6*ns_fine_FC)] = sum_prob_nonex_to_ex_FC

m1FC_coef = zeros(nume_simu, 2); m1FC_se = zeros(nume_simu, 2); m1FC_r2 = zeros(nume_simu)
m2FC_coef = zeros(nume_simu, 3); m2FC_se = zeros(nume_simu, 3); m2FC_r2 = zeros(nume_simu)

for ii = 1:nume_simu
    Atmp_FC = [exportstatus_FC export_share_FC log_value_added_FC future_export_status_FC]
    indx = sample(axes(Atmp_FC, 1), weights_FC, N_sample_FC)
    sd = Atmp_FC[indx, :]
    df = DataFrame(log_value_added=sd[:,3], exportstatus=sd[:,1],
                   export_share=sd[:,2], future_export_status=sd[:,4])
    mdl = lm(@formula(future_export_status ~ 1 + log_value_added), df)
    m1FC_coef[ii,:] = coef(mdl); m1FC_se[ii,:] = stderror(mdl)
    n = nobs(mdl); m1FC_r2[ii] = 1 - (1 - r2(mdl))*(n-1)/(n-3)
end
for ii = 1:nume_simu
    Atmp_FC = [exportstatus_FC export_share_FC log_value_added_FC future_export_status_FC]
    indx2 = sample(axes(Atmp_FC, 1), weights_FC, N_sample2_FC)
    sd2 = Atmp_FC[indx2, :]
    df2 = DataFrame(log_value_added=sd2[:,3], exportstatus=sd2[:,1],
                    export_share=sd2[:,2], future_export_status=sd2[:,4])
    mdl2 = lm(@formula(future_export_status ~ 1 + log_value_added + exportstatus), df2)
    m2FC_coef[ii,:] = coef(mdl2); m2FC_se[ii,:] = stderror(mdl2)
    n2 = nobs(mdl2); m2FC_r2[ii] = 1 - (1 - r2(mdl2))*(n2-1)/(n2-4)
end

println("[tab:Table_regs] GLM regressions done (baseline BL and FC).")

# --- Save Tables/Table_regs.tex ---
let
    # Average across 100 simulations
    c1_va_bl   = mean(m1BL_coef[:,2]);  s1_va_bl   = mean(m1BL_se[:,2])
    c1_r2_bl   = mean(m1BL_r2)
    c2_va_bl   = mean(m2BL_coef[:,2]);  s2_va_bl   = mean(m2BL_se[:,2])
    c2_ex_bl   = mean(m2BL_coef[:,3]);  s2_ex_bl   = mean(m2BL_se[:,3])
    c2_r2_bl   = mean(m2BL_r2)

    c1_va_fc   = mean(m1FC_coef[:,2]);  s1_va_fc   = mean(m1FC_se[:,2])
    c1_r2_fc   = mean(m1FC_r2)
    c2_va_fc   = mean(m2FC_coef[:,2]);  s2_va_fc   = mean(m2FC_se[:,2])
    c2_ex_fc   = mean(m2FC_coef[:,3]);  s2_ex_fc   = mean(m2FC_se[:,3])
    c2_r2_fc   = mean(m2FC_r2)

    f2(x) = @sprintf("%.2f", round(x, digits=2, RoundNearestTiesAway))
    f3(x) = @sprintf("%.3f", round(x, digits=3, RoundNearestTiesAway))
    f4(x) = @sprintf("%.4f", round(x, digits=4, RoundNearestTiesAway))

    tex = """\\begin{table}[!ht]
\\centering
\\caption{Exporting dynamics in the model vs. data}
\\label{table:Table_regs}
\\scalebox{0.80}{
\\begin{tabular}{lcccccccc}
\\\\[-1.8ex]\\hline
\\hline \\\\[-1.8ex]
 & \\multicolumn{2}{c}{Entry cost} & \\multicolumn{2}{c}{Per-period cost} & \\multicolumn{4}{c}{Data} \\\\
\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-9}
 & (1) & (2) & (3) & (4) & (1) & (2) & (3) & (4) \\\\
\\hline \\\\[-1.8ex]
Log value added & $(f2(c1_va_bl)) & $(f2(c2_va_bl)) & $(f2(c1_va_fc)) & $(f2(c2_va_fc)) & $(f2(hfd["tr_va1"])) & $(f2(hfd["tr_va2"])) & $(f2(hfd["tr_va3"])) & $(f2(hfd["tr_va4"])) \\\\
 & ($(f3(s1_va_bl))) & ($(f3(s2_va_bl))) & ($(f3(s1_va_fc))) & ($(f3(s2_va_fc))) & ($(f4(hfd["tr_se1"]))) & ($(f4(hfd["tr_se2"]))) & ($(f4(hfd["tr_se3"]))) & ($(f4(hfd["tr_se4"]))) \\\\
Previous export status & - & $(f2(c2_ex_bl)) & - & $(f2(c2_ex_fc)) & - & $(f2(hfd["tr_ld2"])) & $(f2(hfd["tr_ld3"])) & $(f2(hfd["tr_ld4"])) \\\\
 & & ($(f3(s2_ex_bl))) & & ($(f3(s2_ex_fc))) & & ($(f4(hfd["tr_ld_se2"]))) & ($(f4(hfd["tr_ld_se3"]))) & ($(f4(hfd["tr_ld_se4"]))) \\\\
Previous export share & - & - & - & - & - & - & $(f2(hfd["tr_ls3"])) & $(f2(hfd["tr_ls4"])) \\\\
 & & & & & & & ($(f4(hfd["tr_ls_se3"]))) & ($(f4(hfd["tr_ls_se4"]))) \\\\
\\hline \\\\[-1.8ex]
N & $(N_sample_BL) & $(N_sample2_BL) & $(N_sample_FC) & $(N_sample2_FC) & $(round(Int, hfd["tr_N1"])) & $(round(Int, hfd["tr_N2"])) & $(round(Int, hfd["tr_N3"])) & $(round(Int, hfd["tr_N4"])) \\\\
Adj.\\ R\\textsuperscript{2} & $(f2(c1_r2_bl)) & $(f2(c2_r2_bl)) & $(f2(c1_r2_fc)) & $(f2(c2_r2_fc)) & $(f2(hfd["tr_r2_1"])) & $(f2(hfd["tr_r2_2"])) & $(f2(hfd["tr_r2_3"])) & $(f2(hfd["tr_r2_4"])) \\\\
\\hline
\\hline
\\end{tabular}
}
\\parbox{1.0\\textwidth}{\\caption*{\\scriptsize \\textit{Note:} The coefficient estimates for the linear probability model on simulated data from the quantitative model (with entry cost or with per-period cost) or the Hungarian firm data. Model \\textbf{1} includes only log value-added as a regressor, whereas model \\textbf{2} includes previous export status as well. Because the export share is constant in the model, regressions \\textbf{3} \\& \\textbf{4} cannot be replicated, but they do not change the coefficient estimates of the first two coefficients. The treatment of the regression for the Hungarian data matches that of \\cite{annurev:/content/journals/10.1146/annurev-economics-090919-025159}, with industry dummies for models \\textbf{1}-\\textbf{3}, not for model \\textbf{4}, and year fixed effects for all of them. Standard errors, clustered at the industry level, are reported in parentheses. }}
\\end{table}
"""
    open("../Tables/Table_regs.tex", "w") do io
        write(io, tex)
    end
    println("  Saved → ../Tables/Table_regs.tex")
end

# =============================================================================
# [tab:cap_market_basis]  Capital-market integration starting from trade-only SS
# main.tex: \label{tab:cap_market_basis}  (Appendix)
# Source: plot.jl lines 212–307
# Prerequisites: Steps 1-2 (welfare_change_oCM_from_otrade_trans needed).
# Computes table_9b_results (20×2), table_5b_results (14×2), table_6b_results (4×3)
# Outputs: ../Tables/Table_cap_market_basis.tex
# =============================================================================

# table_9b_results: aggregates with trade-only SS as baseline (cols: trade-only | trade+CM)
table_9b_results = zeros(20, 2)
table_9b_results[1,:]  = [TFP_closed_CM_open_trade[1]/TFP_closed_CM_open_trade[1],
    TFP_open_CM_open_trade[1]/TFP_closed_CM_open_trade[1]]
table_9b_results[2,:]  = [sd_MRPK_closed_CM_open_trade[1], sd_MRPK_open_CM_open_trade[1]]
table_9b_results[3,:]  = [GDP_closed_CM_open_trade[1]/GDP_closed_CM_open_trade[1],
    GDP_open_CM_open_trade[1]/GDP_closed_CM_open_trade[1]]
table_9b_results[4,:]  = [GDP_closed_CM_open_trade[2]/GDP_closed_CM_open_trade[2],
    GDP_open_CM_open_trade[2]/GDP_closed_CM_open_trade[2]]
table_9b_results[5,:]  = [total_consumption_closed_CM_open_trade[1]/total_consumption_closed_CM_open_trade[1],
    total_consumption_open_CM_open_trade[1]/total_consumption_closed_CM_open_trade[1]]
table_9b_results[6,:]  = [capital_demand_closed_CM_open_trade[1]/capital_demand_closed_CM_open_trade[1],
    capital_demand_open_CM_open_trade[1]/capital_demand_closed_CM_open_trade[1]]
table_9b_results[8,:]  = [0, welfare_change_oCM_from_otrade_trans]
table_9b_results[9,:]  = [p90_wealth_closed_CM_open_trade[1], p90_wealth_open_CM_open_trade[1]]
table_9b_results[10,:] = [prices_closed_CM_open_trade[1]/prices_closed_CM_open_trade[5]/
    (prices_closed_CM_open_trade[1]/prices_closed_CM_open_trade[5]),
    prices_open_CM_open_trade[1]/prices_open_CM_open_trade[5]/
    (prices_closed_CM_open_trade[1]/prices_closed_CM_open_trade[5])]
table_9b_results[12,:] = [import_share_closed_CM_open_trade[1], import_share_open_CM_open_trade[1]]
table_9b_results[13,:] = [export_value_closed_CM_open_trade[1]/export_value_closed_CM_open_trade[1],
    export_value_open_CM_open_trade[1]/export_value_closed_CM_open_trade[1]]
table_9b_results[14,:] = [(domestic_pop_closed_CM_open_trade[1]+exporter_pop_closed_CM_open_trade[1]),
    (domestic_pop_open_CM_open_trade[1]+exporter_pop_open_CM_open_trade[1])]
table_9b_results[15,:] = [exporter_pop_closed_CM_open_trade[1]/
    (domestic_pop_closed_CM_open_trade[1]+exporter_pop_closed_CM_open_trade[1]),
    exporter_pop_open_CM_open_trade[1]/
    (domestic_pop_open_CM_open_trade[1]+exporter_pop_open_CM_open_trade[1])]
table_9b_results[16,:] = [prices_closed_CM_open_trade[5], prices_open_CM_open_trade[5]]
table_9b_results[17,:] = [total_credit_closed_CM_open_trade, total_credit_open_CM_open_trade]
table_9b_results[18,:] = [foreign_credit_share_closed_CM_open_trade, foreign_credit_share_open_CM_open_trade]

# table_5b_results: firm-level, trade-only → trade+CM
table_5b_results = zeros(14, 2)
table_5b_results[1,:]  = [sd_MRPK_d_closed_CM_open_trade[1], sd_MRPK_d_open_CM_open_trade[1]]
table_5b_results[2,:]  = [sd_MRPK_x_closed_CM_open_trade[1], sd_MRPK_x_open_CM_open_trade[1]]
table_5b_results[3,:]  = [Misallocation_within_d_closed_CM_open_trade[1], Misallocation_within_d_open_CM_open_trade[1]]
table_5b_results[4,:]  = [Misallocation_within_x_closed_CM_open_trade[1], Misallocation_within_x_open_CM_open_trade[1]]
table_5b_results[5,:]  = [1, domestic_pop_open_CM_open_trade[1]/domestic_pop_closed_CM_open_trade[1]]
table_5b_results[6,:]  = [1, exporter_pop_open_CM_open_trade[1]/exporter_pop_closed_CM_open_trade[1]]
table_5b_results[7,1]  = K_x_closed_CM_open_trade[1]/(K_x_closed_CM_open_trade[1]+K_d_closed_CM_open_trade[1])
table_5b_results[7,2]  = K_x_open_CM_open_trade[1]/(K_x_open_CM_open_trade[1]+K_d_open_CM_open_trade[1])
table_5b_results[8,1]  = L_x_closed_CM_open_trade[1]/(L_x_closed_CM_open_trade[1]+L_d_closed_CM_open_trade[1])
table_5b_results[8,2]  = L_x_open_CM_open_trade[1]/(L_x_open_CM_open_trade[1]+L_d_open_CM_open_trade[1])
table_5b_results[9,:]  = [exporter_pop_closed_CM_open_trade[1]/entry_share_to_exporter_closed_CM_open_trade[1],
    exporter_pop_open_CM_open_trade[1]/entry_share_to_exporter_open_CM_open_trade[1]]
table_5b_results[10,:] = [1,
    sum(k_d_fine_open_CM_open_trade[:,1].*current_domestic_open_CM_open_trade)/domestic_pop_open_CM_open_trade[1] *
    domestic_pop_closed_CM_open_trade[1]/sum(k_d_fine_closed_CM_open_trade[:,1].*current_domestic_closed_CM_open_trade)]
table_5b_results[11,:] = [1,
    sum(k_x_fine_open_CM_open_trade[:,1].*current_exporter_open_CM_open_trade)/exporter_pop_open_CM_open_trade[1] *
    exporter_pop_closed_CM_open_trade[1]/sum(k_x_fine_closed_CM_open_trade[:,1].*current_exporter_closed_CM_open_trade)]
table_5b_results[12,:] = [1,
    sum(k_x_fine_open_CM_open_trade[:,1].*current_exporter_open_CM_open_trade +
        k_d_fine_open_CM_open_trade[:,1].*current_domestic_open_CM_open_trade) /
    (domestic_pop_open_CM_open_trade[1]+exporter_pop_open_CM_open_trade[1]) *
    (domestic_pop_closed_CM_open_trade[1]+exporter_pop_closed_CM_open_trade[1]) /
    sum(k_x_fine_closed_CM_open_trade[:,1].*current_exporter_closed_CM_open_trade +
        k_d_fine_closed_CM_open_trade[:,1].*current_domestic_closed_CM_open_trade)]
table_5b_results[13,:] = [mean_leverage_closed_CM_open_trade[1], mean_leverage_open_CM_open_trade[1]]
table_5b_results[14,:] = [mean_leverage_x_closed_CM_open_trade[1], mean_leverage_x_open_CM_open_trade[1]]

# table_6b_results: exporter distribution (trade-only vs trade+CM; col 1 unused)
y_start_6b = 5; x_start_6b = 18
table_6b_results = zeros(4, 3)
table_6b_results[4,2] = sum(current_exporter_closed_CM_open_trade_mat[(x_start_6b+1):end,(y_start_6b+1):end])/sum(current_exporter_closed_CM_open_trade)
table_6b_results[3,2] = sum(current_exporter_closed_CM_open_trade_mat[1:x_start_6b,(y_start_6b+1):end])/sum(current_exporter_closed_CM_open_trade)
table_6b_results[2,2] = sum(current_exporter_closed_CM_open_trade_mat[(x_start_6b+1):end,1:y_start_6b])/sum(current_exporter_closed_CM_open_trade)
table_6b_results[1,2] = sum(current_exporter_closed_CM_open_trade_mat[1:x_start_6b,1:y_start_6b])/sum(current_exporter_closed_CM_open_trade)
table_6b_results[4,3] = sum(current_exporter_open_CM_open_trade_mat[(x_start_6b+1):end,(y_start_6b+1):end])/sum(current_exporter_open_CM_open_trade)
table_6b_results[3,3] = sum(current_exporter_open_CM_open_trade_mat[1:x_start_6b,(y_start_6b+1):end])/sum(current_exporter_open_CM_open_trade)
table_6b_results[2,3] = sum(current_exporter_open_CM_open_trade_mat[(x_start_6b+1):end,1:y_start_6b])/sum(current_exporter_open_CM_open_trade)
table_6b_results[1,3] = sum(current_exporter_open_CM_open_trade_mat[1:x_start_6b,1:y_start_6b])/sum(current_exporter_open_CM_open_trade)

println("[tab:cap_market_basis] table_9b, table_5b, table_6b computed.")

# --- Save Tables/Table_cap_market_basis.tex ---
let
    r9 = table_9b_results; r5 = table_5b_results; r6 = table_6b_results
    f2(x) = @sprintf("%.2f", round(x, digits=2, RoundNearestTiesAway))
    f1(x) = @sprintf("%.1f", round(x, digits=1, RoundNearestTiesAway))

    tex = """\\begin{table}[!ht]
\\centering
\\caption{Capital market integration starting from trade-open economy}
\\label{table:cap_market_basis}
\\scalebox{0.6}{
\\begin{tabular}{lcc@{\\hspace{3em}}lcc}
\\\\[-1.8ex]\\hline
\\hline \\\\[-1.8ex]
Capital Market & Closed & Integrated & & Closed & Integrated \\\\
\\hline \\\\[-1.8ex]
\\textbf{Productivity} & & & \\textbf{Extensive margin} & & \\\\
TFP & $(fpct(r9[1,1])) & $(fpct(r9[1,2])) & Non-exporting firms & $(fpct(r5[5,1])) & $(fpct(r5[5,2])) \\\\
Standard deviation of ARPK & $(f2(r9[2,1])) & $(f2(r9[2,2])) & Exporting firms & $(fpct(r5[6,1])) & $(fpct(r5[6,2])) \\\\[1ex]
\\textbf{Aggregates} & & & \\textbf{Intensive margin} & & \\\\
GDP & $(fpct(r9[3,1])) & $(fpct(r9[3,2])) & \$ \\%\$ of capital used by exporters & $(fpct(r5[7,1])) & $(fpct(r5[7,2])) \\\\
GDP* & $(fpct(r9[4,1])) & $(fpct(r9[4,2])) & \$ \\%\$ of labor used by exporters & $(fpct(r5[8,1])) & $(fpct(r5[8,2])) \\\\
Consumption & $(fpct(r9[5,1])) & $(f1(100*r9[5,2])) & Avg. duration (years) of export status & $(f1(r5[9,1])) & $(f1(r5[9,2])) \\\\
Capital & $(fpct(r9[6,1])) & $(fpct(r9[6,2])) & Average capital size of non-exporters & $(fpct(r5[10,1])) & $(fpct(r5[10,2])) \\\\[1ex]
\\textbf{Welfare and Inequality} & & & Average capital size of exporters & $(fpct(r5[11,1])) & $(fpct(r5[11,2])) \\\\
Consumption equivalent welfare & 0 & $(f1(100*r9[8,2])) & Average capital size & $(fpct(r5[12,1])) & $(fpct(r5[12,2])) \\\\
Top \$10\\%\$ wealth share & $(fpct(r9[9,1])) & $(fpct(r9[9,2])) & Mean leverage, all firms & $(fpct(r5[13,1])) & $(fpct(r5[13,2])) \\\\[1ex]
\\textbf{Factor prices} & & & Mean leverage, exporters & $(fpct(r5[14,1])) & $(fpct(r5[14,2])) \\\\[1ex]
Real wage & $(fpct(r9[10,1])) & $(fpct(r9[10,2])) & \\textbf{Standard deviation of ARPK} & & \\\\
Interest rate premium \$r - r^*\$ & $(fpct(prices_closed_CM_open_trade[6] - prices_closed_CM_open_trade[7])) & 0 & Non-exporting & $(f2(r5[1,1])) & $(f2(r5[1,2])) \\\\[1ex]
\\textbf{Trade} & & & Exporter & $(f2(r5[2,1])) & $(f2(r5[2,2])) \\\\[1ex]
\$\\frac{\\text{Import}}{\\text{GDP}}\$ & $(fpct(r9[12,1])) & $(fpct(r9[12,2])) & \\textbf{Within type productivity loss} & & \\\\
Export & $(fpct(r9[13,1])) & $(fpct(r9[13,2])) & Non-exporting & $(f1(100*r5[3,1])) & $(f1(100*r5[3,2])) \\\\
Entrepreneurship rate & $(fpct(r9[14,1])) & $(fpct(r9[14,2])) & Exporter & $(f1(100*r5[4,1])) & $(f1(100*r5[4,2])) \\\\[1ex]
Share of exporters & $(fpct(r9[15,1])) & $(fpct(r9[15,2])) & \\textbf{Distribution of exporters (\\%)} & & \\\\
CPI & $(fpct(r9[16,1])) & $(fpct(r9[16,2])) & Low wealth and low productivity & $(fpct(r6[1,2])) & $(fpct(r6[1,3])) \\\\
\$\\frac{\\text{Credit}}{\\text{GDP}}\$ & $(fpct(r9[17,1])) & $(fpct(r9[17,2])) & Low wealth and high productivity & $(fpct(r6[2,2])) & $(fpct(r6[2,3])) \\\\
\$\\frac{\\text{Foreign Credit}}{\\text{Credit}}\$ & 0 & $(fpct(r9[18,2])) & High wealth and low productivity & $(fpct(r6[3,2])) & $(fpct(r6[3,3])) \\\\
 & & & High wealth and high productivity & $(fpct(r6[4,2])) & $(fpct(r6[4,3])) \\\\
\\hline
\\hline
\\end{tabular}
}
\\parbox{1.5\\textwidth}{\\caption*{\\scriptsize \\textit{Note:} The quantities without a data counterpart are normalized to the closed capital market steady state --- the integrated capital markets steady state is compared to the data in Table \\ref{table:Table4}. Asterisks indicate Core variables. Standard deviation of ARPK is calculated within firm status \$\\{d, ex\\}\$. Productivity loss is defined as \$\\frac{TFP_i^{e} - TFP_{i}}{TFP_i^{e}}\$ with \$TFP_i^{e}\$ defined as the productivity of firm status \$\\{d, ex\\}\$ that would occur if all firms would be able to obtain their unconstrained input choice. The productivity threshold is the average entrepreneurial productivity; the wealth threshold is the median wealth of all entrepreneurs in the closed capital market steady state. Numbers indicate percentages.}}
\\end{table}
"""
    open("../Tables/Table_cap_market_basis.tex", "w") do io
        write(io, tex)
    end
    println("  Saved → ../Tables/Table_cap_market_basis.tex")
end

# =============================================================================
# [tab:calibration_fixed_cost]  FC model calibrated parameters vs moments
# main.tex: \label{tab:calibration_fixed_cost}  (Appendix)
# Source: plot_FC.jl lines 33–47
# Prerequisites: Steps 1 and 3 (both BL and FC SS).
# table_1_results_FC (8×4): col1=param, col2=data target, col3=FC model, col4=BL model
# Output: ../Tables/Table_calibration_FC.tex
# =============================================================================

# FC welfare scalars (needed first)
V_saved_initial_FC        = reshape(V_saved_fine_initial_FC[:,:,1], size(V_saved_fine_initial)[1]*3)
V_saved_closed_CM_FC      = reshape(V_saved_fine_closed_CM_open_trade_FC[:,:,1], size(V_saved_fine_initial)[1]*3)
V_saved_open_CM_FC        = reshape(V_saved_fine_open_CM_open_trade_FC[:,:,1], size(V_saved_fine_initial)[1]*3)

welfare_init_closed_CM_stst_FC =
    sum(current_distr_store_initial_FC[:,1] .*
        (exp.((V_saved_closed_CM_FC - V_saved_initial_FC) * (1.0 - Baseline_parameter.β[1])) .- 1))
welfare_init_open_CM_stst_FC =
    sum(current_distr_store_initial_FC[:,1] .*
        (exp.((V_saved_open_CM_FC - V_saved_initial_FC) * (1.0 - Baseline_parameter.β[1])) .- 1))

V_trans_closed_CM_FC =
    reshape(V_saved_store_closed_CM_open_trade_FC[:,:,1,2], size(V_saved_fine_initial)[1]*3)
V_trans_open_CM_FC =
    reshape(V_saved_store_open_CM_open_trade_FC[:,:,1,2], size(V_saved_fine_initial)[1]*3)
V_trans_open_CMdelayed_FC =
    reshape(V_saved_store_open_CMdelayed_open_trade_FC[:,:,1,2], size(V_saved_fine_initial)[1]*3)

welfare_change_clCM_otrade_trans_FC =
    sum(current_distr_store_initial_FC[:,1] .*
        (exp.((V_trans_closed_CM_FC - V_saved_initial_FC) * (1.0 - Baseline_parameter.β[1])) .- 1))
welfare_change_oCM_otrade_trans_FC =
    sum(current_distr_store_initial_FC[:,1] .*
        (exp.((V_trans_open_CM_FC - V_saved_initial_FC) * (1.0 - Baseline_parameter.β[1])) .- 1))
welfare_change_oCM_delayed_trans_FC =
    sum(current_distr_store_initial_FC[:,1] .*
        (exp.((V_trans_open_CMdelayed_FC - V_saved_initial_FC) * (1.0 - Baseline_parameter.β[1])) .- 1))

# FC credit statistics (recomputed cleanly; already computed above but alias here for clarity)
total_credit_open_CM_FC =
    -(domestic_firm_debt_open_CM_open_trade_FC[1] + exporter_firm_debt_open_CM_open_trade_FC[1]) /
    nomGDP_open_CM_open_trade_FC[1]
domestic_credit_open_CM_FC =
    (worker_bond_holding_open_CM_open_trade_FC[1] + domestic_bond_holding_open_CM_open_trade_FC[1] +
     exporter_bond_holding_open_CM_open_trade_FC[1]) / nomGDP_open_CM_open_trade_FC[1]
foreign_credit_open_CM_FC     = total_credit_open_CM_FC - domestic_credit_open_CM_FC
foreign_credit_share_open_CM_FC = foreign_credit_open_CM_FC / total_credit_open_CM_FC

table_1_results_FC = zeros(8, 4)
table_1_results_FC[1,:] = [open_CM_open_trade_parameter_FC.θ[1], 0.62,
    total_credit_open_CM_FC, total_credit_open_CM_open_trade]
table_1_results_FC[2,:] = [open_CM_open_trade_parameter_FC.β[2], 0.05,
    prices_open_CM_open_trade_FC[6], prices_open_CM_open_trade[6]]
table_1_results_FC[3,:] = [open_CM_open_trade_parameter_FC.β[1], 0.53,
    foreign_credit_share_open_CM_FC, foreign_credit_share_open_CM_open_trade]
table_1_results_FC[4,:] = [open_CM_open_trade_parameter_FC.τ[1], 0.21,
    import_share_initial_FC[1], import_share_initial[1]]
table_1_results_FC[5,:] = [open_CM_open_trade_parameter_FC.τ[1], 0.42,
    import_share_open_CM_open_trade_FC[1], import_share_open_CM_open_trade[1]]
table_1_results_FC[6,:] = [open_CM_open_trade_parameter_FC.Fₓ[1], 0.27,
    entry_share_to_exporter_open_CM_open_trade_FC[1]/exporter_pop_open_CM_open_trade_FC[1],
    entry_share_to_exporter_open_CM_open_trade[1]/exporter_pop_open_CM_open_trade[1]]
table_1_results_FC[7,:] = [open_CM_open_trade_parameter_FC.σₛ[1], 0.82,
    sd_growth_rev_open_CM_open_trade_FC[1], sd_growth_rev_open_CM_open_trade[1]]
table_1_results_FC[8,:] = [open_CM_open_trade_parameter_FC.ρ[1], 0.43,
    autocorr_rev_open_CM_open_trade_FC[1], autocorr_rev_open_CM_open_trade[1]]

println("[tab:calibration_fixed_cost] table_1_results_FC (8×4) computed.")

# --- Save Tables/Table_calibration_FC.tex ---
let
    r = table_1_results_FC
    f2(x) = @sprintf("%.2f", round(x, digits=2, RoundNearestTiesAway))

    tex = """\\begin{table}[!ht]
\\centering
\\caption{The calibration without entry cost but with fixed cost compared to baseline}
\\label{tab:calibration_fixed_cost}
\\scalebox{0.55}{
\\begin{tabular}{lccc}
\\\\[-1.8ex]\\hline
\\hline \\\\[-1.8ex]
Target & Data & \$F_{ex} = 0, FM_{ex} = 1.5\$ & \$F_{ex} = 4.5, FM_{ex} = 0\$ \\\\
\\hline \\\\[-1.8ex]
\\textbf{Financial Development} & & & \\\\[1ex]
\$\\text{Total Credit to nonfinancials}\$, \$\\%\\text{GDP}\$ & 62 & $(fpct(r[1,3])) & $(fpct(r[1,4])) \\\\
Bank lending rate in Germany \$r^*\$ & 5 & $(fpct(r[2,3])) & $(fpct(r[2,4])) \\\\
\$\\text{Foreign Credit to nonfinancials}\$, \$\\%\\text{Total credit}\$ & 53 & $(fpct(r[3,3])) & $(fpct(r[3,4])) \\\\[1ex]
\\textbf{Trade} & & & \\\\[1ex]
Initial \$\\frac{\\text{Import}}{\\text{GDP}}\$ & 21 & $(fpct(r[4,3])) & $(fpct(r[4,4])) \\\\
Final \$\\frac{\\text{Import}}{\\text{GDP}}\$ & 42 & $(fpct(r[5,3])) & $(fpct(r[5,4])) \\\\[1ex]
\\textbf{Firm dynamics} & & & \\\\[1ex]
Entry rate to exports & 20 & $(fpct(r[6,3])) & $(fpct(r[6,4])) \\\\
s.d. value added & 0.82 & $(f2(r[7,3])) & $(f2(r[7,4])) \\\\
Auto-correlation of value-added & 0.43 & $(f2(r[8,3])) & $(f2(r[8,4])) \\\\
\\hline
\\hline
\\end{tabular}
}
\\end{table}
"""
    open("../Tables/Table_calibration_FC.tex", "w") do io
        write(io, tex)
    end
    println("  Saved → ../Tables/Table_calibration_FC.tex")
end

# =============================================================================
# [tab:nontargeted_fixed_cost]  FC model non-targeted moments
# main.tex: \label{tab:nontargeted_fixed_cost}  (Appendix)
# Source: plot_FC.jl lines 49–64
# Prerequisites: Steps 1 and 3.
# table_3_results_FC (14×3): col1=data, col2=FC model, col3=BL model
# Output: ../Tables/Table_nontargeted_FC.tex
# =============================================================================

table_3_results_FC = zeros(14, 3)
table_3_results_FC[1,:]  = [1.36, sd_MRPK_open_CM_open_trade_FC[1], sd_MRPK_open_CM_open_trade[1]]
table_3_results_FC[2,:]  = [0.61, sd_growth_k_open_CM_open_trade_FC[1], sd_growth_k_open_CM_open_trade[1]]
table_3_results_FC[3,:]  = [0.29,
    exporter_pop_open_CM_open_trade_FC[1]/(domestic_pop_open_CM_open_trade_FC[1]+exporter_pop_open_CM_open_trade_FC[1]),
    exporter_pop_open_CM_open_trade[1]/(domestic_pop_open_CM_open_trade[1]+exporter_pop_open_CM_open_trade[1])]
table_3_results_FC[4,:]  = [0.46, mean_leverage_open_CM_open_trade_FC[1], mean_leverage_open_CM_open_trade[1]]
table_3_results_FC[5,:]  = [0.51, mean_leverage_x_open_CM_open_trade_FC[1], mean_leverage_x_open_CM_open_trade[1]]
table_3_results_FC[6,:]  = [0.57,
    exporter_firm_debt_open_CM_open_trade_FC[1]/
        (domestic_firm_debt_open_CM_open_trade_FC[1]+exporter_firm_debt_open_CM_open_trade_FC[1]),
    exporter_firm_debt_open_CM_open_trade[1]/
        (domestic_firm_debt_open_CM_open_trade[1]+exporter_firm_debt_open_CM_open_trade[1])]
table_3_results_FC[7,:]  = [0.64,
    K_x_open_CM_open_trade_FC[1]/(K_x_open_CM_open_trade_FC[1]+K_d_open_CM_open_trade_FC[1]),
    K_x_open_CM_open_trade[1]/(K_x_open_CM_open_trade[1]+K_d_open_CM_open_trade[1])]
table_3_results_FC[8,:]  = [0.55,
    L_x_open_CM_open_trade_FC[1]/(L_x_open_CM_open_trade_FC[1]+L_d_open_CM_open_trade_FC[1]),
    L_x_open_CM_open_trade[1]/(L_x_open_CM_open_trade[1]+L_d_open_CM_open_trade[1])]
table_3_results_FC[9,:]  = [0.34,
    prices_open_CM_open_trade_FC[3]/prices_open_CM_open_trade_FC[4]/
        open_CM_open_trade_parameter_FC.L[1]*open_CM_open_trade_parameter_FC.L[2],
    prices_open_CM_open_trade[3]/prices_open_CM_open_trade[4]/
        open_CM_open_trade_parameter.L[1]*open_CM_open_trade_parameter.L[2]]
table_3_results_FC[10,:] = [0.53, p90_wealth_open_CM_open_trade_FC[1], p90_wealth_open_CM_open_trade[1]]
table_3_results_FC[11,:] = [0.34, p90_income_open_CM_open_trade_FC[1], p90_income_open_CM_open_trade[1]]
table_3_results_FC[12,:] = [0.11, p99_income_open_CM_open_trade_FC[1], p99_income_open_CM_open_trade[1]]
table_3_results_FC[13,:] = [0.24, p90_income_initial_FC[1], p90_income_initial[1]]
table_3_results_FC[14,:] = [0.06, p99_income_initial_FC[1], p99_income_initial[1]]

println("[tab:nontargeted_fixed_cost] table_3_results_FC (14×3) computed.")

# --- Save Tables/Table_nontargeted_FC.tex ---
let
    r = table_3_results_FC
    f2(x) = @sprintf("%.2f", round(x, digits=2, RoundNearestTiesAway))

    tex = """\\begin{table}[!ht]
\\centering
\\caption{Moments of without entry cost but with fixed cost compared to baseline}
\\label{tab:nontargeted_fixed_cost}
\\scalebox{0.7}{
\\begin{tabular}{lcccc}
\\hline
\\hline \\\\[-1.8ex]
Description & Data & \$F_{ex} = 0, FM_{ex} = 1.5\$ & \$F_{ex} = 4.5, FM_{ex} = 0\$ & Source \\& Year \\\\
\\hline \\\\[-1.8ex]
\\textbf{Production} & & & & \\\\
Standard deviation of ARPK & $(f2(r[1,1])) & $(f2(r[1,2])) & $(f2(r[1,3])) & Firm level, Hungary \\\\
Standard deviation of log capital growth & $(f2(r[2,1])) & $(f2(r[2,2])) & $(f2(r[2,3])) & Firm level, Hungary \\\\
\\textbf{Exporters} & & & & \\\\
Fraction of firms that export & $(fpct(r[3,1])) & $(fpct(r[3,2])) & $(fpct(r[3,3])) & Table 1 \\\\
Mean leverage, all firms & $(fpct(r[4,1])) & $(fpct(r[4,2])) & $(fpct(r[4,3])) & Table 1 \\\\
Mean leverage, exporters & $(fpct(r[5,1])) & $(fpct(r[5,2])) & $(fpct(r[5,3])) & Table 1 \\\\
Fraction of total debt credited to exporters & $(fpct(r[6,1])) & $(fpct(r[6,2])) & $(fpct(r[6,3])) & Firm level, Hungary \\\\
Fraction of total capital used by exporters & $(fpct(r[7,1])) & $(fpct(r[7,2])) & $(fpct(r[7,3])) & Firm level, Hungary \\\\
Fraction of total employment used by exporters & $(fpct(r[8,1])) & $(fpct(r[8,2])) & $(fpct(r[8,3])) & Firm level, Hungary \\\\
\\textbf{Inequality} & & & & \\\\
GDP per capita, Hungary vs. Germany & $(fpct(r[9,1])) & $(fpct(r[9,2])) & $(fpct(r[9,3])) & WB, 2008 \\\\
Top \$10\\%\$ wealth share & $(fpct(r[10,1])) & $(fpct(r[10,2])) & $(fpct(r[10,3])) & HSO 2014 \\\\
Top \$10\\%\$ income share & $(fpct(r[11,1])) & $(fpct(r[11,2])) & $(fpct(r[11,3])) & WID 2008 \\\\
Top \$1\\%\$ income share & $(fpct(r[12,1])) & $(fpct(r[12,2])) & $(fpct(r[12,3])) & WID 2008 \\\\
Top \$10\\%\$ income share & $(fpct(r[13,1])) & $(fpct(r[13,2])) & $(fpct(r[13,3])) & WID 1991 \\\\
Top \$1\\%\$ income share & $(fpct(r[14,1])) & $(fpct(r[14,2])) & $(fpct(r[14,3])) & WID 1991 \\\\
\\hline
\\hline
\\end{tabular}
}
\\end{table}
"""
    open("../Tables/Table_nontargeted_FC.tex", "w") do io
        write(io, tex)
    end
    println("  Saved → ../Tables/Table_nontargeted_FC.tex")
end

# =============================================================================
# [tab:results_fixed_cost]  FC model trade-liberalization aggregates
# main.tex: \label{tab:results_fixed_cost}  (Appendix)
# Source: plot_FC.jl lines 86–125
# Prerequisites: Steps 1, 3, 4.
# table_4_results_FC (20×6): same structure as table_4_results
# Output: ../Tables/Table_results_FC.tex
# =============================================================================

current_exporter_closed_CM_FC = current_distr_store_closed_CM_open_trade_FC[
    convert(Int64, 2/3*size(current_distr_store_closed_CM_open_trade_FC)[1]+1):end, 1]
current_exporter_open_CM_FC   = current_distr_store_open_CM_open_trade_FC[
    convert(Int64, 2/3*size(current_distr_store_closed_CM_open_trade_FC)[1]+1):end, 1]
current_exporter_initial_FC   = current_distr_store_initial_FC[
    convert(Int64, 2/3*size(current_distr_store_closed_CM_open_trade_FC)[1]+1):end, 1]
current_domestic_open_CM_FC   = current_distr_store_open_CM_open_trade_FC[7201:14400, 1]
current_domestic_closed_CM_FC = current_distr_store_closed_CM_open_trade_FC[7201:14400, 1]
current_domestic_initial_FC   = current_distr_store_initial_FC[7201:14400, 1]

price_union_wide_initial_FC =
    prices_initial_FC[5] * prices_initial_FC[5] * GDP_initial_FC[1] /
    (prices_initial_FC[5]*GDP_initial_FC[1] + GDP_initial_FC[2]) +
    1.0 * GDP_initial_FC[2] / (prices_initial_FC[5]*GDP_initial_FC[1] + GDP_initial_FC[2])
price_union_wide_closed_CM_FC =
    prices_closed_CM_open_trade_FC[5] * prices_closed_CM_open_trade_FC[5] * GDP_closed_CM_open_trade_FC[1] /
    (prices_closed_CM_open_trade_FC[5]*GDP_closed_CM_open_trade_FC[1] + GDP_closed_CM_open_trade_FC[2]) +
    1.0 * GDP_closed_CM_open_trade_FC[2] /
    (prices_closed_CM_open_trade_FC[5]*GDP_closed_CM_open_trade_FC[1] + GDP_closed_CM_open_trade_FC[2])
price_union_wide_open_CM_FC =
    prices_open_CM_open_trade_FC[5] * prices_open_CM_open_trade_FC[5] * GDP_open_CM_open_trade_FC[1] /
    (prices_open_CM_open_trade_FC[5]*GDP_open_CM_open_trade_FC[1] + GDP_open_CM_open_trade_FC[2]) +
    1.0 * GDP_open_CM_open_trade_FC[2] /
    (prices_open_CM_open_trade_FC[5]*GDP_open_CM_open_trade_FC[1] + GDP_open_CM_open_trade_FC[2])

table_4_results_FC = zeros(20, 6)
table_4_results_FC[1,:]  = [TFP_initial_FC[1]/TFP_initial_FC[1], TFP_closed_CM_open_trade_FC[1]/TFP_initial_FC[1],
    TFP_open_CM_open_trade_FC[1]/TFP_initial_FC[1], 100, 116, 132]
table_4_results_FC[2,:]  = [sd_MRPK_initial_FC[1], sd_MRPK_closed_CM_open_trade_FC[1],
    sd_MRPK_open_CM_open_trade_FC[1], 0, 1.27, 1.51]
table_4_results_FC[3,:]  = [GDP_initial_FC[1]/GDP_initial_FC[1], GDP_closed_CM_open_trade_FC[1]/GDP_initial_FC[1],
    GDP_open_CM_open_trade_FC[1]/GDP_initial_FC[1], 100, 158, 249]
table_4_results_FC[4,:]  = [
    prices_initial_FC[5]*GDP_initial_FC[1]/(prices_initial_FC[5]*GDP_initial_FC[1]+GDP_initial_FC[2]) *
        price_union_wide_initial_FC / open_CM_open_trade_parameter_FC.L[1] *
        (open_CM_open_trade_parameter_FC.L[1]+open_CM_open_trade_parameter_FC.L[2]),
    prices_closed_CM_open_trade_FC[5]*GDP_closed_CM_open_trade_FC[1]/
        (prices_closed_CM_open_trade_FC[5]*GDP_closed_CM_open_trade_FC[1]+GDP_closed_CM_open_trade_FC[2]) *
        price_union_wide_closed_CM_FC / open_CM_open_trade_parameter_FC.L[1] *
        (open_CM_open_trade_parameter_FC.L[1]+open_CM_open_trade_parameter_FC.L[2]),
    prices_open_CM_open_trade_FC[5]*GDP_open_CM_open_trade_FC[1]/
        (prices_open_CM_open_trade_FC[5]*GDP_open_CM_open_trade_FC[1]+GDP_open_CM_open_trade_FC[2]) *
        price_union_wide_open_CM_FC / open_CM_open_trade_parameter_FC.L[1] *
        (open_CM_open_trade_parameter_FC.L[1]+open_CM_open_trade_parameter_FC.L[2]),
    0.54, 0.57, 0.64]
table_4_results_FC[5,:]  = [total_consumption_initial_FC[1]/total_consumption_initial_FC[1],
    total_consumption_closed_CM_open_trade_FC[1]/total_consumption_initial_FC[1],
    total_consumption_open_CM_open_trade_FC[1]/total_consumption_initial_FC[1], 100, 116, 133]
table_4_results_FC[6,:]  = [capital_demand_initial_FC[1]/capital_demand_initial_FC[1],
    capital_demand_closed_CM_open_trade_FC[1]/capital_demand_initial_FC[1],
    capital_demand_open_CM_open_trade_FC[1]/capital_demand_initial_FC[1], 100, 124, 154]
table_4_results_FC[8,:]  = [0, welfare_change_clCM_otrade_trans_FC, welfare_change_oCM_otrade_trans_FC, 0, 0, 0]
table_4_results_FC[9,:]  = [p90_income_initial_FC[1], p90_income_closed_CM_open_trade_FC[1],
    p90_income_open_CM_open_trade_FC[1], 0.24, 0.3, 0.34]
table_4_results_FC[10,:] = [prices_initial_FC[1]/prices_initial_FC[5]/(prices_initial_FC[1]/prices_initial_FC[5]),
    prices_closed_CM_open_trade_FC[1]/prices_closed_CM_open_trade_FC[5]/(prices_initial_FC[1]/prices_initial_FC[5]),
    prices_open_CM_open_trade_FC[1]/prices_open_CM_open_trade_FC[5]/(prices_initial_FC[1]/prices_initial_FC[5]),
    0, 100, 144]
table_4_results_FC[11,:] = [prices_initial_FC[6]-prices_initial_FC[7],
    prices_closed_CM_open_trade_FC[6]-prices_closed_CM_open_trade_FC[7],
    prices_open_CM_open_trade_FC[6]-prices_open_CM_open_trade_FC[6], 0.21, 0.05, 0.03]
table_4_results_FC[12,:] = [import_share_initial_FC[1], import_share_closed_CM_open_trade_FC[1],
    import_share_open_CM_open_trade_FC[1], 0.21, 0.42, 0.42]
table_4_results_FC[13,:] = [0, 0, NFA_open_CM_open_trade_FC[1]/nomGDP_open_CM_open_trade_FC[1],
    0.01, -0.05, -0.06]
table_4_results_FC[14,:] = [(domestic_pop_initial_FC[1]+exporter_pop_initial_FC[1]),
    (domestic_pop_closed_CM_open_trade_FC[1]+exporter_pop_closed_CM_open_trade_FC[1]),
    (domestic_pop_open_CM_open_trade_FC[1]+exporter_pop_open_CM_open_trade_FC[1]), 0, 0, 0.07]
table_4_results_FC[15,:] = [exporter_pop_initial_FC[1]/(domestic_pop_initial_FC[1]+exporter_pop_initial_FC[1]),
    exporter_pop_closed_CM_open_trade_FC[1]/(domestic_pop_closed_CM_open_trade_FC[1]+exporter_pop_closed_CM_open_trade_FC[1]),
    exporter_pop_open_CM_open_trade_FC[1]/(domestic_pop_open_CM_open_trade_FC[1]+exporter_pop_open_CM_open_trade_FC[1]),
    0, 0.22, 0.21]
table_4_results_FC[16,:] = [prices_initial_FC[5], prices_closed_CM_open_trade_FC[5],
    prices_open_CM_open_trade_FC[5], 0, prices_closed_CM_open_trade_FC[5],
    prices_closed_CM_open_trade_FC[5]*1.41]

total_credit_initial_FC =
    -(domestic_firm_debt_initial_FC[1]+exporter_firm_debt_initial_FC[1]) / nomGDP_initial_FC[1]
total_credit_closed_CM_FC =
    -(domestic_firm_debt_closed_CM_open_trade_FC[1]+exporter_firm_debt_closed_CM_open_trade_FC[1]) /
    nomGDP_closed_CM_open_trade_FC[1]
total_credit_open_CM_FC2 =
    -(domestic_firm_debt_open_CM_open_trade_FC[1]+exporter_firm_debt_open_CM_open_trade_FC[1]) /
    nomGDP_open_CM_open_trade_FC[1]
table_4_results_FC[17,:] = [total_credit_initial_FC, total_credit_closed_CM_FC,
    total_credit_open_CM_FC2, 43, 54, 62]

println("[tab:results_fixed_cost] table_4_results_FC (20×6) computed.")

# --- Save Tables/Table_results_FC.tex ---
let
    r = table_4_results_FC
    b = table_4_results
    f2(x) = @sprintf("%.2f", round(x, digits=2, RoundNearestTiesAway))
    f1(x) = @sprintf("%.1f", round(x, digits=1, RoundNearestTiesAway))

    tex = """\\begin{table}[!ht]
\\centering
\\caption{Comparing the impact of reforms}
\\label{tab:results_fixed_cost}
\\scalebox{0.7}{
\\begin{tabular}{l ccc | ccc | ccc}
\\\\[-1.8ex]\\hline
\\hline \\\\[-1.8ex]
Entry cost & \\multicolumn{3}{c|}{Per-period} & \\multicolumn{3}{c|}{Sunk} & \\multicolumn{3}{c}{Data} \\\\
Integration & None & Trade & Trade and capital & None & Trade & Trade and capital & 1991 & 2001 & 2008 \\\\
\\hline \\\\[-1.8ex]
\\textbf{Productivity} & & & & & & & & & \\\\
TFP & $(fpct(r[1,1])) & $(fpct(r[1,2])) & $(fpct(r[1,3])) & $(fpct(b[1,1])) & $(fpct(b[1,2])) & $(fpct(b[1,3])) & 100 & 116 & 132 \\\\
s.d of ARPK & $(f2(r[2,1])) & $(f2(r[2,2])) & $(f2(r[2,3])) & $(f2(b[2,1])) & $(f2(b[2,2])) & $(f2(b[2,3])) & - & 1.27 & 1.51 \\\\[1ex]
\\textbf{Aggregates} & & & & & & & & & \\\\
Output & $(fpct(r[3,1])) & $(fpct(r[3,2])) & $(fpct(r[3,3])) & $(fpct(b[3,1])) & $(fpct(b[3,2])) & $(fpct(b[3,3])) & 100 & 158 & 249 \\\\
Relative Output pc. & $(fpct(r[4,1])) & $(fpct(r[4,2])) & $(fpct(r[4,3])) & $(fpct(b[4,1])) & $(fpct(b[4,2])) & $(fpct(b[4,3])) & 54 & 57 & 64 \\\\
Consumption & $(fpct(r[5,1])) & $(fpct(r[5,2])) & $(fpct(r[5,3])) & $(fpct(b[5,1])) & $(fpct(b[5,2])) & $(fpct(b[5,3])) & 100 & 116 & 133 \\\\
Capital & $(fpct(r[6,1])) & $(fpct(r[6,2])) & $(fpct(r[6,3])) & $(fpct(b[6,1])) & $(fpct(b[6,2])) & $(fpct(b[6,3])) & 100 & 124 & 154 \\\\[1ex]
\\textbf{Welfare and Inequality} & & & & & & & & & \\\\
CE Welfare change & 0 & $(f1(100*r[8,2])) & $(f1(100*r[8,3])) & 0 & $(f1(100*b[8,2])) & $(f1(100*b[8,3]))* & & & \\\\
Top \$10\\%\$ income share & $(fpct(r[9,1])) & $(fpct(r[9,2])) & $(fpct(r[9,3])) & $(fpct(b[9,1])) & $(fpct(b[9,2])) & $(fpct(b[9,3])) & 24 & 30 & 34 \\\\[1ex]
\\textbf{Factor prices} & & & & & & & & & \\\\
Real wage & $(fpct(r[10,1])) & $(fpct(r[10,2])) & $(fpct(r[10,3])) & $(fpct(b[10,1])) & $(fpct(b[10,2])) & $(fpct(b[10,3])) & - & 100 & 144 \\\\
Interest rate premium & $(fpct(r[11,1])) & $(fpct(r[11,2])) & $(fpct(r[11,3])) & $(fpct(b[11,1])) & $(fpct(b[11,2])) & $(fpct(b[11,3])) & 21 & 5 & 3 \\\\[1ex]
\\textbf{Trade} & & & & & & & & & \\\\
\$\\frac{\\text{Import}}{\\text{GDP}}\$ & $(fpct(r[12,1])) & $(fpct(r[12,2])) & $(fpct(r[12,3])) & $(fpct(b[12,1])) & $(fpct(b[12,2])) & $(fpct(b[12,3])) & 21 & 42 & 42 \\\\
\$\\frac{\\text{NFA}}{\\text{GDP}}\$ & $(fpct(r[13,1])) & $(fpct(r[13,2])) & $(fpct(r[13,3])) & $(fpct(b[13,1])) & $(fpct(b[13,2])) & $(fpct(b[13,3])) & 1 & -5 & -6 \\\\
Entrepreneurship rate & $(fpct(r[14,1])) & $(fpct(r[14,2])) & $(fpct(r[14,3])) & $(fpct(b[14,1])) & $(fpct(b[14,2])) & $(fpct(b[14,3])) & - & - & 7 \\\\
Share of exporters & $(fpct(r[15,1])) & $(fpct(r[15,2])) & $(fpct(r[15,3])) & $(fpct(b[15,1])) & $(fpct(b[15,2])) & $(fpct(b[15,3])) & - & 22 & 21 \\\\
REER & $(fpct(r[16,1])) & $(fpct(r[16,2])) & $(fpct(r[16,3])) & $(fpct(b[16,1])) & $(fpct(b[16,2])) & $(fpct(b[16,3])) & - & 133 & 187 \\\\
\$\\frac{\\text{Credit}}{\\text{GDP}}\$ & $(fpct(r[17,1])) & $(fpct(r[17,2])) & $(fpct(r[17,3])) & $(fpct(b[17,1])) & $(fpct(b[17,2])) & $(fpct(b[17,3])) & 43 & 54 & 62 \\\\
\\hline
\\hline
\\end{tabular}
}
\\end{table}
"""
    open("../Tables/Table_results_FC.tex", "w") do io
        write(io, tex)
    end
    println("  Saved → ../Tables/Table_results_FC.tex")
end

# =============================================================================
# [tab:results_fixed_cost_micro]  FC model firm-level aggregates
# main.tex: \label{tab:results_fixed_cost_micro}  (Appendix)
# Source: plot_FC.jl lines 127–221
# Prerequisites: Steps 1 and 3.
# table_5_results_FC (22×5): same structure as table_5_results
# Output: ../Tables/Table_results_FC_micro.tex
# =============================================================================

table_5_results_FC = zeros(22, 5)
table_5_results_FC[1,:] = [1,
    domestic_pop_closed_CM_open_trade_FC[1]/domestic_pop_initial_FC[1],
    domestic_pop_open_CM_open_trade_FC[1]/domestic_pop_initial_FC[1], 0.75,
    42781/29958 * domestic_pop_closed_CM_open_trade[1]/domestic_pop_initial[1]]
table_5_results_FC[2,:] = [1,
    exporter_pop_closed_CM_open_trade_FC[1]/exporter_pop_initial_FC[1],
    exporter_pop_open_CM_open_trade_FC[1]/exporter_pop_initial_FC[1], 1.34,
    11228/8414 * exporter_pop_closed_CM_open_trade[1]/exporter_pop_initial[1]]
table_5_results_FC[3,1] = K_x_initial_FC[1]/(K_x_initial_FC[1]+K_d_initial_FC[1])
table_5_results_FC[3,2] = K_x_closed_CM_open_trade_FC[1]/(K_x_closed_CM_open_trade_FC[1]+K_d_closed_CM_open_trade_FC[1])
table_5_results_FC[3,3] = K_x_open_CM_open_trade_FC[1]/(K_x_open_CM_open_trade_FC[1]+K_d_open_CM_open_trade_FC[1])
table_5_results_FC[3,4] = 0.64; table_5_results_FC[3,5] = 0.65
table_5_results_FC[4,1] = L_x_initial_FC[1]/(L_x_initial_FC[1]+L_d_initial_FC[1])
table_5_results_FC[4,2] = L_x_closed_CM_open_trade_FC[1]/(L_x_closed_CM_open_trade_FC[1]+L_d_closed_CM_open_trade_FC[1])
table_5_results_FC[4,3] = L_x_open_CM_open_trade_FC[1]/(L_x_open_CM_open_trade_FC[1]+L_d_open_CM_open_trade_FC[1])
table_5_results_FC[4,4] = 0.57; table_5_results_FC[4,5] = 0.56
table_5_results_FC[5,:] = [exporter_pop_initial_FC[1]/entry_share_to_exporter_initial_FC[1],
    exporter_pop_closed_CM_open_trade_FC[1]/entry_share_to_exporter_closed_CM_open_trade_FC[1],
    exporter_pop_open_CM_open_trade_FC[1]/entry_share_to_exporter_open_CM_open_trade_FC[1], 1.8, 4.8]
table_5_results_FC[6,1:3] = [1,
    sum(k_d_fine_closed_CM_open_trade_FC[:,1].*current_domestic_closed_CM_FC)/
        domestic_pop_closed_CM_open_trade_FC[1] * domestic_pop_initial_FC[1] /
        sum(k_d_fine_initial_FC[:,1].*current_domestic_initial_FC),
    sum(k_d_fine_open_CM_open_trade_FC[:,1].*current_domestic_open_CM_FC)/
        domestic_pop_open_CM_open_trade_FC[1] * domestic_pop_initial_FC[1] /
        sum(k_d_fine_initial_FC[:,1].*current_domestic_initial_FC)]
table_5_results_FC[6,4] = table_5_results_FC[6,2]; table_5_results_FC[6,5] = 1.09*table_5_results_FC[6,4]
table_5_results_FC[7,1:3] = [1,
    sum(k_x_fine_closed_CM_open_trade_FC[:,1].*current_exporter_closed_CM_FC)/
        exporter_pop_closed_CM_open_trade_FC[1] * exporter_pop_initial_FC[1] /
        sum(k_x_fine_initial_FC[:,1].*current_exporter_initial_FC),
    sum(k_x_fine_open_CM_open_trade_FC[:,1].*current_exporter_open_CM_FC)/
        exporter_pop_open_CM_open_trade_FC[1] * exporter_pop_initial_FC[1] /
        sum(k_x_fine_initial_FC[:,1].*current_exporter_initial_FC)]
table_5_results_FC[7,4] = table_5_results_FC[7,2]; table_5_results_FC[7,5] = 1.42*table_5_results_FC[7,4]
table_5_results_FC[8,1:3] = [1,
    sum(k_x_fine_closed_CM_open_trade_FC[:,1].*current_exporter_closed_CM_FC +
        k_d_fine_closed_CM_open_trade_FC[:,1].*current_domestic_closed_CM_FC) /
        (domestic_pop_closed_CM_open_trade_FC[1]+exporter_pop_closed_CM_open_trade_FC[1]) *
        (domestic_pop_initial_FC[1]+exporter_pop_initial_FC[1]) /
        sum(k_x_fine_initial_FC[:,1].*current_exporter_initial_FC +
            k_d_fine_initial_FC[:,1].*current_domestic_initial_FC),
    sum(k_x_fine_open_CM_open_trade_FC[:,1].*current_exporter_open_CM_FC +
        k_d_fine_open_CM_open_trade_FC[:,1].*current_domestic_open_CM_FC) /
        (domestic_pop_open_CM_open_trade_FC[1]+exporter_pop_open_CM_open_trade_FC[1]) *
        (domestic_pop_initial_FC[1]+exporter_pop_initial_FC[1]) /
        sum(k_x_fine_initial_FC[:,1].*current_exporter_initial_FC +
            k_d_fine_initial_FC[:,1].*current_domestic_initial_FC)]
table_5_results_FC[8,4] = table_5_results_FC[8,2]; table_5_results_FC[8,5] = 1.13*table_5_results_FC[8,4]
table_5_results_FC[9,:]  = [sd_MRPK_d_initial_FC[1], sd_MRPK_d_closed_CM_open_trade_FC[1],
    sd_MRPK_d_open_CM_open_trade_FC[1], 1.3, 1.55]
table_5_results_FC[10,:] = [sd_MRPK_x_initial_FC[1], sd_MRPK_x_closed_CM_open_trade_FC[1],
    sd_MRPK_x_open_CM_open_trade_FC[1], 1.21, 1.38]

export_rev_init_FC   = sum(rev_xx_fine_store_initial_FC[:,1] .* current_exporter_initial_FC[:,1])
dom_rev_init_FC      = sum(rev_dx_fine_store_initial_FC[:,1].*current_exporter_initial_FC[:,1] +
    rev_d_fine_store_initial_FC[:,1].*current_domestic_initial_FC[:,1])
export_sales_ratio_init_FC = export_rev_init_FC / (export_rev_init_FC + dom_rev_init_FC)

export_rev_cl_FC  = sum(rev_xx_fine_store_closed_CM_open_trade_FC[:,1].*current_exporter_closed_CM_FC[:,1])
dom_rev_cl_FC     = sum(rev_dx_fine_store_closed_CM_open_trade_FC[:,1].*current_exporter_closed_CM_FC[:,1] +
    rev_d_fine_store_closed_CM_open_trade_FC[:,1].*current_domestic_closed_CM_FC[:,1])
export_sales_ratio_cl_FC = export_rev_cl_FC / (export_rev_cl_FC + dom_rev_cl_FC)

export_rev_op_FC  = sum(rev_xx_fine_store_open_CM_open_trade_FC[:,1].*current_exporter_open_CM_FC[:,1])
dom_rev_op_FC     = sum(rev_dx_fine_store_open_CM_open_trade_FC[:,1].*current_exporter_open_CM_FC[:,1] +
    rev_d_fine_store_open_CM_open_trade_FC[:,1].*current_domestic_open_CM_FC[:,1])
export_sales_ratio_op_FC = export_rev_op_FC / (export_rev_op_FC + dom_rev_op_FC)

table_5_results_FC[11,:] = [export_sales_ratio_init_FC, export_sales_ratio_cl_FC,
    export_sales_ratio_op_FC, 0.06, 0.05]
table_5_results_FC[12,1:3] = table_4_results_FC[15,1:3]
table_5_results_FC[12,4:5] = table_4_results_FC[15,5:6]

exs_init_FC = rev_xx_fine_store_initial_FC[:,1] ./
    (rev_xx_fine_store_initial_FC[:,1] + rev_dx_fine_store_initial_FC[:,1])
tot_rev_exp_init_FC = rev_xx_fine_store_initial_FC[:,1] + rev_dx_fine_store_initial_FC[:,1]
int_margin_init_FC  = sum(exs_init_FC .* tot_rev_exp_init_FC) / sum(tot_rev_exp_init_FC)

exs_cl_FC = rev_xx_fine_store_closed_CM_open_trade_FC[:,1] ./
    (rev_xx_fine_store_closed_CM_open_trade_FC[:,1] + rev_dx_fine_store_closed_CM_open_trade_FC[:,1])
tot_rev_exp_cl_FC = rev_xx_fine_store_closed_CM_open_trade_FC[:,1] + rev_dx_fine_store_closed_CM_open_trade_FC[:,1]
int_margin_cl_FC  = sum(exs_cl_FC .* tot_rev_exp_cl_FC) / sum(tot_rev_exp_cl_FC)

exs_op_FC = rev_xx_fine_store_open_CM_open_trade_FC[:,1] ./
    (rev_xx_fine_store_open_CM_open_trade_FC[:,1] + rev_dx_fine_store_open_CM_open_trade_FC[:,1])
tot_rev_exp_op_FC = rev_xx_fine_store_open_CM_open_trade_FC[:,1] + rev_dx_fine_store_open_CM_open_trade_FC[:,1]
int_margin_op_FC  = sum(exs_op_FC .* tot_rev_exp_op_FC) / sum(tot_rev_exp_op_FC)

table_5_results_FC[13,:] = [int_margin_init_FC, int_margin_cl_FC, int_margin_op_FC, 0.29, 0.26]
table_5_results_FC[14,:] = [
    sum(tot_rev_exp_init_FC)/exporter_pop_initial_FC[1]/sum(rev_d_fine_store_initial_FC[:,1])*domestic_pop_initial_FC[1],
    sum(tot_rev_exp_cl_FC)/exporter_pop_closed_CM_open_trade_FC[1]/sum(rev_d_fine_store_closed_CM_open_trade_FC[:,1])*domestic_pop_closed_CM_open_trade_FC[1],
    sum(tot_rev_exp_op_FC)/exporter_pop_open_CM_open_trade_FC[1]/sum(rev_d_fine_store_open_CM_open_trade_FC[:,1])*domestic_pop_open_CM_open_trade_FC[1],
    2.04, 2.54]
table_5_results_FC[15,:] = [
    entry_share_to_exporter_initial_FC[1]/(entry_share_to_exporter_initial_FC[1]+domestic_pop_initial_FC[1]),
    entry_share_to_exporter_closed_CM_open_trade_FC[1]/
        (exporter_pop_closed_CM_open_trade_FC[1]+domestic_pop_closed_CM_open_trade_FC[1]),
    entry_share_to_exporter_open_CM_open_trade_FC[1]/
        (exporter_pop_open_CM_open_trade_FC[1]+domestic_pop_open_CM_open_trade_FC[1]),
    0.056, 0.058]
table_5_results_FC[16,:] = [
    exit_exporting_to_work_sum_initial_FC[1]/exporter_pop_initial_FC[1],
    exit_exporting_to_work_sum_closed_CM_open_trade_FC[1]/exporter_pop_closed_CM_open_trade_FC[1],
    exit_exporting_to_work_sum_open_CM_open_trade_FC[1]/exporter_pop_open_CM_open_trade_FC[1],
    0.048, 0.048]
table_5_results_FC[17,:] = [Misallocation_within_d_initial_FC[1], Misallocation_within_d_closed_CM_open_trade_FC[1],
    Misallocation_within_d_open_CM_open_trade_FC[1], 0, 0]
table_5_results_FC[18,:] = [Misallocation_within_x_initial_FC[1], Misallocation_within_x_closed_CM_open_trade_FC[1],
    Misallocation_within_x_open_CM_open_trade_FC[1], 0, 0]

# Compute table_6_results_FC: FC exporter distribution by wealth/productivity quadrant
y_start_fc6 = 5; x_start_fc6 = 18
mat_fc6_init = mat_creator(current_exporter_initial_FC)
mat_fc6_cl   = mat_creator(current_exporter_closed_CM_FC)
mat_fc6_op   = mat_creator(current_exporter_open_CM_FC)
table_6_results_FC = zeros(4, 3)
table_6_results_FC[4,1] = sum(mat_fc6_init[(x_start_fc6+1):end,(y_start_fc6+1):end])/sum(current_exporter_initial_FC)
table_6_results_FC[3,1] = sum(mat_fc6_init[1:x_start_fc6,(y_start_fc6+1):end])/sum(current_exporter_initial_FC)
table_6_results_FC[2,1] = sum(mat_fc6_init[(x_start_fc6+1):end,1:y_start_fc6])/sum(current_exporter_initial_FC)
table_6_results_FC[1,1] = sum(mat_fc6_init[1:x_start_fc6,1:y_start_fc6])/sum(current_exporter_initial_FC)
table_6_results_FC[4,2] = sum(mat_fc6_cl[(x_start_fc6+1):end,(y_start_fc6+1):end])/sum(current_exporter_closed_CM_FC)
table_6_results_FC[3,2] = sum(mat_fc6_cl[1:x_start_fc6,(y_start_fc6+1):end])/sum(current_exporter_closed_CM_FC)
table_6_results_FC[2,2] = sum(mat_fc6_cl[(x_start_fc6+1):end,1:y_start_fc6])/sum(current_exporter_closed_CM_FC)
table_6_results_FC[1,2] = sum(mat_fc6_cl[1:x_start_fc6,1:y_start_fc6])/sum(current_exporter_closed_CM_FC)
table_6_results_FC[4,3] = sum(mat_fc6_op[(x_start_fc6+1):end,(y_start_fc6+1):end])/sum(current_exporter_open_CM_FC)
table_6_results_FC[3,3] = sum(mat_fc6_op[1:x_start_fc6,(y_start_fc6+1):end])/sum(current_exporter_open_CM_FC)
table_6_results_FC[2,3] = sum(mat_fc6_op[(x_start_fc6+1):end,1:y_start_fc6])/sum(current_exporter_open_CM_FC)
table_6_results_FC[1,3] = sum(mat_fc6_op[1:x_start_fc6,1:y_start_fc6])/sum(current_exporter_open_CM_FC)

println("[tab:results_fixed_cost_micro] table_5_results_FC and table_6_results_FC computed.")

# --- Save Tables/Table_results_FC_micro.tex ---
let
    r  = table_5_results_FC
    b  = table_5_results
    r6 = table_6_results_FC
    b6 = table_6_results
    f2(x) = @sprintf("%.2f", round(x, digits=2, RoundNearestTiesAway))
    f1(x) = @sprintf("%.1f", round(x, digits=1, RoundNearestTiesAway))
    f1s(x) = let s = f1(x); endswith(s, ".0") ? s[1:end-2] : s; end

    tex = """\\begin{table}[!ht]
\\centering
\\caption{Comparing the impact of reforms}
\\label{tab:results_fixed_cost_micro}
\\scalebox{0.7}{
\\begin{tabular}{l ccc | ccc | ccc}
\\\\[-1.8ex]\\hline
\\hline \\\\[-1.8ex]
Entry cost & \\multicolumn{3}{c|}{Per-period} & \\multicolumn{3}{c|}{Sunk} & \\multicolumn{3}{c}{Data} \\\\
Integration & None & Trade & Trade and capital & None & Trade & Trade and capital & 2001 & 2008 \\\\
\\hline \\\\[-1.8ex]
\\textbf{Extensive margin} & & & & & & & & & \\\\
Non-exporting firms & $(fpct(r[1,1])) & $(fpct(r[1,2])) & $(fpct(r[1,3])) & $(fpct(b[1,1])) & $(fpct(b[1,2])) & $(fpct(b[1,3])) & $(fpct(b[1,4])) & $(fpct(b[1,5])) \\\\
Exporting firms & $(fpct(r[2,1])) & $(fpct(r[2,2])) & $(fpct(r[2,3])) & $(fpct(b[2,1])) & $(fpct(b[2,2])) & $(fpct(b[2,3])) & $(fpct(b[2,4])) & $(fpct(b[2,5])) \\\\[1ex]
\\textbf{Intensive margin} & & & & & & & & & \\\\
\$ \\%\$ of capital used by exporters & $(fpct(r[3,1])) & $(fpct(r[3,2])) & $(fpct(r[3,3])) & $(fpct(b[3,1])) & $(fpct(b[3,2])) & $(fpct(b[3,3])) & $(fpct(b[3,4])) & $(fpct(b[3,5])) \\\\
\$ \\%\$ of labor used by exporters & $(fpct(r[4,1])) & $(fpct(r[4,2])) & $(fpct(r[4,3])) & $(fpct(b[4,1])) & $(fpct(b[4,2])) & $(fpct(b[4,3])) & $(fpct(b[4,4])) & $(fpct(b[4,5])) \\\\
Avg. duration (years) of export status & $(f1(r[5,1])) & $(f1(r[5,2])) & $(f1(r[5,3])) & $(f1(b[5,1])) & $(f1(b[5,2])) & $(f1(b[5,3])) & $(f1(b[5,4])) & $(f1(b[5,5])) \\\\
Average capital size of non-exporters & $(fpct(r[6,1])) & $(fpct(r[6,2])) & $(fpct(r[6,3])) & $(fpct(b[6,1])) & $(fpct(b[6,2])) & $(fpct(b[6,3])) & $(fpct(b[6,4])) & $(fpct(b[6,5])) \\\\
Average capital size of exporters & $(fpct(r[7,1])) & $(fpct(r[7,2])) & $(fpct(r[7,3])) & $(fpct(b[7,1])) & $(fpct(b[7,2])) & $(fpct(b[7,3])) & $(fpct(b[7,4])) & $(fpct(b[7,5])) \\\\
Average capital size & $(fpct(r[8,1])) & $(fpct(r[8,2])) & $(fpct(r[8,3])) & $(fpct(b[8,1])) & $(fpct(b[8,2])) & $(fpct(b[8,3])) & $(fpct(b[8,4])) & $(fpct(b[8,5])) \\\\[1ex]
\\textbf{Standard deviation of ARPK} & & & & & & & & & \\\\
 Non-exporting & $(f2(r[9,1])) & $(f2(r[9,2])) & $(f2(r[9,3])) & $(f2(b[9,1])) & $(f2(b[9,2])) & $(f2(b[9,3])) & $(@sprintf("%.1f", b[9,4])) & $(f2(b[9,5])) \\\\
 Exporter & $(f2(r[10,1])) & $(f2(r[10,2])) & $(f2(r[10,3])) & $(f2(b[10,1])) & $(f2(b[10,2])) & $(f2(b[10,3])) & $(f2(b[10,4])) & $(f2(b[10,5])) \\\\[1ex]
\\textbf{Export sales decomposition} & & & & & & & & & \\\\
Export-sales ratio & $(fpct(r[11,1])) & $(fpct(r[11,2])) & $(fpct(r[11,3])) & $(fpct(b[11,1])) & $(fpct(b[11,2])) & $(fpct(b[11,3])) & $(fpct(b[11,4])) & $(fpct(b[11,5])) \\\\
Extensive margin & $(fpct(r[12,1])) & $(fpct(r[12,2])) & $(fpct(r[12,3])) & $(fpct(b[12,1])) & $(fpct(b[12,2])) & $(fpct(b[12,3])) & $(fpct(b[12,4])) & $(fpct(b[12,5])) \\\\
Intensive margin & $(fpct(r[13,1])) & $(fpct(r[13,2])) & $(fpct(r[13,3])) & $(fpct(b[13,1])) & $(fpct(b[13,2])) & $(fpct(b[13,3])) & $(fpct(b[13,4])) & $(fpct(b[13,5])) \\\\
Exporter premium & $(fpct(r[14,1])) & $(fpct(r[14,2])) & $(fpct(r[14,3])) & $(fpct(b[14,1])) & $(fpct(b[14,2])) & $(fpct(b[14,3])) & $(fpct(b[14,4])) & $(fpct(b[14,5])) \\\\
Starter rate & $(fpct(r[15,1])) & $(fpct(r[15,2])) & $(fpct(r[15,3])) & $(fpct(b[15,1])) & $(fpct(b[15,2])) & $(fpct(b[15,3])) & $(fpct(b[15,4])) & $(fpct(b[15,5])) \\\\
Stopper rate & $(fpct(r[16,1])) & $(fpct(r[16,2])) & $(fpct(r[16,3])) & $(fpct(b[16,1])) & $(fpct(b[16,2])) & $(fpct(b[16,3])) & $(fpct(b[16,4])) & $(fpct(b[16,5])) \\\\
\\textbf{Productivity loss} & & & & & & & & & \\\\
 Non-exporting & $(f1s(100*r[17,1])) & $(f1s(100*r[17,2])) & $(f1s(100*r[17,3])) & $(f1s(100*b[17,1])) & $(f1s(100*b[17,2])) & $(f1s(100*b[17,3])) \\\\
 Exporter & $(f1s(100*r[18,1])) & $(f1s(100*r[18,2])) & $(f1s(100*r[18,3])) & $(f1s(100*b[18,1])) & $(f1s(100*b[18,2])) & $(f1s(100*b[18,3])) \\\\
\\textbf{Distribution of exporters} & & & & & & & & & \\\\
 Low wealth and low productivity & $(fpct(r6[1,1])) & $(fpct(r6[1,2])) & $(fpct(r6[1,3])) & $(fpct(b6[1,1])) & $(fpct(b6[1,2])) & $(fpct(b6[1,3])) & \\multirow{2}{*}{13} & \\multirow{2}{*}{5} \\\\
Low wealth and high productivity & $(fpct(r6[2,1])) & $(fpct(r6[2,2])) & $(fpct(r6[2,3])) & $(fpct(b6[2,1])) & $(fpct(b6[2,2])) & $(fpct(b6[2,3])) \\\\
High wealth and low productivity & $(fpct(r6[3,1])) & $(fpct(r6[3,2])) & $(fpct(r6[3,3])) & $(fpct(b6[3,1])) & $(fpct(b6[3,2])) & $(fpct(b6[3,3])) & \\multirow{2}{*}{87} & \\multirow{2}{*}{95} \\\\
High wealth and high productivity & $(fpct(r6[4,1])) & $(fpct(r6[4,2])) & $(fpct(r6[4,3])) & $(fpct(b6[4,1])) & $(fpct(b6[4,2])) & $(fpct(b6[4,3])) \\\\
\\hline
\\hline
\\end{tabular}
}
\\end{table}
"""
    open("../Tables/Table_results_FC_micro.tex", "w") do io
        write(io, tex)
    end
    println("  Saved → ../Tables/Table_results_FC_micro.tex")
end

# =============================================================================
# [fig:combined_changes]  Capital, profit, exit maps — trade+CM vs trade-only
# main.tex: \label{fig:combined_changes} (appendix)
# Figures: Figure5a (capital ∆ integrated CM, Blues)
#          Figure5b (capital ∆ trade-only,    Oranges)
#          Figure6a (profit ∆ integrated CM,  Blues)
#          Figure6b (profit ∆ trade-only,     Oranges)
#          Figure7a (exit decision ∆ open CM vs initial, Blues)
#          Figure7b (exit decision ∆ closed CM vs initial, Oranges)
# Source: plot.jl lines 664–712
# =============================================================================
let
    # Figure5a — capital per exporter: open_CM_open_trade / initial (Blues)
    y_value = k_x_fine_open_CM_open_trade[:,1] ./ k_x_fine_initial[:,1]
    y_mat   = mat_creator(y_value)
    y_mat_applied = 100.0 .* y_mat[10:36, 1:100] .- 100
    contourf(y_mat_applied', levels=5, c=cgrad(:Blues_6, categorical=true, rev=true),
            xlabel="Productivity", ylabel="Wealth", axis=nothing, lw=0,
            tickfontsize=14, xguidefontsize=18, yguidefontsize=18,
            legendfontsize=18)
    savefig("../Figures/Figure5a.pdf")
    println("[fig:combined_changes] Saved Figure5a.pdf")

    # Figure5b — capital per exporter: closed_CM_open_trade / initial (Oranges)
    y_value = k_x_fine_closed_CM_open_trade[:,1] ./ k_x_fine_initial[:,1]
    y_mat   = mat_creator(y_value)
    y_mat_applied = 100.0 .* y_mat[10:36, 1:100] .- 100
    contourf(y_mat_applied', levels=5, c=cgrad(:Oranges_6, categorical=true, rev=true),
            xlabel="Productivity", ylabel="Wealth", axis=nothing, lw=0,
            tickfontsize=14, xguidefontsize=18, yguidefontsize=18,
            legendfontsize=18)
    savefig("../Figures/Figure5b.pdf")
    println("[fig:combined_changes] Saved Figure5b.pdf")

    # Figure6a — profit: open_CM_open_trade / initial (Blues)
    y_value = Profit_fine_x_open_CM_open_trade[:,1] ./ Profit_fine_x_initial[:,1]
    y_mat   = mat_creator(y_value)
    y_mat_applied = 100.0 .* y_mat[10:36, 1:100] .- 100
    contourf(y_mat_applied', levels=5, c=cgrad(:Blues_6, categorical=true, rev=true),
            xlabel="Productivity", ylabel="Wealth", axis=nothing, lw=0,
            tickfontsize=14, xguidefontsize=18, yguidefontsize=18,
            legendfontsize=18)
    savefig("../Figures/Figure6a.pdf")
    println("[fig:combined_changes] Saved Figure6a.pdf")

    # Figure6b — profit: closed_CM_open_trade / initial (Oranges)
    y_value = Profit_fine_x_closed_CM_open_trade[:,1] ./ Profit_fine_x_initial[:,1]
    y_mat   = mat_creator(y_value)
    y_mat_applied = 100.0 .* y_mat[10:36, 1:100] .- 100
    contourf(y_mat_applied', levels=5, c=cgrad(:Oranges_6, categorical=true, rev=true),
            xlabel="Productivity", ylabel="Wealth", axis=nothing, lw=0,
            tickfontsize=14, xguidefontsize=18, yguidefontsize=18,
            legendfontsize=18)
    savefig("../Figures/Figure6b.pdf")
    println("[fig:combined_changes] Saved Figure6b.pdf")

    # Figure7a — exit decision: open_CM vs initial (Blues 3-level)
    y_value = 2 .* (future_occupation_fine_open_CM_open_trade[:,3,1,1] .- 1) .+
              (future_occupation_fine_initial[:,3,1,1] .- 1)
    y_mat   = mat_creator(y_value)
    y_mat_applied = y_mat[1:36, 1:30]
    contourf(y_mat_applied', levels=3, c=cgrad(:Blues_3, categorical=true, rev=false),
            xlabel="Productivity", ylabel="Wealth", axis=nothing, colorbar=nothing, lw=0,
            tickfontsize=14, xguidefontsize=18, yguidefontsize=18,
            legendfontsize=18)
    savefig("../Figures/Figure7a.pdf")
    println("[fig:combined_changes] Saved Figure7a.pdf")

    # Figure7b — exit decision: closed_CM vs initial (Oranges 3-level)
    y_value = 2 .* (future_occupation_fine_closed_CM_open_trade[:,3,1,1] .- 1) .+
              (future_occupation_fine_initial[:,3,1,1] .- 1)
    y_mat   = mat_creator(y_value)
    y_mat_applied = y_mat[1:36, 1:30]
    contourf(y_mat_applied', levels=3, c=cgrad(:Oranges_3, categorical=true, rev=false),
            xlabel="Productivity", ylabel="Wealth", axis=nothing, colorbar=nothing, lw=0,
            tickfontsize=14, xguidefontsize=18, yguidefontsize=18,
            legendfontsize=18)
    savefig("../Figures/Figure7b.pdf")
    println("[fig:combined_changes] Saved Figure7b.pdf")
end

# =============================================================================
# [fig:Figure_CMI_only]  Capital-market-integration marginal-effect maps
# main.tex: \label{fig:Figure_CMI_only}
# Figures: Figure5c (capital ∆ open_CM / closed_CM, Blues)
#          Figure6c (profit ∆ open_CM / closed_CM, Blues)
#          Figure7c (exit decision open_CM vs closed_CM, Blues)
# Source: plot.jl lines 273–278, 695–700, 714–718
# =============================================================================
let
    # Figure5c — capital: open_CM_open_trade / closed_CM_open_trade (Blues)
    y_value = k_x_fine_open_CM_open_trade[:,1] ./ k_x_fine_closed_CM_open_trade[:,1]
    y_mat   = mat_creator(y_value)
    y_mat_applied = 100.0 .* y_mat[10:36, 1:100] .- 100
    contourf(y_mat_applied', levels=5, c=cgrad(:Blues_6, categorical=true, rev=true),
            xlabel="Productivity", ylabel="Wealth", axis=nothing, lw=0,
            tickfontsize=14, xguidefontsize=18, yguidefontsize=18,
            legendfontsize=18)
    savefig("../Figures/Figure5c.pdf")
    println("[fig:Figure_CMI_only] Saved Figure5c.pdf")

    # Figure6c — profit: open_CM_open_trade / closed_CM_open_trade (Blues)
    y_value = Profit_fine_x_open_CM_open_trade[:,1] ./ Profit_fine_x_closed_CM_open_trade[:,1]
    y_mat   = mat_creator(y_value)
    y_mat_applied = 100.0 .* y_mat[10:36, 1:100] .- 100
    contourf(y_mat_applied', levels=5, c=cgrad(:Blues_6, categorical=true, rev=true),
            xlabel="Productivity", ylabel="Wealth", axis=nothing, lw=0,
            tickfontsize=14, xguidefontsize=18, yguidefontsize=18,
            legendfontsize=18)
    savefig("../Figures/Figure6c.pdf")
    println("[fig:Figure_CMI_only] Saved Figure6c.pdf")

    # Figure7c — exit decision: open_CM vs closed_CM (Blues 3-level)
    y_value = 2 .* (future_occupation_fine_open_CM_open_trade[:,3,1,1] .- 1) .+
              (future_occupation_fine_closed_CM_open_trade[:,3,1,1] .- 1)
    y_mat   = mat_creator(y_value)
    y_mat_applied = y_mat[1:36, 1:30]
    contourf(y_mat_applied', levels=3, c=cgrad(:Blues_3, categorical=true, rev=false),
            xlabel="Productivity", ylabel="Wealth", axis=nothing, colorbar=nothing, lw=0,
            tickfontsize=14, xguidefontsize=18, yguidefontsize=18,
            legendfontsize=18)
    savefig("../Figures/Figure7c.pdf")
    println("[fig:Figure_CMI_only] Saved Figure7c.pdf")
end

# =============================================================================
# [fig:Figure12]  Transition dynamics after trade liberalization
# main.tex: \label{fig:Figure12}
# Figures: Figure12a (TFP), 12b (capital), 12c (consumption), 12d (exporters),
#          12e (non-exporters), 12f (exporter share), 12g (K/exporter),
#          12h (K/non-exporter), 12i (capital share exporters),
#          12k (s.d. ARPK all), 12l (s.d. ARPK domestic), 12m (s.d. ARPK exporter)
#          [12j skipped — commented out in main.tex]
# Source: plot.jl lines 730–882
# =============================================================================
let
    # --- TFP (Figure12a) ---
    TFP_corr_open    = TFP_open_CM_open_trade[1] / TFP_store_open_CM_open_trade[1,end]
    TFP_corr_closed  = TFP_closed_CM_open_trade[1] / TFP_store_closed_CM_open_trade[1,end-1]
    TFP_corr_delayed = TFP_open_CM_open_trade[1] / TFP_store_open_CMdelayed_open_trade[1,end-1]
    TFP_smooth_open    = movmean(100 .* (TFP_corr_open    .* TFP_store_open_CM_open_trade[1,:] ./ TFP_initial[1] .- 1), 5)
    TFP_smooth_closed  = movmean(100 .* (TFP_corr_closed  .* TFP_store_closed_CM_open_trade[1,:] ./ TFP_initial[1] .- 1), 5)
    TFP_smooth_delayed = movmean(100 .* (TFP_corr_delayed .* TFP_store_open_CMdelayed_open_trade[1,:] ./ TFP_initial[1] .- 1), 5)
    plot([range(1991,2020)], [TFP_smooth_open, TFP_smooth_closed, TFP_smooth_delayed],
         legend=(0.5,0.7), linewidth=5, linestyle=[:solid :dash :dot], grid=false,
         tickfontsize=14, ylims=[4,35],
         label=["Immediate full integration" "Historical" "Only trade"], legendfont=font(14))
    savefig("../Figures/Figure12a.pdf")
    println("[fig:Figure12] Saved Figure12a.pdf")

    # --- Capital (Figure12b) ---
    cap_corr_open    = capital_demand_open_CM_open_trade[1] / capital_supply_future_open_CMdelayed_open_trade[1,end]
    cap_corr_closed  = capital_demand_closed_CM_open_trade[1] / capital_supply_future_closed_CM_open_trade[1,end]
    cap_corr_delayed = capital_demand_open_CM_open_trade[1] / capital_supply_future_open_CMdelayed_open_trade[1,end-1]
    cap_smooth_open    = movmean(100 .* (cap_corr_open    .* capital_supply_future_open_CM_open_trade[1,:] ./ capital_demand_initial[1] .- 1), 5)
    cap_smooth_closed  = movmean(100 .* (cap_corr_closed  .* capital_supply_future_closed_CM_open_trade[1,:] ./ capital_demand_initial[1] .- 1), 5)
    cap_smooth_delayed = movmean(100 .* (cap_corr_delayed .* capital_supply_future_open_CMdelayed_open_trade[1,:] ./ capital_demand_initial[1] .- 1), 5)
    plot([range(1991,2022)], [cap_smooth_open, cap_smooth_closed, cap_smooth_delayed],
         legend=nothing, linewidth=5, linestyle=[:solid :dash :dot], grid=false,
         tickfontsize=14, ylims=[-10,40])
    savefig("../Figures/Figure12b.pdf")
    println("[fig:Figure12] Saved Figure12b.pdf")

    # --- Consumption (Figure12c) ---
    con_corr_open    = total_consumption_open_CM_open_trade[1] / total_consumption_store_open_CM_open_trade[1,end]
    con_corr_closed  = total_consumption_closed_CM_open_trade[1] / total_consumption_store_closed_CM_open_trade[1,end]
    con_corr_delayed = total_consumption_open_CM_open_trade[1] / total_consumption_store_open_CMdelayed_open_trade[1,end-1]
    con_smooth_open    = movmean(100 .* (con_corr_open    .* total_consumption_store_open_CM_open_trade[1,:] ./ total_consumption_initial[1] .- 1), 5)
    con_smooth_closed  = movmean(100 .* (con_corr_closed  .* total_consumption_store_closed_CM_open_trade[1,:] ./ total_consumption_initial[1] .- 1), 5)
    con_smooth_delayed = movmean(100 .* (con_corr_delayed .* total_consumption_store_open_CMdelayed_open_trade[1,:] ./ total_consumption_initial[1] .- 1), 5)
    plot([range(1991,2020)], [con_smooth_open, con_smooth_closed, con_smooth_delayed],
         legend=nothing, linewidth=5, linestyle=[:solid :dash :dot], grid=false,
         tickfontsize=14, ylims=[0,10])
    savefig("../Figures/Figure12c.pdf")
    println("[fig:Figure12] Saved Figure12c.pdf")

    # --- Exporters pop (Figure12d) ---
    exp_corr_open    = exporter_pop_open_CM_open_trade[1] / exporter_pop_store_open_CM_open_trade[1,end]
    exp_corr_closed  = exporter_pop_closed_CM_open_trade[1] / exporter_pop_store_closed_CM_open_trade[1,end]
    exp_corr_delayed = exporter_pop_open_CM_open_trade[1] / exporter_pop_store_open_CMdelayed_open_trade[1,end-1]
    exp_smooth_open    = movmean(100 .* (exp_corr_open    .* exporter_pop_store_open_CM_open_trade[1,:] ./ exporter_pop_initial[1] .- 1), 5)
    exp_smooth_closed  = movmean(100 .* (exp_corr_closed  .* exporter_pop_store_closed_CM_open_trade[1,:] ./ exporter_pop_initial[1] .- 1), 5)
    exp_smooth_delayed = movmean(100 .* (exp_corr_delayed .* exporter_pop_store_open_CMdelayed_open_trade[1,:] ./ exporter_pop_initial[1] .- 1), 5)
    plot([range(1991,2020)], [exp_smooth_open, exp_smooth_closed, exp_smooth_delayed],
         legend=nothing, linewidth=3, linestyle=[:solid :dash :dot], grid=false,
         tickfontsize=14, ylims=[0.0,45])
    savefig("../Figures/Figure12d.pdf")
    println("[fig:Figure12] Saved Figure12d.pdf")

    # --- Domestic pop (Figure12e) ---
    dom_corr_open    = domestic_pop_open_CM_open_trade[1] / domestic_pop_store_open_CM_open_trade[1,end]
    dom_corr_closed  = domestic_pop_closed_CM_open_trade[1] / domestic_pop_store_closed_CM_open_trade[1,end]
    dom_corr_delayed = domestic_pop_open_CM_open_trade[1] / domestic_pop_store_open_CMdelayed_open_trade[1,end-1]
    dom_smooth_open    = movmean(100 .* (dom_corr_open    .* domestic_pop_store_open_CM_open_trade[1,:] ./ domestic_pop_initial[1] .- 1), 5)
    dom_smooth_closed  = movmean(100 .* (dom_corr_closed  .* domestic_pop_store_closed_CM_open_trade[1,:] ./ domestic_pop_initial[1] .- 1), 5)
    dom_smooth_delayed = movmean(100 .* (dom_corr_delayed .* domestic_pop_store_open_CMdelayed_open_trade[1,:] ./ domestic_pop_initial[1] .- 1), 5)
    plot([range(1991,2020)], [dom_smooth_open, dom_smooth_closed, dom_smooth_delayed],
         legend=nothing, linewidth=3, linestyle=[:solid :dash :dot], grid=false,
         tickfontsize=14, ylims=[-30.0,10])
    savefig("../Figures/Figure12e.pdf")
    println("[fig:Figure12] Saved Figure12e.pdf")

    # --- Exporter share (Figure12f) ---
    share_open    = movmean(exporter_pop_store_open_CM_open_trade[1,:] ./ (exporter_pop_store_open_CM_open_trade[1,:] .+ domestic_pop_store_open_CM_open_trade[1,:]), 5)
    share_closed  = movmean(exporter_pop_store_closed_CM_open_trade[1,:] ./ (exporter_pop_store_closed_CM_open_trade[1,:] .+ domestic_pop_store_closed_CM_open_trade[1,:]), 5)
    share_delayed = movmean(exporter_pop_store_open_CMdelayed_open_trade[1,:] ./ (exporter_pop_store_open_CMdelayed_open_trade[1,:] .+ domestic_pop_store_open_CMdelayed_open_trade[1,:]), 5)
    plot([range(1991,2020)], 100 .* [share_open, share_closed, share_delayed],
         legend=nothing, linewidth=3, linestyle=[:solid :dash :dot], grid=false,
         tickfontsize=14, ylims=[30,50])
    savefig("../Figures/Figure12f.pdf")
    println("[fig:Figure12] Saved Figure12f.pdf")

    # --- Capital per exporter (Figure12g) ---
    Kx_init_per_exp = K_x_initial[1] / exporter_pop_initial[1]
    Kxcorr_open    = K_x_open_CM_open_trade[1]/exporter_pop_open_CM_open_trade[1] * exporter_pop_store_open_CM_open_trade[1,end-1] / (K_x_ratio_store_open_CM_open_trade[1,end] * capital_supply_future_open_CM_open_trade[1,end])
    Kxcorr_closed  = K_x_closed_CM_open_trade[1]/exporter_pop_closed_CM_open_trade[1] * exporter_pop_store_closed_CM_open_trade[1,end-1] / (K_x_ratio_store_closed_CM_open_trade[1,end] * capital_supply_future_closed_CM_open_trade[1,end])
    Kxcorr_delayed = K_x_open_CM_open_trade[1]/exporter_pop_open_CM_open_trade[1] * exporter_pop_store_open_CMdelayed_open_trade[1,end-1] / (K_x_ratio_store_open_CMdelayed_open_trade[1,end-1] * capital_supply_future_open_CMdelayed_open_trade[1,end])
    Kx_smooth_open    = movmean(100 .* ((Kxcorr_open    .* K_x_ratio_store_open_CM_open_trade[1,:] .* capital_supply_future_open_CM_open_trade[1,2:(end-1)]) ./ exporter_pop_store_open_CM_open_trade[1,:] ./ Kx_init_per_exp .- 1), 5)
    Kx_smooth_closed  = movmean(100 .* ((Kxcorr_closed  .* K_x_ratio_store_closed_CM_open_trade[1,:] .* capital_supply_future_closed_CM_open_trade[1,2:(end-1)]) ./ exporter_pop_store_closed_CM_open_trade[1,:] ./ Kx_init_per_exp .- 1), 5)
    Kx_smooth_delayed = movmean(100 .* ((Kxcorr_delayed .* K_x_ratio_store_open_CMdelayed_open_trade[1,:] .* capital_supply_future_open_CMdelayed_open_trade[1,2:(end-1)]) ./ exporter_pop_store_open_CMdelayed_open_trade[1,:] ./ Kx_init_per_exp .- 1), 5)
    plot([range(1991,2020)], [Kx_smooth_open, Kx_smooth_closed, Kx_smooth_delayed],
         legend=nothing, linewidth=3, linestyle=[:solid :dash :dot], grid=false,
         tickfontsize=14, ylims=[-20.0,40])
    savefig("../Figures/Figure12g.pdf")
    println("[fig:Figure12] Saved Figure12g.pdf")

    # --- Capital per non-exporter (Figure12h) ---
    Kd_init_per_dom = K_d_initial[1] / domestic_pop_initial[1]
    Kdcorr_open    = K_d_open_CM_open_trade[1]/domestic_pop_open_CM_open_trade[1] * domestic_pop_store_open_CM_open_trade[1,end-1] / (K_d_ratio_store_open_CM_open_trade[1,end] * capital_supply_future_open_CM_open_trade[1,end])
    Kdcorr_closed  = K_d_closed_CM_open_trade[1]/domestic_pop_closed_CM_open_trade[1] * domestic_pop_store_closed_CM_open_trade[1,end-1] / (K_d_ratio_store_closed_CM_open_trade[1,end] * capital_supply_future_closed_CM_open_trade[1,end])
    Kdcorr_delayed = K_d_open_CM_open_trade[1]/domestic_pop_open_CM_open_trade[1] * domestic_pop_store_open_CMdelayed_open_trade[1,end-1] / (K_d_ratio_store_open_CMdelayed_open_trade[1,end-1] * capital_supply_future_open_CMdelayed_open_trade[1,end])
    Kd_smooth_open    = movmean(100 .* ((Kdcorr_open    .* K_d_ratio_store_open_CM_open_trade[1,:] .* capital_supply_future_open_CM_open_trade[1,2:(end-1)]) ./ domestic_pop_store_open_CM_open_trade[1,:] ./ Kd_init_per_dom .- 1), 5)
    Kd_smooth_closed  = movmean(100 .* ((Kdcorr_closed  .* K_d_ratio_store_closed_CM_open_trade[1,:] .* capital_supply_future_closed_CM_open_trade[1,2:(end-1)]) ./ domestic_pop_store_closed_CM_open_trade[1,:] ./ Kd_init_per_dom .- 1), 5)
    Kd_smooth_delayed = movmean(100 .* ((Kdcorr_delayed .* K_d_ratio_store_open_CMdelayed_open_trade[1,:] .* capital_supply_future_open_CMdelayed_open_trade[1,2:(end-1)]) ./ domestic_pop_store_open_CMdelayed_open_trade[1,:] ./ Kd_init_per_dom .- 1), 5)
    plot([range(1991,2020)], [Kd_smooth_open, Kd_smooth_closed, Kd_smooth_delayed],
         legend=nothing, linewidth=3, linestyle=[:solid :dash :dot], grid=false,
         tickfontsize=14, ylims=[-20.0,20])
    savefig("../Figures/Figure12h.pdf")
    println("[fig:Figure12] Saved Figure12h.pdf")

    # --- Capital share of exporters (Figure12i) ---
    Kshare_open    = movmean(100 .* K_x_ratio_store_open_CM_open_trade[1,1:end], 5)
    Kshare_closed  = movmean(100 .* K_x_ratio_store_closed_CM_open_trade[1,1:end], 5)
    Kshare_delayed = movmean(100 .* K_x_ratio_store_open_CMdelayed_open_trade[1,1:end], 5)
    plot([range(1991,2020)], [Kshare_open, Kshare_closed, Kshare_delayed],
         legend=nothing, linewidth=3, linestyle=[:solid :dash :dot], grid=false,
         tickfontsize=14, ylims=[45,80])
    savefig("../Figures/Figure12i.pdf")
    println("[fig:Figure12] Saved Figure12i.pdf")

    # --- s.d. ARPK all (Figure12k) ---
    sdMRPK_corr_open    = sd_MRPK_open_CM_open_trade[1] / sd_MRPK_store_open_CM_open_trade[1,end]
    sdMRPK_corr_closed  = sd_MRPK_closed_CM_open_trade[1] / sd_MRPK_store_closed_CM_open_trade[1,end]
    sdMRPK_corr_delayed = sd_MRPK_open_CM_open_trade[1] / sd_MRPK_store_open_CMdelayed_open_trade[1,end-1]
    sdMRPK_smooth_open    = movmean(sdMRPK_corr_open    .* sd_MRPK_store_open_CM_open_trade[1,:], 5)
    sdMRPK_smooth_closed  = movmean(sdMRPK_corr_closed  .* sd_MRPK_store_closed_CM_open_trade[1,:], 5)
    sdMRPK_smooth_delayed = movmean(sdMRPK_corr_delayed .* sd_MRPK_store_open_CMdelayed_open_trade[1,:], 5)
    plot([range(1991,2020)], [sdMRPK_smooth_open, sdMRPK_smooth_closed, sdMRPK_smooth_delayed],
         legend=nothing, linewidth=3, linestyle=[:solid :dash :dot], grid=false,
         tickfontsize=14, ylims=[0.3,0.55])
    savefig("../Figures/Figure12k.pdf")
    println("[fig:Figure12] Saved Figure12k.pdf")

    # --- s.d. ARPK domestic (Figure12l) ---
    sdMRPKd_corr_open    = sd_MRPK_d_open_CM_open_trade[1] / sd_MRPK_d_store_open_CM_open_trade[1,end]
    sdMRPKd_corr_closed  = sd_MRPK_d_closed_CM_open_trade[1] / sd_MRPK_d_store_closed_CM_open_trade[1,end]
    sdMRPKd_corr_delayed = sd_MRPK_d_open_CM_open_trade[1] / sd_MRPK_d_store_open_CMdelayed_open_trade[1,end-1]
    sdMRPKd_smooth_open    = movmean(sdMRPKd_corr_open    .* sd_MRPK_d_store_open_CM_open_trade[1,:], 5)
    sdMRPKd_smooth_closed  = movmean(sdMRPKd_corr_closed  .* sd_MRPK_d_store_closed_CM_open_trade[1,:], 5)
    sdMRPKd_smooth_delayed = movmean(sdMRPKd_corr_delayed .* sd_MRPK_d_store_open_CMdelayed_open_trade[1,:], 5)
    plot([range(1991,2020)], [sdMRPKd_smooth_open, sdMRPKd_smooth_closed, sdMRPKd_smooth_delayed],
         legend=nothing, linewidth=5, linestyle=[:solid :dash :dot], grid=false,
         tickfontsize=14, ylims=[0.32,0.55])
    savefig("../Figures/Figure12l.pdf")
    println("[fig:Figure12] Saved Figure12l.pdf")

    # --- s.d. ARPK exporter (Figure12m) ---
    sdMRPKx_corr_open    = sd_MRPK_x_open_CM_open_trade[1] / sd_MRPK_x_store_open_CM_open_trade[1,end]
    sdMRPKx_corr_closed  = sd_MRPK_x_closed_CM_open_trade[1] / sd_MRPK_x_store_closed_CM_open_trade[1,end]
    sdMRPKx_corr_delayed = sd_MRPK_x_open_CM_open_trade[1] / sd_MRPK_x_store_open_CMdelayed_open_trade[1,end-1]
    sdMRPKx_smooth_open    = movmean(sdMRPKx_corr_open    .* sd_MRPK_x_store_open_CM_open_trade[1,:], 5)
    sdMRPKx_smooth_closed  = movmean(sdMRPKx_corr_closed  .* sd_MRPK_x_store_closed_CM_open_trade[1,:], 5)
    sdMRPKx_smooth_delayed = movmean(sdMRPKx_corr_delayed .* sd_MRPK_x_store_open_CMdelayed_open_trade[1,:], 5)
    # Note: plot.jl uses sdMRPKd_smooth_delayed (domestic) for the delayed line in Figure12m
    plot([range(1991,2020)], [sdMRPKx_smooth_open, sdMRPKx_smooth_closed, sdMRPKd_smooth_delayed],
         legend=nothing, linewidth=5, linestyle=[:solid :dash :dot], grid=false,
         tickfontsize=14, ylims=[0.25,0.55])
    savefig("../Figures/Figure12m.pdf")
    println("[fig:Figure12] Saved Figure12m.pdf")
end

# =============================================================================
# [fig:Figure12c]  DiD simulation — incumbent exporters by wealth quartile
# main.tex: \label{fig:Figure12c}
# Figures: Figure12q (above/below median wealth, immediate integration)
#          Figure12r (above/below median wealth, delayed integration)
# Source: plot.jl lines 288–392
# Note: requires current_*_mat distributions and K_x_ratio_store_open_CM_open_trade_initial
#       from transition case 0 (open_CM_open_trade_initial).
# =============================================================================
let
    y_start = 5   # median wealth cut-off index (as in plot.jl line 292)
    T_cap   = 10  # delay length (set by transition case 4)

    # Matrix representations of distributions in closed_CM_open_trade (starting state)
    k_x_cl_mat = mat_creator(k_x_fine_closed_CM_open_trade[:,1])
    k_d_cl_mat = mat_creator(k_d_fine_closed_CM_open_trade[:,1])
    k_x_op_mat = mat_creator(k_x_fine_open_CM_open_trade[:,1])
    k_d_op_mat = mat_creator(k_d_fine_open_CM_open_trade[:,1])

    occ_cl_val = (future_occupation_fine_closed_CM_open_trade[:,3,1,1] .== 3.0)
    occ_op_val = (future_occupation_fine_open_CM_open_trade[:,3,1,1] .== 3.0)
    occ_cl_mat = mat_creator(occ_cl_val)
    occ_op_mat = mat_creator(occ_op_val)

    past_exp_cl = distr_current_closed_CM_open_trade[(Baseline_parameter.n_fine[1]*Baseline_parameter.n_fine[2]*2+1):(Baseline_parameter.n_fine[1]*Baseline_parameter.n_fine[2]*3),1]
    past_exp_op = distr_current_open_CM_open_trade[(Baseline_parameter.n_fine[1]*Baseline_parameter.n_fine[2]*2+1):(Baseline_parameter.n_fine[1]*Baseline_parameter.n_fine[2]*3),1]
    past_exp_cl_mat = mat_creator(past_exp_cl)
    past_exp_op_mat = mat_creator(past_exp_op)

    denom_cl = sum(current_exporter_closed_CM_open_trade_mat .* k_x_cl_mat .+ current_domestic_closed_CM_open_trade_mat .* k_d_cl_mat)
    denom_op = sum(current_exporter_open_CM_open_trade_mat   .* k_x_op_mat .+ current_domestic_open_CM_open_trade_mat   .* k_d_op_mat)

    # Above-median wealth share of incumbent exporter capital
    cap_share_cl_med = sum(past_exp_cl_mat[:,(y_start+1):end] .* k_x_cl_mat[:,(y_start+1):end] .* occ_cl_mat[:,(y_start+1):end]) / denom_cl
    cap_share_op_med = sum(past_exp_op_mat[:,(y_start+1):end] .* k_x_op_mat[:,(y_start+1):end] .* occ_op_mat[:,(y_start+1):end]) / denom_op
    cap_share_median = range(cap_share_cl_med, cap_share_op_med, length=30)
    cap_inc_trans_med = cap_share_median .* K_x_ratio_store_open_CM_open_trade_initial[1,:] .* capital_supply_future_open_CM_open_trade_initial[1,2:(end-1)]
    K_x_med = movmean(100 .* (cap_inc_trans_med ./ (cap_share_cl_med * K_x_closed_CM_open_trade[1]) .- 1), 5)

    # Below-median wealth share
    cap_share_cl_bmed = sum(past_exp_cl_mat[:,1:(y_start+1)] .* k_x_cl_mat[:,1:(y_start+1)] .* occ_cl_mat[:,1:(y_start+1)]) / denom_cl
    cap_share_op_bmed = sum(past_exp_op_mat[:,1:(y_start+1)] .* k_x_op_mat[:,1:(y_start+1)] .* occ_op_mat[:,1:(y_start+1)]) / denom_op
    cap_share_bmedian = range(cap_share_cl_bmed, cap_share_op_bmed, length=30)
    cap_inc_trans_bmed = cap_share_bmedian .* K_x_ratio_store_open_CM_open_trade_initial[1,:] .* capital_supply_future_open_CM_open_trade_initial[1,2:(end-1)]
    K_x_bmed = movmean(100 .* (cap_inc_trans_bmed ./ (cap_share_cl_bmed * K_x_closed_CM_open_trade[1]) .- 1), 5)

    plot([range(2001,2008)], [K_x_med[1:8], K_x_bmed[1:8]],
         legend=(0.1,0.9), label=["Above median wealth" "Below median wealth"],
         legendfont=font(14), linewidth=3, linestyle=[:solid :dash],
         grid=false, tickfontsize=14, ylims=[0.0,50])
    savefig("../Figures/Figure12q.pdf")
    println("[fig:Figure12c] Saved Figure12q.pdf")

    # Delayed integration version (Figure12r)
    cap_inc_trans_med_del  = cap_share_median  .* K_x_ratio_store_open_CMdelayed_open_trade[1,:] .* capital_supply_future_open_CMdelayed_open_trade[1,2:(end-1)]
    cap_inc_trans_bmed_del = cap_share_bmedian .* K_x_ratio_store_open_CMdelayed_open_trade[1,:] .* capital_supply_future_open_CMdelayed_open_trade[1,2:(end-1)]
    K_x_med_del  = 100 .* (cap_inc_trans_med_del  ./ cap_inc_trans_med_del[T_cap]  .- 1)
    K_x_bmed_del = 100 .* (cap_inc_trans_bmed_del ./ cap_inc_trans_bmed_del[T_cap] .- 1)
    plot([range(2001,2008)], [K_x_med_del[T_cap:(T_cap+7)], K_x_bmed_del[T_cap:(T_cap+7)]],
         legend=nothing, linewidth=3, linestyle=[:solid :dash],
         grid=false, tickfontsize=14, ylims=[-20.0,25])
    savefig("../Figures/Figure12r.pdf")
    println("[fig:Figure12c] Saved Figure12r.pdf")
end

# =============================================================================
# [fig:Figure12b]  Transition dynamics — top 10% distributional shares
# main.tex: \label{fig:Figure12b}
# Figures: Figure12n (top 10% wealth share)
#          Figure12o (top 10% income share)
#          Figure12p (top 10% consumption share)
# Source: plot.jl lines 912–944
# =============================================================================
let
    # Figure12n — p90 wealth share
    pw_corr_open    = p90_wealth_open_CM_open_trade[1] / p90_wealth_store_open_CM_open_trade[1,end]
    pw_corr_closed  = p90_wealth_closed_CM_open_trade[1] / p90_wealth_store_closed_CM_open_trade[1,end]
    pw_corr_delayed = p90_wealth_open_CM_open_trade[1] / p90_wealth_store_open_CMdelayed_open_trade[1,end-1]
    pw_smooth_open    = movmean(pw_corr_open    .* p90_wealth_store_open_CM_open_trade[1,:], 5)
    pw_smooth_closed  = movmean(pw_corr_closed  .* p90_wealth_store_closed_CM_open_trade[1,:], 5)
    pw_smooth_delayed = movmean(pw_corr_delayed .* p90_wealth_store_open_CMdelayed_open_trade[1,:], 5)
    plot(100 .* [pw_smooth_open, pw_smooth_closed, pw_smooth_delayed],
         legend=nothing, linewidth=3, linestyle=[:solid :dash :dot],
         grid=false, tickfontsize=14, ylims=[40,58])
    savefig("../Figures/Figure12n.pdf")
    println("[fig:Figure12b] Saved Figure12n.pdf")

    # Figure12o — p90 income share
    pi_corr_open    = p90_income_open_CM_open_trade[1] / p90_income_store_open_CM_open_trade[1,end]
    pi_corr_closed  = p90_income_closed_CM_open_trade[1] / p90_income_store_closed_CM_open_trade[1,end]
    pi_corr_delayed = p90_income_open_CM_open_trade[1] / p90_income_store_open_CMdelayed_open_trade[1,end-1]
    pi_smooth_open    = movmean(pi_corr_open    .* p90_income_store_open_CM_open_trade[1,:], 5)
    pi_smooth_closed  = movmean(pi_corr_closed  .* p90_income_store_closed_CM_open_trade[1,:], 5)
    pi_smooth_delayed = movmean(pi_corr_delayed .* p90_income_store_open_CMdelayed_open_trade[1,:], 5)
    plot(100 .* [pi_smooth_open, pi_smooth_closed, pi_smooth_delayed],
         legend=nothing, linewidth=5, linestyle=[:solid :dash :dot],
         grid=false, tickfontsize=14, ylims=[23,30])
    savefig("../Figures/Figure12o.pdf")
    println("[fig:Figure12b] Saved Figure12o.pdf")

    # Figure12p — p90 consumption share
    pc_corr_open    = p90_cons_open_CM_open_trade[1] / p90_cons_store_open_CM_open_trade[1,end]
    pc_corr_closed  = p90_cons_closed_CM_open_trade[1] / p90_cons_store_closed_CM_open_trade[1,end]
    pc_corr_delayed = p90_cons_open_CM_open_trade[1] / p90_cons_store_open_CMdelayed_open_trade[1,end-1]
    pc_smooth_open    = movmean(pc_corr_open    .* p90_cons_store_open_CM_open_trade[1,:], 5)
    pc_smooth_closed  = movmean(pc_corr_closed  .* p90_cons_store_closed_CM_open_trade[1,:], 5)
    pc_smooth_delayed = movmean(pc_corr_delayed .* p90_cons_store_open_CMdelayed_open_trade[1,:], 5)
    plot(100 .* [pc_smooth_open, pc_smooth_closed, pc_smooth_delayed],
         legend=nothing, linewidth=5, linestyle=[:solid :dash :dot],
         grid=false, tickfontsize=14, ylims=[15,25])
    savefig("../Figures/Figure12p.pdf")
    println("[fig:Figure12b] Saved Figure12p.pdf")
end

# =============================================================================
# [fig:welfare]  Consumption-equivalent welfare contour maps
# main.tex: \label{fig:welfare}
# Figures: Figure9a  (domestic producers, open CM, Blues)   — labeled "Workers" in paper
#          Figure9b  (domestic producers, closed CM, Oranges)
#          Figure10a (exporters, open CM, Blues)
#          Figure10b (exporters, closed CM, Oranges)
#          Figure11a (workers, open-closed CM difference, Blues)
#          Figure11c (exporters, open-closed CM difference, Blues)
# Source: plot.jl lines 947–996
# Note: Figure8a/8b (workers absolute welfare) not in main.tex — omitted.
#       Figure11b (domestic producers difference) not in main.tex — omitted.
# =============================================================================
let
    β = Baseline_parameter.β[1]

    # --- open_CM_open_trade vs initial ---
    y_pop_open = exp.((V_trans_open_CM_open_trade .- V_saved_initial) .* (1.0 - β)) .- 1

    # Figure9a — domestic producers, open CM (labeled "Workers" in paper)
    y_mat = mat_creator(100.0 .* y_pop_open[7201:14400])
    contourf(y_mat', levels=5, c=cgrad(:Blues_6, categorical=true, rev=true),
            xlabel="Productivity", ylabel="Wealth", axis=nothing, lw=0,
            tickfontsize=18, xguidefontsize=18, yguidefontsize=18,
            legendfontsize=18)
    savefig("../Figures/Figure9a.pdf")
    println("[fig:welfare] Saved Figure9a.pdf")

    # Figure10a — exporters, open CM
    y_mat = mat_creator(100.0 .* y_pop_open[14401:end])
    contourf(y_mat', levels=5, c=cgrad(:Blues_6, categorical=true, rev=true),
            xlabel="Productivity", ylabel="Wealth", axis=nothing, lw=0,
            tickfontsize=18, xguidefontsize=18, yguidefontsize=18,
            legendfontsize=18)
    savefig("../Figures/Figure10a.pdf")
    println("[fig:welfare] Saved Figure10a.pdf")

    # --- closed_CM_open_trade vs initial ---
    y_pop_closed = exp.((V_trans_closed_CM_open_trade .- V_saved_initial) .* (1.0 - β)) .- 1

    # Figure9b — domestic producers, closed CM
    y_mat = mat_creator(100.0 .* y_pop_closed[7201:14400])
    contourf(y_mat', levels=5, c=cgrad(:Oranges_6, categorical=true, rev=true),
            xlabel="Productivity", ylabel="Wealth", axis=nothing, lw=0,
            tickfontsize=18, xguidefontsize=18, yguidefontsize=18,
            legendfontsize=18)
    savefig("../Figures/Figure9b.pdf")
    println("[fig:welfare] Saved Figure9b.pdf")

    # Figure10b — exporters, closed CM
    y_mat = mat_creator(100.0 .* y_pop_closed[14401:end])
    contourf(y_mat', levels=5, c=cgrad(:Oranges_6, categorical=true, rev=true),
            xlabel="Productivity", ylabel="Wealth", axis=nothing, lw=0,
            tickfontsize=18, xguidefontsize=18, yguidefontsize=18,
            legendfontsize=18)
    savefig("../Figures/Figure10b.pdf")
    println("[fig:welfare] Saved Figure10b.pdf")

    # --- open_CM vs closed_CM difference ---
    y_pop_diff = exp.((V_trans_open_CM_open_trade .- V_trans_closed_CM_open_trade) .* (1.0 - β)) .- 1

    # Figure11a — workers (1:7200), open-vs-closed difference
    y_mat = mat_creator(100.0 .* y_pop_diff[1:7200])
    contourf(y_mat', levels=5, c=cgrad(:Blues_6, categorical=true, rev=true),
            xlabel="Productivity", ylabel="Wealth", axis=nothing, lw=0,
            tickfontsize=18, xguidefontsize=18, yguidefontsize=18,
            legendfontsize=18)
    savefig("../Figures/Figure11a.pdf")
    println("[fig:welfare] Saved Figure11a.pdf")

    # Figure11c — exporters, open-vs-closed difference
    y_mat = mat_creator(100.0 .* y_pop_diff[14401:end])
    contourf(y_mat', levels=5, c=cgrad(:Blues_6, categorical=true, rev=true),
            xlabel="Productivity", ylabel="Wealth", axis=nothing, lw=0,
            colorbar_tickfontsize=40, colorbar_titlefontsize=40,
            guidefontsize=18, legendfontsize=18, titlefontsize=30,
            fontfamily="Computer Modern")
    savefig("../Figures/Figure11c.pdf")
    println("[fig:welfare] Saved Figure11c.pdf")
end

# =============================================================================
# [voting]  Majority voting on capital market reforms — Table_voting.tex
# Three scenarios, all from the household's perspective at the start of each SS:
#   (1) Full integration  (open CM + trade)  vs initial SS  — signal: V_trans_open_CM - V_init
#   (2) Trade only        (closed CM + trade) vs initial SS  — signal: V_trans_closed_CM - V_init
#   (3) CM integration    (given trade done)  vs trade-only  — signal: V_trans_CM_initial - V_trade
# Distributions: (1)+(2) use mu_init; (3) uses mu_0 (trade-only SS)
# For each scenario × occupation: Pop%, Frac.yes (=M^yes within occ.), Avg. CEV
# =============================================================================
let
    β      = Baseline_parameter.β[1]
    μ_init = current_distr_store_initial[:,1]
    μ₀     = current_distr_store_closed_CM_open_trade[:,1]

    # helper: returns (pops[3], frac_yes[3], avg_cev[3], pop_tot, frac_yes_tot, avg_cev_tot)
    function vstats(signal, μ)
        idxs   = [1:7200, 7201:14400, 14401:length(μ)]
        cev    = exp.(signal .* (1 - β)) .- 1
        pops   = [sum(μ[i])                         for i in idxs]
        my     = [sum((signal[i] .> 0) .* μ[i])     for i in idxs]
        ac     = [sum(cev[i] .* μ[i]) / pops[k]     for (k,i) in enumerate(idxs)]
        pt     = sum(pops)
        return pops, my ./ pops, ac, pt, sum(my)/pt, sum(cev .* μ)/pt
    end

    scenarios = [
        ("Full integration",  V_trans_open_CM_open_trade    .- V_saved_initial,                  μ_init),
        ("Trade only",        V_trans_closed_CM_open_trade  .- V_saved_initial,                  μ_init),
        ("CM integration",    V_trans_open_CM_open_trade_initial .- V_saved_closed_CM_open_trade, μ₀),
    ]

    # collect results: results[s] = (pops, frac_yes, avg_cev, pop_tot, fy_tot, ac_tot)
    results = [vstats(sig, μ) for (_, sig, μ) in scenarios]

    occ = ["Workers", "Domestic prod.", "Exporters", "\\hline\nTotal"]
    f1(x) = @sprintf("%.1f", round(round(100*x, digits=2, RoundNearestTiesAway), digits=1, RoundNearestTiesAway))

    # console summary
    for (k, (lbl, _, _)) in enumerate(scenarios)
        pops, fy, ac, pt, fy_t, ac_t = results[k]
        println("\n=== $lbl ===")
        println("Occupation       | Pop%  | Frac.yes | Avg.CEV")
        for (i, nm) in enumerate(["Workers","Domestic","Exporters"])
            @printf("  %-14s | %5.1f | %8.1f | %7.1f\n", nm, 100pops[i], 100fy[i], 100ac[i])
        end
        @printf("  %-14s | %5.1f | %8.1f | %7.1f\n", "Total", 100pt, 100fy_t, 100ac_t)
    end

    # build LaTeX data rows
    # columns: occ | pop_init | fy1 | ac1 | fy2 | ac2 | pop_mu0 | fy3 | ac3
    p_i = results[1][1];  p_0 = results[3][1]
    p_i_tot = results[1][4]; p_0_tot = results[3][4]

    body = join([
        "  $(["Workers","Domestic prod.","Exporters"][i]) & $(f1(p_i[i])) & $(f1(results[1][2][i])) & $(f1(results[1][3][i])) & $(f1(results[2][2][i])) & $(f1(results[2][3][i])) & $(f1(p_0[i])) & $(f1(results[3][2][i])) & $(f1(results[3][3][i])) \\\\"
        for i in 1:3
    ], "\n")
    tot = "  Total & $(f1(p_i_tot)) & $(f1(results[1][5])) & $(f1(results[1][6])) & $(f1(results[2][5])) & $(f1(results[2][6])) & $(f1(p_0_tot)) & $(f1(results[3][5])) & $(f1(results[3][6])) \\\\"

    tex = """\\begin{table}[!ht]
\\centering
\\caption{Majority voting on capital market reforms}
\\label{table:Table_voting}
\\scalebox{0.75}{
\\begin{tabular}{l c cc cc | c cc}
\\\\[-1.8ex]\\hline
\\hline \\\\[-1.8ex]
 & \\multicolumn{5}{c|}{From pre-reform (initial) SS} & \\multicolumn{3}{c}{From trade-only SS} \\\\
\\cmidrule(lr){2-6}\\cmidrule(lr){7-9}
 & Pop. & \\multicolumn{2}{c}{Full integration} & \\multicolumn{2}{c|}{Trade only} & Pop. & \\multicolumn{2}{c}{CM integration} \\\\
\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){8-9}
Occupation & (\\%) & \$M^{yes}\$ & Avg.\\ CEV & \$M^{yes}\$ & Avg.\\ CEV & (\\%) & \$M^{yes}\$ & Avg.\\ CEV \\\\
\\hline \\\\[-1.8ex]
$body
\\hline
$tot
\\hline
\\hline
\\end{tabular}
}
\\parbox{0.95\\textwidth}{\\caption*{\\scriptsize \\textit{Note:} \$M^{yes}\$ is the fraction of households within each occupation whose transition-path value function is strictly higher under the reform than under the status quo: \$M^{yes}_g = \\int_g \\mathbf{1}[V^{\\text{reform}}(a,z,s) > V^{\\text{status quo}}(a,z,s)]\\,d\\mu / \\int_g d\\mu\$. The status quo for columns 1--2 is the pre-reform (initial) steady state; for column 3 it is the trade-only (closed capital markets) steady state. Avg.\\ CEV is the population-weighted average consumption-equivalent welfare change within each occupation, in percent. Cross-check: aggregate Avg.\\ CEV in the Total row equals 13.1\\%, 4.6\\%, and 6.2\\% for columns 1--3 respectively.}}
\\end{table}
"""
    open("../Tables/Table_voting.tex", "w") do io
        write(io, tex)
    end
    println("\n  Saved → ../Tables/Table_voting.tex")
end

# =============================================================================
# [table:Table8]  Welfare change following external reforms  (inline → file)
# main.tex: \label{table:Table8}
# Source: plot.jl lines 888–906 (table_8_results 17×3)
# All welfare values are percentages (×100 of decimal CEV).
# Delayed-integration values appear in square brackets in column 3.
# TFP decomposition rows are hardcoded (post-processing results).
# =============================================================================
let
    r = zeros(17, 3)
    r[1,:]  = [1.0, GDP_closed_CM_open_trade[1]/GDP_initial[1], GDP_open_CM_open_trade[1]/GDP_initial[1]]
    r[2,:]  = [1.0, total_consumption_closed_CM_open_trade[1]/total_consumption_initial[1],
               total_consumption_open_CM_open_trade[1]/total_consumption_initial[1]]
    r[3,:]  = [1.0, capital_demand_closed_CM_open_trade[1]/capital_demand_initial[1],
               capital_demand_open_CM_open_trade[1]/capital_demand_initial[1]]
    r[4,:]  = [0.0, welfare_init_closed_CM_open_trade_stst, welfare_init_open_CM_open_trade_stst]
    r[5,:]  = [0.0, welfare_change_clCM_otrade_trans, welfare_change_oCM_otrade_trans]
    r[6,:]  = [0.0, welfare_init_closed_CM_open_trade_stst_F, welfare_init_open_CM_open_trade_stst_F]
    r[7,:]  = [0.0, welfare_change_clCM_otrade_trans_F, welfare_change_oCM_otrade_trans_F]
    r[8,:]  = [0.0, welfare_clCM_otrade_totalEU, welfare_oCM_otrade_totalEU]
    r[9,:]  = [p90_wealth_initial[1], p90_wealth_closed_CM_open_trade[1], p90_wealth_open_CM_open_trade[1]]
    r[10,:] = [p90_income_initial[1], p90_income_closed_CM_open_trade[1], p90_income_open_CM_open_trade[1]]
    r[11,:] = [p90_cons_initial[1], p90_cons_closed_CM_open_trade[1], p90_cons_open_CM_open_trade[1]]
    r[12,:] = [wealth_of_exporters_initial[1], wealth_of_exporters_closed_CM_open_trade[1], wealth_of_exporters_open_CM_open_trade[1]]
    r[13,:] = [wealth_of_workers_initial[1], wealth_of_workers_closed_CM_open_trade[1], wealth_of_workers_open_CM_open_trade[1]]
    r[14,:] = [cons_of_exporters_initial[1], cons_of_exporters_closed_CM_open_trade[1], cons_of_exporters_open_CM_open_trade[1]]
    r[15,:] = [cons_of_workers_initial[1], cons_of_workers_closed_CM_open_trade[1], cons_of_workers_open_CM_open_trade[1]]
    r[16,:] = [1.0,
               prices_closed_CM_open_trade[1]/prices_closed_CM_open_trade[5] / (prices_initial[1]/prices_initial[5]),
               prices_open_CM_open_trade[1]/prices_open_CM_open_trade[5] / (prices_initial[1]/prices_initial[5])]
    r[17,:] = [prices_initial[6]-prices_initial[7],
               prices_closed_CM_open_trade[6]-prices_closed_CM_open_trade[7],
               0.0]

    f0(x) = @sprintf("%.0f", round(round(round(x * 100, digits=2, RoundNearestTiesAway), digits=1, RoundNearestTiesAway), digits=0, RoundNearestTiesAway))
    f1(x) = @sprintf("%.1f", round(round(x * 100, digits=2, RoundNearestTiesAway), digits=1, RoundNearestTiesAway))
    f2(x) = @sprintf("%.2f", round(x * 100, digits=2, RoundNearestTiesAway))

    tex = """\\begin{table}[!t]
\\centering
\\caption{Welfare change following external reforms}
\\label{table:Table8}
  \\scalebox{0.80}{
\\begin{tabular}{l ccc ccc ccc}
\\\\[-1.8ex]\\hline
\\hline \\\\[-1.8ex]
Integration & None & Trade & Trade and capital \\\\
\\hline \\\\[-1.8ex]
\\textbf{Aggregates}& &   &   \\\\
Income & $(f0(r[1,1]))& $(f0(r[1,2]))& $(f0(r[1,3])) \\\\
Consumption & $(f0(r[2,1])) & $(f0(r[2,2])) & $(f1(r[2,3])) \\\\
Capital & $(f0(r[3,1]))& $(f0(r[3,2]))& $(f0(r[3,3])) \\\\[1ex]
\\textbf{Welfare change}& &   &   \\\\
Steady-state only  & 0 & $(f1(r[4,2])) & $(f1(r[4,3])) \\\\
Transition dynamics  & 0 & $(f1(r[5,2])) & $(f1(r[5,3]))[$(f1(welfare_change_oCM_otrade_delayed10_trans))] \\\\
Steady-state only*  & 0 & $(f1(r[6,2])) & $(f1(r[6,3])) \\\\
Transition dynamics*  & 0 & $(f1(r[7,2])) & $(f1(r[7,3]))[$(f1(welfare_change_oCM_otrade_delayed10_trans_F))] \\\\
EU welfare  & 0 & $(f1(r[8,2])) & $(f1(r[8,3]))[$(f1(welfare_delayed_totalEU))] \\\\[1ex]
\\textbf{Inequality}& &   &   \\\\
Top \$10\\%\$ wealth share &  $(f0(r[9,1])) & $(f0(r[9,2])) & $(f0(r[9,3])) \\\\
Top \$10\\%\$ income share & $(f0(r[10,1])) & $(f0(r[10,2])) & $(f0(r[10,3])) \\\\
Top \$10\\%\$ consumption share & $(f0(r[11,1])) & $(f0(r[11,2])) & $(f0(r[11,3])) \\\\
Consumption share of exporters & $(f0(r[14,1])) & $(f0(r[14,2])) & $(f0(r[14,3])) \\\\
Consumption share of workers & $(f0(r[15,1])) & $(f0(r[15,2])) & $(f0(r[15,3])) \\\\
Wealth owned by exporters & $(f0(r[12,1])) & $(f0(r[12,2])) & $(f0(r[12,3])) \\\\
Wealth owned by workers & $(f0(r[13,1])) & $(f0(r[13,2])) & $(f0(r[13,3])) \\\\[1ex]
\\textbf{Factor prices} & &   &\\\\
Real wage & $(f0(r[16,1])) & $(f0(r[16,2])) & $(f0(r[16,3])) \\\\
Interest rate premium \$r - r^*\$ & $(f0(r[17,1])) & $(f0(r[17,2])) & $(f0(r[17,3])) \\\\
\\hline
\\textbf{TFP loss decomposition} & &   &   \\\\
\\textit{Factors}& \$\\mathbf{-3}\$& \$\\mathbf{40}\$& \$\\mathbf{100}\$  \\\\
Within &  13 & 27 &  36\\\\
Between & -16 & 13 &  64 \\\\
External Reforms & 103 & 60 & 0  \\\\
\\hline
\\hline
\\textit{Total TFP loss \$\\%\$ }& \$\\mathbf{19}\$& \$\\mathbf{11}\$& \$\\mathbf{16}\$  \\\\
\\hline
\\end{tabular}
}
\\parbox{0.65\\textwidth}{\\caption*{\\scriptsize \\textit{Note:}  Compared to the pre-reform economy, different welfare measures change across different steady states. Utilitarian welfare change is the average consumption bundle required for individual households to achieve the same utility in the initial steady state. Asterisks denote Core moments, and square brackets indicate a ten-year delay in implementing capital markets. The EU-wide welfare change accounts for the transition path and refers to joint, population-weighted, NMS and Core welfare.}}
\\end{table}
"""
    open("../Tables/Table8.tex", "w") do io
        write(io, tex)
    end
    println("[table:Table8] Saved → ../Tables/Table8.tex")
end

# =============================================================================
# DEV ECONOMY WELFARE CALCULATIONS
# Needed for Table_combined (Table10/11) and Table4b (Core economy).
# Source: plot_dev.jl lines 1–47, 1055–1079; plot.jl lines 1145–1211
# Prerequisites: Step 5 (developed_steady_states.jl) + Step 2 (cases 6-7).
# =============================================================================

# Dev NMS value functions (from developed_steady_states.jl)
V_saved_initial_dev_NMS      = reshape(V_saved_fine_initial_dev[:,:,1],
                                        size(V_saved_fine_initial_dev)[1]*3)
V_saved_closed_CM_open_trade_dev_NMS = reshape(V_saved_fine_closed_CM_open_trade_dev[:,:,1],
                                                size(V_saved_fine_initial_dev)[1]*3)
V_saved_open_CM_open_trade_dev_NMS   = reshape(V_saved_fine_open_CM_open_trade_dev[:,:,1],
                                                size(V_saved_fine_initial_dev)[1]*3)

# SS welfare: dev scenarios vs dev initial (NMS)
welfare_init_closed_CM_open_trade_stst_dev =
    sum(current_distr_store_initial_dev[:,1] .*
        (exp.((V_saved_closed_CM_open_trade_dev_NMS - V_saved_initial_dev_NMS) *
              (1.0 - Developed_parameter.β[1])) .- 1))
welfare_init_open_CM_open_trade_stst_dev =
    sum(current_distr_store_initial_dev[:,1] .*
        (exp.((V_saved_open_CM_open_trade_dev_NMS - V_saved_initial_dev_NMS) *
              (1.0 - Developed_parameter.β[1])) .- 1))

# Transition welfare: dev transitions (from Step 2, cases 6-7)
V_trans_closed_CM_open_trade_dev_NMS =
    reshape(V_saved_store_closed_CM_open_trade_dev[:,:,1,2], size(V_saved_fine_initial_dev)[1]*3)
V_trans_open_CM_open_trade_dev_NMS =
    reshape(V_saved_store_open_CM_open_trade_dev[:,:,1,2],   size(V_saved_fine_initial_dev)[1]*3)
welfare_change_clCM_otrade_trans_dev =
    sum(current_distr_store_initial_dev[:,1] .*
        (exp.((V_trans_closed_CM_open_trade_dev_NMS - V_saved_initial_dev_NMS) *
              (1.0 - Developed_parameter.β[1])) .- 1))
welfare_change_oCM_otrade_trans_dev =
    sum(current_distr_store_initial_dev[:,1] .*
        (exp.((V_trans_open_CM_open_trade_dev_NMS - V_saved_initial_dev_NMS) *
              (1.0 - Developed_parameter.β[1])) .- 1))

# Cross-model welfare: developed initial vs baseline initial (uses both β values)
# V_saved_initial is already in scope from the shared welfare section
welfare_init_developed_stst =
    sum(current_distr_store_initial[:,1] .*
        (exp.(V_saved_initial_dev_NMS .* (1.0 - Developed_parameter.β[1]) .-
              V_saved_initial        .* (1.0 - Baseline_parameter.β[1])) .- 1))

# Dev FOREIGN (Core) value functions
V_saved_initial_dev_F_vec      = reshape(V_saved_fine_initial_dev[:,:,2],
                                          size(V_saved_fine_initial_dev)[1]*3)
V_saved_closed_CM_open_trade_dev_F_vec = reshape(V_saved_fine_closed_CM_open_trade_dev[:,:,2],
                                                   size(V_saved_fine_initial_dev)[1]*3)
V_saved_open_CM_open_trade_dev_F_vec   = reshape(V_saved_fine_open_CM_open_trade_dev[:,:,2],
                                                   size(V_saved_fine_initial_dev)[1]*3)

# V_saved_initial_F is already in scope from the shared welfare section
welfare_init_developed_stst_F =
    sum(current_distr_store_initial[:,2] ./ Developed_parameter.L[2] .*
        (exp.(V_saved_initial_dev_F_vec .* (1.0 - Developed_parameter.β[2]) .-
              V_saved_initial_F         .* (1.0 - Baseline_parameter.β[2])) .- 1))

# Dev foreign SS welfare
welfare_init_closed_CM_open_trade_stst_dev_F =
    sum(current_distr_store_initial_dev[:,2] ./ Baseline_parameter.L[2] .*
        (exp.((V_saved_closed_CM_open_trade_dev_F_vec - V_saved_initial_dev_F_vec) *
              (1.0 - Developed_parameter.β[1])) .- 1))
welfare_init_open_CM_open_trade_stst_dev_F =
    sum(current_distr_store_initial_dev[:,2] ./ Baseline_parameter.L[2] .*
        (exp.((V_saved_open_CM_open_trade_dev_F_vec - V_saved_initial_dev_F_vec) *
              (1.0 - Developed_parameter.β[1])) .- 1))

# Dev foreign transition welfare (from Step 2 cases 6-7)
V_trans_closed_CM_open_trade_dev_F_vec =
    reshape(V_saved_store_closed_CM_open_trade_dev[:,:,2,2], size(V_saved_fine_initial_dev)[1]*3)
V_trans_open_CM_open_trade_dev_F_vec =
    reshape(V_saved_store_open_CM_open_trade_dev[:,:,2,2],   size(V_saved_fine_initial_dev)[1]*3)
welfare_change_clCM_otrade_trans_dev_F =
    sum(current_distr_store_initial_dev[:,2] ./ Developed_parameter.L[2] .*
        (exp.((V_trans_closed_CM_open_trade_dev_F_vec - V_saved_initial_dev_F_vec) *
              (1.0 - Developed_parameter.β[2])) .- 1))
welfare_change_oCM_otrade_trans_dev_F =
    sum(current_distr_store_initial_dev[:,2] ./ Developed_parameter.L[2] .*
        (exp.((V_trans_open_CM_open_trade_dev_F_vec - V_saved_initial_dev_F_vec) *
              (1.0 - Developed_parameter.β[2])) .- 1))

# Dev credit variables: NMS (plot_dev.jl lines 44-47, 120-133)
total_credit_initial_dev =
    -(domestic_firm_debt_initial_dev[1] + exporter_firm_debt_initial_dev[1]) /
    nomGDP_initial_dev[1]
domestic_credit_initial_dev =
    (worker_bond_holding_initial_dev[1] + domestic_bond_holding_initial_dev[1] +
     exporter_bond_holding_initial_dev[1]) ./ nomGDP_initial_dev[1]
foreign_credit_initial_dev       = total_credit_initial_dev - domestic_credit_initial_dev
foreign_credit_share_initial_dev = foreign_credit_initial_dev / total_credit_initial_dev

total_credit_closed_CM_open_trade_dev =
    -(domestic_firm_debt_closed_CM_open_trade_dev[1] +
      exporter_firm_debt_closed_CM_open_trade_dev[1]) / nomGDP_closed_CM_open_trade_dev[1]
domestic_credit_closed_CM_open_trade_dev =
    (worker_bond_holding_closed_CM_open_trade_dev[1] +
     domestic_bond_holding_closed_CM_open_trade_dev[1] +
     exporter_bond_holding_closed_CM_open_trade_dev[1]) ./ nomGDP_closed_CM_open_trade_dev[1]
foreign_credit_closed_CM_open_trade_dev =
    total_credit_closed_CM_open_trade_dev - domestic_credit_closed_CM_open_trade_dev
foreign_credit_share_closed_CM_open_trade_dev =
    foreign_credit_closed_CM_open_trade_dev / total_credit_closed_CM_open_trade_dev

total_credit_open_CM_open_trade_dev =
    -(domestic_firm_debt_open_CM_open_trade_dev[1] +
      exporter_firm_debt_open_CM_open_trade_dev[1]) / nomGDP_open_CM_open_trade_dev[1]
domestic_credit_open_CM_open_trade_dev =
    (worker_bond_holding_open_CM_open_trade_dev[1] +
     domestic_bond_holding_open_CM_open_trade_dev[1] +
     exporter_bond_holding_open_CM_open_trade_dev[1]) ./ nomGDP_open_CM_open_trade_dev[1]
foreign_credit_open_CM_open_trade_dev =
    total_credit_open_CM_open_trade_dev - domestic_credit_open_CM_open_trade_dev
foreign_credit_share_open_CM_open_trade_dev =
    foreign_credit_open_CM_open_trade_dev / total_credit_open_CM_open_trade_dev

# Core credit variables (plot_dev.jl lines 296-329)
total_credit_initial_F =
    -(domestic_firm_debt_initial[2] + exporter_firm_debt_initial[2]) / nomGDP_initial[2]
domestic_credit_initial_F =
    (worker_bond_holding_initial[2] + domestic_bond_holding_initial[2] +
     exporter_bond_holding_initial[2]) ./ nomGDP_initial[2]
foreign_credit_initial_F       = total_credit_initial_F - domestic_credit_initial_F
foreign_credit_share_initial_F = foreign_credit_initial_F / total_credit_initial_F

total_credit_closed_CM_open_trade_F =
    -(domestic_firm_debt_closed_CM_open_trade[2] +
      exporter_firm_debt_closed_CM_open_trade[2]) / nomGDP_closed_CM_open_trade[2]
domestic_credit_closed_CM_open_trade_F =
    (worker_bond_holding_closed_CM_open_trade[2] +
     domestic_bond_holding_closed_CM_open_trade[2] +
     exporter_bond_holding_closed_CM_open_trade[2]) ./ nomGDP_closed_CM_open_trade[2]
foreign_credit_closed_CM_open_trade_F =
    total_credit_closed_CM_open_trade_F - domestic_credit_closed_CM_open_trade_F
foreign_credit_share_closed_CM_open_trade_F =
    foreign_credit_closed_CM_open_trade_F / total_credit_closed_CM_open_trade_F

total_credit_open_CM_open_trade_F =
    -(domestic_firm_debt_open_CM_open_trade[2] +
      exporter_firm_debt_open_CM_open_trade[2]) / nomGDP_open_CM_open_trade[2]
domestic_credit_open_CM_open_trade_F =
    (worker_bond_holding_open_CM_open_trade[2] +
     domestic_bond_holding_open_CM_open_trade[2] +
     exporter_bond_holding_open_CM_open_trade[2]) ./ nomGDP_open_CM_open_trade[2]
foreign_credit_open_CM_open_trade_F =
    total_credit_open_CM_open_trade_F - domestic_credit_open_CM_open_trade_F
foreign_credit_share_open_CM_open_trade_F =
    foreign_credit_open_CM_open_trade_F / total_credit_open_CM_open_trade_F

total_credit_open_CM_closed_trade_F =
    -(domestic_firm_debt_open_CM_closed_trade[2] +
      exporter_firm_debt_open_CM_closed_trade[2]) / nomGDP_open_CM_closed_trade[2]
domestic_credit_open_CM_closed_trade_F =
    (worker_bond_holding_open_CM_closed_trade[2] +
     domestic_bond_holding_open_CM_closed_trade[2] +
     exporter_bond_holding_open_CM_closed_trade[2]) ./ nomGDP_open_CM_closed_trade[2]
foreign_credit_open_CM_closed_trade_F =
    total_credit_open_CM_closed_trade_F - domestic_credit_open_CM_closed_trade_F
foreign_credit_share_open_CM_closed_trade_F =
    foreign_credit_open_CM_closed_trade_F / total_credit_open_CM_closed_trade_F

total_credit_initial_dev_F =
    -(domestic_firm_debt_initial_dev[2] + exporter_firm_debt_initial_dev[2]) /
    nomGDP_initial_dev[2]
domestic_credit_initial_dev_F =
    (worker_bond_holding_initial_dev[2] + domestic_bond_holding_initial_dev[2] +
     exporter_bond_holding_initial_dev[2]) ./ nomGDP_initial_dev[2]
foreign_credit_initial_dev_F       = total_credit_initial_dev_F - domestic_credit_initial_dev_F
foreign_credit_share_initial_dev_F = foreign_credit_initial_dev_F / total_credit_initial_dev_F

total_credit_closed_CM_open_trade_dev_F =
    -(domestic_firm_debt_closed_CM_open_trade_dev[2] +
      exporter_firm_debt_closed_CM_open_trade_dev[2]) / nomGDP_closed_CM_open_trade_dev[2]
domestic_credit_closed_CM_open_trade_dev_F =
    (worker_bond_holding_closed_CM_open_trade_dev[2] +
     domestic_bond_holding_closed_CM_open_trade_dev[2] +
     exporter_bond_holding_closed_CM_open_trade_dev[2]) ./ nomGDP_closed_CM_open_trade_dev[2]
foreign_credit_closed_CM_open_trade_dev_F =
    total_credit_closed_CM_open_trade_dev_F - domestic_credit_closed_CM_open_trade_dev_F
foreign_credit_share_closed_CM_open_trade_dev_F =
    foreign_credit_closed_CM_open_trade_dev_F / total_credit_closed_CM_open_trade_dev_F

total_credit_open_CM_open_trade_dev_F =
    -(domestic_firm_debt_open_CM_open_trade_dev[2] +
      exporter_firm_debt_open_CM_open_trade_dev[2]) / nomGDP_open_CM_open_trade_dev[2]
domestic_credit_open_CM_open_trade_dev_F =
    (worker_bond_holding_open_CM_open_trade_dev[2] +
     domestic_bond_holding_open_CM_open_trade_dev[2] +
     exporter_bond_holding_open_CM_open_trade_dev[2]) ./ nomGDP_open_CM_open_trade_dev[2]
foreign_credit_open_CM_open_trade_dev_F =
    total_credit_open_CM_open_trade_dev_F - domestic_credit_open_CM_open_trade_dev_F
foreign_credit_share_open_CM_open_trade_dev_F =
    foreign_credit_open_CM_open_trade_dev_F / total_credit_open_CM_open_trade_dev_F

# open_CM_closed_trade credit for NMS (Table9 / plot.jl lines 1135-1138)
total_credit_open_CM_closed_trade =
    -(domestic_firm_debt_open_CM_closed_trade[1] +
      exporter_firm_debt_open_CM_closed_trade[1]) / nomGDP_open_CM_closed_trade[1]
domestic_credit_open_CM_closed_trade =
    (worker_bond_holding_open_CM_closed_trade[1] +
     domestic_bond_holding_open_CM_closed_trade[1] +
     exporter_bond_holding_open_CM_closed_trade[1]) ./ nomGDP_open_CM_closed_trade[1]
foreign_credit_open_CM_closed_trade =
    total_credit_open_CM_closed_trade - domestic_credit_open_CM_closed_trade
foreign_credit_share_open_CM_closed_trade =
    foreign_credit_open_CM_closed_trade / total_credit_open_CM_closed_trade

println("Dev welfare and credit calculations done.")

# =============================================================================
# [table:Table9]  Standalone capital market integration  (computed)
# main.tex: \label{table:Table9}
# Source: plot.jl lines 1119–1141
# Columns: None (initial), Capital (open_CM_closed_trade), Trade+cap (open_CM_open_trade)
# =============================================================================
let
    f0(x) = @sprintf("%.0f", round(round(round(x * 100, digits=2, RoundNearestTiesAway), digits=1, RoundNearestTiesAway), digits=0, RoundNearestTiesAway))
    f1(x) = @sprintf("%.1f", round(round(x * 100, digits=2, RoundNearestTiesAway), digits=1, RoundNearestTiesAway))
    f2(x) = @sprintf("%.2f", round(x, digits=2, RoundNearestTiesAway))        # sd_MRPK shown as-is
    fip(x) = @sprintf("%.0f", round(round(round(x * 100, digits=2, RoundNearestTiesAway), digits=1, RoundNearestTiesAway), digits=0, RoundNearestTiesAway)) # interest rate premium

    r = zeros(20, 3)
    r[1,:] = [1.0,
              TFP_open_CM_closed_trade[1] / TFP_initial[1],
              TFP_open_CM_open_trade[1]   / TFP_initial[1]]
    r[2,:] = [sd_MRPK_initial[1], sd_MRPK_open_CM_closed_trade[1], sd_MRPK_open_CM_open_trade[1]]
    r[3,:] = [1.0,
              GDP_open_CM_closed_trade[1] / GDP_initial[1],
              GDP_open_CM_open_trade[1]   / GDP_initial[1]]
    r[4,:] = [1.0,
              prices_open_CM_closed_trade[3] / prices_initial[3],
              prices_open_CM_open_trade[3]   / prices_initial[3]]
    r[5,:] = [1.0,
              total_consumption_open_CM_closed_trade[1] / total_consumption_initial[1],
              total_consumption_open_CM_open_trade[1]   / total_consumption_initial[1]]
    r[6,:] = [1.0,
              capital_demand_open_CM_closed_trade[1] / capital_demand_initial[1],
              capital_demand_open_CM_open_trade[1]   / capital_demand_initial[1]]
    r[8,:] = [0.0, welfare_change_oCM_cltrade_trans, welfare_change_oCM_otrade_trans]
    r[9,:] = [p90_wealth_initial[1], p90_wealth_open_CM_closed_trade[1],
              p90_wealth_open_CM_open_trade[1]]
    r[10,:] = [1.0,
               prices_open_CM_closed_trade[1]/prices_open_CM_closed_trade[5] /
                   (prices_initial[1]/prices_initial[5]),
               prices_open_CM_open_trade[1]/prices_open_CM_open_trade[5] /
                   (prices_initial[1]/prices_initial[5])]
    # Interest rate premium: r - r*.  CM-open → premium = 0.
    r[11,:] = [prices_initial[6] - prices_initial[7], 0.0, 0.0]
    r[12,:] = [import_share_initial[1], import_share_open_CM_closed_trade[1],
               import_share_open_CM_open_trade[1]]
    r[13,:] = [export_value_initial[1]       / nomGDP_initial[2],
               export_value_open_CM_closed_trade[1] / nomGDP_open_CM_closed_trade[2],
               export_value_open_CM_open_trade[1]   / nomGDP_open_CM_open_trade[2]]
    r[14,:] = [domestic_pop_initial[1] + exporter_pop_initial[1],
               domestic_pop_open_CM_closed_trade[1] + exporter_pop_open_CM_closed_trade[1],
               domestic_pop_open_CM_open_trade[1]   + exporter_pop_open_CM_open_trade[1]]
    r[15,:] = [exporter_pop_initial[1] / (domestic_pop_initial[1] + exporter_pop_initial[1]),
               exporter_pop_open_CM_closed_trade[1] /
                   (domestic_pop_open_CM_closed_trade[1] + exporter_pop_open_CM_closed_trade[1]),
               exporter_pop_open_CM_open_trade[1] /
                   (domestic_pop_open_CM_open_trade[1] + exporter_pop_open_CM_open_trade[1])]
    r[16,:] = [prices_initial[5], prices_open_CM_closed_trade[5], prices_open_CM_open_trade[5]]
    r[17,:] = [total_credit_initial, total_credit_open_CM_closed_trade,
               total_credit_open_CM_open_trade]
    r[18,:] = [foreign_credit_share_initial, foreign_credit_share_open_CM_closed_trade,
               foreign_credit_share_open_CM_open_trade]

    println("[table:Table9] table_9_results computed.")
    println(round.(100 .* r, digits=1))

    tex = """\\begin{table}[!t]
\\centering
\\caption{Standalone capital market integration}
\\label{table:Table9}
  \\scalebox{0.80}{
\\begin{tabular}{l ccc ccc ccc}
\\\\[-1.8ex]\\hline
\\hline \\\\[-1.8ex]
Integration & None & Capital & Trade and capital \\\\
\\hline \\\\[-1.8ex]
\\textbf{Productivity}& &   &   \\\\
TFP & $(f0(r[1,1]))  & $(f0(r[1,2]))& $(f0(r[1,3])) \\\\
s.d. \$mrpk\$ & $(f2(r[2,1]))  & $(f2(r[2,2]))& $(f2(r[2,3])) \\\\[1ex]
\\textbf{Aggregates}& &   &   \\\\
Output & $(f0(r[3,1]))& $(f0(r[3,2]))& $(f0(r[3,3])) \\\\
Income & $(f0(r[4,1]))& $(f0(r[4,2]))& $(f0(r[4,3])) \\\\
Consumption & $(f0(r[5,1])) & $(f0(r[5,2])) & $(f1(r[5,3])) \\\\
Capital & $(f0(r[6,1]))& $(f0(r[6,2]))& $(f0(r[6,3])) \\\\[1ex]
\\textbf{Welfare and Inequality}& &   &   \\\\
Transition dynamics  & 0 & $(f0(r[8,2])) & $(f0(r[8,3])) \\\\
Top \$10\\%\$ wealth share &  $(f0(r[9,1])) & $(f0(r[9,2])) & $(f0(r[9,3])) \\\\[1ex]
\\textbf{Factor prices} & &   &\\\\
Real wage & $(f0(r[10,1])) & $(f0(r[10,2])) & $(f0(r[10,3]))\\\\
Interest rate premium \$r - r^*\$ & $(fip(r[11,1])) & 0 & 0 \\\\[1ex]
\\textbf{Trade}& &   &\\\\
\$\\frac{\\text{Import}}{\\text{GDP}}\$ & $(f0(r[12,1])) & $(f0(r[12,2])) & $(f0(r[12,3])) \\\\
\$\\frac{\\text{Export}}{\\text{GDP}^*}\$ & $(f0(r[13,1])) & $(f0(r[13,2])) & $(f0(r[13,3])) \\\\
Entrepreneurship rate & $(f0(r[14,1])) & $(f0(r[14,2])) & $(f0(r[14,3])) \\\\
Share of exporters & $(f0(r[15,1])) & $(f0(r[15,2])) & $(f0(r[15,3])) \\\\
CPI & $(f0(r[16,1])) & $(f0(r[16,2])) & $(f0(r[16,3])) \\\\
\$\\frac{\\text{Credit}}{\\text{GDP}}\$ & $(f0(r[17,1])) & $(f0(r[17,2])) & $(f0(r[17,3])) \\\\
\$\\frac{\\text{Foreign Credit}}{\\text{Credit}}\$ & $(f0(r[18,1])) & $(f0(r[18,2])) & $(f0(r[18,3]))\\\\
\\hline
\\end{tabular}}
\\end{table}
"""
    open("../Tables/Table9.tex", "w") do io
        write(io, tex)
    end
    println("[table:Table9] Saved → ../Tables/Table9.tex")
end

# =============================================================================
# [table:combined]  Financial development and external reforms  (computed)
# main.tex: \label{table:combined}
# Source: plot.jl lines 1145–1211; plot_dev.jl tables
# Table10: NMS Low (initial) vs High (initial_dev) financial development
# Table11: Reforms in the high-financial-development economy (3 scenarios)
# =============================================================================
let
    f0(x)  = @sprintf("%.0f", round(round(round(x * 100, digits=2, RoundNearestTiesAway), digits=1, RoundNearestTiesAway), digits=0, RoundNearestTiesAway))
    f1(x)  = @sprintf("%.1f", round(round(x * 100, digits=2, RoundNearestTiesAway), digits=1, RoundNearestTiesAway))
    f2(x)  = @sprintf("%.2f", round(x, digits=2, RoundNearestTiesAway))        # sd_MRPK as-is
    fip(x) = @sprintf("%.0f", round(round(round(x * 100, digits=2, RoundNearestTiesAway), digits=1, RoundNearestTiesAway), digits=0, RoundNearestTiesAway))  # interest rate

    # ── Table10 ────────────────────────────────────────────────────────────────
    t10 = zeros(20, 2)
    t10[1,:]  = [1.0, TFP_initial_dev[1] / TFP_initial[1]]
    t10[2,:]  = [sd_MRPK_initial[1], sd_MRPK_initial_dev[1]]
    t10[3,:]  = [1.0, GDP_initial_dev[1] / GDP_initial[1]]
    t10[4,:]  = [1.0, prices_initial_dev[3] / prices_initial[3]]
    t10[5,:]  = [1.0, total_consumption_initial_dev[1] / total_consumption_initial[1]]
    t10[6,:]  = [1.0, capital_demand_initial_dev[1] / capital_demand_initial[1]]
    t10[7,:]  = [0.0, welfare_init_developed_stst]
    t10[8,:]  = [p90_wealth_initial[1], p90_wealth_initial_dev[1]]
    t10[9,:]  = [p90_income_initial[1], p90_income_initial_dev[1]]
    t10[10,:] = [p90_cons_initial[1],   p90_cons_initial_dev[1]]
    t10[11,:] = [1.0,
                 prices_initial_dev[1]/prices_initial_dev[5] /
                     (prices_initial[1]/prices_initial[5])]
    t10[12,:] = [prices_initial[6] - prices_initial[7],
                 prices_initial_dev[6] - prices_initial_dev[7]]
    t10[13,:] = [import_share_initial[1], import_share_initial_dev[1]]
    t10[14,:] = [export_value_initial[1] / nomGDP_initial[2],
                 export_value_initial_dev[1] / nomGDP_initial_dev[2]]
    t10[15,:] = [domestic_pop_initial[1] + exporter_pop_initial[1],
                 domestic_pop_initial_dev[1] + exporter_pop_initial_dev[1]]
    t10[16,:] = [exporter_pop_initial[1] /
                     (domestic_pop_initial[1] + exporter_pop_initial[1]),
                 exporter_pop_initial_dev[1] /
                     (domestic_pop_initial_dev[1] + exporter_pop_initial_dev[1])]
    t10[17,:] = [prices_initial[5], prices_initial_dev[5]]
    t10[18,:] = [total_credit_initial, total_credit_initial_dev]
    t10[19,:] = [foreign_credit_share_initial, foreign_credit_share_initial_dev]

    # ── Table11 ────────────────────────────────────────────────────────────────
    t11 = zeros(21, 3)
    t11[1,:]  = [1.0,
                 TFP_closed_CM_open_trade_dev[1] / TFP_initial_dev[1],
                 TFP_open_CM_open_trade_dev[1]   / TFP_initial_dev[1]]
    t11[2,:]  = [sd_MRPK_initial_dev[1],
                 sd_MRPK_closed_CM_open_trade_dev[1],
                 sd_MRPK_open_CM_open_trade_dev[1]]
    t11[3,:]  = [1.0,
                 GDP_closed_CM_open_trade_dev[1] / GDP_initial_dev[1],
                 GDP_open_CM_open_trade_dev[1]   / GDP_initial_dev[1]]
    t11[4,:]  = [1.0,
                 prices_closed_CM_open_trade_dev[3] / prices_initial_dev[3],
                 prices_open_CM_open_trade_dev[3]   / prices_initial_dev[3]]
    t11[5,:]  = [1.0,
                 total_consumption_closed_CM_open_trade_dev[1] / total_consumption_initial_dev[1],
                 total_consumption_open_CM_open_trade_dev[1]   / total_consumption_initial_dev[1]]
    t11[6,:]  = [1.0,
                 capital_demand_closed_CM_open_trade_dev[1] / capital_demand_initial_dev[1],
                 capital_demand_open_CM_open_trade_dev[1]   / capital_demand_initial_dev[1]]
    t11[7,:]  = [0.0, welfare_init_closed_CM_open_trade_stst_dev,
                 welfare_init_open_CM_open_trade_stst_dev]
    t11[8,:]  = [0.0, welfare_change_clCM_otrade_trans_dev,
                 welfare_change_oCM_otrade_trans_dev]
    t11[9,:]  = [p90_wealth_initial_dev[1], p90_wealth_closed_CM_open_trade_dev[1],
                 p90_wealth_open_CM_open_trade_dev[1]]
    t11[10,:] = [p90_income_initial_dev[1], p90_income_closed_CM_open_trade_dev[1],
                 p90_income_open_CM_open_trade_dev[1]]
    t11[11,:] = [p90_cons_initial_dev[1], p90_cons_closed_CM_open_trade_dev[1],
                 p90_cons_open_CM_open_trade_dev[1]]
    t11[12,:] = [1.0,
                 prices_closed_CM_open_trade_dev[1]/prices_closed_CM_open_trade_dev[5] /
                     (prices_initial_dev[1]/prices_initial_dev[5]),
                 prices_open_CM_open_trade_dev[1]/prices_open_CM_open_trade_dev[5] /
                     (prices_initial_dev[1]/prices_initial_dev[5])]
    t11[14,:] = [prices_initial_dev[6] - prices_initial_dev[7],
                 prices_closed_CM_open_trade_dev[6] - prices_closed_CM_open_trade_dev[7],
                 prices_open_CM_open_trade_dev[6]   - prices_open_CM_open_trade_dev[6]]
    t11[15,:] = [import_share_initial_dev[1], import_share_closed_CM_open_trade_dev[1],
                 import_share_open_CM_open_trade_dev[1]]
    t11[16,:] = [export_value_initial_dev[1] / nomGDP_initial_dev[2],
                 export_value_closed_CM_open_trade_dev[1] / nomGDP_closed_CM_open_trade_dev[2],
                 export_value_open_CM_open_trade_dev[1]   / nomGDP_open_CM_open_trade_dev[2]]
    t11[17,:] = [domestic_pop_initial_dev[1] + exporter_pop_initial_dev[1],
                 domestic_pop_closed_CM_open_trade_dev[1] + exporter_pop_closed_CM_open_trade_dev[1],
                 domestic_pop_open_CM_open_trade_dev[1]   + exporter_pop_open_CM_open_trade_dev[1]]
    t11[18,:] = [exporter_pop_initial_dev[1] /
                     (domestic_pop_initial_dev[1] + exporter_pop_initial_dev[1]),
                 exporter_pop_closed_CM_open_trade_dev[1] /
                     (domestic_pop_closed_CM_open_trade_dev[1] + exporter_pop_closed_CM_open_trade_dev[1]),
                 exporter_pop_open_CM_open_trade_dev[1] /
                     (domestic_pop_open_CM_open_trade_dev[1] + exporter_pop_open_CM_open_trade_dev[1])]
    t11[19,:] = [prices_initial_dev[5], prices_closed_CM_open_trade_dev[5],
                 prices_open_CM_open_trade_dev[5]]
    t11[20,:] = [total_credit_initial_dev, total_credit_closed_CM_open_trade_dev,
                 total_credit_open_CM_open_trade_dev]
    t11[21,:] = [foreign_credit_share_initial_dev, foreign_credit_share_closed_CM_open_trade_dev,
                 foreign_credit_share_open_CM_open_trade_dev]

    println("[table:combined] table_10_results and table_11_results computed.")

    tex = """\\begin{table}[!t]
\\centering
\\caption{Financial development and external reforms}
\\label{table:combined}
\\begin{minipage}[t]{0.48\\textwidth}
\\centering
\\caption*{Effects of financial development}
\\label{table:Table10}
\\vspace{0pt}
\\scalebox{0.70}{
\\begin{tabular}{l cc}
\\hline
\\hline
Financial development & Low  & High \\\\
\\hline
\\textbf{Productivity}& &    \\\\
TFP & $(f0(t10[1,1]))  & $(f0(t10[1,2])) \\\\
Standard deviation of ARPK & $(f2(t10[2,1]))  & $(f2(t10[2,2])) \\\\[1ex]
\\textbf{Aggregates}& &     \\\\
Output & $(f0(t10[3,1]))& $(f0(t10[3,2])) \\\\
Income & $(f0(t10[4,1]))& $(f0(t10[4,2])) \\\\
Consumption & $(f0(t10[5,1])) & $(f0(t10[5,2])) \\\\
Capital & $(f0(t10[6,1]))& $(f0(t10[6,2])) \\\\[1ex]
\\textbf{Welfare and Inequality}& &    \\\\
Steady state  & - & $(f0(t10[7,2])) \\\\
Top \$10\\%\$ wealth share &  $(f0(t10[8,1])) & $(f0(t10[8,2])) \\\\
Top \$10\\%\$ income share &  $(f0(t10[9,1])) & $(f0(t10[9,2])) \\\\
Top \$10\\%\$ consumption share &  $(f0(t10[10,1])) & $(f0(t10[10,2])) \\\\[1ex]
\\textbf{Factor prices} & &  \\\\
Real wage & $(f0(t10[11,1])) & $(f0(t10[11,2]))\\\\
\$r - r^*\$ & $(fip(t10[12,1])) & $(fip(t10[12,2])) \\\\[1ex]
\\textbf{Trade}& &  \\\\
\$\\frac{\\text{Import}}{\\text{GDP}}\$ & $(f0(t10[13,1])) & $(f0(t10[13,2])) \\\\
\$\\frac{\\text{Export}}{\\text{GDP}^*}\$ & $(f0(t10[14,1])) &  $(f0(t10[14,2])) \\\\
Entrepreneurship rate & $(f0(t10[15,1])) & $(f0(t10[15,2])) \\\\
Share of exporters & $(f0(t10[16,1])) & $(f0(t10[16,2])) \\\\
RER & $(f0(t10[17,1])) & $(f0(t10[17,2])) \\\\
\$\\frac{\\text{Credit}}{\\text{GDP}}\$& $(f0(t10[18,1])) & $(f0(t10[18,2])) \\\\
\$\\frac{\\text{Foreign Credit}}{\\text{Credit}}\$ & 0 & 0 \\\\
\\hline
\\end{tabular}}
\\end{minipage}%
\\hfill
\\begin{minipage}[t]{0.48\\textwidth}
\\centering
\\caption*{Reforms \\& high financial development}
\\label{table:Table11}
\\vspace{0pt}
\\scalebox{0.70}{
\\begin{tabular}{l ccc}
\\hline
\\hline
Integration & None  & Trade & Trade + cap. \\\\
\\hline
\\textbf{Productivity}& &   &   \\\\
TFP & $(f0(t11[1,1]))  & $(f0(t11[1,2]))& $(f0(t11[1,3])) \\\\
Std. dev. of ARPK & $(f2(t11[2,1]))  & $(f2(t11[2,2]))& $(f2(t11[2,3])) \\\\[1ex]
\\textbf{Aggregates}& &   &   \\\\
Output & $(f0(t11[3,1]))& $(f0(t11[3,2]))& $(f0(t11[3,3])) \\\\
Income & $(f0(t11[4,1]))& $(f0(t11[4,2]))& $(f0(t11[4,3])) \\\\
Consumption & $(f0(t11[5,1])) & $(f0(t11[5,2]))& $(f0(t11[5,3])) \\\\
Capital & $(f0(t11[6,1]))& $(f0(t11[6,2]))& $(f0(t11[6,3])) \\\\[1ex]
\\textbf{Welfare change}& &   &   \\\\
Steady state only  & 0 & $(f1(t11[7,2])) & $(f1(t11[7,3])) \\\\
Transition dynamics  & 0 & $(f1(t11[8,2])) & $(f1(t11[8,3])) \\\\[1ex]
\\textbf{Inequality}& &   &   \\\\
Top \$10\\%\$ wealth share &  $(f0(t11[9,1])) & $(f0(t11[9,2])) & $(f0(t11[9,3])) \\\\
Top \$10\\%\$ income share & $(f0(t11[10,1])) & $(f0(t11[10,2])) & $(f0(t11[10,3])) \\\\
Top \$10\\%\$ cons. share & $(f0(t11[11,1])) & $(f0(t11[11,2])) & $(f0(t11[11,3])) \\\\[1ex]
\\textbf{Factor prices} & &   &\\\\
Real wage & $(f0(t11[12,1])) & $(f0(t11[12,2])) & $(f0(t11[12,3])) \\\\
\$r - r^*\$ & $(fip(t11[14,1])) & $(fip(t11[14,2])) & $(fip(t11[14,3])) \\\\[1ex]
\\textbf{Trade}& &   &\\\\
\$\\frac{\\text{Import}}{\\text{GDP}}\$ & $(f0(t11[15,1])) & $(f0(t11[15,2])) & $(f0(t11[15,3])) \\\\
\$\\frac{\\text{Export}}{\\text{GDP}^*}\$ & $(f0(t11[16,1])) & $(f0(t11[16,2])) & $(f0(t11[16,3])) \\\\
Entrepreneurship rate & $(f0(t11[17,1])) & $(f0(t11[17,2])) & $(f0(t11[17,3])) \\\\
Share of exporters & $(f0(t11[18,1])) & $(f0(t11[18,2])) & $(f0(t11[18,3])) \\\\
RER & $(f0(t11[19,1])) & $(f0(t11[19,2])) & $(f0(t11[19,3])) \\\\
\$\\frac{\\text{Credit}}{\\text{GDP}}\$ & $(f0(t11[20,1])) & $(f0(t11[20,2])) & $(f0(t11[20,3])) \\\\
\$\\frac{\\text{Foreign Credit}}{\\text{Credit}}\$ & 0 & 0 & $(f0(t11[21,3]))\\\\
\\hline
\\end{tabular} }
\\end{minipage}
\\end{table}
"""
    open("../Tables/Table_combined.tex", "w") do io
        write(io, tex)
    end
    println("[table:combined] Saved → ../Tables/Table_combined.tex")
end

# =============================================================================
# [table:Table4b]  Effect of reforms on the Core economy  (computed)
# main.tex: \label{table:Table4b}
# Source: plot_dev.jl lines 234–345  (table_12_results, 26×8)
# Columns: None, Trade, Trade+cap, Capital, Development, Trade(dev), Trade+cap(dev)
# All refer to the Core (foreign, country 2) economy.
# =============================================================================
let
    f0(x)  = @sprintf("%.0f", round(round(round(x * 100, digits=2, RoundNearestTiesAway), digits=1, RoundNearestTiesAway), digits=0, RoundNearestTiesAway))
    f1(x)  = @sprintf("%.1f", round(round(x * 100, digits=2, RoundNearestTiesAway), digits=1, RoundNearestTiesAway))
    f2(x)  = @sprintf("%.2f", round(x, digits=2, RoundNearestTiesAway))         # sd_MRPK as-is
    fip(x) = @sprintf("%.1f", round(round(x * 100, digits=2, RoundNearestTiesAway), digits=1, RoundNearestTiesAway))   # interest rate r*
    fimp(x) = @sprintf("%.1f", round(round(x * 100, digits=2, RoundNearestTiesAway), digits=1, RoundNearestTiesAway))  # import share (1 decimal)

    # Build table_12_results (7 data columns = cols 2:8 of plot_dev.jl's 8-col matrix)
    # Rows 1-25 map to plot_dev.jl rows 1-25 (col 1 = row index discarded here)
    r = zeros(26, 7)

    # Row 1: TFP
    r[1,:] = [TFP_initial[2] / TFP_initial[2],
              TFP_closed_CM_open_trade[2] / TFP_initial[2],
              TFP_open_CM_open_trade[2]   / TFP_initial[2],
              TFP_open_CM_closed_trade[2] / TFP_initial[2],
              TFP_initial_dev[2]          / TFP_initial[2],
              TFP_closed_CM_open_trade_dev[2] / TFP_initial_dev[2],
              TFP_open_CM_open_trade_dev[2]   / TFP_initial_dev[2]]
    # Row 2: sd_MRPK (Core)
    r[2,:] = [sd_MRPK_initial[2], sd_MRPK_closed_CM_open_trade[2],
              sd_MRPK_open_CM_open_trade[2], sd_MRPK_open_CM_closed_trade[2],
              sd_MRPK_initial_dev[2], sd_MRPK_closed_CM_open_trade_dev[2],
              sd_MRPK_open_CM_open_trade_dev[2]]
    # Row 3: GDP (Output)
    r[3,:] = [1.0,
              GDP_closed_CM_open_trade[2] / GDP_initial[2],
              GDP_open_CM_open_trade[2]   / GDP_initial[2],
              GDP_open_CM_closed_trade[2] / GDP_initial[2],
              GDP_initial_dev[2]          / GDP_initial[2],
              GDP_closed_CM_open_trade_dev[2] / GDP_initial_dev[2],
              GDP_open_CM_open_trade_dev[2]   / GDP_initial_dev[2]]
    # Row 4: Income (foreign wage = prices[4])
    r[4,:] = [1.0,
              prices_closed_CM_open_trade[4] / prices_initial[4],
              prices_open_CM_open_trade[4]   / prices_initial[4],
              prices_open_CM_closed_trade[4] / prices_initial[4],
              prices_initial_dev[4]          / prices_initial[4],
              prices_closed_CM_open_trade_dev[4] / prices_initial_dev[4],
              prices_open_CM_open_trade_dev[4]   / prices_initial_dev[4]]
    # Row 5: Consumption
    r[5,:] = [1.0,
              total_consumption_closed_CM_open_trade[2] / total_consumption_initial[2],
              total_consumption_open_CM_open_trade[2]   / total_consumption_initial[2],
              total_consumption_open_CM_closed_trade[2] / total_consumption_initial[2],
              total_consumption_initial_dev[2]          / total_consumption_initial[2],
              total_consumption_closed_CM_open_trade_dev[2] / total_consumption_initial_dev[2],
              total_consumption_open_CM_open_trade_dev[2]   / total_consumption_initial_dev[2]]
    # Row 6: Capital
    r[6,:] = [1.0,
              capital_demand_closed_CM_open_trade[2] / capital_demand_initial[2],
              capital_demand_open_CM_open_trade[2]   / capital_demand_initial[2],
              capital_demand_open_CM_closed_trade[2] / capital_demand_initial[2],
              capital_demand_initial_dev[2]          / capital_demand_initial[2],
              capital_demand_closed_CM_open_trade_dev[2] / capital_demand_initial_dev[2],
              capital_demand_open_CM_open_trade_dev[2]   / capital_demand_initial_dev[2]]
    # Row 7: CE welfare (Core transition dynamics)
    r[7,:] = [0.0, welfare_change_clCM_otrade_trans_F, welfare_change_oCM_otrade_trans_F,
              welfare_change_oCM_cltrade_trans_F, welfare_init_developed_stst_F,
              welfare_change_clCM_otrade_trans_dev_F, welfare_change_oCM_otrade_trans_dev_F]
    # Row 8: EU welfare (population-weighted NMS+Core)
    r[8,:] = [0.0,
              (welfare_change_clCM_otrade_trans + Baseline_parameter.L[2] * welfare_change_clCM_otrade_trans_F) / (Baseline_parameter.L[2] + 1),
              (welfare_change_oCM_otrade_trans  + Baseline_parameter.L[2] * welfare_change_oCM_otrade_trans_F)  / (Baseline_parameter.L[2] + 1),
              (welfare_change_oCM_cltrade_trans + Baseline_parameter.L[2] * welfare_change_oCM_cltrade_trans_F) / (Baseline_parameter.L[2] + 1),
              (welfare_init_developed_stst      + Baseline_parameter.L[2] * welfare_init_developed_stst_F)      / (Baseline_parameter.L[2] + 1),
              (welfare_change_clCM_otrade_trans_dev + Baseline_parameter.L[2] * welfare_change_clCM_otrade_trans_dev_F) / (Baseline_parameter.L[2] + 1),
              (welfare_change_oCM_otrade_trans_dev  + Baseline_parameter.L[2] * welfare_change_oCM_otrade_trans_dev_F)  / (Baseline_parameter.L[2] + 1)]
    # Row 9: Top 10% wealth share (Core)
    r[9,:] = [p90_wealth_initial[2], p90_wealth_closed_CM_open_trade[2],
              p90_wealth_open_CM_open_trade[2], p90_wealth_open_CM_closed_trade[2],
              p90_wealth_initial_dev[2], p90_wealth_closed_CM_open_trade_dev[2],
              p90_wealth_open_CM_open_trade_dev[2]]
    # Row 10: Top 10% income share
    r[10,:] = [p90_income_initial[2], p90_income_closed_CM_open_trade[2],
               p90_income_open_CM_open_trade[2], p90_income_open_CM_closed_trade[2],
               p90_income_initial_dev[2], p90_income_closed_CM_open_trade_dev[2],
               p90_income_open_CM_open_trade_dev[2]]
    # Row 11: Top 10% consumption share
    r[11,:] = [p90_cons_initial[2], p90_cons_closed_CM_open_trade[2],
               p90_cons_open_CM_open_trade[2], p90_cons_open_CM_closed_trade[2],
               p90_cons_initial_dev[2], p90_cons_closed_CM_open_trade_dev[2],
               p90_cons_open_CM_open_trade_dev[2]]
    # Row 12: Wealth of exporters
    r[12,:] = [wealth_of_exporters_initial[2], wealth_of_exporters_closed_CM_open_trade[2],
               wealth_of_exporters_open_CM_open_trade[2], wealth_of_exporters_open_CM_closed_trade[2],
               wealth_of_exporters_initial_dev[2], wealth_of_exporters_closed_CM_open_trade_dev[2],
               wealth_of_exporters_open_CM_open_trade_dev[2]]
    # Row 13: Real wage (Core, prices[2] = Core wage normalized)
    r[13,:] = [1.0,
               prices_closed_CM_open_trade[2] / prices_initial[2],
               prices_open_CM_open_trade[2]   / prices_initial[2],
               prices_open_CM_closed_trade[2] / prices_initial[2],
               prices_initial_dev[2]          / prices_initial[2],
               prices_closed_CM_open_trade_dev[2] / prices_initial_dev[2],
               prices_open_CM_open_trade_dev[2]   / prices_initial_dev[2]]
    # Row 14: r* (Core borrowing rate; use prices[6] when CM open, prices[7] otherwise)
    r[14,:] = [prices_initial[7],
               prices_closed_CM_open_trade[7],
               prices_open_CM_open_trade[6],
               prices_open_CM_closed_trade[6],
               prices_initial_dev[7],
               prices_closed_CM_open_trade_dev[7],
               prices_open_CM_open_trade_dev[6]]
    # Row 15: Import/GDP* (Core imports / Core GDP)
    r[15,:] = [import_share_initial[2], import_share_closed_CM_open_trade[2],
               import_share_open_CM_open_trade[2], import_share_open_CM_closed_trade[2],
               import_share_initial_dev[2], import_share_closed_CM_open_trade_dev[2],
               import_share_open_CM_open_trade_dev[2]]
    # Row 16: Export/GDP (Core exports / NMS GDP)
    r[16,:] = [export_value_initial[2]       / nomGDP_initial[1],
               export_value_closed_CM_open_trade[2] / nomGDP_closed_CM_open_trade[1],
               export_value_open_CM_open_trade[2]   / nomGDP_open_CM_open_trade[1],
               export_value_open_CM_closed_trade[2] / nomGDP_open_CM_closed_trade[1],
               export_value_initial_dev[2]          / nomGDP_initial_dev[1],
               export_value_closed_CM_open_trade_dev[2] / nomGDP_closed_CM_open_trade_dev[1],
               export_value_open_CM_open_trade_dev[2]   / nomGDP_open_CM_open_trade_dev[1]]
    # Row 17: Entrepreneurship rate (Core)
    r[17,:] = [(domestic_pop_initial[2] + exporter_pop_initial[2]) / Baseline_parameter.L[2],
               (domestic_pop_closed_CM_open_trade[2] + exporter_pop_closed_CM_open_trade[2]) / Baseline_parameter.L[2],
               (domestic_pop_open_CM_open_trade[2]   + exporter_pop_open_CM_open_trade[2])   / Baseline_parameter.L[2],
               (domestic_pop_open_CM_closed_trade[2] + exporter_pop_open_CM_closed_trade[2]) / Baseline_parameter.L[2],
               (domestic_pop_initial_dev[2]          + exporter_pop_initial_dev[2])           / Baseline_parameter.L[2],
               (domestic_pop_closed_CM_open_trade_dev[2] + exporter_pop_closed_CM_open_trade_dev[2]) / Baseline_parameter.L[2],
               (domestic_pop_open_CM_open_trade_dev[2]   + exporter_pop_open_CM_open_trade_dev[2])   / Baseline_parameter.L[2]]
    # Row 18: Share of exporters (Core)
    r[18,:] = [exporter_pop_initial[2] /
                   (domestic_pop_initial[2] + exporter_pop_initial[2]),
               exporter_pop_closed_CM_open_trade[2] /
                   (domestic_pop_closed_CM_open_trade[2] + exporter_pop_closed_CM_open_trade[2]),
               exporter_pop_open_CM_open_trade[2] /
                   (domestic_pop_open_CM_open_trade[2] + exporter_pop_open_CM_open_trade[2]),
               exporter_pop_open_CM_closed_trade[2] /
                   (domestic_pop_open_CM_closed_trade[2] + exporter_pop_open_CM_closed_trade[2]),
               exporter_pop_initial_dev[2] /
                   (domestic_pop_initial_dev[2] + exporter_pop_initial_dev[2]),
               exporter_pop_closed_CM_open_trade_dev[2] /
                   (domestic_pop_closed_CM_open_trade_dev[2] + exporter_pop_closed_CM_open_trade_dev[2]),
               exporter_pop_open_CM_open_trade_dev[2] /
                   (domestic_pop_open_CM_open_trade_dev[2] + exporter_pop_open_CM_open_trade_dev[2])]
    # Row 19: 1/P* (inverse of NMS CPI/RER)
    r[19,:] = [1/prices_initial[5], 1/prices_closed_CM_open_trade[5],
               1/prices_open_CM_open_trade[5], 1/prices_open_CM_closed_trade[5],
               1/prices_initial_dev[5], 1/prices_closed_CM_open_trade_dev[5],
               1/prices_open_CM_open_trade_dev[5]]
    # Row 20: Credit/GDP (Core)
    r[20,:] = [total_credit_initial_F, total_credit_closed_CM_open_trade_F,
               total_credit_open_CM_open_trade_F, total_credit_open_CM_closed_trade_F,
               total_credit_initial_dev_F, total_credit_closed_CM_open_trade_dev_F,
               total_credit_open_CM_open_trade_dev_F]
    # Row 21: Credit Abroad / Credit (Core)
    r[21,:] = [foreign_credit_share_initial_F, foreign_credit_share_closed_CM_open_trade_F,
               foreign_credit_share_open_CM_open_trade_F, foreign_credit_share_open_CM_closed_trade_F,
               foreign_credit_share_initial_dev_F, foreign_credit_share_closed_CM_open_trade_dev_F,
               foreign_credit_share_open_CM_open_trade_dev_F]
    # Row 22: Sector sd_MRPK domestic (Core)
    r[22,:] = [sd_MRPK_d_initial[2], sd_MRPK_d_closed_CM_open_trade[2],
               sd_MRPK_d_open_CM_open_trade[2], sd_MRPK_d_open_CM_closed_trade[2],
               sd_MRPK_d_initial_dev[2], sd_MRPK_d_closed_CM_open_trade_dev[2],
               sd_MRPK_d_open_CM_open_trade_dev[2]]
    # Row 23: Sector sd_MRPK exporter (Core)
    r[23,:] = [sd_MRPK_x_initial[2], sd_MRPK_x_closed_CM_open_trade[2],
               sd_MRPK_x_open_CM_open_trade[2], sd_MRPK_x_open_CM_closed_trade[2],
               sd_MRPK_x_initial_dev[2], sd_MRPK_x_closed_CM_open_trade_dev[2],
               sd_MRPK_x_open_CM_open_trade_dev[2]]
    # Row 24: Within-sector productivity loss — domestic (Core)
    r[24,:] = [Misallocation_within_d_initial[2], Misallocation_within_d_closed_CM_open_trade[2],
               Misallocation_within_d_open_CM_open_trade[2], Misallocation_within_d_open_CM_closed_trade[2],
               Misallocation_within_d_initial_dev[2], Misallocation_within_d_closed_CM_open_trade_dev[2],
               Misallocation_within_d_open_CM_open_trade_dev[2]]
    # Row 25: Within-sector productivity loss — exporter (Core)
    r[25,:] = [Misallocation_within_x_initial[2], Misallocation_within_x_closed_CM_open_trade[2],
               Misallocation_within_x_open_CM_open_trade[2], Misallocation_within_x_open_CM_closed_trade[2],
               Misallocation_within_x_initial_dev[2], Misallocation_within_x_closed_CM_open_trade_dev[2],
               Misallocation_within_x_open_CM_open_trade_dev[2]]

    println("[table:Table4b] table_12_results computed.")
    println(round.(100 .* r, digits=1))

    # CE welfare col 3 gets asterisk (delayed-10yr footnote in original)
    ce3 = @sprintf("%.1f", round(r[7,3]*100, digits=1, RoundNearestTiesAway)) * "\\textsuperscript{*}"

    tex = """\\begin{table}[!t]
\\centering
\\caption{Effect of reforms on the Core economy}
\\label{table:Table4b}
  \\scalebox{0.60}{
\\begin{tabular}{l ccc ccc ccc ccc ccc ccc ccc}
\\\\[-1.8ex]\\hline
\\hline \\\\[-1.8ex]
Integration & None  & Trade & Trade and capital & Capital & Development  & Trade & Trade and capital\\\\
\\hline \\\\[-1.8ex]
\\textbf{Productivity}& &   &   \\\\
TFP & $(f1(r[1,1]))  & $(f1(r[1,2]))& $(f1(r[1,3]))& $(f1(r[1,4])) & $(f1(r[1,5])) & $(f1(r[1,6])) & $(f1(r[1,7])) \\\\
Standard deviation of ARPK & $(f2(r[2,1]))  & $(f2(r[2,2]))& $(f2(r[2,3])) & $(f2(r[2,4]))& $(f2(r[2,5])) & $(f2(r[2,6])) &$(f2(r[2,7]))\\\\[1ex]
\\textbf{Aggregates}& &   &   \\\\
Output & $(f0(r[3,1]))& $(f0(r[3,2]))& $(f0(r[3,3])) &$(f0(r[3,4])) & $(f0(r[3,5])) & $(f0(r[3,6])) &$(f0(r[3,7]))\\\\
Income & $(f0(r[4,1]))& $(f0(r[4,2]))& $(f0(r[4,3]))&$(f0(r[4,4]))&$(f0(r[4,5]))&$(f0(r[4,6]))&$(f0(r[4,7])) \\\\
Consumption & $(f0(r[5,1])) & $(f0(r[5,2])) & $(f0(r[5,3])) &$(f0(r[5,4])) &$(f0(r[5,5])) &$(f0(r[5,6])) &$(f0(r[5,7])) \\\\
Capital & $(f0(r[6,1]))& $(f0(r[6,2]))& $(f0(r[6,3])) & $(f0(r[6,4])) & $(f0(r[6,5])) & $(f0(r[6,6])) & $(f0(r[6,7])) \\\\[1ex]
\\textbf{Welfare and Inequality}& &   &   \\\\
Consumption equivalent welfare  & 0 & $(f1(r[7,2]))& $(ce3) & $(f1(r[7,4])) & $(f1(r[7,5])) & $(f1(r[7,6])) & $(f1(r[7,7]))  \\\\
EU welfare & 0 & $(f1(r[8,2]))& $(f1(r[8,3])) & $(f1(r[8,4])) & $(f1(r[8,5])) & $(f1(r[8,6])) & $(f1(r[8,7]))  \\\\
Top \$10\\%\$ wealth share &  $(f0(r[9,1])) & $(f0(r[9,2])) & $(f0(r[9,3]))& $(f0(r[9,4])) & $(f0(r[9,5]))& $(f0(r[9,6]))& $(f0(r[9,7]))  \\\\
Top \$10\\%\$ income share &  $(f0(r[10,1])) & $(f0(r[10,2])) & $(f0(r[10,3]))& $(f0(r[10,4])) & $(f0(r[10,5]))& $(f0(r[10,6]))& $(f0(r[10,7]))  \\\\
Top \$10\\%\$ consumption share &  $(f0(r[11,1])) & $(f0(r[11,2])) & $(f0(r[11,3]))& $(f0(r[11,4])) & $(f0(r[11,5]))& $(f0(r[11,6]))& $(f0(r[11,7]))  \\\\
Wealth owned by exporters & $(f0(r[12,1])) & $(f0(r[12,2])) & $(f0(r[12,3]))& $(f0(r[12,4])) & $(f0(r[12,5]))& $(f0(r[12,6])) & $(f0(r[12,7]))  \\\\[1ex]
\\textbf{Factor prices} & &   &\\\\
Real wage & $(f0(r[13,1])) & $(f0(r[13,2])) & $(f0(r[13,3])) & $(f0(r[13,4])) & $(f0(r[13,5])) & $(f0(r[13,6])) & $(f0(r[13,7]))  \\\\
Interest rate \$\\%\$: \$r^*\$ & $(fip(r[14,1])) & $(fip(r[14,2])) & $(fip(r[14,3])) & $(fip(r[14,4])) & $(fip(r[14,5])) & $(fip(r[14,6])) & $(fip(r[14,7])) \\\\[1ex]
\\textbf{Trade}& &   &\\\\
\$\\frac{\\text{Import}}{\\text{GDP}^*}\$ & $(fimp(r[15,1])) & $(fimp(r[15,2])) & $(fimp(r[15,3])) & $(fimp(r[15,4])) &$(fimp(r[15,5])) & $(fimp(r[15,6])) & $(fimp(r[15,7]))  \\\\
Entrepreneurship rate & $(f0(r[17,1])) & $(f0(r[17,2])) & $(f0(r[17,3])) & $(f0(r[17,4])) & $(f0(r[17,5])) & $(f0(r[17,6]))& $(f0(r[17,7])) \\\\
Share of exporters & $(f0(r[18,1])) & $(f0(r[18,2])) & $(f0(r[18,3])) & $(f0(r[18,4])) & $(f0(r[18,5])) & $(f0(r[18,6])) & $(f0(r[18,7])) \\\\
\$\\frac{1}{P^*}\$ & $(f0(r[19,1])) & $(f0(r[19,2])) & $(f0(r[19,3])) & $(f0(r[19,4])) & $(f0(r[19,5])) & $(f0(r[19,6])) & $(f0(r[19,7])) \\\\
\$\\frac{\\text{Credit}}{\\text{GDP}}\$ & $(f0(r[20,1])) & $(f0(r[20,2])) & $(f0(r[20,3])) & $(f0(r[20,4])) & $(f0(r[20,5])) & $(f0(r[20,6])) & $(f0(r[20,7]))\\\\
\$\\frac{\\text{Credit Abroad}}{\\text{Credit}}\$ & $(f1(r[21,1])) & $(f1(r[21,2])) & $(f1(r[21,3])) & $(f1(r[21,4])) & $(f1(r[21,5]))& $(f1(r[21,6])) & $(f1(r[21,7]))\\\\[1ex]
\\textbf{Sector s.d. \$mrpk\$} & & & \\\\
 Domestic & $(f2(r[22,1])) & $(f2(r[22,2])) & $(f2(r[22,3])) &$(f2(r[22,4])) & $(f2(r[22,5])) & $(f2(r[22,6])) & $(f2(r[22,7])) \\\\
 Exporter& $(f2(r[23,1])) & $(f2(r[23,2])) & $(f2(r[23,3])) &$(f2(r[23,4])) & $(f2(r[23,5])) & $(f2(r[23,6])) & $(f2(r[23,7])) \\\\[1ex]
\\textbf{Within Sector productivity loss} & & & \\\\
 Domestic  & $(f1(r[24,1])) & $(f1(r[24,2])) & $(f1(r[24,3])) & $(f1(r[24,4])) & $(f1(r[24,5])) & $(f1(r[24,6])) & $(f1(r[24,7])) \\\\
 Exporter & $(f1(r[25,1])) & $(f1(r[25,2])) & $(f1(r[25,3])) & $(f1(r[25,4])) & $(f1(r[25,5])) & $(f1(r[25,6])) & $(f1(r[25,7])) \\\\
\\hline
\\end{tabular} }
\\parbox{0.85\\textwidth}{\\caption*{\\scriptsize \\textit{Note:}  Postponing by 10 years increases gains by \$0.1 \\%\$. The last two columns are relative to the ``Development'' steady state. Except for EU welfare, all variables refer to the Core economy.}}
\\end{table}
"""
    open("../Tables/Table4b.tex", "w") do io
        write(io, tex)
    end
    println("[table:Table4b] Saved → ../Tables/Table4b.tex")
end

# =============================================================================
# [table:Table_asym]  Symmetric vs asymmetric trade liberalization welfare (computed)
# main.tex: \label{table:Table_asym}
# Source: plot_asym.jl lines 1–58
# Columns: Closed (trade only), Integrated (full), Delayed integration
# =============================================================================
let
    β1 = Baseline_parameter.β[1]
    β2 = Baseline_parameter.β[2]

    # NMS asymmetric welfare (from Step 2c: V_saved_store_*_asym)
    V_asym_cl  = reshape(V_saved_store_closed_CM_open_trade_asym[:,:,1,2],
                          size(V_saved_fine_initial)[1]*3)
    V_asym_op  = reshape(V_saved_store_open_CM_open_trade_asym[:,:,1,2],
                          size(V_saved_fine_initial)[1]*3)
    V_asym_del = reshape(V_saved_store_open_CMdelayed_open_trade_asym[:,:,1,2],
                          size(V_saved_fine_initial)[1]*3)

    welf_NMS_asym_cl  = sum(current_distr_store_initial[:,1] .*
                             (exp.((V_asym_cl  - V_saved_initial) .* (1.0 - β1)) .- 1))
    welf_NMS_asym_op  = sum(current_distr_store_initial[:,1] .*
                             (exp.((V_asym_op  - V_saved_initial) .* (1.0 - β1)) .- 1))
    welf_NMS_asym_del = sum(current_distr_store_initial[:,1] .*
                             (exp.((V_asym_del - V_saved_initial) .* (1.0 - β1)) .- 1))

    # Core asymmetric welfare
    # V_saved_initial_F is already in scope from the shared welfare section
    V_asym_cl_F  = reshape(V_saved_store_closed_CM_open_trade_asym[:,:,2,2],
                            size(V_saved_fine_initial)[1]*3)
    V_asym_op_F  = reshape(V_saved_store_open_CM_open_trade_asym[:,:,2,2],
                            size(V_saved_fine_initial)[1]*3)
    V_asym_del_F = reshape(V_saved_store_open_CMdelayed_open_trade_asym[:,:,2,2],
                            size(V_saved_fine_initial)[1]*3)

    welf_Core_asym_cl  = sum(current_distr_store_initial[:,2] ./ Baseline_parameter.L[2] .*
                              (exp.((V_asym_cl_F  - V_saved_initial_F) .* (1.0 - β2)) .- 1))
    welf_Core_asym_op  = sum(current_distr_store_initial[:,2] ./ Baseline_parameter.L[2] .*
                              (exp.((V_asym_op_F  - V_saved_initial_F) .* (1.0 - β2)) .- 1))
    welf_Core_asym_del = sum(current_distr_store_initial[:,2] ./ Baseline_parameter.L[2] .*
                              (exp.((V_asym_del_F - V_saved_initial_F) .* (1.0 - β2)) .- 1))

    # Symmetric welfare already computed in the shared welfare section:
    #   NMS: welfare_change_clCM_otrade_trans, welfare_change_oCM_otrade_trans,
    #        welfare_change_oCM_otrade_delayed10_trans
    #  Core: welfare_change_clCM_otrade_trans_F, welfare_change_oCM_otrade_trans_F,
    #        welfare_change_oCM_otrade_delayed10_trans_F

    f1(x) = @sprintf("%.1f", round(round(x * 100, digits=2, RoundNearestTiesAway), digits=1, RoundNearestTiesAway))

    println("[table:Table_asym] Welfare values:")
    println("  NMS sym:  ", round.([welfare_change_clCM_otrade_trans, welfare_change_oCM_otrade_trans, welfare_change_oCM_otrade_delayed10_trans] .* 100, digits=1))
    println("  NMS asym: ", round.([welf_NMS_asym_cl, welf_NMS_asym_op, welf_NMS_asym_del] .* 100, digits=1))
    println("  Core sym: ", round.([welfare_change_clCM_otrade_trans_F, welfare_change_oCM_otrade_trans_F, welfare_change_oCM_otrade_delayed10_trans_F] .* 100, digits=1))
    println("  Core asym:", round.([welf_Core_asym_cl, welf_Core_asym_op, welf_Core_asym_del] .* 100, digits=1))

    tex = """\\begin{table}[h]
 \\caption{Comparing welfare gains following a symmetric vs asymmetric trade liberalization }
    \\label{table:Table_asym}
    \\centering
    \\scalebox{0.8}{%
    \\begin{tabular}{l ccc ccc ccc}
        \\hline
        \\hline \\\\[-1.8ex]
        Capital Market & Closed & Integrated & Delayed integration \\\\
        \\hline
        NMS CE Welfare - symmetric  & $(f1(welfare_change_clCM_otrade_trans)) & $(f1(welfare_change_oCM_otrade_trans)) & $(f1(welfare_change_oCM_otrade_delayed10_trans))  \\\\
        NMS CE Welfare - asymmetric  & $(f1(welf_NMS_asym_cl))  & $(f1(welf_NMS_asym_op)) & $(f1(welf_NMS_asym_del)) \\\\
        Core CE Welfare - symmetric  &  $(f1(welfare_change_clCM_otrade_trans_F))& $(f1(welfare_change_oCM_otrade_trans_F)) & $(f1(welfare_change_oCM_otrade_delayed10_trans_F))  \\\\
        Core CE Welfare - asymmetric  & $(f1(welf_Core_asym_cl))& $(f1(welf_Core_asym_op)) & $(f1(welf_Core_asym_del)) \\\\
        \\hline
    \\end{tabular}
}
\\end{table}
"""
    open("../Tables/Table_asym.tex", "w") do io
        write(io, tex)
    end
    println("[table:Table_asym] Saved → ../Tables/Table_asym.tex")
end

# =============================================================================
# [fig:Figure12_cons]  Transition dynamics of consumption — symmetric + asymmetric
# main.tex: \label{fig:Figure12_cons}
# Figures: Figure12c.pdf  (already saved in [fig:Figure12] above)
#          Figure12c_asym.pdf  (from asymmetric liberalization, Step 2c)
# Source: plot_asym.jl lines 106–116
# =============================================================================
let
    # Figure12c already saved above ([fig:Figure12]).
    # Figure12c_asym:
    con_corr_open_a    = total_consumption_open_CM_open_trade[1] / total_consumption_store_open_CM_open_trade_asym[1,end]
    con_corr_closed_a  = total_consumption_closed_CM_open_trade[1] / total_consumption_store_closed_CM_open_trade_asym[1,end]
    con_corr_delayed_a = total_consumption_open_CM_open_trade[1] / total_consumption_store_open_CMdelayed_open_trade_asym[1,end-1]
    con_smooth_open_a    = movmean(100 .* (con_corr_open_a    .* total_consumption_store_open_CM_open_trade_asym[1,:] ./ total_consumption_initial[1] .- 1), 5)
    con_smooth_closed_a  = movmean(100 .* (con_corr_closed_a  .* total_consumption_store_closed_CM_open_trade_asym[1,:] ./ total_consumption_initial[1] .- 1), 5)
    con_smooth_delayed_a = movmean(100 .* (con_corr_delayed_a .* total_consumption_store_open_CMdelayed_open_trade_asym[1,:] ./ total_consumption_initial[1] .- 1), 5)
    plot([range(1991,2020)], [con_smooth_open_a, con_smooth_closed_a, con_smooth_delayed_a],
         legend=nothing, linewidth=5, linestyle=[:solid :dash :dot], grid=false,
         tickfontsize=14, ylims=[0,10])
    savefig("../Figures/Figure12c_asym.pdf")
    println("[fig:Figure12_cons] Saved Figure12c_asym.pdf")
end

# =============================================================================
# [fig:Figure12_TFP]  Transition dynamics of TFP — symmetric + asymmetric
# main.tex: \label{fig:Figure12_TFP}
# Figures: Figure12a.pdf  (already saved in [fig:Figure12] above)
#          Figure12a_asym.pdf  (from asymmetric liberalization, Step 2c)
# Source: plot_asym.jl lines 78–89
# =============================================================================
let
    # Figure12a already saved above.
    # Figure12a_asym:
    TFP_corr_open_a    = TFP_open_CM_open_trade[1] / TFP_store_open_CM_open_trade_asym[1,end]
    TFP_corr_closed_a  = TFP_closed_CM_open_trade[1] / TFP_store_closed_CM_open_trade_asym[1,end-1]
    TFP_corr_delayed_a = TFP_open_CM_open_trade[1] / TFP_store_open_CMdelayed_open_trade_asym[1,end-1]
    TFP_smooth_open_a    = movmean(100 .* (TFP_corr_open_a    .* TFP_store_open_CM_open_trade_asym[1,:] ./ TFP_initial[1] .- 1), 5)
    TFP_smooth_closed_a  = movmean(100 .* (TFP_corr_closed_a  .* TFP_store_closed_CM_open_trade_asym[1,:] ./ TFP_initial[1] .- 1), 5)
    TFP_smooth_delayed_a = movmean(100 .* (TFP_corr_delayed_a .* TFP_store_open_CMdelayed_open_trade_asym[1,:] ./ TFP_initial[1] .- 1), 5)
    plot([range(1991,2020)], [TFP_smooth_open_a, TFP_smooth_closed_a, TFP_smooth_delayed_a],
         legend=nothing, linewidth=5, linestyle=[:solid :dash :dot], grid=false,
         tickfontsize=14, ylims=[0,35])
    savefig("../Figures/Figure12a_asym.pdf")
    println("[fig:Figure12_TFP] Saved Figure12a_asym.pdf")
end

# =============================================================================
# [fig:Figure12_ARPK]  Transition dynamics of s.d. ARPK — symmetric + asymmetric
# main.tex: \label{fig:Figure12_ARPK}
# Figures: Figure12k.pdf  (already saved in [fig:Figure12] above)
#          Figure12k_asym.pdf  (from asymmetric liberalization, Step 2c)
# Source: plot_asym.jl lines 196–206
# =============================================================================
let
    # Figure12k already saved above.
    # Figure12k_asym:
    sdMRPK_corr_open_a    = sd_MRPK_open_CM_open_trade[1] / sd_MRPK_store_open_CM_open_trade_asym[1,end]
    sdMRPK_corr_closed_a  = sd_MRPK_closed_CM_open_trade[1] / sd_MRPK_store_closed_CM_open_trade_asym[1,end]
    sdMRPK_corr_delayed_a = sd_MRPK_open_CM_open_trade[1] / sd_MRPK_store_open_CMdelayed_open_trade_asym[1,end-1]
    sdMRPK_smooth_open_a    = movmean(sdMRPK_corr_open_a    .* sd_MRPK_store_open_CM_open_trade_asym[1,:], 5)
    sdMRPK_smooth_closed_a  = movmean(sdMRPK_corr_closed_a  .* sd_MRPK_store_closed_CM_open_trade_asym[1,:], 5)
    sdMRPK_smooth_delayed_a = movmean(sdMRPK_corr_delayed_a .* sd_MRPK_store_open_CMdelayed_open_trade_asym[1,:], 5)
    plot([range(1991,2020)], [sdMRPK_smooth_open_a, sdMRPK_smooth_closed_a, sdMRPK_smooth_delayed_a],
         legend=nothing, linewidth=5, linestyle=[:solid :dash :dot], grid=false,
         tickfontsize=14, ylims=[0.3,0.55])
    savefig("../Figures/Figure12k_asym.pdf")
    println("[fig:Figure12_ARPK] Saved Figure12k_asym.pdf")
end

# =============================================================================
# [fig:Figure12_asym]  Full set of transition dynamics under asymmetric liberalization
# main.tex: \label{fig:Figure12_asym}
# Figures: Figure12a_asym (TFP)           — already saved above
#          Figure12b_asym (capital)
#          Figure12c_asym (consumption)   — already saved above
#          Figure12d_asym (exporters)
#          Figure12e_asym (non-exporters)
#          Figure12f_asym (exporter share)
#          Figure12g_asym (K/exporter)
#          Figure12h_asym (K/non-exporter)
#          Figure12i_asym (capital share exporters)
#          Figure12k_asym (s.d. ARPK all) — already saved above
#          Figure12l_asym (s.d. ARPK domestic)
#          Figure12m_asym (s.d. ARPK exporter)
# Source: plot_asym.jl lines 91–228
# =============================================================================
let
    # Figure12b_asym — capital
    cap_corr_open_a    = capital_demand_open_CM_open_trade[1] / capital_supply_future_open_CM_open_trade_asym[1,end]
    cap_corr_closed_a  = capital_demand_closed_CM_open_trade[1] / capital_supply_future_closed_CM_open_trade_asym[1,end]
    cap_corr_delayed_a = capital_demand_open_CM_open_trade[1] / capital_supply_future_open_CMdelayed_open_trade_asym[1,end-1]
    cap_smooth_open_a    = movmean(100 .* (cap_corr_open_a    .* capital_supply_future_open_CM_open_trade_asym[1,:] ./ capital_demand_initial[1] .- 1), 5)
    cap_smooth_closed_a  = movmean(100 .* (cap_corr_closed_a  .* capital_supply_future_closed_CM_open_trade_asym[1,:] ./ capital_demand_initial[1] .- 1), 5)
    cap_smooth_delayed_a = movmean(100 .* (cap_corr_delayed_a .* capital_supply_future_open_CMdelayed_open_trade_asym[1,:] ./ capital_demand_initial[1] .- 1), 5)
    plot([range(1991,2022)], [cap_smooth_open_a, cap_smooth_closed_a, cap_smooth_delayed_a],
         legend=nothing, linewidth=5, linestyle=[:solid :dash :dot], grid=false,
         tickfontsize=14, ylims=[-10,40])
    savefig("../Figures/Figure12b_asym.pdf")
    println("[fig:Figure12_asym] Saved Figure12b_asym.pdf")

    # Figure12d_asym — exporters
    exp_corr_open_a    = exporter_pop_open_CM_open_trade[1] / exporter_pop_store_open_CM_open_trade_asym[1,end]
    exp_corr_closed_a  = exporter_pop_closed_CM_open_trade[1] / exporter_pop_store_closed_CM_open_trade_asym[1,end]
    exp_corr_delayed_a = exporter_pop_open_CM_open_trade[1] / exporter_pop_store_open_CMdelayed_open_trade_asym[1,end-1]
    exp_smooth_open_a    = movmean(100 .* (exp_corr_open_a    .* exporter_pop_store_open_CM_open_trade_asym[1,:] ./ exporter_pop_initial[1] .- 1), 5)
    exp_smooth_closed_a  = movmean(100 .* (exp_corr_closed_a  .* exporter_pop_store_closed_CM_open_trade_asym[1,:] ./ exporter_pop_initial[1] .- 1), 5)
    exp_smooth_delayed_a = movmean(100 .* (exp_corr_delayed_a .* exporter_pop_store_open_CMdelayed_open_trade_asym[1,:] ./ exporter_pop_initial[1] .- 1), 5)
    plot([range(1991,2020)], [exp_smooth_open_a, exp_smooth_closed_a, exp_smooth_delayed_a],
         legend=nothing, linewidth=5, linestyle=[:solid :dash :dot], grid=false,
         tickfontsize=14, ylims=[0.0,45])
    savefig("../Figures/Figure12d_asym.pdf")
    println("[fig:Figure12_asym] Saved Figure12d_asym.pdf")

    # Figure12e_asym — non-exporters
    dom_corr_open_a    = domestic_pop_open_CM_open_trade[1] / domestic_pop_store_open_CM_open_trade_asym[1,end]
    dom_corr_closed_a  = domestic_pop_closed_CM_open_trade[1] / domestic_pop_store_closed_CM_open_trade_asym[1,end]
    dom_corr_delayed_a = domestic_pop_open_CM_open_trade[1] / domestic_pop_store_open_CMdelayed_open_trade_asym[1,end-1]
    dom_smooth_open_a    = movmean(100 .* (dom_corr_open_a    .* domestic_pop_store_open_CM_open_trade_asym[1,:] ./ domestic_pop_initial[1] .- 1), 5)
    dom_smooth_closed_a  = movmean(100 .* (dom_corr_closed_a  .* domestic_pop_store_closed_CM_open_trade_asym[1,:] ./ domestic_pop_initial[1] .- 1), 5)
    dom_smooth_delayed_a = movmean(100 .* (dom_corr_delayed_a .* domestic_pop_store_open_CMdelayed_open_trade_asym[1,:] ./ domestic_pop_initial[1] .- 1), 5)
    plot([range(1991,2020)], [dom_smooth_open_a, dom_smooth_closed_a, dom_smooth_delayed_a],
         legend=nothing, linewidth=5, linestyle=[:solid :dash :dot], grid=false,
         tickfontsize=14, ylims=[-30.0,30.0])
    savefig("../Figures/Figure12e_asym.pdf")
    println("[fig:Figure12_asym] Saved Figure12e_asym.pdf")

    # Figure12f_asym — exporter share
    share_open_a    = movmean(exporter_pop_store_open_CM_open_trade_asym[1,:] ./ (exporter_pop_store_open_CM_open_trade_asym[1,:] .+ domestic_pop_store_open_CM_open_trade_asym[1,:]), 5)
    share_closed_a  = movmean(exporter_pop_store_closed_CM_open_trade_asym[1,:] ./ (exporter_pop_store_closed_CM_open_trade_asym[1,:] .+ domestic_pop_store_closed_CM_open_trade_asym[1,:]), 5)
    share_delayed_a = movmean(exporter_pop_store_open_CMdelayed_open_trade_asym[1,:] ./ (exporter_pop_store_open_CMdelayed_open_trade_asym[1,:] .+ domestic_pop_store_open_CMdelayed_open_trade_asym[1,:]), 5)
    plot([range(1991,2020)], 100 .* [share_open_a, share_closed_a, share_delayed_a],
         legend=nothing, linewidth=5, linestyle=[:solid :dash :dot], grid=false,
         tickfontsize=14, ylims=[30,50])
    savefig("../Figures/Figure12f_asym.pdf")
    println("[fig:Figure12_asym] Saved Figure12f_asym.pdf")

    # Figure12g_asym — capital per exporter
    Kx_init_per_exp = K_x_initial[1] / exporter_pop_initial[1]
    Kxcorr_open_a    = K_x_open_CM_open_trade[1]/exporter_pop_open_CM_open_trade[1] * exporter_pop_store_open_CM_open_trade_asym[1,end-1] / (K_x_ratio_store_open_CM_open_trade[1,end] * capital_supply_future_open_CM_open_trade[1,end])
    Kxcorr_closed_a  = K_x_closed_CM_open_trade[1]/exporter_pop_closed_CM_open_trade[1] * exporter_pop_store_closed_CM_open_trade_asym[1,end-1] / (K_x_ratio_store_closed_CM_open_trade[1,end] * capital_supply_future_closed_CM_open_trade[1,end])
    Kxcorr_delayed_a = K_x_open_CM_open_trade[1]/exporter_pop_open_CM_open_trade[1] * exporter_pop_store_open_CMdelayed_open_trade_asym[1,end-1] / (K_x_ratio_store_open_CMdelayed_open_trade[1,end-1] * capital_supply_future_open_CMdelayed_open_trade[1,end])
    Kx_smooth_open_a    = movmean(100 .* ((Kxcorr_open_a    .* K_x_ratio_store_open_CM_open_trade_asym[1,:] .* capital_supply_future_open_CM_open_trade_asym[1,2:(end-1)]) ./ exporter_pop_store_open_CM_open_trade_asym[1,:] ./ Kx_init_per_exp .- 1), 5)
    Kx_smooth_closed_a  = movmean(100 .* ((Kxcorr_closed_a  .* K_x_ratio_store_closed_CM_open_trade_asym[1,:] .* capital_supply_future_closed_CM_open_trade_asym[1,2:(end-1)]) ./ exporter_pop_store_closed_CM_open_trade_asym[1,:] ./ Kx_init_per_exp .- 1), 5)
    Kx_smooth_delayed_a = movmean(100 .* ((Kxcorr_delayed_a .* K_x_ratio_store_open_CMdelayed_open_trade_asym[1,:] .* capital_supply_future_open_CMdelayed_open_trade_asym[1,2:(end-1)]) ./ exporter_pop_store_open_CMdelayed_open_trade_asym[1,:] ./ Kx_init_per_exp .- 1), 5)
    plot([range(1991,2020)], [Kx_smooth_open_a, Kx_smooth_closed_a, Kx_smooth_delayed_a],
         legend=nothing, linewidth=5, linestyle=[:solid :dash :dot], grid=false,
         tickfontsize=14, ylims=[-20.0,40])
    savefig("../Figures/Figure12g_asym.pdf")
    println("[fig:Figure12_asym] Saved Figure12g_asym.pdf")

    # Figure12h_asym — capital per non-exporter
    Kd_init_per_dom = K_d_initial[1] / domestic_pop_initial[1]
    Kdcorr_open_a    = K_d_open_CM_open_trade[1]/domestic_pop_open_CM_open_trade[1] * domestic_pop_store_open_CM_open_trade_asym[1,end-1] / (K_d_ratio_store_open_CM_open_trade_asym[1,end] * capital_supply_future_open_CM_open_trade_asym[1,end])
    Kdcorr_closed_a  = K_d_closed_CM_open_trade[1]/domestic_pop_closed_CM_open_trade[1] * domestic_pop_store_closed_CM_open_trade_asym[1,end-1] / (K_d_ratio_store_closed_CM_open_trade_asym[1,end] * capital_supply_future_closed_CM_open_trade_asym[1,end])
    Kdcorr_delayed_a = K_d_open_CM_open_trade[1]/domestic_pop_open_CM_open_trade[1] * domestic_pop_store_open_CMdelayed_open_trade_asym[1,end-1] / (K_d_ratio_store_open_CMdelayed_open_trade_asym[1,end-1] * capital_supply_future_open_CMdelayed_open_trade_asym[1,end])
    Kd_smooth_open_a    = movmean(100 .* ((Kdcorr_open_a    .* K_d_ratio_store_open_CM_open_trade_asym[1,:] .* capital_supply_future_open_CM_open_trade_asym[1,2:(end-1)]) ./ domestic_pop_store_open_CM_open_trade_asym[1,:] ./ Kd_init_per_dom .- 1), 5)
    Kd_smooth_closed_a  = movmean(100 .* ((Kdcorr_closed_a  .* K_d_ratio_store_closed_CM_open_trade_asym[1,:] .* capital_supply_future_closed_CM_open_trade_asym[1,2:(end-1)]) ./ domestic_pop_store_closed_CM_open_trade_asym[1,:] ./ Kd_init_per_dom .- 1), 5)
    Kd_smooth_delayed_a = movmean(100 .* ((Kdcorr_delayed_a .* K_d_ratio_store_open_CMdelayed_open_trade_asym[1,:] .* capital_supply_future_open_CMdelayed_open_trade_asym[1,2:(end-1)]) ./ domestic_pop_store_open_CMdelayed_open_trade_asym[1,:] ./ Kd_init_per_dom .- 1), 5)
    plot([range(1991,2020)], [Kd_smooth_open_a, Kd_smooth_closed_a, Kd_smooth_delayed_a],
         legend=nothing, linewidth=5, linestyle=[:solid :dash :dot], grid=false,
         tickfontsize=14, ylims=[-20.0,20])
    savefig("../Figures/Figure12h_asym.pdf")
    println("[fig:Figure12_asym] Saved Figure12h_asym.pdf")

    # Figure12i_asym — capital share of exporters
    Kshare_open_a    = movmean(100 .* K_x_ratio_store_open_CM_open_trade_asym[1,1:end], 5)
    Kshare_closed_a  = movmean(100 .* K_x_ratio_store_closed_CM_open_trade_asym[1,1:end], 5)
    Kshare_delayed_a = movmean(100 .* K_x_ratio_store_open_CMdelayed_open_trade_asym[1,1:end], 5)
    plot([range(1991,2020)], [Kshare_open_a, Kshare_closed_a, Kshare_delayed_a],
         legend=nothing, linewidth=5, linestyle=[:solid :dash :dot], grid=false,
         tickfontsize=14, ylims=[30,80])
    savefig("../Figures/Figure12i_asym.pdf")
    println("[fig:Figure12_asym] Saved Figure12i_asym.pdf")

    # Figure12l_asym — s.d. ARPK domestic
    sdMRPKd_corr_open_a    = sd_MRPK_d_open_CM_open_trade[1] / sd_MRPK_d_store_open_CM_open_trade_asym[1,end]
    sdMRPKd_corr_closed_a  = sd_MRPK_d_closed_CM_open_trade[1] / sd_MRPK_d_store_closed_CM_open_trade_asym[1,end]
    sdMRPKd_corr_delayed_a = sd_MRPK_d_open_CM_open_trade[1] / sd_MRPK_d_store_open_CMdelayed_open_trade_asym[1,end-1]
    sdMRPKd_smooth_open_a    = movmean(sdMRPKd_corr_open_a    .* sd_MRPK_d_store_open_CM_open_trade_asym[1,:], 5)
    sdMRPKd_smooth_closed_a  = movmean(sdMRPKd_corr_closed_a  .* sd_MRPK_d_store_closed_CM_open_trade_asym[1,:], 5)
    sdMRPKd_smooth_delayed_a = movmean(sdMRPKd_corr_delayed_a .* sd_MRPK_d_store_open_CMdelayed_open_trade_asym[1,:], 5)
    plot([range(1991,2020)], [sdMRPKd_smooth_open_a, sdMRPKd_smooth_closed_a, sdMRPKd_smooth_delayed_a],
         legend=nothing, linewidth=5, linestyle=[:solid :dash :dot], grid=false,
         tickfontsize=14, ylims=[0.32,0.55])
    savefig("../Figures/Figure12l_asym.pdf")
    println("[fig:Figure12_asym] Saved Figure12l_asym.pdf")

    # Figure12m_asym — s.d. ARPK exporter
    sdMRPKx_corr_open_a    = sd_MRPK_x_open_CM_open_trade[1] / sd_MRPK_x_store_open_CM_open_trade_asym[1,end]
    sdMRPKx_corr_closed_a  = sd_MRPK_x_closed_CM_open_trade[1] / sd_MRPK_x_store_closed_CM_open_trade_asym[1,end]
    sdMRPKx_corr_delayed_a = sd_MRPK_x_open_CM_open_trade[1] / sd_MRPK_x_store_open_CMdelayed_open_trade_asym[1,end-1]
    sdMRPKx_smooth_open_a    = movmean(sdMRPKx_corr_open_a    .* sd_MRPK_x_store_open_CM_open_trade_asym[1,:], 5)
    sdMRPKx_smooth_closed_a  = movmean(sdMRPKx_corr_closed_a  .* sd_MRPK_x_store_closed_CM_open_trade_asym[1,:], 5)
    sdMRPKx_smooth_delayed_a = movmean(sdMRPKx_corr_delayed_a .* sd_MRPK_x_store_open_CMdelayed_open_trade_asym[1,:], 5)
    # Note: delayed line uses domestic delayed (matching plot_asym.jl line 227)
    plot([range(1991,2020)], [sdMRPKx_smooth_open_a, sdMRPKx_smooth_closed_a, sdMRPKd_smooth_delayed_a],
         legend=nothing, linewidth=5, linestyle=[:solid :dash :dot], grid=false,
         tickfontsize=14, ylims=[0.25,0.55])
    savefig("../Figures/Figure12m_asym.pdf")
    println("[fig:Figure12_asym] Saved Figure12m_asym.pdf")
end

println("\n=== figures_tables.jl: completed sections ===")
println("  [fig:Figure1]                  Figure1.tikz (foreign credit, imports, tariffs)")
println("  [table:Table1]                 table_1_results  (8×3)")
println("  [table:Table2]                 table_2_results  (6×2)")
println("  [fig:Figure4]                  Figure4a.pdf, Figure4b.pdf")
println("  [table:Table4]                 table_4_results  (20×6)")
println("  [table:Table5]                 table_5_results  (22×5) + table_6_results (4×3)")
println("  [table:Table_regs]             BL + FC GLM regressions")
println("  [table:cap_market_basis]       table_9b (20×2), table_5b (14×2), table_6b (4×3)")
println("  [tab:calibration_fixed_cost]   table_1_results_FC (8×4)")
println("  [tab:nontargeted_fixed_cost]   table_3_results_FC (14×3)")
println("  [tab:results_fixed_cost]       table_4_results_FC (20×6)")
println("  [tab:results_fixed_cost_micro] table_5_results_FC (22×5)")
println("  [fig:combined_changes]         Figure5a, 5b, 6a, 6b, 7a, 7b")
println("  [fig:Figure_CMI_only]          Figure5c, 6c, 7c")
println("  [fig:Figure12]                 Figure12a–m (skip j)")
println("  [fig:Figure12c]                Figure12q, Figure12r  (DiD simulation)")
println("  [fig:Figure12b]                Figure12n, Figure12o, Figure12p")
println("  [fig:welfare]                  Figure9a, 9b, 10a, 10b, 11a, 11c")
println("  [table:Table8]                 Welfare change (17-row computed table)")
println("  [table:Table9]                 Standalone CM integration (computed)")
println("  [table:combined]               Financial development Table10 + Table11 (computed)")
println("  [table:Table4b]                Core economy effects (computed)")
println("  [table:Table_asym]             Asymmetric vs symmetric welfare (computed)")
println("  [fig:Figure12_cons]            Figure12c + Figure12c_asym")
println("  [fig:Figure12_TFP]             Figure12a + Figure12a_asym")
println("  [fig:Figure12_ARPK]            Figure12k + Figure12k_asym")
println("  [fig:Figure12_asym]            Figure12b/d/e/f/g/h/i/l/m_asym")