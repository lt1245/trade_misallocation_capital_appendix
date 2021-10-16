#Plotting and tables
#include("main_script.jl")
# Shut down workers
#rmprocs(2,3)
# Welfare comparison - between steady states:

Phi_fine_aug = kron(Matrix(1.0I, 3, 3),row_kron(Baseline_parameter.Country_spec_p[1][6],funbase(Baseline_parameter.fspace_a, Baseline_parameter.Country_spec_p[1][3][:,1])));
welfare_init_init =  (Phi_fine_aug * vcat(vcat(coeff_final_initial[:,1,1],coeff_final_initial[:,2,1]),coeff_final_initial[:,3,1]));
welfare_init_closed_CM_open_trade =  (Phi_fine_aug * vcat(vcat(coeff_final_closed_CM_open_trade[:,1,1],coeff_final_closed_CM_open_trade[:,2,1]),coeff_final_closed_CM_open_trade[:,3,1]));
welfare_init_open_CM_open_trade =  (Phi_fine_aug * vcat(vcat(coeff_final_open_CM_open_trade[:,1,1],coeff_final_open_CM_open_trade[:,2,1]),coeff_final_open_CM_open_trade[:,3,1]));
welfare_init_closed_CM_open_trade_stst = sum(current_distr_store_initial[:,1] .* (exp.((welfare_init_closed_CM_open_trade - welfare_init_init) * (1.0 - Baseline_parameter.β[1])) - ones(Baseline_parameter.n_fine[1] * Baseline_parameter.n_fine[2] * 3)));
welfare_init_open_CM_open_trade_stst = sum(current_distr_store_initial[:,1] .* (exp.((welfare_init_open_CM_open_trade - welfare_init_init) * (1.0 - Baseline_parameter.β[1])) - ones(Baseline_parameter.n_fine[1] * Baseline_parameter.n_fine[2] * 3)));

#load transition dynamics - for now, leave out case_final = 5,8:
#case_valid = [1,2,3,4,6,7]
#For now, simply evaluate transition.jl at different case_final values
# Transition dynamics: assumes that the transition.jl is run
Phi_fine_aug = kron(Matrix(1.0I, 3, 3),row_kron(Baseline_parameter.Country_spec_p[1][6],funbase(Baseline_parameter.fspace_a, Baseline_parameter.Country_spec_p[1][3][:,1])));
welfare_init_init =  (Phi_fine_aug * vcat(vcat(coeff_final_initial[:,1,1],coeff_final_initial[:,2,1]),coeff_final_initial[:,3,1]));
welfare_init_closed_CM_open_trade_trans =  (Phi_fine_aug * vcat(vcat(coeff_store_closed_CM_open_trade[:,1,1,2],coeff_store_closed_CM_open_trade[:,2,1,2]),coeff_store_closed_CM_open_trade[:,3,1,2]));
welfare_init_open_CM_open_trade_trans =  (Phi_fine_aug * vcat(vcat(coeff_store_open_CM_open_trade[:,1,1,2],coeff_store_open_CM_open_trade[:,2,1,2]),coeff_store_open_CM_open_trade[:,3,1,2]));
welfare_init_open_CM_closed_trade_trans =  (Phi_fine_aug * vcat(vcat(coeff_store_open_CM_closed_trade[:,1,1,2],coeff_store_open_CM_closed_trade[:,2,1,2]),coeff_store_open_CM_closed_trade[:,3,1,2]));
welfare_init_open_CMdelayed_open_trade_trans =  (Phi_fine_aug * vcat(vcat(coeff_store_open_CMdelayed_open_trade[:,1,1,2],coeff_store_open_CMdelayed_open_trade[:,2,1,2]),coeff_store_open_CMdelayed_open_trade[:,3,1,2]));
welfare_10yr_closed_CM_open_trade_trans =  (Phi_fine_aug * vcat(vcat(coeff_store_closed_CM_open_trade[:,1,1,11],coeff_store_closed_CM_open_trade[:,2,1,11]),coeff_store_closed_CM_open_trade[:,3,1,11]));
welfare_10yr_open_CM_open_trade_trans =  (Phi_fine_aug * vcat(vcat(coeff_store_open_CM_open_trade[:,1,1,11],coeff_store_open_CM_open_trade[:,2,1,11]),coeff_store_open_CM_open_trade[:,3,1,11]));
welfare_change_clCM_otrade_trans = sum(current_distr_store_initial[:,1] .* (exp.((welfare_init_closed_CM_open_trade_trans - welfare_init_init) * (1.0 - Baseline_parameter.β[1]))
- ones(Baseline_parameter.n_fine[1] * Baseline_parameter.n_fine[2] * 3)));
welfare_change_oCM_otrade_trans = sum(current_distr_store_initial[:,1] .* (exp.((welfare_init_open_CM_open_trade_trans - welfare_init_init) * (1.0 - Baseline_parameter.β[1]))
- ones(Baseline_parameter.n_fine[1] * Baseline_parameter.n_fine[2] * 3)));
welfare_change_oCM_cltrade_trans = sum(current_distr_store_initial[:,1] .* (exp.((welfare_init_open_CM_closed_trade_trans - welfare_init_init) * (1.0 - Baseline_parameter.β[1]))
- ones(Baseline_parameter.n_fine[1] * Baseline_parameter.n_fine[2] * 3)));
welfare_change_oCM_otrade_delayed10_trans = sum(current_distr_store_initial[:,1] .* (exp.((welfare_init_open_CMdelayed_open_trade_trans - welfare_init_init) * (1.0 - Baseline_parameter.β[1]))
- ones(Baseline_parameter.n_fine[1] * Baseline_parameter.n_fine[2] * 3)));

#Table 1
total_credit_open_CM_open_trade = -(domestic_firm_debt_open_CM_open_trade[1] + exporter_firm_debt_open_CM_open_trade[1])/nomGDP_open_CM_open_trade[1];
domestic_credit_open_CM_open_trade = (worker_bond_holding_open_CM_open_trade[1]+domestic_bond_holding_open_CM_open_trade[1]+exporter_bond_holding_open_CM_open_trade[1])./nomGDP_open_CM_open_trade[1]; # domestic_credit to NFC/GDP
foreign_credit_open_CM_open_trade= total_credit_open_CM_open_trade - domestic_credit_open_CM_open_trade;
foreign_credit_share_open_CM_open_trade =foreign_credit_open_CM_open_trade/total_credit_open_CM_open_trade;
table_1_results = zeros(8,3);
table_1_results[1,:] = [open_CM_open_trade_parameter.θ[1],0.62,total_credit_open_CM_open_trade]
table_1_results[2,:] = [open_CM_open_trade_parameter.β[2],0.05,prices_open_CM_open_trade[6]]
table_1_results[3,:] = [open_CM_open_trade_parameter.β[1],0.53,foreign_credit_share_open_CM_open_trade]
table_1_results[4,:] = [Baseline_parameter.τ[1],0.21,import_share_initial[1]]
table_1_results[5,:] = [open_CM_open_trade_parameter.τ[1],0.42,import_share_open_CM_open_trade[1]]
table_1_results[6,:] = [open_CM_open_trade_parameter.Fₓ[1],0.27,entry_share_to_exporter_open_CM_open_trade[1]/exporter_pop_open_CM_open_trade[1]]
table_1_results[7,:] = [open_CM_open_trade_parameter.σₛ[1],0.86,sd_growth_rev_open_CM_open_trade[1]]
table_1_results[8,:] = [open_CM_open_trade_parameter.ρ[1],0.4,autocorr_rev_open_CM_open_trade[1]]

table_2_results = zeros(6,2);
table_2_results[1,:] = [open_CM_open_trade_parameter.L[1],0]
table_2_results[2,:] = [open_CM_open_trade_parameter.L[2],0]
table_2_results[3,:] = [open_CM_open_trade_parameter.σ[1],0]
table_2_results[4,:] = [open_CM_open_trade_parameter.θ[2],0]
table_2_results[5,:] = [open_CM_open_trade_parameter.δ[1],0]
table_2_results[6,:] = [open_CM_open_trade_parameter.Fₓ[2],prices_open_CM_open_trade[4]^(1/4)/prices_open_CM_open_trade[3]^(1/4)/prices_open_CM_open_trade[1]*prices_open_CM_open_trade[2] * prices_open_CM_open_trade[5]]

table_3_results = zeros(13,2);
table_3_results[1,:] = [1.06,sd_MRPK_open_CM_open_trade[1]]
table_3_results[2,:] = [0.72,sd_growth_k_open_CM_open_trade[1]]
table_3_results[3,:] = [0.38,exporter_pop_open_CM_open_trade[1]/(domestic_pop_open_CM_open_trade[1] + exporter_pop_open_CM_open_trade[1])]
table_3_results[4,:] = [0.39,exporter_firm_debt_open_CM_open_trade[1]/(domestic_firm_debt_open_CM_open_trade[1] + exporter_firm_debt_open_CM_open_trade[1])]
table_3_results[5,:] = [0.67,mean_leverage_open_CM_open_trade[1]]
table_3_results[6,:] = [0.56,mean_leverage_x_open_CM_open_trade[1]]
table_3_results[7,:] = [0.02,fraction_zombie_exporter_open_CM_open_trade[1]]
table_3_results[8,:] = [0.2,prices_open_CM_open_trade[3]/prices_open_CM_open_trade[4]/ open_CM_open_trade_parameter.L[1]*open_CM_open_trade_parameter.L[2]]
table_3_results[9,:] = [0.53,p90_wealth_open_CM_open_trade[1]]
table_3_results[10,:] = [0.34,p90_income_open_CM_open_trade[1]]
table_3_results[11,:] = [0.11,p99_income_open_CM_open_trade[1]]
table_3_results[12,:] = [0.24,p90_income_initial[1]]
table_3_results[13,:] = [0.06,p99_income_initial[1]]
#Table 4
table_4_results = zeros(20,3);
table_4_results[1,:] = [TFP_initial[1]./TFP_initial[1],TFP_closed_CM_open_trade[1]./TFP_initial[1],TFP_open_CM_open_trade[1]./TFP_initial[1]];
table_4_results[2,:] = [sd_MRPK_initial[1],sd_MRPK_closed_CM_open_trade[1],sd_MRPK_open_CM_open_trade[1]];
table_4_results[3,:] = [GDP_initial[1]./GDP_initial[1],GDP_closed_CM_open_trade[1]./GDP_initial[1],GDP_open_CM_open_trade[1]./GDP_initial[1]];
table_4_results[4,:] = [prices_initial[3]./prices_initial[3],prices_closed_CM_open_trade[3]./prices_initial[3],prices_open_CM_open_trade[3]./prices_initial[3]];
table_4_results[5,:] = [total_consumption_initial[1]./total_consumption_initial[1],total_consumption_closed_CM_open_trade[1]./total_consumption_initial[1],total_consumption_open_CM_open_trade[1]./total_consumption_initial[1]];
table_4_results[6,:] = [capital_demand_initial[1]./capital_demand_initial[1],capital_demand_closed_CM_open_trade[1]./capital_demand_initial[1],capital_demand_open_CM_open_trade[1]./capital_demand_initial[1]];
table_4_results[8,:] = [0,welfare_change_clCM_otrade_trans,welfare_change_oCM_otrade_trans];
table_4_results[9,:] = [p90_wealth_initial[1],p90_wealth_closed_CM_open_trade[1],p90_wealth_open_CM_open_trade[1]];
table_4_results[10,:] = [prices_initial[1]/prices_initial[5]/(prices_initial[1]/prices_initial[5]),prices_closed_CM_open_trade[1]/prices_closed_CM_open_trade[5]/(prices_initial[1]/prices_initial[5]),prices_open_CM_open_trade[1]/prices_open_CM_open_trade[5]/(prices_initial[1]/prices_initial[5])];
table_4_results[11,:] = [prices_initial[6]-prices_initial[7],prices_closed_CM_open_trade[6]-prices_closed_CM_open_trade[7],prices_open_CM_open_trade[6]-prices_open_CM_open_trade[6]];
table_4_results[12,:] = [import_share_initial[1],import_share_closed_CM_open_trade[1],import_share_open_CM_open_trade[1]];
table_4_results[13,:] = [export_value_initial[1]./nomGDP_initial[2],export_value_closed_CM_open_trade[1]./nomGDP_closed_CM_open_trade[2],export_value_open_CM_open_trade[1]./nomGDP_open_CM_open_trade[2]];
table_4_results[14,:] = [(domestic_pop_initial[1]+ exporter_pop_initial[1]),(domestic_pop_closed_CM_open_trade[1]+ exporter_pop_closed_CM_open_trade[1]),(domestic_pop_open_CM_open_trade[1]+ exporter_pop_open_CM_open_trade[1])];
table_4_results[15,:] = [exporter_pop_initial[1]./(domestic_pop_initial[1]+ exporter_pop_initial[1]),exporter_pop_closed_CM_open_trade[1]./(domestic_pop_closed_CM_open_trade[1]+ exporter_pop_closed_CM_open_trade[1]),exporter_pop_open_CM_open_trade[1]./(domestic_pop_open_CM_open_trade[1]+ exporter_pop_open_CM_open_trade[1])];
table_4_results[16,:] = [prices_initial[5],prices_closed_CM_open_trade[5],prices_open_CM_open_trade[5]];

total_credit_initial = -(domestic_firm_debt_initial[1] + exporter_firm_debt_initial[1])/nomGDP_initial[1];
domestic_credit_initial = (worker_bond_holding_initial[1]+domestic_bond_holding_initial[1]+exporter_bond_holding_initial[1])./nomGDP_initial[1]; # domestic_credit to NFC/GDP
foreign_credit_initial= total_credit_initial - domestic_credit_initial;
foreign_credit_share_initial =foreign_credit_initial/total_credit_initial;

total_credit_closed_CM_open_trade = -(domestic_firm_debt_closed_CM_open_trade[1] + exporter_firm_debt_closed_CM_open_trade[1])/nomGDP_closed_CM_open_trade[1];
domestic_credit_closed_CM_open_trade = (worker_bond_holding_closed_CM_open_trade[1]+domestic_bond_holding_closed_CM_open_trade[1]+exporter_bond_holding_closed_CM_open_trade[1])./nomGDP_closed_CM_open_trade[1]; # domestic_credit to NFC/GDP
foreign_credit_closed_CM_open_trade= total_credit_closed_CM_open_trade - domestic_credit_closed_CM_open_trade;
foreign_credit_share_closed_CM_open_trade =foreign_credit_closed_CM_open_trade/total_credit_closed_CM_open_trade;
table_4_results[17,:] = [total_credit_initial,total_credit_closed_CM_open_trade,total_credit_open_CM_open_trade];
table_4_results[18,:] = [foreign_credit_share_initial,foreign_credit_share_closed_CM_open_trade,foreign_credit_share_open_CM_open_trade];
# for now, only steady state welfare. Update it later

table_5_results = zeros(7,3);
table_5_results[1,:] = [sd_MRPK_d_initial[1],sd_MRPK_d_closed_CM_open_trade[1],sd_MRPK_d_open_CM_open_trade[1]];
table_5_results[2,:] = [sd_MRPK_x_initial[1],sd_MRPK_x_closed_CM_open_trade[1],sd_MRPK_x_open_CM_open_trade[1]];
table_5_results[3,:] = [Misallocation_within_d_initial[1],Misallocation_within_d_closed_CM_open_trade[1],Misallocation_within_d_open_CM_open_trade[1]];
table_5_results[4,:] = [Misallocation_within_x_initial[1],Misallocation_within_x_closed_CM_open_trade[1],Misallocation_within_x_open_CM_open_trade[1]];
table_5_results[5,:] = [exporter_pop_initial[1]./(domestic_pop_initial[1]+ exporter_pop_initial[1]),exporter_pop_closed_CM_open_trade[1]./(domestic_pop_closed_CM_open_trade[1]+ exporter_pop_closed_CM_open_trade[1]),exporter_pop_open_CM_open_trade[1]./(domestic_pop_open_CM_open_trade[1]+ exporter_pop_open_CM_open_trade[1])];
table_5_results[6,:] = [fraction_zombie_exporter_initial[1],fraction_zombie_exporter_closed_CM_open_trade[1],fraction_zombie_exporter_open_CM_open_trade[1]];
table_5_results[7,:] = [exporter_pop_initial[1]/entry_share_to_exporter_initial[1],exporter_pop_closed_CM_open_trade[1]/entry_share_to_exporter_closed_CM_open_trade[1],exporter_pop_open_CM_open_trade[1]/entry_share_to_exporter_open_CM_open_trade[1]];
# Selecting unproductive and wealthy exporters:
table_6_results = zeros(4,3);
current_exporter_initial = current_distr_store_initial[convert(Int64,2/3*size(current_distr_store_initial)[1]+1):end,1];
current_exporter_closed_CM_open_trade = current_distr_store_closed_CM_open_trade[convert(Int64,2/3*size(current_distr_store_initial)[1]+1):end,1];
current_exporter_open_CM_open_trade = current_distr_store_open_CM_open_trade[convert(Int64,2/3*size(current_distr_store_initial)[1]+1):end,1];

function mat_creator(y_value)
    y_mat = zeros(36,200);
    for x_index = 1:36
        for y_index=1:200
            y_mat[x_index,y_index] = y_value[200*(x_index-1) + y_index];
        end
    end
    return y_mat
end
current_exporter_initial_mat = mat_creator(current_exporter_initial );
x_start = 18;
y_start = 8;
table_6_results[4,1] = sum(current_exporter_initial_mat[(x_start+1):end,(y_start+1):end])/sum(current_exporter_initial);
table_6_results[3,1] = sum(current_exporter_initial_mat[1:x_start,(y_start+1):end])/sum(current_exporter_initial);
table_6_results[2,1] = sum(current_exporter_initial_mat[(x_start+1):end,1:y_start])/sum(current_exporter_initial);
table_6_results[1,1] =sum(current_exporter_initial_mat[1:x_start,1:y_start])/sum(current_exporter_initial);
current_exporter_closed_CM_open_trade_mat = mat_creator(current_exporter_closed_CM_open_trade );
table_6_results[4,2] = sum(current_exporter_closed_CM_open_trade_mat[(x_start+1):end,(y_start+1):end])/sum(current_exporter_closed_CM_open_trade);
table_6_results[3,2] = sum(current_exporter_closed_CM_open_trade_mat[1:x_start,(y_start+1):end])/sum(current_exporter_closed_CM_open_trade);
table_6_results[2,2] = sum(current_exporter_closed_CM_open_trade_mat[(x_start+1):end,1:y_start])/sum(current_exporter_closed_CM_open_trade);
table_6_results[1,2] =sum(current_exporter_closed_CM_open_trade_mat[1:x_start,1:y_start])/sum(current_exporter_closed_CM_open_trade);
current_exporter_open_CM_open_trade_mat = mat_creator(current_exporter_open_CM_open_trade );
table_6_results[4,3] = sum(current_exporter_open_CM_open_trade_mat[(x_start+1):end,(y_start+1):end])/sum(current_exporter_open_CM_open_trade);
table_6_results[3,3] = sum(current_exporter_open_CM_open_trade_mat[1:x_start,(y_start+1):end])/sum(current_exporter_open_CM_open_trade);
table_6_results[2,3] = sum(current_exporter_open_CM_open_trade_mat[(x_start+1):end,1:y_start])/sum(current_exporter_open_CM_open_trade);
table_6_results[1,3] =sum(current_exporter_open_CM_open_trade_mat[1:x_start,1:y_start])/sum(current_exporter_open_CM_open_trade);
# Experimental productivity decomposition
# COnstruct the counterfactual with frictionless average productivity:

#TFP_second_best_initial = (TFP_d_closed_CM_open_trade_fl .* domestic_pop_initial ./domestic_pop_closed_CM_open_trade_fl .* K_d_ratio_eff_closed_CM_open_trade.^Baseline_parameter.α₁_eff .* L_d_ratio_eff_closed_CM_open_trade.^Baseline_parameter.α₂_eff
#    + TFP_x_closed_CM_open_trade_fl .* exporter_pop_initial ./exporter_pop_closed_CM_open_trade_fl.* K_x_ratio_eff_closed_CM_open_trade.^Baseline_parameter.α₁_eff .* L_x_ratio_eff_closed_CM_open_trade.^Baseline_parameter.α₂_eff).^(Baseline_parameter.σ/(Baseline_parameter.σ-1));
#TFP_second_best_open_CM_open_trade = (TFP_d_closed_CM_open_trade_fl .* domestic_pop_open_CM_open_trade ./domestic_pop_closed_CM_open_trade_fl .* K_d_ratio_eff_closed_CM_open_trade.^Baseline_parameter.α₁_eff .* L_d_ratio_eff_closed_CM_open_trade.^Baseline_parameter.α₂_eff
#        + TFP_x_closed_CM_open_trade_fl .* exporter_pop_open_CM_open_trade ./exporter_pop_closed_CM_open_trade_fl .* K_x_ratio_eff_closed_CM_open_trade.^Baseline_parameter.α₁_eff .* L_x_ratio_eff_closed_CM_open_trade.^Baseline_parameter.α₂_eff).^(Baseline_parameter.σ/(Baseline_parameter.σ-1));
#TFP_second_best_closed_CM_open_trade = (TFP_d_closed_CM_open_trade_fl .* domestic_pop_closed_CM_open_trade ./domestic_pop_closed_CM_open_trade_fl .* K_d_ratio_eff_closed_CM_open_trade.^Baseline_parameter.α₁_eff .* L_d_ratio_eff_closed_CM_open_trade.^Baseline_parameter.α₂_eff
#                + TFP_x_closed_CM_open_trade_fl .* exporter_pop_closed_CM_open_trade ./exporter_pop_closed_CM_open_trade_fl .* K_x_ratio_eff_closed_CM_open_trade.^Baseline_parameter.α₁_eff .* L_x_ratio_eff_closed_CM_open_trade.^Baseline_parameter.α₂_eff).^(Baseline_parameter.σ/(Baseline_parameter.σ-1));
table_7_results = zeros(4,3);
total_loss_initial = (1- TFP_initial[1]/TFP_second_best_open_CM_open_trade_dev[1]);
total_loss_openCM_open_trade = (1- TFP_open_CM_open_trade[1]/TFP_second_best_open_CM_open_trade_dev[1]);
total_loss_closedCM_open_trade = (1- TFP_closed_CM_open_trade[1]/TFP_second_best_open_CM_open_trade_dev[1]);
within_loss_initial = (TFP_within_initial[1] - TFP_initial[1])/TFP_second_best_open_CM_open_trade_dev[1];
within_loss_openCM_open_trade = (TFP_within_open_CM_open_trade[1] - TFP_open_CM_open_trade[1])/TFP_second_best_open_CM_open_trade_dev[1];
within_loss_closedCM_open_trade = (TFP_within_closed_CM_open_trade[1] - TFP_closed_CM_open_trade[1])/TFP_second_best_open_CM_open_trade_dev[1];
factor_loss_initial = (TFP_second_best_initial_dev[1]- TFP_within_initial[1])/TFP_second_best_open_CM_open_trade_dev[1];
factor_loss_openCM_open_trade = (TFP_second_best_open_CM_open_trade_dev[1]- TFP_within_open_CM_open_trade[1])/TFP_second_best_open_CM_open_trade_dev[1];
factor_loss_closedCM_open_trade = (TFP_second_best_closed_CM_open_trade_dev[1] - TFP_within_closed_CM_open_trade[1])/TFP_second_best_open_CM_open_trade_dev[1];
reform_loss_initial = (TFP_second_best_open_CM_open_trade_dev[1]- TFP_second_best_initial_dev[1])/TFP_second_best_open_CM_open_trade_dev[1];
reform_loss_openCM_open_trade = (TFP_second_best_open_CM_open_trade_dev[1]- TFP_second_best_open_CM_open_trade_dev[1])/TFP_second_best_open_CM_open_trade_dev[1];
reform_loss_closedCM_open_trade = (TFP_second_best_open_CM_open_trade_dev[1] - TFP_second_best_closed_CM_open_trade_dev[1])/TFP_second_best_open_CM_open_trade_dev[1];

table_7_results[1,1] = within_loss_initial/total_loss_initial;
table_7_results[1,2] = within_loss_closedCM_open_trade/total_loss_closedCM_open_trade;
table_7_results[1,3] = within_loss_openCM_open_trade/total_loss_openCM_open_trade;
table_7_results[2,1] = factor_loss_initial/total_loss_initial;
table_7_results[2,2] = factor_loss_closedCM_open_trade/total_loss_closedCM_open_trade;
table_7_results[2,3] = factor_loss_openCM_open_trade/total_loss_openCM_open_trade;
table_7_results[3,1] = reform_loss_initial/total_loss_initial;
table_7_results[3,2] = reform_loss_closedCM_open_trade/total_loss_closedCM_open_trade;
table_7_results[3,3] = reform_loss_openCM_open_trade/total_loss_openCM_open_trade;
table_7_results[4,1] = total_loss_initial;
table_7_results[4,2] = total_loss_closedCM_open_trade;
table_7_results[4,3] = total_loss_openCM_open_trade;
table_7_results = 100 * table_7_results;
# Figures
# Figure 4 Decisions
using Plots; plotlyjs()
y_value =  k_x_fine_open_CM_open_trade[:,1]./k_opt_x_fine_open_CM_open_trade[:,1];
y_mat = mat_creator(y_value);
y_mat_applied = 100.0 * y_mat[10:36,1:100];
contour(y_mat_applied', fill=true,levels = 5, c =cgrad(:Blues_6, rev = true),xlabel = "Productivity",ylabel = "Wealth",axis=nothing,tickfontsize = 14,xguidefontsize=18,yguidefontsize=18,legendfontsize=18)
#;colorbar_cticks  = [0,20,40,60,80])
savefig("Figure4a.pdf")

y_value = future_occupation_fine_open_CM_open_trade[:,3,1,1];
y_mat = mat_creator(y_value);
y_mat_applied = y_mat[1:36,1:30];
contour(y_mat_applied', fill=true,levels = 1, c =cgrad(:Blues_3, rev = false),xlabel = "Productivity",ylabel = "Wealth",axis=nothing,colorbar=nothing,tickfontsize = 14,xguidefontsize=18,yguidefontsize=18,legendfontsize=18)
savefig("Figure4b.pdf")

y_value =  k_x_fine_open_CM_open_trade[:,1]./k_x_fine_initial[:,1];
y_mat = mat_creator(y_value);
y_mat_applied = 100.0 * y_mat[10:36,1:100] .-100 ;
contour(y_mat_applied', fill=true,levels = 5, c =cgrad(:Blues_6, rev = true),xlabel = "Productivity",ylabel = "Wealth",axis=nothing,tickfontsize = 14,xguidefontsize=18,yguidefontsize=18,legendfontsize=18)
#;colorbar_cticks  = [0,20,40,60,80])
savefig("Figure5a.pdf")

y_value =  k_x_fine_closed_CM_open_trade[:,1]./k_x_fine_initial[:,1];
y_mat = mat_creator(y_value);
y_mat_applied = 100.0 * y_mat[10:36,1:100] .-100 ;
contour(y_mat_applied', fill=true,levels = 5, c =cgrad(:Blues_6, rev = true),xlabel = "Productivity",ylabel = "Wealth",axis=nothing,tickfontsize = 14,xguidefontsize=18,yguidefontsize=18,legendfontsize=18)
#;colorbar_cticks  = [0,20,40,60,80])
savefig("Figure5b.pdf")
y_value =  Profit_fine_x_open_CM_open_trade[:,1]./Profit_fine_x_initial[:,1];
y_mat = mat_creator(y_value);
y_mat_applied = 100.0 * y_mat[10:36,1:100] .-100 ;
contour(y_mat_applied', fill=true,levels = 5, c =cgrad(:Blues_6, rev = true),xlabel = "Productivity",ylabel = "Wealth",axis=nothing,tickfontsize = 14,xguidefontsize=18,yguidefontsize=18,legendfontsize=18)
#;colorbar_cticks  = [0,20,40,60,80])
savefig("Figure6a.pdf")

y_value =  Profit_fine_x_closed_CM_open_trade[:,1]./Profit_fine_x_initial[:,1];
y_mat = mat_creator(y_value);
y_mat_applied = 100.0 * y_mat[10:36,1:100] .-100 ;
contour(y_mat_applied', fill=true,levels = 5, c =cgrad(:Blues_6, rev = true),xlabel = "Productivity",ylabel = "Wealth",axis=nothing,tickfontsize = 14,xguidefontsize=18,yguidefontsize=18,legendfontsize=18)
#;colorbar_cticks  = [0,20,40,60,80])
savefig("Figure6b.pdf")

y_value = 2*(future_occupation_fine_open_CM_open_trade[:,3,1,1].-1) + (future_occupation_fine_initial[:,3,1,1] .-1) ;
y_mat = mat_creator(y_value);
y_mat_applied = y_mat[1:36,1:30];
contour(y_mat_applied', fill=true,levels = 3, c =cgrad(:Blues_3, rev = false),xlabel = "Productivity",ylabel = "Wealth",axis=nothing,colorbar=nothing,tickfontsize = 14,xguidefontsize=18,yguidefontsize=18,legendfontsize=18)
savefig("Figure7a.pdf")

y_value = 2*(future_occupation_fine_closed_CM_open_trade[:,3,1,1].-1) + (future_occupation_fine_initial[:,3,1,1] .-1) ;
y_mat = mat_creator(y_value);
y_mat_applied = y_mat[1:36,1:30];
contour(y_mat_applied', fill=true,levels = 3, c =cgrad(:Blues_3, rev = false),xlabel = "Productivity",ylabel = "Wealth",axis=nothing,colorbar=nothing,tickfontsize = 14,xguidefontsize=18,yguidefontsize=18,legendfontsize=18)
savefig("Figure7b.pdf")

# Welfare plots
y_value_pop = (exp.((welfare_init_open_CM_open_trade - welfare_init_init) * (1.0 - Baseline_parameter.β[1])) - ones(Baseline_parameter.n_fine[1] * Baseline_parameter.n_fine[2] * 3));
y_value = y_value_pop[1:7200]
y_mat = mat_creator(y_value); #Workers
y_mat_applied = 100.0 * y_mat;
contour(y_mat_applied', fill=true,levels = 5, c =cgrad(:Blues_6, rev = true),xlabel = "Productivity",ylabel = "Wealth",axis=nothing,tickfontsize = 18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18)
savefig("Figure8a.pdf")
y_value = y_value_pop[7201:14400]
y_mat = mat_creator(y_value); #Domestic producers
y_mat_applied = 100.0 * y_mat;
contour(y_mat_applied', fill=true,levels = 5, c =cgrad(:Blues_6, rev = true),xlabel = "Productivity",ylabel = "Wealth",axis=nothing,tickfontsize = 18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18)
savefig("Figure9a.pdf")
y_value = y_value_pop[14401:end]
y_mat = mat_creator(y_value); #Exporters
y_mat_applied = 100.0 * y_mat;
contour(y_mat_applied', fill=true,levels = 5, c =cgrad(:Blues_6, rev = true),xlabel = "Productivity",ylabel = "Wealth",axis=nothing,tickfontsize = 18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18)
savefig("Figure10a.pdf")

y_value_pop = (exp.((welfare_init_closed_CM_open_trade - welfare_init_init) * (1.0 - Baseline_parameter.β[1])) - ones(Baseline_parameter.n_fine[1] * Baseline_parameter.n_fine[2] * 3));
y_value = y_value_pop[1:7200]
y_mat = mat_creator(y_value); #Workers
y_mat_applied = 100.0 * y_mat;
contour(y_mat_applied', fill=true,levels = 5, c =cgrad(:Blues_6, rev = true),xlabel = "Productivity",ylabel = "Wealth",axis=nothing,tickfontsize = 18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18)
savefig("Figure8b.pdf")
y_value = y_value_pop[7201:14400]
y_mat = mat_creator(y_value); #Domestic producers
y_mat_applied = 100.0 * y_mat;
contour(y_mat_applied', fill=true,levels = 5, c =cgrad(:Blues_6, rev = true),xlabel = "Productivity",ylabel = "Wealth",axis=nothing,tickfontsize = 18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18)
savefig("Figure9b.pdf")
y_value = y_value_pop[14401:end]
y_mat = mat_creator(y_value); #Exporters
y_mat_applied = 100.0 * y_mat;
contour(y_mat_applied', fill=true,levels = 5, c =cgrad(:Blues_6, rev = true),xlabel = "Productivity",ylabel = "Wealth",axis=nothing,tickfontsize = 18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18)
savefig("Figure10b.pdf")
y_value_pop = (exp.((welfare_init_open_CM_open_trade - welfare_init_closed_CM_open_trade) * (1.0 - Baseline_parameter.β[1])) - ones(Baseline_parameter.n_fine[1] * Baseline_parameter.n_fine[2] * 3));
y_value = y_value_pop[1:7200]
y_mat = mat_creator(y_value); #Workers
y_mat_applied = 100.0 * y_mat;
contour(y_mat_applied', fill=true,levels = 5, c =cgrad(:Blues_6, rev = true),xlabel = "Productivity",ylabel = "Wealth",axis=nothing,tickfontsize = 18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18)
savefig("Figure11a.pdf")
y_value = y_value_pop[7201:14400]
y_mat = mat_creator(y_value); #Domestic producers
y_mat_applied = 100.0 * y_mat;
contour(y_mat_applied', fill=true,levels = 5, c =cgrad(:Blues_6, rev = true),xlabel = "Productivity",ylabel = "Wealth",axis=nothing,tickfontsize = 18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18)
savefig("Figure11b.pdf")
y_value = y_value_pop[14401:end]
y_mat = mat_creator(y_value); #Exporters
y_mat_applied = 100.0 * y_mat;
contour(y_mat_applied', fill=true,levels = 5, c =cgrad(:Blues_6, rev = true),xlabel = "Productivity",ylabel = "Wealth",axis=nothing,colorbar_tickfontsize = 40,colorbar_titlefontsize= 40,guidefontsize=18,legendfontsize=18, titlefontsize = 30)
savefig("Figure11c.pdf")
# Welfare table 8
table_8_results = zeros(11,3);
table_8_results[1,:] = [GDP_initial[1]./GDP_initial[1],GDP_closed_CM_open_trade[1]./GDP_initial[1],GDP_open_CM_open_trade[1]./GDP_initial[1]];
table_8_results[2,:] = [total_consumption_initial[1]./total_consumption_initial[1],total_consumption_closed_CM_open_trade[1]./total_consumption_initial[1],total_consumption_open_CM_open_trade[1]./total_consumption_initial[1]];
table_8_results[3,:] = [capital_demand_initial[1]./capital_demand_initial[1],capital_demand_closed_CM_open_trade[1]./capital_demand_initial[1],capital_demand_open_CM_open_trade[1]./capital_demand_initial[1]];
table_8_results[4,:] = [0,welfare_init_closed_CM_open_trade_stst,welfare_init_open_CM_open_trade_stst];
table_8_results[5,:] = [0,welfare_change_clCM_otrade_trans,welfare_change_oCM_otrade_trans];
table_8_results[6,:] = [p90_wealth_initial[1],p90_wealth_closed_CM_open_trade[1],p90_wealth_open_CM_open_trade[1]];
table_8_results[7,:] = [p90_income_initial[1],p90_income_closed_CM_open_trade[1],p90_income_open_CM_open_trade[1]];
table_8_results[8,:] = [p90_cons_initial[1],p90_cons_closed_CM_open_trade[1],p90_cons_open_CM_open_trade[1]];
table_8_results[9,:] = [prices_initial[1]/prices_initial[5]/(prices_initial[1]/prices_initial[5]),prices_closed_CM_open_trade[1]/prices_closed_CM_open_trade[5]/(prices_initial[1]/prices_initial[5]),prices_open_CM_open_trade[1]/prices_open_CM_open_trade[5]/(prices_initial[1]/prices_initial[5])];
table_8_results[10,:] = [prices_initial[6]-prices_initial[7],prices_closed_CM_open_trade[6]-prices_closed_CM_open_trade[7],prices_open_CM_open_trade[6]-prices_open_CM_open_trade[6]];
table_8_results[11,:] = [wealth_of_exporters_initial[1],wealth_of_exporters_closed_CM_open_trade[1],wealth_of_exporters_open_CM_open_trade[1]];


#welfare_change_02_10yr_init_distr_trans = sum(current_distr_store_initial[:,1] .* (exp.( (1.0 - Baseline_parameter.β[1])* (welfare_10yr_open_CM_open_trade_trans * Baseline_parameter.β[1]^10
#- welfare_init_init  + welfare_init_closed_CM_open_trade_trans - welfare_10yr_closed_CM_open_trade_trans* Baseline_parameter.β[1]^10 )) .-1));
# Now true delayed
# Now the figures - smooth them for plotting
function movmean(array::Array{Float64,1},window::Int64)
    array_size = size(array)[1];
    array_smooth = zeros(array_size);
    for i = 1:array_size
        if i<window
            array_smooth[i] = sum(array[1:i])/i
        else
            array_smooth[i] = sum(array[(i-window+1):i])/window
        end
    end
    return array_smooth
end
TFP_open_CM_open_trade_plotcorrection = TFP_open_CM_open_trade[1]./ TFP_store_open_CM_open_trade[1,end];
TFP_closed_CM_open_trade_plotcorrection = TFP_closed_CM_open_trade[1]./ TFP_store_closed_CM_open_trade[1,end-1];
TFP_open_CMdelayed_open_trade_plotcorrection = TFP_open_CM_open_trade[1]./ TFP_store_open_CMdelayed_open_trade[1,end-1];

TFP_smooth_open =100 * (TFP_open_CM_open_trade_plotcorrection * TFP_store_open_CM_open_trade[1,:]/TFP_initial[1].-1);
TFP_smooth_closed = 100 * (TFP_closed_CM_open_trade_plotcorrection * TFP_store_closed_CM_open_trade[1,:]/TFP_initial[1].-1);
TFP_smooth_open_delayed =100 * (TFP_open_CMdelayed_open_trade_plotcorrection * TFP_store_open_CMdelayed_open_trade[1,:]/TFP_initial[1].-1);
TFP_smooth_open_delayed= movmean(TFP_smooth_open_delayed,5);
TFP_smooth_open = movmean(TFP_smooth_open,5);
TFP_smooth_closed = movmean(TFP_smooth_closed,5);
plot([TFP_smooth_open,TFP_smooth_closed,TFP_smooth_open_delayed],legend = nothing, linewidth = 3,linestyle = [:solid :dash :dot],grid = false,tickfontsize = 14,ylims = [0,35])
savefig("Figure12a.pdf")
total_production_open_CM_open_trade_plotcorrection = total_production_open_CM_open_trade[1]./ total_production_store_open_CM_open_trade[1,end];
total_production_closed_CM_open_trade_plotcorrection = total_production_closed_CM_open_trade[1]./ total_production_store_closed_CM_open_trade[1,end];
total_production_open_CMdelayed_open_trade_plotcorrection = total_production_open_CM_open_trade[1]./ total_production_store_open_CMdelayed_open_trade[1,end-1];
total_production_smooth_open_delayed =100 * (total_production_open_CMdelayed_open_trade_plotcorrection * total_production_store_open_CMdelayed_open_trade[1,:]/total_production_initial[1].-1);
total_production_smooth_open_delayed= movmean(total_production_smooth_open_delayed,5);
total_production_smooth_open = 100 * (total_production_open_CM_open_trade_plotcorrection*total_production_store_open_CM_open_trade[1,:]/total_production_initial[1].-1);
total_production_smooth_closed = 100 * (total_production_closed_CM_open_trade_plotcorrection*total_production_store_closed_CM_open_trade[1,:]/total_production_initial[1].-1);
total_production_smooth_open = movmean(total_production_smooth_open,5);
total_production_smooth_closed = movmean(total_production_smooth_closed,5);
plot([total_production_smooth_open,total_production_smooth_closed,total_production_smooth_open_delayed],legend = nothing, linewidth = 3,linestyle = [:solid :dash :dot],grid = false,tickfontsize = 14,ylims = [0,45])
savefig("Figure12b.pdf")
fraction_zombie_exporter_open_CMdelayed_open_trade_plotcorrection = fraction_zombie_exporter_open_CM_open_trade[1]./ fraction_zombie_exporter_store_open_CMdelayed_open_trade[1,end-1];
fraction_zombie_exporter_smooth_open_delayed =100 * (fraction_zombie_exporter_open_CMdelayed_open_trade_plotcorrection * fraction_zombie_exporter_store_open_CMdelayed_open_trade[1,:]);
fraction_zombie_exporter_smooth_open_delayed= movmean(fraction_zombie_exporter_smooth_open_delayed,5);
fraction_zombie_exporter_open_CM_open_trade_plotcorrection = fraction_zombie_exporter_open_CM_open_trade[1]./ fraction_zombie_exporter_store_open_CM_open_trade[1,end];
fraction_zombie_exporter_closed_CM_open_trade_plotcorrection = fraction_zombie_exporter_closed_CM_open_trade[1]./ fraction_zombie_exporter_store_closed_CM_open_trade[1,end];
fraction_zombie_exporter_smooth_open = 100 * (fraction_zombie_exporter_open_CM_open_trade_plotcorrection*fraction_zombie_exporter_store_open_CM_open_trade[1,:]);
fraction_zombie_exporter_smooth_closed = 100 * (fraction_zombie_exporter_closed_CM_open_trade_plotcorrection*fraction_zombie_exporter_store_closed_CM_open_trade[1,:]);
fraction_zombie_exporter_smooth_open = movmean(fraction_zombie_exporter_smooth_open,5);
fraction_zombie_exporter_smooth_closed = movmean(fraction_zombie_exporter_smooth_closed,5);
plot([fraction_zombie_exporter_smooth_open,fraction_zombie_exporter_smooth_closed,fraction_zombie_exporter_smooth_open_delayed],legend = nothing, linewidth = 3,linestyle = [:solid :dash :dot],grid = false,tickfontsize = 14,ylims = [0,16])
savefig("Figure12c.pdf")
sd_MRPK_open_CMdelayed_open_trade_plotcorrection = sd_MRPK_open_CM_open_trade[1]./ sd_MRPK_store_open_CMdelayed_open_trade[1,end-1];
sd_MRPK_smooth_open_delayed = (sd_MRPK_open_CMdelayed_open_trade_plotcorrection * sd_MRPK_store_open_CMdelayed_open_trade[1,:]);
sd_MRPK_smooth_open_delayed= movmean(sd_MRPK_smooth_open_delayed,5);
sd_MRPK_open_CM_open_trade_plotcorrection = sd_MRPK_open_CM_open_trade[1]./ sd_MRPK_store_open_CM_open_trade[1,end];
sd_MRPK_closed_CM_open_trade_plotcorrection = sd_MRPK_closed_CM_open_trade[1]./ sd_MRPK_store_closed_CM_open_trade[1,end];
sd_MRPK_smooth_open = (sd_MRPK_open_CM_open_trade_plotcorrection * sd_MRPK_store_open_CM_open_trade[1,:]);
sd_MRPK_smooth_closed =  (sd_MRPK_closed_CM_open_trade_plotcorrection * sd_MRPK_store_closed_CM_open_trade[1,:]);
sd_MRPK_smooth_open = movmean(sd_MRPK_smooth_open,5);
sd_MRPK_smooth_closed = movmean(sd_MRPK_smooth_closed,5);
plot([sd_MRPK_smooth_open,sd_MRPK_smooth_closed,sd_MRPK_smooth_open_delayed],legend = nothing, linewidth =3,linestyle = [:solid :dash :dot],grid = false,tickfontsize = 14,ylims = [0,0.55])
savefig("Figure12d.pdf")
sd_MRPK_d_open_CMdelayed_open_trade_plotcorrection = sd_MRPK_d_open_CM_open_trade[1]./ sd_MRPK_d_store_open_CMdelayed_open_trade[1,end-1];
sd_MRPK_d_smooth_open_delayed = (sd_MRPK_d_open_CMdelayed_open_trade_plotcorrection * sd_MRPK_d_store_open_CMdelayed_open_trade[1,:]);
sd_MRPK_d_smooth_open_delayed= movmean(sd_MRPK_d_smooth_open_delayed,5);
sd_MRPK_d_open_CM_open_trade_plotcorrection = sd_MRPK_d_open_CM_open_trade[1]./ sd_MRPK_d_store_open_CM_open_trade[1,end];
sd_MRPK_d_closed_CM_open_trade_plotcorrection = sd_MRPK_d_closed_CM_open_trade[1]./ sd_MRPK_d_store_closed_CM_open_trade[1,end];
sd_MRPK_d_smooth_open = (sd_MRPK_d_open_CM_open_trade_plotcorrection*sd_MRPK_d_store_open_CM_open_trade[1,:]);
sd_MRPK_d_smooth_closed =  (sd_MRPK_d_closed_CM_open_trade_plotcorrection*sd_MRPK_d_store_closed_CM_open_trade[1,:]);
sd_MRPK_d_smooth_open = movmean(sd_MRPK_d_smooth_open,5);
sd_MRPK_d_smooth_closed = movmean(sd_MRPK_d_smooth_closed,5);
plot([sd_MRPK_d_smooth_open,sd_MRPK_d_smooth_closed,sd_MRPK_d_smooth_open_delayed],legend = nothing, linewidth = 3,linestyle = [:solid :dash :dot],grid = false,tickfontsize = 14,ylims = [0,0.55])
savefig("Figure12e.pdf")
sd_MRPK_x_open_CMdelayed_open_trade_plotcorrection = sd_MRPK_x_open_CM_open_trade[1]./ sd_MRPK_x_store_open_CMdelayed_open_trade[1,end-1];
sd_MRPK_x_smooth_open_delayed = (sd_MRPK_x_open_CMdelayed_open_trade_plotcorrection * sd_MRPK_x_store_open_CMdelayed_open_trade[1,:]);
sd_MRPK_x_smooth_open_delayed= movmean(sd_MRPK_x_smooth_open_delayed,5);
sd_MRPK_x_open_CM_open_trade_plotcorrection = sd_MRPK_x_open_CM_open_trade[1]./ sd_MRPK_x_store_open_CM_open_trade[1,end];
sd_MRPK_x_closed_CM_open_trade_plotcorrection = sd_MRPK_x_closed_CM_open_trade[1]./ sd_MRPK_x_store_closed_CM_open_trade[1,end];
sd_MRPK_x_smooth_open = (sd_MRPK_x_open_CM_open_trade_plotcorrection*sd_MRPK_x_store_open_CM_open_trade[1,:]);
sd_MRPK_x_smooth_closed =  (sd_MRPK_x_closed_CM_open_trade_plotcorrection*sd_MRPK_x_store_closed_CM_open_trade[1,:]);
sd_MRPK_x_smooth_open = movmean(sd_MRPK_x_smooth_open,5);
sd_MRPK_x_smooth_closed = movmean(sd_MRPK_x_smooth_closed,5);
plot([sd_MRPK_x_smooth_open,sd_MRPK_x_smooth_closed,sd_MRPK_d_smooth_open_delayed],legend = nothing, linewidth = 3,linestyle = [:solid :dash :dot],grid = false,tickfontsize = 14,ylims = [0,0.55])
savefig("Figure12f.pdf")
p90_wealth_open_CMdelayed_open_trade_plotcorrection = p90_wealth_open_CM_open_trade[1]./ p90_wealth_store_open_CMdelayed_open_trade[1,end-1];
p90_wealth_smooth_open_delayed = (p90_wealth_open_CMdelayed_open_trade_plotcorrection * p90_wealth_store_open_CMdelayed_open_trade[1,:]);
p90_wealth_smooth_open_delayed= movmean(p90_wealth_smooth_open_delayed,5);
p90_wealth_open_CM_open_trade_plotcorrection = p90_wealth_open_CM_open_trade[1]./ p90_wealth_store_open_CM_open_trade[1,end];
p90_wealth_closed_CM_open_trade_plotcorrection = p90_wealth_closed_CM_open_trade[1]./ p90_wealth_store_closed_CM_open_trade[1,end];
p90_wealth_smooth_open = (p90_wealth_open_CM_open_trade_plotcorrection * p90_wealth_store_open_CM_open_trade[1,:]);
p90_wealth_smooth_closed =  (p90_wealth_closed_CM_open_trade_plotcorrection * p90_wealth_store_closed_CM_open_trade[1,:]);
p90_wealth_smooth_open = movmean(p90_wealth_smooth_open,5);
p90_wealth_smooth_closed = movmean(p90_wealth_smooth_closed,5);
plot([p90_wealth_smooth_open,p90_wealth_smooth_closed,p90_wealth_smooth_open_delayed],legend = nothing, linewidth =3,linestyle = [:solid :dash :dot],grid = false,tickfontsize = 14,ylims = [0.1,0.58])
savefig("Figure12g.pdf")
p90_income_open_CMdelayed_open_trade_plotcorrection = p90_income_open_CM_open_trade[1]./ p90_income_store_open_CMdelayed_open_trade[1,end-1];
p90_income_smooth_open_delayed = (p90_income_open_CMdelayed_open_trade_plotcorrection * p90_income_store_open_CMdelayed_open_trade[1,:]);
p90_income_smooth_open_delayed= movmean(p90_income_smooth_open_delayed,5);
p90_income_open_CM_open_trade_plotcorrection = p90_income_open_CM_open_trade[1]./ p90_income_store_open_CM_open_trade[1,end];
p90_income_closed_CM_open_trade_plotcorrection = p90_income_closed_CM_open_trade[1]./ p90_income_store_closed_CM_open_trade[1,end];
p90_income_smooth_open = (p90_income_open_CM_open_trade_plotcorrection*p90_income_store_open_CM_open_trade[1,:]);
p90_income_smooth_closed =  (p90_income_closed_CM_open_trade_plotcorrection*p90_income_store_closed_CM_open_trade[1,:]);
p90_income_smooth_open = movmean(p90_income_smooth_open,5);
p90_income_smooth_closed = movmean(p90_income_smooth_closed,5);
plot([p90_income_smooth_open,p90_income_smooth_closed,p90_income_smooth_open_delayed],legend = nothing, linewidth = 3,linestyle = [:solid :dash :dot],grid = false,tickfontsize = 14,ylims = [0.1,0.3])
savefig("Figure12h.pdf")
p90_cons_open_CMdelayed_open_trade_plotcorrection = p90_cons_open_CM_open_trade[1]./ p90_cons_store_open_CMdelayed_open_trade[1,end-1];
p90_cons_smooth_open_delayed = (p90_cons_open_CMdelayed_open_trade_plotcorrection * p90_cons_store_open_CMdelayed_open_trade[1,:]);
p90_cons_smooth_open_delayed= movmean(p90_cons_smooth_open_delayed,5);
p90_cons_open_CM_open_trade_plotcorrection = p90_cons_open_CM_open_trade[1]./ p90_cons_store_open_CM_open_trade[1,end];
p90_cons_closed_CM_open_trade_plotcorrection = p90_cons_closed_CM_open_trade[1]./ p90_cons_store_closed_CM_open_trade[1,end];
p90_cons_smooth_open = (p90_cons_open_CM_open_trade_plotcorrection*p90_cons_store_open_CM_open_trade[1,:]);
p90_cons_smooth_closed =  (p90_cons_closed_CM_open_trade_plotcorrection*p90_cons_store_closed_CM_open_trade[1,:]);
p90_cons_smooth_open = movmean(p90_cons_smooth_open,5);
p90_cons_smooth_closed = movmean(p90_cons_smooth_closed,5);
plot([p90_cons_smooth_open,p90_cons_smooth_closed,p90_cons_smooth_open_delayed],legend = nothing, linewidth = 3,linestyle = [:solid :dash :dot],grid = false,tickfontsize = 14,ylims = [0.1,0.25])
savefig("Figure12i.pdf")
#Only capital market integration
table_9_results = zeros(20,3);
table_9_results[1,:] = [TFP_initial[1]./TFP_initial[1],TFP_open_CM_closed_trade[1]./TFP_initial[1],TFP_open_CM_open_trade[1]./TFP_initial[1]];
table_9_results[2,:] = [sd_MRPK_initial[1],sd_MRPK_open_CM_closed_trade[1],sd_MRPK_open_CM_open_trade[1]];
table_9_results[3,:] = [GDP_initial[1]./GDP_initial[1],GDP_open_CM_closed_trade[1]./GDP_initial[1],GDP_open_CM_open_trade[1]./GDP_initial[1]];
table_9_results[4,:] = [prices_initial[3]./prices_initial[3],prices_open_CM_closed_trade[3]./prices_initial[3],prices_open_CM_open_trade[3]./prices_initial[3]];
table_9_results[5,:] = [total_consumption_initial[1]./total_consumption_initial[1],total_consumption_open_CM_closed_trade[1]./total_consumption_initial[1],total_consumption_open_CM_open_trade[1]./total_consumption_initial[1]];
table_9_results[6,:] = [capital_demand_initial[1]./capital_demand_initial[1],capital_demand_open_CM_closed_trade[1]./capital_demand_initial[1],capital_demand_open_CM_open_trade[1]./capital_demand_initial[1]];
table_9_results[8,:] = [0,welfare_change_oCM_cltrade_trans,welfare_change_oCM_otrade_trans];
table_9_results[9,:] = [p90_wealth_initial[1],p90_wealth_open_CM_closed_trade[1],p90_wealth_open_CM_open_trade[1]];
table_9_results[10,:] = [prices_initial[1]/prices_initial[5]/(prices_initial[1]/prices_initial[5]),prices_open_CM_closed_trade[1]/prices_open_CM_closed_trade[5]/(prices_initial[1]/prices_initial[5]),prices_open_CM_open_trade[1]/prices_open_CM_open_trade[5]/(prices_initial[1]/prices_initial[5])];
table_9_results[12,:] = [import_share_initial[1],import_share_open_CM_closed_trade[1],import_share_open_CM_open_trade[1]];
table_9_results[13,:] = [export_value_initial[1]./nomGDP_initial[2],export_value_open_CM_closed_trade[1]./nomGDP_open_CM_closed_trade[2],export_value_open_CM_open_trade[1]./nomGDP_open_CM_open_trade[2]];
table_9_results[14,:] = [(domestic_pop_initial[1]+ exporter_pop_initial[1]),(domestic_pop_open_CM_closed_trade[1]+ exporter_pop_open_CM_closed_trade[1]),(domestic_pop_open_CM_open_trade[1]+ exporter_pop_open_CM_open_trade[1])];
table_9_results[15,:] = [exporter_pop_initial[1]./(domestic_pop_initial[1]+ exporter_pop_initial[1]),exporter_pop_open_CM_closed_trade[1]./(domestic_pop_open_CM_closed_trade[1]+ exporter_pop_open_CM_closed_trade[1]),exporter_pop_open_CM_open_trade[1]./(domestic_pop_open_CM_open_trade[1]+ exporter_pop_open_CM_open_trade[1])];
table_9_results[16,:] = [prices_initial[5],prices_open_CM_closed_trade[5],prices_open_CM_open_trade[5]];

total_credit_open_CM_closed_trade = -(domestic_firm_debt_open_CM_closed_trade[1] + exporter_firm_debt_open_CM_closed_trade[1])/nomGDP_open_CM_closed_trade[1];
domestic_credit_open_CM_closed_trade = (worker_bond_holding_open_CM_closed_trade[1]+domestic_bond_holding_open_CM_closed_trade[1]+exporter_bond_holding_open_CM_closed_trade[1])./nomGDP_open_CM_closed_trade[1]; # domestic_credit to NFC/GDP
foreign_credit_open_CM_closed_trade= total_credit_open_CM_closed_trade - domestic_credit_open_CM_closed_trade;
foreign_credit_share_open_CM_closed_trade =foreign_credit_open_CM_closed_trade/total_credit_open_CM_closed_trade;

table_9_results[17,:] = [total_credit_initial,total_credit_open_CM_closed_trade,total_credit_open_CM_open_trade];
table_9_results[18,:] = [foreign_credit_share_initial,foreign_credit_share_open_CM_closed_trade,foreign_credit_share_open_CM_open_trade];

# Economy with higher financial development
# Steady state welfare for development:
welfare_init_developed =  (Phi_fine_aug * vcat(vcat(coeff_final_initial_dev[:,1,1],coeff_final_initial_dev[:,2,1]),coeff_final_initial_dev[:,3,1]));
welfare_init_developed_stst = sum(current_distr_store_initial[:,1] .* (exp.((welfare_init_developed* (1.0 - Developed_parameter.β[1]) - welfare_init_init * (1.0 - Baseline_parameter.β[1]))) - ones(Baseline_parameter.n_fine[1] * Baseline_parameter.n_fine[2] * 3)));
# First, compare the effects of financial development:
table_10_results = zeros(20,2);
table_10_results[1,:] = [TFP_initial[1]./TFP_initial[1],TFP_initial_dev[1]./TFP_initial[1]];
table_10_results[2,:] = [sd_MRPK_initial[1],sd_MRPK_initial_dev[1]];
table_10_results[3,:] = [GDP_initial[1]./GDP_initial[1],GDP_initial_dev[1]./GDP_initial[1]];
table_10_results[4,:] = [prices_initial[3]./prices_initial[3],prices_initial_dev[3]./prices_initial[3]];
table_10_results[5,:] = [total_consumption_initial[1]./total_consumption_initial[1],total_consumption_initial_dev[1]./total_consumption_initial[1]];
table_10_results[6,:] = [capital_demand_initial[1]./capital_demand_initial[1],capital_demand_initial_dev[1]./capital_demand_initial[1]];
table_10_results[7,:] = [0,welfare_init_developed_stst];
table_10_results[8,:] = [p90_wealth_initial[1],p90_wealth_initial_dev[1]];
table_10_results[9,:] = [p90_income_initial[1],p90_income_initial_dev[1]];
table_10_results[10,:] = [p90_cons_initial[1],p90_cons_initial_dev[1]];
table_10_results[11,:] = [prices_initial[1]/prices_initial[5]/(prices_initial[1]/prices_initial[5]),prices_initial_dev[1]/prices_initial_dev[5]/(prices_initial[1]/prices_initial[5])];
table_10_results[12,:] = [prices_initial[6]-prices_initial[7],prices_initial_dev[6]-prices_initial_dev[7]];
table_10_results[13,:] = [import_share_initial[1],import_share_initial_dev[1]];
table_10_results[14,:] = [export_value_initial[1]./nomGDP_initial[2],export_value_initial_dev[1]./nomGDP_initial_dev[2]];
table_10_results[15,:] = [(domestic_pop_initial[1]+ exporter_pop_initial[1]),(domestic_pop_initial_dev[1]+ exporter_pop_initial_dev[1])];
table_10_results[16,:] = [exporter_pop_initial[1]./(domestic_pop_initial[1]+ exporter_pop_initial[1]),exporter_pop_initial_dev[1]./(domestic_pop_initial_dev[1]+ exporter_pop_initial_dev[1])];
table_10_results[17,:] = [prices_initial[5],prices_initial_dev[5]];

total_credit_initial_dev = -(domestic_firm_debt_initial_dev[1] + exporter_firm_debt_initial_dev[1])/nomGDP_initial_dev[1];
domestic_credit_initial_dev = (worker_bond_holding_initial_dev[1]+domestic_bond_holding_initial_dev[1]+exporter_bond_holding_initial_dev[1])./nomGDP_initial_dev[1]; # domestic_credit to NFC/GDP
foreign_credit_initial_dev= total_credit_initial_dev - domestic_credit_initial_dev;
foreign_credit_share_initial_dev =foreign_credit_initial_dev/total_credit_initial_dev;

table_10_results[18,:] = [total_credit_initial,total_credit_initial_dev];
table_10_results[19,:] = [foreign_credit_share_initial,foreign_credit_share_initial_dev];

# Welfare comparison - between steady states:

Phi_fine_aug = kron(Matrix(1.0I, 3, 3),row_kron(Developed_parameter.Country_spec_p[1][6],funbase(Developed_parameter.fspace_a, Developed_parameter.Country_spec_p[1][3][:,1])));
welfare_init_init_dev =  (Phi_fine_aug * vcat(vcat(coeff_final_initial_dev[:,1,1],coeff_final_initial_dev[:,2,1]),coeff_final_initial_dev[:,3,1]));
welfare_init_closed_CM_open_trade_dev =  (Phi_fine_aug * vcat(vcat(coeff_final_closed_CM_open_trade_dev[:,1,1],coeff_final_closed_CM_open_trade_dev[:,2,1]),coeff_final_closed_CM_open_trade_dev[:,3,1]));
welfare_init_open_CM_open_trade_dev =  (Phi_fine_aug * vcat(vcat(coeff_final_open_CM_open_trade_dev[:,1,1],coeff_final_open_CM_open_trade_dev[:,2,1]),coeff_final_open_CM_open_trade_dev[:,3,1]));
welfare_init_closed_CM_open_trade_stst_dev = sum(current_distr_store_initial_dev[:,1] .* (exp.((welfare_init_closed_CM_open_trade_dev - welfare_init_init_dev) * (1.0 - Developed_parameter.β[1])) - ones(Developed_parameter.n_fine[1] * Developed_parameter.n_fine[2] * 3)));
welfare_init_open_CM_open_trade_stst_dev = sum(current_distr_store_initial_dev[:,1] .* (exp.((welfare_init_open_CM_open_trade_dev - welfare_init_init_dev) * (1.0 - Developed_parameter.β[1])) - ones(Developed_parameter.n_fine[1] * Developed_parameter.n_fine[2] * 3)));
#Transition dynamics
welfare_init_closed_CM_open_trade_trans_dev =  (Phi_fine_aug * vcat(vcat(coeff_store_closed_CM_open_trade_dev[:,1,1,2],coeff_store_closed_CM_open_trade_dev[:,2,1,2]),coeff_store_closed_CM_open_trade_dev[:,3,1,2]));
welfare_init_open_CM_open_trade_trans_dev =  (Phi_fine_aug * vcat(vcat(coeff_store_open_CM_open_trade_dev[:,1,1,2],coeff_store_open_CM_open_trade_dev[:,2,1,2]),coeff_store_open_CM_open_trade_dev[:,3,1,2]));
#welfare_init_open_CM_closed_trade_trans =  (Phi_fine_aug * vcat(vcat(coeff_store_open_CM_closed_trade[:,1,1,2],coeff_store_open_CM_closed_trade[:,2,1,2]),coeff_store_open_CM_closed_trade[:,3,1,2]));
welfare_change_clCM_otrade_trans_dev = sum(current_distr_store_initial_dev[:,1] .* (exp.((welfare_init_closed_CM_open_trade_trans_dev - welfare_init_init_dev) * (1.0 - Developed_parameter.β[1]))
- ones(Baseline_parameter.n_fine[1] * Baseline_parameter.n_fine[2] * 3)));
welfare_change_oCM_otrade_trans_dev = sum(current_distr_store_initial_dev[:,1] .* (exp.((welfare_init_open_CM_open_trade_trans_dev - welfare_init_init_dev) * (1.0 - Developed_parameter.β[1]))
- ones(Baseline_parameter.n_fine[1] * Baseline_parameter.n_fine[2] * 3)));
#Table 11 - the effects of reforms in a more advanced Economy
table_11_results = zeros(21,3);
table_11_results[1,:] = [TFP_initial_dev[1]./TFP_initial_dev[1],TFP_closed_CM_open_trade_dev[1]./TFP_initial_dev[1],TFP_open_CM_open_trade_dev[1]./TFP_initial_dev[1]];
table_11_results[2,:] = [sd_MRPK_initial_dev[1],sd_MRPK_closed_CM_open_trade_dev[1],sd_MRPK_open_CM_open_trade_dev[1]];
table_11_results[3,:] = [GDP_initial_dev[1]./GDP_initial_dev[1],GDP_closed_CM_open_trade_dev[1]./GDP_initial_dev[1],GDP_open_CM_open_trade_dev[1]./GDP_initial_dev[1]];
table_11_results[4,:] = [prices_initial_dev[3]./prices_initial_dev[3],prices_closed_CM_open_trade_dev[3]./prices_initial_dev[3],prices_open_CM_open_trade_dev[3]./prices_initial_dev[3]];
table_11_results[5,:] = [total_consumption_initial_dev[1]./total_consumption_initial_dev[1],total_consumption_closed_CM_open_trade_dev[1]./total_consumption_initial_dev[1],total_consumption_open_CM_open_trade_dev[1]./total_consumption_initial_dev[1]];
table_11_results[6,:] = [capital_demand_initial_dev[1]./capital_demand_initial_dev[1],capital_demand_closed_CM_open_trade_dev[1]./capital_demand_initial_dev[1],capital_demand_open_CM_open_trade_dev[1]./capital_demand_initial_dev[1]];
table_11_results[7,:] = [0,welfare_init_closed_CM_open_trade_stst_dev,welfare_init_open_CM_open_trade_stst_dev];
table_11_results[8,:] = [0,welfare_change_clCM_otrade_trans_dev,welfare_change_oCM_otrade_trans_dev];
table_11_results[9,:] = [p90_wealth_initial_dev[1],p90_wealth_closed_CM_open_trade_dev[1],p90_wealth_open_CM_open_trade_dev[1]];
table_11_results[10,:] = [p90_income_initial_dev[1],p90_income_closed_CM_open_trade_dev[1],p90_income_open_CM_open_trade_dev[1]];
table_11_results[11,:] = [p90_cons_initial_dev[1],p90_cons_closed_CM_open_trade_dev[1],p90_cons_open_CM_open_trade_dev[1]];
table_11_results[12,:] = [prices_initial_dev[1]/prices_initial_dev[5]/(prices_initial_dev[1]/prices_initial_dev[5]),prices_closed_CM_open_trade_dev[1]/prices_closed_CM_open_trade_dev[5]/(prices_initial_dev[1]/prices_initial_dev[5]),prices_open_CM_open_trade_dev[1]/prices_open_CM_open_trade_dev[5]/(prices_initial_dev[1]/prices_initial_dev[5])];
table_11_results[14,:] = [prices_initial_dev[6]-prices_initial_dev[7],prices_closed_CM_open_trade_dev[6]-prices_closed_CM_open_trade_dev[7],prices_open_CM_open_trade_dev[6]-prices_open_CM_open_trade_dev[6]];
table_11_results[15,:] = [import_share_initial_dev[1],import_share_closed_CM_open_trade_dev[1],import_share_open_CM_open_trade_dev[1]];
table_11_results[16,:] = [export_value_initial_dev[1]./nomGDP_initial_dev[2],export_value_closed_CM_open_trade_dev[1]./nomGDP_closed_CM_open_trade_dev[2],export_value_open_CM_open_trade_dev[1]./nomGDP_open_CM_open_trade_dev[2]];
table_11_results[17,:] = [(domestic_pop_initial_dev[1]+ exporter_pop_initial_dev[1]),(domestic_pop_closed_CM_open_trade_dev[1]+ exporter_pop_closed_CM_open_trade_dev[1]),(domestic_pop_open_CM_open_trade_dev[1]+ exporter_pop_open_CM_open_trade_dev[1])];
table_11_results[18,:] = [exporter_pop_initial_dev[1]./(domestic_pop_initial_dev[1]+ exporter_pop_initial_dev[1]),exporter_pop_closed_CM_open_trade_dev[1]./(domestic_pop_closed_CM_open_trade_dev[1]+ exporter_pop_closed_CM_open_trade_dev[1]),exporter_pop_open_CM_open_trade_dev[1]./(domestic_pop_open_CM_open_trade_dev[1]+ exporter_pop_open_CM_open_trade_dev[1])];
table_11_results[19,:] = [prices_initial_dev[5],prices_closed_CM_open_trade_dev[5],prices_open_CM_open_trade_dev[5]];

total_credit_closed_CM_open_trade_dev = -(domestic_firm_debt_closed_CM_open_trade_dev[1] + exporter_firm_debt_closed_CM_open_trade_dev[1])/nomGDP_closed_CM_open_trade_dev[1];
domestic_credit_closed_CM_open_trade_dev = (worker_bond_holding_closed_CM_open_trade_dev[1]+domestic_bond_holding_closed_CM_open_trade_dev[1]+exporter_bond_holding_closed_CM_open_trade_dev[1])./nomGDP_closed_CM_open_trade_dev[1]; # domestic_credit to NFC/GDP
foreign_credit_closed_CM_open_trade_dev= total_credit_closed_CM_open_trade_dev - domestic_credit_closed_CM_open_trade_dev;
foreign_credit_share_closed_CM_open_trade_dev =foreign_credit_closed_CM_open_trade_dev/total_credit_closed_CM_open_trade_dev;

total_credit_open_CM_open_trade_dev = -(domestic_firm_debt_open_CM_open_trade_dev[1] + exporter_firm_debt_open_CM_open_trade_dev[1])/nomGDP_open_CM_open_trade_dev[1];
domestic_credit_open_CM_open_trade_dev = (worker_bond_holding_open_CM_open_trade_dev[1]+domestic_bond_holding_open_CM_open_trade_dev[1]+exporter_bond_holding_open_CM_open_trade_dev[1])./nomGDP_open_CM_open_trade_dev[1]; # domestic_credit to NFC/GDP
foreign_credit_open_CM_open_trade_dev= total_credit_open_CM_open_trade_dev - domestic_credit_open_CM_open_trade_dev;
foreign_credit_share_open_CM_open_trade_dev =foreign_credit_open_CM_open_trade_dev/total_credit_open_CM_open_trade_dev;

domestic_credit_closed_CM_open_trade_dev = (worker_bond_holding_closed_CM_open_trade_dev[1]+domestic_bond_holding_closed_CM_open_trade_dev[1]+exporter_bond_holding_closed_CM_open_trade_dev[1])./nomGDP_closed_CM_open_trade_dev[1]; # domestic_credit to NFC/GDP
foreign_credit_closed_CM_open_trade_dev= -((domestic_firm_debt_closed_CM_open_trade_dev[1] + exporter_firm_debt_closed_CM_open_trade_dev[1])/nomGDP_closed_CM_open_trade_dev[1] + domestic_credit_closed_CM_open_trade_dev);
domestic_credit_open_CM_open_trade_dev = (worker_bond_holding_open_CM_open_trade_dev[1]+domestic_bond_holding_open_CM_open_trade_dev[1]+exporter_bond_holding_open_CM_open_trade_dev[1])./nomGDP_open_CM_open_trade_dev[1]; # domestic_credit to NFC/GDP
foreign_credit_open_CM_open_trade_dev= -((domestic_firm_debt_open_CM_open_trade_dev[1] + exporter_firm_debt_open_CM_open_trade_dev[1])/nomGDP_open_CM_open_trade_dev[1] + domestic_credit_open_CM_open_trade_dev);
table_11_results[20,:] = [total_credit_initial_dev,total_credit_closed_CM_open_trade_dev,total_credit_open_CM_open_trade_dev];
table_11_results[21,:] = [foreign_credit_share_initial_dev,foreign_credit_share_closed_CM_open_trade_dev,foreign_credit_share_open_CM_open_trade_dev];
