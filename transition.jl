#cd("C:\\Users\\laszl\\Dropbox\\Julia")
#using SharedArrays, Distributed

#include("main_script.jl")
#@everywhere begin
case_final =7; # Put to the main script
load_solution = 1; # 1 if load saved solutions
save_solution = 0; # 1 to save solution at the end
if case_final <5
    # Jacobian calculations:
#    function f(prices)
#        residual = Residual_stst(prices, Baseline_parameter);
#        residual_nonwalras = zeros(7);
#        residual_nonwalras[1:3] = residual[1:3];
#        residual_nonwalras[4:7] = residual[5:8];
#        return residual_nonwalras
#    end
#    function findiff_jacob_f_restricted(x)
#        G = FiniteDiff.finite_difference_jacobian(f,x);
#        return G
#    end
#    Jac_init = FiniteDiff.finite_difference_jacobian(f, prices_initial)
    Jac_init = 1.0*Matrix(I, 7, 7);
else
#    function f(prices)
#        residual = Residual_stst(prices, Developed_parameter);
#        residual_nonwalras = zeros(7);
#        residual_nonwalras[1:3] = residual[1:3];
#        residual_nonwalras[4:7] = residual[5:8];
#        return residual_nonwalras
#    end
#    function findiff_jacob_f_restricted(x)
#        G = FiniteDiff.finite_difference_jacobian(f,x);
#        return G
#    end
#    Jac_init = FiniteDiff.finite_difference_jacobian(f, prices_initial_dev)
    Jac_init = 1.0*Matrix(I, 7, 7);
end
if case_final == 1
    #Opening up with closed capital markets
    price_start = copy(prices_initial);
    coeff_begin = copy(coeff_final_initial);
    init_distr = copy(distr_current_initial);
    init_capital = copy(capital_demand_initial);
    price_finish = copy(prices_open_CM_closed_trade);
    coeff_end = copy(coeff_final_open_CM_closed_trade);
    fini_distr = copy(distr_current_open_CM_closed_trade);
    fini_capital = copy(capital_demand_open_CM_closed_trade);
    T = 30; #Initial guess for the transition length
    T_end = 15; # Initial guess - convergence after
    parameter_end = copy(open_CM_closed_trade_parameter);
    if load_solution == 1
        dataframe_prices = CSV.read("vec_price_openCM_notrade.csv", DataFrame,delim = ',',header=false);
        vec_price_trans = dataframe_prices.Column1;
    end
elseif case_final == 2
    #Opening up with closed capital markets
    price_start = copy(prices_initial);
    coeff_begin = copy(coeff_final_initial);
    init_distr = copy(distr_current_initial);
    init_capital = copy(capital_demand_initial);
    price_finish = copy(prices_closed_CM_open_trade);
    coeff_end = copy(coeff_final_closed_CM_open_trade);
    fini_distr = copy(distr_current_closed_CM_open_trade);
    fini_capital = copy(capital_demand_closed_CM_open_trade);
    T = 30; #Initial guess for the transition length
    T_end = 15; # Initial guess - convergence after
    parameter_end = copy(closed_CM_open_trade_parameter);
    if load_solution == 1
        dataframe_prices = CSV.read("vec_price_clCM.csv", DataFrame,delim = ',',header=false);
        vec_price_trans = dataframe_prices.Column1;
    end
elseif case_final == 3
    #Opening up with open capital markets
    price_start = copy(prices_initial);
    coeff_begin = copy(coeff_final_initial);
    init_distr = copy(distr_current_initial);
    init_capital = copy(capital_demand_initial);
    price_finish = copy(prices_open_CM_open_trade);
    coeff_end = copy(coeff_final_open_CM_open_trade);
    fini_distr = copy(distr_current_open_CM_open_trade);
    fini_capital = copy(capital_demand_open_CM_open_trade);
    T = 30; #Initial guess for the transition length
    T_end = 15; # Initial guess - convergence after
    parameter_end = copy(open_CM_open_trade_parameter);
    if load_solution == 1
        dataframe_prices = CSV.read("vec_price_openCM.csv", DataFrame,delim = ',',header=false);
        vec_price_trans = dataframe_prices.Column1;
    end
elseif case_final == 4
    #Opening up with closed capital markets - but
    #Open capital markets at T_cap
    # Assume solution for case_final 2 is available already
    # This is the expected plan
    T_cap = 10;
    T_end = 5;
    price_start = copy(prices_initial);
    coeff_begin = copy(coeff_final_initial);
    init_distr = copy(distr_current_initial);
    init_capital = copy(capital_demand_initial);
    price_finish = copy(prices_open_CM_open_trade);
    coeff_end = copy(coeff_final_open_CM_open_trade);
    fini_distr = copy(distr_current_open_CM_open_trade);
    fini_capital = copy(capital_demand_open_CM_open_trade);
    T = 30; #Initial guess for the transition length
    #T_end = 15 # Initial guess - convergence after
    parameter_end = copy(open_CM_open_trade_parameter);
    parameter_mid = copy(closed_CM_open_trade_parameter);
    price_mid = copy(prices_closed_CM_open_trade);
    # Load both solution of prices:
    if load_solution == 1
        dataframe_prices = CSV.read("vec_price_openCM_delayed.csv", DataFrame,delim = ',',header=false);
        vec_price_trans = dataframe_prices.Column1;
    end
elseif case_final == 5
    #Opening up with closed capital markets
    price_start = copy(prices_initial_dev);
    coeff_begin = copy(coeff_final_initial_dev);
    init_distr = copy(distr_current_initial_dev);
    init_capital = copy(capital_demand_initial_dev);
    price_finish = copy(prices_open_CM_closed_trade_dev);
    coeff_end = copy(coeff_final_open_CM_closed_trade_dev);
    fini_distr = copy(distr_current_open_CM_closed_trade_dev);
    fini_capital = copy(capital_demand_open_CM_closed_trade_dev);
    T = 30; #Initial guess for the transition length
    T_end = 15; # Initial guess - convergence after
    parameter_end = copy(open_CM_closed_trade_dev_parameter);
    if load_solution == 1
        dataframe_prices = CSV.read("vec_price_openCM_notrade_dev.csv", DataFrame,delim = ',',header=false);
        vec_price_trans = dataframe_prices.Column1;
    end
elseif case_final == 6
    #Opening up with closed capital markets
    price_start = copy(prices_initial_dev);
    coeff_begin = copy(coeff_final_initial_dev);
    init_distr = copy(distr_current_initial_dev);
    init_capital = copy(capital_demand_initial_dev);
    price_finish = copy(prices_closed_CM_open_trade_dev);
    coeff_end = copy(coeff_final_closed_CM_open_trade_dev);
    fini_distr = copy(distr_current_closed_CM_open_trade_dev);
    fini_capital = copy(capital_demand_closed_CM_open_trade_dev);
    T = 30; #Initial guess for the transition length
    T_end = 15; # Initial guess - convergence after
    parameter_end = copy(closed_CM_open_trade_dev_parameter);
    if load_solution == 1
        dataframe_prices = CSV.read("vec_price_clCM_dev.csv", DataFrame,delim = ',',header=false);
        vec_price_trans = dataframe_prices.Column1;
    end
elseif case_final == 7
    #Opening up with open capital markets
    price_start = copy(prices_initial_dev);
    coeff_begin = copy(coeff_final_initial_dev);
    init_distr = copy(distr_current_initial_dev);
    init_capital = copy(capital_demand_initial_dev);
    price_finish = copy(prices_open_CM_open_trade_dev);
    coeff_end = copy(coeff_final_open_CM_open_trade_dev);
    fini_distr = copy(distr_current_open_CM_open_trade_dev);
    fini_capital = copy(capital_demand_open_CM_open_trade_dev);
    T = 30; #Initial guess for the transition length
    T_end = 15; # Initial guess - convergence after
    parameter_end = copy(open_CM_open_trade_dev_parameter);
    if load_solution == 1
        dataframe_prices = CSV.read("vec_price_openCM_dev.csv", DataFrame,delim = ',',header=false);
        vec_price_trans = dataframe_prices.Column1;
    end
elseif case_final == 8
    #Opening up with closed capital markets - but
    #Open capital markets at T_cap
    # Assume solution for case_final 2 is available already
    # This is the expected plan
    T_cap = 10;
    T_end = 5;
    price_start = copy(prices_initial_dev);
    coeff_begin = copy(coeff_final_initial_dev);
    init_distr = copy(distr_current_initial_dev);
    init_capital = copy(capital_demand_initial_dev);
    price_finish = copy(prices_open_CM_open_trade_dev);
    coeff_end = copy(coeff_final_open_CM_open_trade_dev);
    fini_distr = copy(distr_current_open_CM_open_trade_dev);
    fini_capital = copy(capital_demand_open_CM_open_trade_dev);
    T = 30; #Initial guess for the transition length
    #T_end = 15 # Initial guess - convergence after
    parameter_end = copy(open_CM_open_trade_dev_parameter);
    parameter_mid = copy(closed_CM_open_trade_dev_parameter);
    price_mid = copy(prices_closed_CM_open_trade_dev);
    # Load both solution of prices:
    if load_solution == 1
        dataframe_prices = CSV.read("vec_price_openCM_delayed_dev.csv", DataFrame,delim = ',',header=false);
        vec_price_trans = dataframe_prices.Column1;
    end
end
#Final jacobian calculations:
if case_final <5
#    if parameter_end.openness ==0
#        function f(prices)
#            residual = Residual_stst(prices, parameter_end);
#            residual_nonwalras = zeros(7);
#            residual_nonwalras[1:3] = residual[1:3];
#            residual_nonwalras[4:7] = residual[5:8];
#            return residual_nonwalras
#        end
#    elseif parameter_end.openness ==1
#        function f(prices)
#            residual = Residual_stst(prices, parameter_end);
#            residual_nonwalras = zeros(6);
#            residual_nonwalras[1:3] = residual[1:3];
#            residual_nonwalras[4:6] = residual[5:7];
#            return residual_nonwalras
#        end
#    end
else
    if parameter_end.openness ==0
        Jac_fini  = 1.0*Matrix(I, 7, 7);
#        function f(prices)
#            residual = Residual_stst_frictionless(prices, parameter_end);
#            residual_nonwalras = zeros(7);
#            residual_nonwalras[1:3] = residual[1:3];
#            residual_nonwalras[4:7] = residual[5:8];
#            return residual_nonwalras
#        end
    elseif parameter_end.openness ==1
        Jac_fini  = 1.0*Matrix(I, 6, 6);
#        function f(prices)
#            residual = Residual_stst_frictionless(prices, parameter_end);
#            residual_nonwalras = zeros(6);
#            residual_nonwalras[1:3] = residual[1:3];
#            residual_nonwalras[4:6] = residual[5:7];
#            return residual_nonwalras
#        end
    end
end
#function findiff_jacob_f_restricted(x)
#    G = FiniteDiff.finite_difference_jacobian(f,x);
#    return G
#end
#Jac_fini = FiniteDiff.finite_difference_jacobian(f, price_finish)

#Jac_fini =gradient1(f, price_finish)
T_tau_transition = 4; # Slower introduction in the cuts of trade liberalization
τ_trans = ones(2,T+2);
τ_trans[1,:] = parameter_end.τ[1] * τ_trans[1,:];
τ_trans[2,:] = parameter_end.τ[2] * τ_trans[2,:];
τ_trans[1,1:T_tau_transition] = range(Baseline_parameter.τ[1],parameter_end.τ[1],length = T_tau_transition);
τ_trans[2,1:T_tau_transition] = range(Baseline_parameter.τ[2],parameter_end.τ[2],length = T_tau_transition);
if case_final == 1
    openness_transition = ones(T+2);
    openness_transition[1:2] .= 0;
elseif case_final == 2
    openness_transition = 0* ones(T+2);
elseif case_final == 3
    openness_transition = ones(T+2);
    openness_transition[1:2] .= 0;
elseif case_final == 4
    openness_transition = ones(T+2);
    openness_transition[1:(T_cap+2)] .= 0;
elseif case_final == 5
    openness_transition = ones(T+2);
    openness_transition[1:2] .= 0;
elseif case_final == 6
    openness_transition = 0* ones(T+2);
elseif case_final == 7
    openness_transition = ones(T+2);
    openness_transition[1:2] .= 0;
elseif case_final == 8
    openness_transition = ones(T+2);
    openness_transition[1:(T_cap+2)] .= 0;
end
openness_transition  = trunc.(Int, openness_transition);

coeff_store = zeros(size(coeff_end,1),size(coeff_end,2),parameter_end.country_no,T+2);
distr_store = zeros(size(init_distr,1),parameter_end.country_no,T+2);
distr_store[:,:,1] = init_distr;
distr_store[:,:,2] = init_distr;
distr_store[:,:,T + 2] = fini_distr;
if case_final == 1
    parameter_end.bounds[1,2] = 1.5 *maximum([price_start[1:2];price_finish[1:2]]);
    parameter_end.bounds[1,1] = 0.5 *minimum([price_start[1:2];price_finish[1:2]]);
    parameter_end.bounds[2,2] = 1.5 *maximum([price_start[6:7];price_finish[6]]);
    parameter_end.bounds[2,1] = 0.5*minimum([price_start[6:7];price_finish[6]]);
    parameter_end.bounds[3,2] = 1.5 *maximum([price_start[5];price_finish[5]]);
    parameter_end.bounds[3,1] = 0.5;
    parameter_end.bounds[4,2] = 1.5 *maximum([price_start[3:4];price_finish[3:4]]);
    parameter_end.bounds[4,1] = 0.5 *minimum([price_start[3:4];price_finish[3:4]]);
elseif case_final ==2
    parameter_end.bounds[1,2] = 1.5 *maximum([price_start[1:2];price_finish[1:2]]);
    parameter_end.bounds[1,1] = 0.5 *minimum([price_start[1:2];price_finish[1:2]]);
    parameter_end.bounds[2,2] = 1.5 *maximum([price_start[6:7];price_finish[6:7]]);
    parameter_end.bounds[2,1] = 0.5*minimum([price_start[6:7];price_finish[6:7]]);
    parameter_end.bounds[3,2] = 1.5 *maximum([price_start[5];price_finish[5]]);
    parameter_end.bounds[3,1] = 0.5;
    parameter_end.bounds[4,2] = 1.5 *maximum([price_start[3:4];price_finish[3:4]]);
    parameter_end.bounds[4,1] = 0.5 *minimum([price_start[3:4];price_finish[3:4]]);
elseif case_final == 3
    parameter_end.bounds[1,2] = 1.5 *maximum([price_start[1:2];price_finish[1:2]]);
    parameter_end.bounds[1,1] = 0.5 *minimum([price_start[1:2];price_finish[1:2]]);
    parameter_end.bounds[2,2] = 1.5 *maximum([price_start[6:7];price_finish[6]]);
    parameter_end.bounds[2,1] = 0.5*minimum([price_start[6:7];price_finish[6]]);
    parameter_end.bounds[3,2] = 1.5 *maximum([price_start[5];price_finish[5]]);
    parameter_end.bounds[3,1] = 0.5;
    parameter_end.bounds[4,2] = 1.5 *maximum([price_start[3:4];price_finish[3:4]]);
    parameter_end.bounds[4,1] = 0.5 *minimum([price_start[3:4];price_finish[3:4]]);
elseif case_final == 4
    parameter_end.bounds[1,2] = 1.5 *maximum([price_start[1:2];price_finish[1:2]]);
    parameter_end.bounds[1,1] = 0.5 *minimum([price_start[1:2];price_finish[1:2]]);
    parameter_end.bounds[2,2] = 1.5 *maximum([price_start[6:7];price_finish[6]]);
    parameter_end.bounds[2,1] = 0.5*minimum([price_start[6:7];price_finish[6]]);
    parameter_end.bounds[3,2] = 1.5 *maximum([price_start[5];price_finish[5]]);
    parameter_end.bounds[3,1] = 0.5;
    parameter_end.bounds[4,2] = 1.5 *maximum([price_start[3:4];price_finish[3:4]]);
    parameter_end.bounds[4,1] = 0.5 *minimum([price_start[3:4];price_finish[3:4]]);
elseif case_final == 5
    parameter_end.bounds[1,2] = 1.5 *maximum([price_start[1:2];price_finish[1:2]]);
    parameter_end.bounds[1,1] = 0.5 *minimum([price_start[1:2];price_finish[1:2]]);
    parameter_end.bounds[2,2] = 1.5 *maximum([price_start[6:7];price_finish[6]]);
    parameter_end.bounds[2,1] = 0.5*minimum([price_start[6:7];price_finish[6]]);
    parameter_end.bounds[3,2] = 1.5 *maximum([price_start[5];price_finish[5]]);
    parameter_end.bounds[3,1] = 0.5;
    parameter_end.bounds[4,2] = 1.5 *maximum([price_start[3:4];price_finish[3:4]]);
    parameter_end.bounds[4,1] = 0.5 *minimum([price_start[3:4];price_finish[3:4]]);
elseif case_final ==6
    parameter_end.bounds[1,2] = 1.5 *maximum([price_start[1:2];price_finish[1:2]]);
    parameter_end.bounds[1,1] = 0.5 *minimum([price_start[1:2];price_finish[1:2]]);
    parameter_end.bounds[2,2] = 1.5 *maximum([price_start[6:7];price_finish[6:7]]);
    parameter_end.bounds[2,1] = 0.5*minimum([price_start[6:7];price_finish[6:7]]);
    parameter_end.bounds[3,2] = 1.5 *maximum([price_start[5];price_finish[5]]);
    parameter_end.bounds[3,1] = 0.5;
    parameter_end.bounds[4,2] = 1.5 *maximum([price_start[3:4];price_finish[3:4]]);
    parameter_end.bounds[4,1] = 0.5 *minimum([price_start[3:4];price_finish[3:4]]);
elseif case_final == 7
    parameter_end.bounds[1,2] = 1.5 *maximum([price_start[1:2];price_finish[1:2]]);
    parameter_end.bounds[1,1] = 0.5 *minimum([price_start[1:2];price_finish[1:2]]);
    parameter_end.bounds[2,2] = 1.5 *maximum([price_start[6:7];price_finish[6]]);
    parameter_end.bounds[2,1] = 0.5*minimum([price_start[6:7];price_finish[6]]);
    parameter_end.bounds[3,2] = 1.5 *maximum([price_start[5];price_finish[5]]);
    parameter_end.bounds[3,1] = 0.5;
    parameter_end.bounds[4,2] = 1.5 *maximum([price_start[3:4];price_finish[3:4]]);
    parameter_end.bounds[4,1] = 0.5 *minimum([price_start[3:4];price_finish[3:4]]);
elseif case_final == 8
    parameter_end.bounds[1,2] = 1.5 *maximum([price_start[1:2];price_finish[1:2]]);
    parameter_end.bounds[1,1] = 0.5 *minimum([price_start[1:2];price_finish[1:2]]);
    parameter_end.bounds[2,2] = 1.5 *maximum([price_start[6:7];price_finish[6]]);
    parameter_end.bounds[2,1] = 0.5*minimum([price_start[6:7];price_finish[6]]);
    parameter_end.bounds[3,2] = 1.5 *maximum([price_start[5];price_finish[5]]);
    parameter_end.bounds[3,1] = 0.5;
    parameter_end.bounds[4,2] = 1.5 *maximum([price_start[3:4];price_finish[3:4]]);
    parameter_end.bounds[4,1] = 0.5 *minimum([price_start[3:4];price_finish[3:4]]);
end
#T_end = 10; % Initial guess for the length of already converged solution

coeff_store[:,:,:,1] = coeff_begin;
coeff_store[:,:,:,T+2] = coeff_end;

capital_trans = zeros(parameter_end.country_no,T+2);
capital_trans[:,1] = init_capital;
capital_trans[:,2] = init_capital;
capital_trans[:,T+2] = fini_capital;
if case_final == 1
    price_trans = zeros(7,T+2);
    price_trans[1:6,end] = price_finish;
    price_trans[1:7,1] = price_start;
    price_trans[1:7,2] = price_start;
    x = repeat(transpose(repeat(range(0,1,length = T_end),1,1)),6,1);
    price_trans[1:(end-1),3:(T_end+2)] = repeat(price_start[[1,2,3,4,5,7]],1,T_end) +  x.*repeat((price_finish - price_start[[1,2,3,4,5,7]]),1,T_end);
    price_trans[1:(end-1),(T_end+3):(end-1)] = repeat(price_finish,1,T-T_end-1);
    if  load_solution == 1
        price_trans[:,2:(T+1)] = reshape(vec_price_trans,7,T);
    end
elseif case_final ==2
    price_trans = repeat(price_finish,1,T+2);
    price_trans[:,1] = price_start;
    x = repeat(transpose(repeat(range(0,1,length = T-T_end),1,1)),7,1);
    price_trans[:,2:(T-T_end+1)] = repeat(price_start,1,T-T_end) +  x.*repeat((price_finish - price_start),1,T-T_end);
    if  load_solution == 1
        price_trans[:,2:(T+1)] = reshape(vec_price_trans,7,T);
    end
elseif case_final == 3
    price_trans = zeros(7,T+2);
    price_trans[1:6,end] = price_finish;
    price_trans[1:7,1] = price_start;
    price_trans[1:7,2] = price_start;
    x = repeat(transpose(repeat(range(0,1,length = T_end),1,1)),6,1);
    price_trans[1:(end-1),3:(T_end+2)] = repeat(price_start[[1,2,3,4,5,7]],1,T_end) +  x.*repeat((price_finish - price_start[[1,2,3,4,5,7]]),1,T_end);
    price_trans[1:(end-1),(T_end+3):(end-1)] = repeat(price_finish,1,T-T_end-1);
    if  load_solution == 1
        price_trans[:,2:(T+1)] = reshape(vec_price_trans,7,T);
    end
elseif case_final == 4
    price_trans = zeros(7,T+2);
    price_trans[1:6,end] = price_finish;
    price_trans[1:7,1] = price_start;
    x = repeat(transpose(repeat(range(0,1,length = T_cap-T_end),1,1)),7,1);
    price_trans[:,2:(T_cap-T_end+1)] = repeat(price_start,1,T_cap-T_end) +  x.*repeat((price_mid - price_start),1,T_cap-T_end);
    price_trans[:,(T_cap-T_end+2):(T_cap+2)] .= price_mid;

    x = repeat(transpose(repeat(range(0,1,length = T_end),1,1)),6,1);
    price_trans[1:(end-1),3:(T_end+2)] = repeat(price_start[[1,2,3,4,5,7]],1,T_end) +  x.*repeat((price_finish - price_start[[1,2,3,4,5,7]]),1,T_end);
    price_trans[1:(end-1),(T_end+3):(end-1)] = repeat(price_finish,1,T-T_end-1);
    if  load_solution == 1
        price_trans[:,2:(T+1)] = reshape(vec_price_trans,7,T);
    end
elseif case_final == 5
    price_trans = zeros(7,T+2);
    price_trans[1:6,end] = price_finish;
    price_trans[1:7,1] = price_start;
    price_trans[1:7,2] = price_start;
    x = repeat(transpose(repeat(range(0,1,length = T_end),1,1)),6,1);
    price_trans[1:(end-1),3:(T_end+2)] = repeat(price_start[[1,2,3,4,5,7]],1,T_end) +  x.*repeat((price_finish - price_start[[1,2,3,4,5,7]]),1,T_end);
    price_trans[1:(end-1),(T_end+3):(end-1)] = repeat(price_finish,1,T-T_end-1);
    if  load_solution == 1
        price_trans[:,2:(T+1)] = reshape(vec_price_trans,7,T);
    end
elseif case_final ==6
    price_trans = repeat(price_finish,1,T+2);
    price_trans[:,1] = price_start;
    x = repeat(transpose(repeat(range(0,1,length = T-T_end),1,1)),7,1);
    price_trans[:,2:(T-T_end+1)] = repeat(price_start,1,T-T_end) +  x.*repeat((price_finish - price_start),1,T-T_end);
    if  load_solution == 1
        price_trans[:,2:(T+1)] = reshape(vec_price_trans,7,T);
    end
elseif case_final == 7
    price_trans = zeros(7,T+2);
    price_trans[1:6,end] = price_finish;
    price_trans[1:7,1] = price_start;
    price_trans[1:7,2] = price_start;
    x = repeat(transpose(repeat(range(0,1,length = T_end),1,1)),6,1);
    price_trans[1:(end-1),3:(T_end+2)] = repeat(price_start[[1,2,3,4,5,7]],1,T_end) +  x.*repeat((price_finish - price_start[[1,2,3,4,5,7]]),1,T_end);
    price_trans[1:(end-1),(T_end+3):(end-1)] = repeat(price_finish,1,T-T_end-1);
    if  load_solution == 1
        price_trans[:,2:(T+1)] = reshape(vec_price_trans,7,T);
    end
elseif case_final == 8
    price_trans = zeros(7,T+2);
    price_trans[1:6,end] = price_finish;
    price_trans[1:7,1] = price_start;
    x = repeat(transpose(repeat(range(0,1,length = T_cap-T_end),1,1)),7,1);
    price_trans[:,2:(T_cap-T_end+1)] = repeat(price_start,1,T_cap-T_end) +  x.*repeat((price_mid - price_start),1,T_cap-T_end);
    price_trans[:,(T_cap-T_end+2):(T_cap+2)] .= price_mid;

    x = repeat(transpose(repeat(range(0,1,length = T_end),1,1)),6,1);
    price_trans[1:(end-1),3:(T_end+2)] = repeat(price_start[[1,2,3,4,5,7]],1,T_end) +  x.*repeat((price_finish - price_start[[1,2,3,4,5,7]]),1,T_end);
    price_trans[1:(end-1),(T_end+3):(end-1)] = repeat(price_finish,1,T-T_end-1);
    if  load_solution == 1
        price_trans[:,2:(T+1)] = reshape(vec_price_trans,7,T);
    end
end

vec_price_trans = copy(reshape(price_trans[:,2:(end-1)],(T)*7));
@everywhere begin
    function local_var_creator_policy(country_no::Int64,ns_tmp::Int64,ns_tmp_fine::Int64,size_iid_cost_val::Int64,T::Int64)
        a_prime_fine = SharedArray{Float64}(ns_tmp_fine,3,country_no,size_iid_cost_val,T);#Array{Float64}(undef, ns_tmp_fine,3,country_no,size_iid_cost_val,T);#
        future_occupation_fine = SharedArray{Float64}(ns_tmp_fine,3,country_no,size_iid_cost_val,T);#Array{Float64}(undef, ns_tmp_fine,3,country_no,size_iid_cost_val,T);#
        cons_fine = SharedArray{Float64}(ns_tmp_fine,3,country_no,size_iid_cost_val,T);#Array{Float64}(undef, ns_tmp_fine,3,country_no,size_iid_cost_val,T);#
        rev_d_fine_store = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
        rev_dx_fine_store = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
        rev_xx_fine_store = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
        k_x_fine = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
        k_d_fine = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
        l_x_fine = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
        l_d_fine = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
        #Profits:
        Profit_fine_x = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
        Profit_fine_d = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
        #Output
        output_d_fine_store = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
        output_dx_fine_store = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
        output_xx_fine_store = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
        #Prices
        price_d_fine_store = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
        price_dx_fine_store = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
        price_xx_fine_store = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
        #Country aggregates
        labor_excess_demand_store_percent = SharedArray{Float64}(country_no,T);#Array{Float64}(undef,country_no,T);#
        excess_demand_final_good_store_percent = SharedArray{Float64}(country_no,T);#Array{Float64}(undef,country_no,T);#
        total_demand_final_good = SharedArray{Float64}(country_no,T);#Array{Float64}(undef,country_no,T);
        export_price_sum_store = SharedArray{Float64}(country_no,T);#Array{Float64}(undef,country_no,T);#
        domestic_price_sum_store = SharedArray{Float64}(country_no,T);#Array{Float64}(undef,country_no,T);#SharedArray{Float64}(country_no,T);
        NFA_store = SharedArray{Float64}(country_no,T);#Array{Float64}(undef,country_no,T);#SharedArray{Float64}(country_no,T);
        asset_supply_store = SharedArray{Float64}(country_no,T);#Array{Float64}(undef,country_no,T);#SharedArray{Float64}(country_no,T);
        asset_demand_store = SharedArray{Float64}(country_no,T);#Array{Float64}(undef,country_no,T);#SharedArray{Float64}(country_no,T);
        return (a_prime_fine,future_occupation_fine,cons_fine,rev_d_fine_store,
            rev_dx_fine_store,rev_xx_fine_store,k_x_fine,k_d_fine,
            l_x_fine,l_d_fine,Profit_fine_x,
            Profit_fine_d,output_d_fine_store,output_dx_fine_store,output_xx_fine_store,
            price_d_fine_store,price_dx_fine_store,price_xx_fine_store,labor_excess_demand_store_percent,
            excess_demand_final_good_store_percent,total_demand_final_good,export_price_sum_store,domestic_price_sum_store,NFA_store,
            asset_supply_store,asset_demand_store)
    end
function local_var_creator_policy_detailed(country_no::Int64,ns_tmp::Int64,ns_tmp_fine::Int64,size_iid_cost_val::Int64,T::Int64)
    a_prime_fine = SharedArray{Float64}(ns_tmp_fine,3,country_no,size_iid_cost_val,T);#Array{Float64}(undef, ns_tmp_fine,3,country_no,size_iid_cost_val,T);#
    future_occupation_fine = SharedArray{Float64}(ns_tmp_fine,3,country_no,size_iid_cost_val,T);#Array{Float64}(undef, ns_tmp_fine,3,country_no,size_iid_cost_val,T);#
    cons_fine = SharedArray{Float64}(ns_tmp_fine,3,country_no,size_iid_cost_val,T);#Array{Float64}(undef, ns_tmp_fine,3,country_no,size_iid_cost_val,T);#
    rev_d_fine_store = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
    rev_dx_fine_store = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
    rev_xx_fine_store = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
    k_x_fine = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
    k_d_fine = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
    l_x_fine = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
    l_d_fine = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
    #Profits:
    Profit_fine_x = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
    Profit_fine_d = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
    #Output
    output_d_fine_store = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
    output_dx_fine_store = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
    output_xx_fine_store = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
    #Prices
    price_d_fine_store = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
    price_dx_fine_store = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
    price_xx_fine_store = SharedArray{Float64}(ns_tmp_fine,country_no,T);#Array{Float64}(undef, ns_tmp_fine,country_no,T);#
    #Country aggregates
    labor_excess_demand_store_percent = SharedArray{Float64}(country_no,T);#Array{Float64}(undef,country_no,T);#
    excess_demand_final_good_store_percent = SharedArray{Float64}(country_no,T);#Array{Float64}(undef,country_no,T);#
    total_demand_final_good = SharedArray{Float64}(country_no,T);#Array{Float64}(undef,country_no,T);
    export_price_sum_store = SharedArray{Float64}(country_no,T);#Array{Float64}(undef,country_no,T);#
    domestic_price_sum_store = SharedArray{Float64}(country_no,T);#Array{Float64}(undef,country_no,T);#SharedArray{Float64}(country_no,T);
    NFA_store = SharedArray{Float64}(country_no,T);#Array{Float64}(undef,country_no,T);#SharedArray{Float64}(country_no,T);
    asset_supply_store = SharedArray{Float64}(country_no,T);#Array{Float64}(undef,country_no,T);#SharedArray{Float64}(country_no,T);
    asset_demand_store = SharedArray{Float64}(country_no,T);#Array{Float64}(undef,country_no,T);#SharedArray{Float64}(country_no,T);
    # Additional, only required for details
    lambdda_d_fine_store = SharedArray{Float64}(ns_tmp_fine,country_no,T);
    lambdda_x_fine_store = SharedArray{Float64}(ns_tmp_fine,country_no,T);
    p90_wealth_store = SharedArray{Float64}(country_no,T);
    p90_cons_store = SharedArray{Float64}(country_no,T);
    p90_income_store = SharedArray{Float64}(country_no,T);
    fraction_zombie_exporter_store = SharedArray{Float64}(country_no,T);
    domestic_prod_store = SharedArray{Float64}(country_no,T);
    export_prod_store = SharedArray{Float64}(country_no,T);
    export_value_store = SharedArray{Float64}(country_no,T);
    sd_MRPK_d_store = SharedArray{Float64}(country_no,T);
    sd_MRPK_x_store = SharedArray{Float64}(country_no,T);
    sd_MRPK_store = SharedArray{Float64}(country_no,T);
    D_d_denom_store = SharedArray{Float64}(country_no,T);
    D_d_store = SharedArray{Float64}(country_no,T);
    D_x_denom_store = SharedArray{Float64}(country_no,T);
    D_x_store = SharedArray{Float64}(country_no,T);
    D_d_denom_eff_store = SharedArray{Float64}(country_no,T);
    D_d_eff_store = SharedArray{Float64}(country_no,T);
    D_x_denom_eff_store = SharedArray{Float64}(country_no,T);
    D_x_eff_store = SharedArray{Float64}(country_no,T);
    domestic_pop_store = SharedArray{Float64}(country_no,T);
    exporter_pop_store = SharedArray{Float64}(country_no,T);
    return (a_prime_fine,future_occupation_fine,cons_fine,rev_d_fine_store,
        rev_dx_fine_store,rev_xx_fine_store,k_x_fine,k_d_fine,
        l_x_fine,l_d_fine,Profit_fine_x,
        Profit_fine_d,output_d_fine_store,output_dx_fine_store,output_xx_fine_store,
        price_d_fine_store,price_dx_fine_store,price_xx_fine_store,labor_excess_demand_store_percent,
        excess_demand_final_good_store_percent,total_demand_final_good,export_price_sum_store,domestic_price_sum_store,NFA_store,
        asset_supply_store,asset_demand_store,lambdda_d_fine_store,lambdda_x_fine_store,p90_wealth_store,p90_cons_store,p90_income_store,
        fraction_zombie_exporter_store,domestic_prod_store,export_prod_store,export_value_store,sd_MRPK_d_store,sd_MRPK_x_store,sd_MRPK_store,
        D_d_denom_store,D_d_store,D_x_denom_store,D_x_store,D_d_denom_eff_store,D_d_eff_store,D_x_denom_eff_store,D_x_eff_store,domestic_pop_store,
        exporter_pop_store)
end
function Residual_transition_backward(
        country::Int64,price_final_prev::Float64,price_final_current::Float64, R::Float64,W_loc::Float64,r_loc::Float64 ,s_curr::Array{Float64,2},
        constant_x::Float64,constant_d::Float64,constant_z::Array{Float64,1},constant_z_bar::Array{Float64,1},
        β_loc::Float64,θ_loc::Float64,τ_loc::Float64,σ::Float64,α₁_eff::Float64,α₂_eff::Float64,α₁::Float64,α₂::Float64,
        ones_tmp::Array{Float64,1},ns::Int64,FM::Array{Float64,1},FMₓ::Array{Float64,1},V_tmp::Array{Float64,3},
        F::Array{Float64,1}, Fₓ::Array{Float64,1},size_iid_cost_val::Int64,iid_cost_value::Array{Float64,1},iid_cost_prob::Array{Float64,1},Exit::Array{Float64,1},
        Exitₓ::Array{Float64,1},a_min::Float64,a_max::Float64,α::Float64,fspace_a::Dict{Symbol,Any},
        fspace_a_fine::Dict{Symbol,Any},openness::Int64,agrid_fine::Array{Float64,1},coefficients_next_tmp::Array{Float64,2},
        P_kron::SparseMatrixCSC{Float64,Int64},Phi::SparseMatrixCSC{Float64,Int64},Phi_z::SparseMatrixCSC{Float64,Int64},income_mat::Array{Float64,4},x_tmp::Array{Float64,3},
        coeff_next::Array{Float64,2})


    (Profits_d,Profits_x,output_dx,output_xx,output_d,labor_d,labor_x,k_choice_d,k_choice_x,price_dx,price_d,
        price_xx,lambdda_d,lambdda_x,rev_d,rev_dx,rev_xx ) = true_profit(price_final_prev,R,W_loc,s_curr,
        constant_x,constant_d,constant_z,constant_z_bar,θ_loc,τ_loc,σ,α₁_eff,α₂_eff,α₁,α₂,ones_tmp);
    (income_mat,max_x) = income_creator(W_loc,ones_tmp,Profits_d,Profits_x,FM,FMₓ,country,
        r_loc,s_curr,iid_cost_value,F, Fₓ,size_iid_cost_val,income_mat,ns,a_min,a_max);
    (coeff, conv) = Bellman_iteration(coefficients_next_tmp,coeff_next,ns,x_tmp,V_tmp,ones_tmp,size_iid_cost_val,
        iid_cost_prob,income_mat,P_kron,Phi,a_min,Phi_z,β_loc,
        fspace_a,price_final_current,max_x);
    return coeff
end

function Residual_transition_forward(country::Int64,price_final_prev::Float64,price_final_current::Float64, R::Float64,W_loc::Float64,r_loc::Float64 ,s_fine::Array{Float64,2},
constant_x::Float64,constant_d::Float64,constant_z_fine::Array{Float64,1},constant_z_bar_fine::Array{Float64,1},
β_loc::Float64,θ_loc::Float64,τ_loc::Float64,σ::Float64,α₁_eff::Float64,α₂_eff::Float64,α₁::Float64,α₂::Float64,
ones_tmp_fine::Array{Float64,1},ns_fine::Int64,FM::Array{Float64,1},FMₓ::Array{Float64,1},V_tmp::Array{Float64,3},
F::Array{Float64,1}, Fₓ::Array{Float64,1},size_iid_cost_val::Int64,iid_cost_value::Array{Float64,1},iid_cost_prob::Array{Float64,1},Exit::Array{Float64,1},
Exitₓ::Array{Float64,1},a_min::Float64,a_max::Float64,α::Float64,fspace_a::Dict{Symbol,Any},
fspace_a_fine::Dict{Symbol,Any},openness::Int64,agrid_fine::Array{Float64,1},coefficients_next_tmp::Array{Float64,2},
P_kron_fine::SparseMatrixCSC{Float64,Int64},P_kron1::SparseMatrixCSC{Float64,Int64},Phi::SparseMatrixCSC{Float64,Int64},Phi_z_fine::SparseMatrixCSC{Float64,Int64},income_mat_fine::Array{Float64,4},x_tmp::Array{Float64,3},
 cons_fine_local::Array{Float64,3},a_prime_fine_local::Array{Float64,3},future_occupation_fine_local::Array{Float64,3},Q_trans::SparseMatrixCSC{Float64,Int64},distr_prev::Array{Float64,1})
#distr_prev = distr_store_tmp[:,country,t-1];
#coefficients_next_tmp = coeff_store_tmp[:,:,country,t+1];
    (Profits_d_fine,Profits_x_fine,output_dx_fine,output_xx_fine,output_d_fine,labor_d_fine,labor_x_fine
        ,k_choice_d_fine,k_choice_x_fine,price_dx_fine,price_d_fine,price_xx_fine,lambdda_d_fine,
        lambdda_x_fine,rev_d_fine,rev_dx_fine,rev_xx_fine) = true_profit(price_final_prev,R,W_loc,s_fine
        ,constant_x,constant_d,constant_z_fine,constant_z_bar_fine,θ_loc,τ_loc,σ,α₁_eff,α₂_eff,α₁,α₂,ones_tmp_fine);
    (income_mat_fine,max_x_fine) = income_creator(W_loc,ones_tmp_fine,Profits_d_fine,Profits_x_fine,FM,FMₓ,country,
            r_loc,s_fine,iid_cost_value,F, Fₓ,size_iid_cost_val,income_mat_fine,ns_fine,a_min,a_max);
    (Q_trans_prime,cons_fine_local,a_prime_fine_local,
    future_occupation_fine_local) = Q_transition(coefficients_next_tmp,
        ns_fine,ones_tmp_fine,size_iid_cost_val,iid_cost_prob,income_mat_fine,
        P_kron_fine,a_min,Phi_z_fine,β_loc,fspace_a,fspace_a_fine,
        price_final_current,max_x_fine,P_kron1,cons_fine_local,a_prime_fine_local,
        future_occupation_fine_local,Q_trans);
    distr_current = Q_trans_prime*distr_prev;
    return (a_prime_fine_local,future_occupation_fine_local,
    cons_fine_local,rev_d_fine,rev_dx_fine,rev_xx_fine,k_choice_x_fine,k_choice_d_fine,labor_x_fine,labor_d_fine
    ,Profits_x_fine,Profits_d_fine,output_d_fine,output_dx_fine,output_xx_fine,
    price_d_fine,price_dx_fine,price_xx_fine,distr_current,lambdda_d_fine,lambdda_x_fine)
end


function Residual_transition_backward_aggregates(country::Int64,output_final_loc::Float64,price_final_prev::Float64,price_final_current::Float64, R::Float64,W_loc::Float64,r_loc::Float64 ,L_loc::Float64,s_fine::Array{Float64,2},
σ::Float64,ns_fine::Int64,FM::Array{Float64,1},FMₓ::Array{Float64,1},banking_cost_loc::Float64,δ::Float64,
F::Array{Float64,1}, Fₓ::Array{Float64,1},size_iid_cost_val::Int64,iid_cost_value::Array{Float64,1},iid_cost_prob::Array{Float64,1},Exit::Array{Float64,1},
Exitₓ::Array{Float64,1},capital_supply_future_tmp::Float64,
a_prime_fine_local::Array{Float64,2},future_occupation_fine_local::Array{Float64,2},cons_fine_local::Array{Float64,2},
rev_d_fine::Array{Float64,1},rev_dx_fine::Array{Float64,1},rev_xx_fine::Array{Float64,1},k_choice_x_fine::Array{Float64,1},
k_choice_d_fine::Array{Float64,1},labor_x_fine::Array{Float64,1},labor_d_fine::Array{Float64,1},
Profits_x_fine::Array{Float64,1},Profits_d_fine::Array{Float64,1},output_d_fine::Array{Float64,1},
output_dx_fine::Array{Float64,1},output_xx_fine::Array{Float64,1},price_d_fine::Array{Float64,1},
price_dx_fine::Array{Float64,1},price_xx_fine::Array{Float64,1}, distr_prev::Array{Float64,1},capital_supply_past_tmp::Float64 = 0)

    # Calculate aggregates
    # Distribution of current workers:
    worker_past_dist = distr_prev[(ns_fine *0 + 1):(ns_fine *1)];
    domestic_past_dist = distr_prev[(ns_fine *1 + 1):(ns_fine *2)];
    exporter_past_dist = distr_prev[(ns_fine *2 + 1):(ns_fine *3)];

    # Decisions avgs:
    future_occupation_fine_avgs =  Array{SparseMatrixCSC{Float64,Int64}}(undef, 3,3);
    cons_fine_avgs =  Array{SparseMatrixCSC{Float64,Int64}}(undef, 3,1);
    cons_fine_second_moment =  Array{SparseMatrixCSC{Float64,Int64}}(undef, 3,1);
    #here
    for j=1:3
     for jj=1:3
         future_occupation_fine_avgs[j,jj] = spzeros(ns_fine);
         for j_cost = 1:size_iid_cost_val
             future_occupation_fine_avgs[j,jj]=future_occupation_fine_avgs[j,jj] + iid_cost_prob[j_cost]*(
                 future_occupation_fine_local[:,j,j_cost] .== jj);
         end
     end
     cons_fine_avgs[j,1] = spzeros(ns_fine);
     cons_fine_second_moment[j,1] = spzeros(ns_fine);
     for j_cost = 1:size_iid_cost_val
             cons_fine_avgs[j,1]=cons_fine_avgs[j,1] + iid_cost_prob[j_cost
             ]*cons_fine_local[:,j,j_cost];
             cons_fine_second_moment[j,1]=cons_fine_second_moment[j,1] + iid_cost_prob[j_cost
             ]*cons_fine_local[:,j,j_cost].^2.0;
     end
    end
    # Distribution of current workers:
    stay_workers = worker_past_dist.*future_occupation_fine_avgs[1,1];
    exit_domestic = domestic_past_dist.*future_occupation_fine_avgs[2,1];
    exit_exporting_to_work = exporter_past_dist.*future_occupation_fine_avgs[3,1];
    current_workers = stay_workers + exit_domestic + exit_exporting_to_work;

    entrants_domestic_from_workers = worker_past_dist.*future_occupation_fine_avgs[1,2];
    incumbents_domestic = domestic_past_dist.*future_occupation_fine_avgs[2,2];
    exit_exporting_to_domestic = exporter_past_dist.*future_occupation_fine_avgs[3,2];
    current_domestic = entrants_domestic_from_workers + incumbents_domestic + exit_exporting_to_domestic;

    entrants_to_exporters = worker_past_dist.*future_occupation_fine_avgs[1,3
     ] +  domestic_past_dist.*future_occupation_fine_avgs[2,3];
    incumbents_exporter = exporter_past_dist.*future_occupation_fine_avgs[3,3];
    current_exporter = entrants_to_exporters + incumbents_exporter;
    current_distr_store_tmp = copy(distr_prev);
    current_distr_store_tmp[(ns_fine *0 + 1):(ns_fine *1)] = current_workers;
    current_distr_store_tmp[(ns_fine *1 + 1):(ns_fine *2)] = current_domestic;
    current_distr_store_tmp[(ns_fine *2 + 1):(ns_fine *3)] = current_exporter;
    # Fixed cost accounting:
    # entrants for a given cost realization
    entry_costs_payed =  Array{SparseMatrixCSC{Float64,Int64}}(undef, 2,1);;
    entry_costs_payed[1,1] = spzeros(ns_fine);
    entry_costs_payed[2,1]= spzeros(ns_fine);
    for j_cost = 1:size_iid_cost_val
     entry_costs_payed[1,1] = entry_costs_payed[1,1] + iid_cost_prob[j_cost]*iid_cost_value[
         j_cost]*worker_past_dist.*(future_occupation_fine_local[:,1,j_cost] .== 2); # from worker to domestic
     entry_costs_payed[2,1]  = entry_costs_payed[2,1] + iid_cost_prob[j_cost]*iid_cost_value[
         j_cost]*(worker_past_dist.*(future_occupation_fine_local[:,1,j_cost] .== 3) #from worker to exporter
         + domestic_past_dist.*(future_occupation_fine_local[:,2,j_cost] .== 3)); # from domestic to exporter
    end

    domestic_prod_tmp = sum(current_exporter .* output_dx_fine.^((σ - 1)/σ) +
         current_domestic.* output_d_fine.^((σ - 1)/σ));
    export_prod_tmp = sum(current_exporter.* output_xx_fine .^((σ - 1)/σ));

    domestic_price_sum_tmp = sum(current_exporter.* price_dx_fine.^(1 - σ) + current_domestic.* price_d_fine.^(1 - σ));
    export_price_sum_tmp = sum(current_exporter.* price_xx_fine.^(1 - σ)) ;
    export_value_tmp = sum(current_exporter.* price_xx_fine .* output_xx_fine);
    exit_share_to_domestic_tmp = sum(exit_exporting_to_domestic);
    entry_share_to_domestic_tmp =sum(entrants_domestic_from_workers);
    exit_domestic_sum_tmp = sum(exit_domestic);
    exit_exporting_to_work_sum_tmp = sum(exit_exporting_to_work);
    exit_share_to_worker_tmp = exit_domestic_sum_tmp + exit_exporting_to_work_sum_tmp;
    entry_share_to_exporter_tmp = sum(entrants_to_exporters);
    exporter_pop_tmp = sum(current_exporter);
    domestic_pop_tmp = sum(current_domestic);
    worker_pop_tmp =  L_loc - exporter_pop_tmp - domestic_pop_tmp;
    L_x_tmp = sum(current_exporter.*labor_x_fine);
    L_d_tmp =sum(current_domestic.*labor_d_fine);
    total_entry_cost_tmp = sum(entry_costs_payed[1,1])  * F[country,1] + sum(entry_costs_payed[2,1]
     ) * Fₓ[country,1];
    total_incumbents_cost_tmp = exporter_pop_tmp* FMₓ[country,1] + domestic_pop_tmp * FM[country,1];
    total_exit_cost_tmp= exit_domestic_sum_tmp*Exit[country,1] + exit_exporting_to_work_sum_tmp * Exitₓ[country,1];
    labor_demand_tmp =  L_x_tmp + L_d_tmp + total_entry_cost_tmp + total_incumbents_cost_tmp + total_exit_cost_tmp;
    labor_excess_demand_tmp = labor_demand_tmp - worker_pop_tmp;
    K_x_tmp = sum(current_exporter.*k_choice_x_fine);
    K_d_tmp= sum(current_domestic.*k_choice_d_fine);
    capital_demand_tmp = (K_d_tmp + K_x_tmp);

    worker_bond_holding_tmp = s_fine[:,1] .*current_workers;
    worker_bond_holding_sum_tmp = sum(worker_bond_holding_tmp);
    domestic_bond_holding_tmp = (s_fine[:,1] - price_final_prev * k_choice_d_fine) .*current_domestic;
    domestic_bond_holding_sum_tmp = sum(domestic_bond_holding_tmp[domestic_bond_holding_tmp .> 0.0]);
    domestic_firm_debt_tmp = sum(domestic_bond_holding_tmp[domestic_bond_holding_tmp .<0.0]);
    exporter_bond_holding_tmp = (s_fine[:,1]- price_final_prev * k_choice_x_fine) .*current_exporter;
    exporter_bond_holding_sum_tmp = sum(exporter_bond_holding_tmp[exporter_bond_holding_tmp .> 0.0]);
    exporter_firm_debt_tmp = sum(exporter_bond_holding_tmp[exporter_bond_holding_tmp .<0.0]);
    asset_supply_tmp = worker_bond_holding_sum_tmp + domestic_bond_holding_sum_tmp + exporter_bond_holding_sum_tmp;
    asset_demand_tmp = - (domestic_firm_debt_tmp + exporter_firm_debt_tmp);
    NFA_tmp = asset_supply_tmp - (1 + banking_cost_loc) * asset_demand_tmp;
    total_consumption_tmp = sum(worker_past_dist.* cons_fine_avgs[1,1] + domestic_past_dist.* cons_fine_avgs[2,1]
     + exporter_past_dist.* cons_fine_avgs[3,1]);
    Aggr_bank_cost_tmp = ((banking_cost_loc)*asset_demand_tmp)/price_final_current;
    # normalize as it is denominated in foreign final good!
    x_k_tmp = K_x_tmp/capital_demand_tmp; # capital allocated to exports
    x_l_tmp = L_x_tmp/labor_demand_tmp; # labor allocated to exports
    investment_tmp = capital_supply_future_tmp - (1 -δ)* capital_demand_tmp;
    total_demand_final_good_tmp = total_consumption_tmp  + investment_tmp;
    excess_demand_final_good_tmp = total_demand_final_good_tmp - output_final_loc;

    labor_excess_demand_tmp_percent = labor_excess_demand_tmp./L_loc;
    excess_demand_final_good_tmp_percent = excess_demand_final_good_tmp/(output_final_loc + total_demand_final_good_tmp);
    if capital_supply_past_tmp >0
        NFA_tmp = capital_supply_past_tmp-capital_demand_tmp;
    end
    return (labor_excess_demand_tmp_percent,excess_demand_final_good_tmp_percent,
    total_demand_final_good_tmp,export_price_sum_tmp, domestic_price_sum_tmp,NFA_tmp,asset_supply_tmp,
    asset_demand_tmp,capital_demand_tmp)
end

function Residual_transition_backward_aggregates_detailed(country::Int64,output_final_loc::Float64,price_final_prev::Float64,price_final_current::Float64, R::Float64,W_loc::Float64,r_loc::Float64 ,L_loc::Float64,s_fine::Array{Float64,2},
σ::Float64,ns_fine::Int64,FM::Array{Float64,1},FMₓ::Array{Float64,1},banking_cost_loc::Float64,δ::Float64,
F::Array{Float64,1}, Fₓ::Array{Float64,1},size_iid_cost_val::Int64,iid_cost_value::Array{Float64,1},iid_cost_prob::Array{Float64,1},Exit::Array{Float64,1},
Exitₓ::Array{Float64,1},capital_supply_future_tmp::Float64,
a_prime_fine_local::Array{Float64,2},future_occupation_fine_local::Array{Float64,2},cons_fine_local::Array{Float64,2},
rev_d_fine::Array{Float64,1},rev_dx_fine::Array{Float64,1},rev_xx_fine::Array{Float64,1},k_choice_x_fine::Array{Float64,1},
k_choice_d_fine::Array{Float64,1},labor_x_fine::Array{Float64,1},labor_d_fine::Array{Float64,1},
Profits_x_fine::Array{Float64,1},Profits_d_fine::Array{Float64,1},output_d_fine::Array{Float64,1},
output_dx_fine::Array{Float64,1},output_xx_fine::Array{Float64,1},price_d_fine::Array{Float64,1},
price_dx_fine::Array{Float64,1},price_xx_fine::Array{Float64,1}, distr_prev::Array{Float64,1},lambdda_d_fine::Array{Float64,1},lambdda_x_fine::Array{Float64,1},
ones_tmp_fine::Array{Float64,1},z_tilde_fine::Array{Float64,1},agrid_fine::Array{Float64,1},α₁_eff::Float64,α₂_eff::Float64,capital_supply_past_tmp::Float64 = 0)
    # Calculate aggregates
    # Distribution of current workers:

    worker_past_dist = distr_prev[(ns_fine *0 + 1):(ns_fine *1)];
    domestic_past_dist = distr_prev[(ns_fine *1 + 1):(ns_fine *2)];
    exporter_past_dist = distr_prev[(ns_fine *2 + 1):(ns_fine *3)];

    # Decisions avgs:
    future_occupation_fine_avgs =  Array{SparseMatrixCSC{Float64,Int64}}(undef, 3,3);
    cons_fine_avgs =  Array{SparseMatrixCSC{Float64,Int64}}(undef, 3,1);
    cons_fine_second_moment =  Array{SparseMatrixCSC{Float64,Int64}}(undef, 3,1);
    #here
    for j=1:3
     for jj=1:3
         future_occupation_fine_avgs[j,jj] = spzeros(ns_fine);
         for j_cost = 1:size_iid_cost_val
             future_occupation_fine_avgs[j,jj]=future_occupation_fine_avgs[j,jj] + iid_cost_prob[j_cost]*(
                 future_occupation_fine_local[:,j,j_cost] .== jj);
         end
     end
     cons_fine_avgs[j,1] = spzeros(ns_fine);
     cons_fine_second_moment[j,1] = spzeros(ns_fine);
     for j_cost = 1:size_iid_cost_val
             cons_fine_avgs[j,1]=cons_fine_avgs[j,1] + iid_cost_prob[j_cost
             ]*cons_fine_local[:,j,j_cost];
             cons_fine_second_moment[j,1]=cons_fine_second_moment[j,1] + iid_cost_prob[j_cost
             ]*cons_fine_local[:,j,j_cost].^2.0;
     end
    end
    # Distribution of current workers:
    stay_workers = worker_past_dist.*future_occupation_fine_avgs[1,1];
    exit_domestic = domestic_past_dist.*future_occupation_fine_avgs[2,1];
    exit_exporting_to_work = exporter_past_dist.*future_occupation_fine_avgs[3,1];
    current_workers = stay_workers + exit_domestic + exit_exporting_to_work;

    entrants_domestic_from_workers = worker_past_dist.*future_occupation_fine_avgs[1,2];
    incumbents_domestic = domestic_past_dist.*future_occupation_fine_avgs[2,2];
    exit_exporting_to_domestic = exporter_past_dist.*future_occupation_fine_avgs[3,2];
    current_domestic = entrants_domestic_from_workers + incumbents_domestic + exit_exporting_to_domestic;

    entrants_to_exporters = worker_past_dist.*future_occupation_fine_avgs[1,3
     ] +  domestic_past_dist.*future_occupation_fine_avgs[2,3];
    incumbents_exporter = exporter_past_dist.*future_occupation_fine_avgs[3,3];
    current_exporter = entrants_to_exporters + incumbents_exporter;
    current_distr_store_tmp = copy(distr_prev);
    current_distr_store_tmp[(ns_fine *0 + 1):(ns_fine *1)] = current_workers;
    current_distr_store_tmp[(ns_fine *1 + 1):(ns_fine *2)] = current_domestic;
    current_distr_store_tmp[(ns_fine *2 + 1):(ns_fine *3)] = current_exporter;
    # Fixed cost accounting:
    # entrants for a given cost realization
    entry_costs_payed =  Array{SparseMatrixCSC{Float64,Int64}}(undef, 2,1);;
    entry_costs_payed[1,1] = spzeros(ns_fine);
    entry_costs_payed[2,1]= spzeros(ns_fine);
    for j_cost = 1:size_iid_cost_val
     entry_costs_payed[1,1] = entry_costs_payed[1,1] + iid_cost_prob[j_cost]*iid_cost_value[
         j_cost]*worker_past_dist.*(future_occupation_fine_local[:,1,j_cost] .== 2); # from worker to domestic
     entry_costs_payed[2,1]  = entry_costs_payed[2,1] + iid_cost_prob[j_cost]*iid_cost_value[
         j_cost]*(worker_past_dist.*(future_occupation_fine_local[:,1,j_cost] .== 3) #from worker to exporter
         + domestic_past_dist.*(future_occupation_fine_local[:,2,j_cost] .== 3)); # from domestic to exporter
    end

    domestic_prod_tmp = sum(current_exporter .* output_dx_fine.^((σ - 1)/σ) +
         current_domestic.* output_d_fine.^((σ - 1)/σ));
    export_prod_tmp = sum(current_exporter.* output_xx_fine .^((σ - 1)/σ));

    domestic_price_sum_tmp = sum(current_exporter.* price_dx_fine.^(1 - σ) + current_domestic.* price_d_fine.^(1 - σ));
    export_price_sum_tmp = sum(current_exporter.* price_xx_fine.^(1 - σ)) ;
    export_value_tmp = sum(current_exporter.* price_xx_fine .* output_xx_fine);
    exit_share_to_domestic_tmp = sum(exit_exporting_to_domestic);
    entry_share_to_domestic_tmp =sum(entrants_domestic_from_workers);
    exit_domestic_sum_tmp = sum(exit_domestic);
    exit_exporting_to_work_sum_tmp = sum(exit_exporting_to_work);
    exit_share_to_worker_tmp = exit_domestic_sum_tmp + exit_exporting_to_work_sum_tmp;
    entry_share_to_exporter_tmp = sum(entrants_to_exporters);
    exporter_pop_tmp = sum(current_exporter);
    domestic_pop_tmp = sum(current_domestic);
    worker_pop_tmp =  L_loc - exporter_pop_tmp - domestic_pop_tmp;
    L_x_tmp = sum(current_exporter.*labor_x_fine);
    L_d_tmp =sum(current_domestic.*labor_d_fine);
    total_entry_cost_tmp = sum(entry_costs_payed[1,1])  * F[country,1] + sum(entry_costs_payed[2,1]
     ) * Fₓ[country,1];
    total_incumbents_cost_tmp = exporter_pop_tmp* FMₓ[country,1] + domestic_pop_tmp * FM[country,1];
    total_exit_cost_tmp= exit_domestic_sum_tmp*Exit[country,1] + exit_exporting_to_work_sum_tmp * Exitₓ[country,1];
    labor_demand_tmp =  L_x_tmp + L_d_tmp + total_entry_cost_tmp + total_incumbents_cost_tmp + total_exit_cost_tmp;
    labor_excess_demand_tmp = labor_demand_tmp - worker_pop_tmp;
    K_x_tmp = sum(current_exporter.*k_choice_x_fine);
    K_d_tmp= sum(current_domestic.*k_choice_d_fine);
    capital_demand_tmp = (K_d_tmp + K_x_tmp);

    worker_bond_holding_tmp = s_fine[:,1] .*current_workers;
    worker_bond_holding_sum_tmp = sum(worker_bond_holding_tmp);
    domestic_bond_holding_tmp = (s_fine[:,1] - price_final_prev * k_choice_d_fine) .*current_domestic;
    domestic_bond_holding_sum_tmp = sum(domestic_bond_holding_tmp[domestic_bond_holding_tmp .> 0.0]);
    domestic_firm_debt_tmp = sum(domestic_bond_holding_tmp[domestic_bond_holding_tmp .<0.0]);
    exporter_bond_holding_tmp = (s_fine[:,1]- price_final_prev * k_choice_x_fine) .*current_exporter;
    exporter_bond_holding_sum_tmp = sum(exporter_bond_holding_tmp[exporter_bond_holding_tmp .> 0.0]);
    exporter_firm_debt_tmp = sum(exporter_bond_holding_tmp[exporter_bond_holding_tmp .<0.0]);
    asset_supply_tmp = worker_bond_holding_sum_tmp + domestic_bond_holding_sum_tmp + exporter_bond_holding_sum_tmp;
    asset_demand_tmp = - (domestic_firm_debt_tmp + exporter_firm_debt_tmp);
    NFA_tmp = asset_supply_tmp - (1 + banking_cost_loc) * asset_demand_tmp;
    total_consumption_tmp = sum(worker_past_dist.* cons_fine_avgs[1,1] + domestic_past_dist.* cons_fine_avgs[2,1]
     + exporter_past_dist.* cons_fine_avgs[3,1]);
    Aggr_bank_cost_tmp = ((banking_cost_loc)*asset_demand_tmp)/price_final_current;
    # normalize as it is denominated in foreign final good!
    x_k_tmp = K_x_tmp/capital_demand_tmp; # capital allocated to exports
    x_l_tmp = L_x_tmp/labor_demand_tmp; # labor allocated to exports
    investment_tmp = capital_supply_future_tmp - (1 -δ)* capital_demand_tmp;
    total_demand_final_good_tmp = total_consumption_tmp  + investment_tmp;
    excess_demand_final_good_tmp = total_demand_final_good_tmp - output_final_loc;

    labor_excess_demand_tmp_percent = labor_excess_demand_tmp./L_loc;
    excess_demand_final_good_tmp_percent = excess_demand_final_good_tmp/(output_final_loc + total_demand_final_good_tmp);
    if capital_supply_past_tmp >0
        NFA_tmp = capital_supply_past_tmp-capital_demand_tmp;
    end
    # include: TFP GDP Cons sd_mrpk sd_mrpk_d sd_mrpk_x Top 10 % wealth/income/cons
    # try fraction of exporting firms and zombie firms
    # Do additional calculations/statistics
    mean_MRPK_d_tmp=@fastmath sum(current_domestic.*log.((lambdda_d_fine
    + R*ones_tmp_fine))./domestic_pop_tmp);
    mean_MRPK_x_tmp=@fastmath sum(current_exporter.*log.((lambdda_x_fine
    + R*ones_tmp_fine))./exporter_pop_tmp);
    mean_MRPK_tmp = (domestic_pop_tmp* mean_MRPK_d_tmp +
     exporter_pop_tmp *mean_MRPK_x_tmp)/(
     domestic_pop_tmp + exporter_pop_tmp) ;
    sd_MRPK_d_tmp =@fastmath  (sum(current_domestic.*(log.((lambdda_d_fine
     + R*ones_tmp_fine))).^2/domestic_pop_tmp) -
     mean_MRPK_d_tmp^2)^(1/2);
    sd_MRPK_x_tmp =@fastmath  (sum(current_exporter.*(log.((lambdda_x_fine
     + R*ones_tmp_fine))).^2/sum(current_exporter) )- mean_MRPK_x_tmp.^2)^(1/2);
    sd_MRPK_tmp =@fastmath ((sum(current_domestic.*(log.((lambdda_d_fine
     + R*ones_tmp_fine))).^2 )+ sum(current_exporter.*(log.((lambdda_x_fine + R*ones_tmp_fine))).^2))/
     (domestic_pop_tmp + exporter_pop_tmp)
     - mean_MRPK_tmp^2)^(1/2);
    mean_logProd_d_tmp = @fastmath sum((current_domestic.*log.(z_tilde_fine))/domestic_pop_tmp);
    mean_logProd_x_tmp = @fastmath sum((current_exporter.*log.(z_tilde_fine))/exporter_pop_tmp);
    mean_logProd_tmp =@fastmath sum((current_exporter + current_domestic).*log.(z_tilde_fine))/(exporter_pop_tmp + domestic_pop_tmp);

    sd_logProd_d_tmp =@fastmath (sum((current_domestic.*(log.(z_tilde_fine)).^2)/domestic_pop_tmp) - mean_logProd_d_tmp^2)^(1/2);
    sd_logProd_x_tmp =@fastmath (sum((current_exporter.*(log.(z_tilde_fine)).^2)/exporter_pop_tmp) - mean_logProd_x_tmp^2)^(1/2);
    sd_logProd_tmp =@fastmath (sum((current_exporter + current_domestic).*(log.(z_tilde_fine)).^2)/(exporter_pop_tmp + domestic_pop_tmp)
    - mean_logProd_tmp^2)^(1/2);

    cov_TFPR_z_d_tmp = @fastmath sum(current_domestic.*(α₁_eff* log.(z_tilde_fine).*log.((lambdda_d_fine +
    R*ones_tmp_fine))))/domestic_pop_tmp- α₁_eff*mean_MRPK_d_tmp * mean_logProd_d_tmp;
    cov_TFPR_z_x_tmp = @fastmath sum(current_exporter.*(α₁_eff* log.(z_tilde_fine).*log.(lambdda_x_fine +
    R*ones_tmp_fine)))/exporter_pop_tmp - α₁_eff*mean_MRPK_x_tmp * mean_logProd_x_tmp;
    cov_TFPR_z_tmp =@fastmath ((sum( current_exporter.*(log.((lambdda_x_fine + R*ones_tmp_fine)) .* α₁_eff.* log.(z_tilde_fine) )
    ) + sum( current_domestic.*(log.((lambdda_d_fine + R*ones_tmp_fine)).* α₁_eff.* log.(z_tilde_fine) )))
    )/ (domestic_pop_tmp + exporter_pop_tmp) - α₁_eff*mean_MRPK_tmp * mean_logProd_tmp;
    corr_TFPR_z_d_tmp = cov_TFPR_z_d_tmp / (α₁_eff* sd_MRPK_d_tmp* sd_logProd_d_tmp );
    corr_TFPR_z_x_tmp = cov_TFPR_z_x_tmp / (α₁_eff* sd_MRPK_x_tmp* sd_logProd_x_tmp );
    corr_TFPR_z_tmp = cov_TFPR_z_tmp / (α₁_eff* sd_MRPK_tmp* sd_logProd_tmp );

    # Average profits:
    avg_Pi_x_tmp = sum(current_exporter.*Profits_x_fine)/exporter_pop_tmp;
    avg_Pi_d_tmp = sum(current_domestic.*Profits_d_fine)/domestic_pop_tmp;
    # Save values for the construction of TFP.
    #D_d_denom_tmp =  sum(current_domestic.*(z_tilde_fine.^(σ-1) .*(lambdda_d_fine + R*ones_tmp_fine).^(-α*(σ-1))));
    #D_x_denom_tmp =  sum(current_exporter.*(z_tilde_fine.^(σ-1) .*(lambdda_x_fine + R*ones_tmp_fine).^(-α*(σ-1))));
    #D_d_tmp = sum(current_domestic.*(z_tilde_fine.^σ .*(lambdda_d_fine + R*ones_tmp_fine).^(-α*(σ-1))));
    #D_x_tmp = sum(current_exporter.*(z_tilde_fine.^σ .*(lambdda_x_fine + R*ones_tmp_fine).^(-α*(σ-1))));
    #D_d_denom_eff_tmp =  sum(current_domestic.*(z_tilde_fine.^(σ-1) .*(R).^(-α*(σ-1))));
    #D_x_denom_eff_tmp =  sum(current_exporter.*(z_tilde_fine.^(σ-1) .*(R).^(-α*(σ-1))));
    #D_d_eff_tmp = sum(current_domestic.*(z_tilde_fine.^σ .*(R).^(-α*(σ-1))));
    #D_x_eff_tmp = sum(current_exporter.*(z_tilde_fine.^σ .*(R).^(-α*(σ-1))));
    D_d_denom_tmp =  sum(current_domestic.*(z_tilde_fine.^σ .*(lambdda_d_fine + R*ones_tmp_fine).^((α₂_eff-1)*σ)));
    D_d_tmp =  sum(current_domestic.*(z_tilde_fine.^σ .*(lambdda_d_fine + R*ones_tmp_fine).^(-α₁_eff*σ)));
    #D_d_denom_tmp = D_d_K_tmp.^α₁_eff*D_d_tmp^α₂_eff;
    D_x_denom_tmp =  sum(current_exporter.*(z_tilde_fine.^σ .*(lambdda_x_fine + R*ones_tmp_fine).^((α₂_eff-1)*σ)));
    D_x_tmp =  sum(current_exporter.*(z_tilde_fine.^σ .*(lambdda_x_fine + R*ones_tmp_fine).^(-α₁_eff*σ)));
    #D_x_denom_tmp = D_x_K_tmp.^α₁_eff*D_x_tmp^α₂_eff;
    D_d_denom_eff_tmp =  sum(current_domestic.*(z_tilde_fine.^σ .*R^((α₂_eff-1)*σ)));
    D_d_eff_tmp =  sum(current_domestic.*(z_tilde_fine.^σ*R^(-α₁_eff*σ)));
    #D_d_denom_eff_tmp = D_d_K_eff_tmp.^α₁_eff*D_d_eff_tmp^α₂_eff;
    D_x_denom_eff_tmp =  sum(current_exporter.*(z_tilde_fine.^σ .*R^((α₂_eff-1)*σ)));
    D_x_eff_tmp =  sum(current_exporter.*(z_tilde_fine.^σ .*R^(-α₁_eff*σ)));
    #D_x_denom_eff_tmp = D_x_K_eff_tmp.^α₁_eff*D_x_eff_tmp^α₂_eff;
    #D_d_denom_tmp =  sum(current_domestic.*(z_tilde_fine.^σ .*(lambdda_d_fine + R*ones_tmp_fine).^((α₂_eff-1)*σ)));
    #D_x_denom_tmp =  sum(current_exporter.*(z_tilde_fine.^σ .*(lambdda_x_fine + R*ones_tmp_fine).^((α₂_eff-1)*σ)));
    #D_x_tmp = sum(current_exporter.*(z_tilde_fine.^σ .*(lambdda_x_fine + R*ones_tmp_fine).^(-α₁_eff*σ)));
    #D_d_tmp = sum(current_domestic.*(z_tilde_fine.^σ .*(lambdda_d_fine + R*ones_tmp_fine).^(-α₁_eff*σ)));
    #D_d_denom_eff_tmp =  sum(current_domestic.*(z_tilde_fine.^σ .*(R).^((α₂_eff-1)*σ)));
    #D_x_denom_eff_tmp =  sum(current_exporter'*(z_tilde_fine.^σ .*(R).^((α₂_eff-1)*σ)));
    #D_x_eff_tmp = sum(current_exporter'*(z_tilde_fine.^σ .*(R).^(-α₁_eff*σ)));
    #D_d_eff_tmp = sum(current_domestic'*(z_tilde_fine.^σ .*(R).^(-α₁_eff*σ)));
    #Inequality

    income_1_fine = W_loc* ones_tmp_fine;
    income_2_fine = Profits_d_fine - W_loc * FM[country,1] * ones_tmp_fine;
    income_3_fine = Profits_x_fine - W_loc * FMₓ[country,1]* ones_tmp_fine;

    mean_cons_tmp = total_consumption_tmp/ L_loc;
    sd_cons_tmp = (sum(worker_past_dist.*cons_fine_second_moment[1,1] + domestic_past_dist.*cons_fine_second_moment[2,1] +
    exporter_past_dist.*cons_fine_second_moment[3,1])/L_loc - mean_cons_tmp^2)^(1/2);
    mean_income_tmp = sum(current_workers.*income_1_fine + current_domestic.*income_2_fine +
    current_exporter.*income_3_fine)/L_loc;
    sd_income_tmp = (sum(current_workers.*income_1_fine.^2 + current_domestic.*income_2_fine.^2 +
    current_exporter.*income_3_fine.^2)/L_loc - mean_income_tmp^2)^(1/2);
    mean_wealth_tmp = sum(worker_past_dist.*s_fine[:,1] + domestic_past_dist.*s_fine[:,1] +
    exporter_past_dist.*s_fine[:,1])/L_loc;
    sd_wealth_tmp = (sum(worker_past_dist.*s_fine[:,1].^2 + domestic_past_dist.*s_fine[:,1].^2 +
    exporter_past_dist.*s_fine[:,1].^2)/L_loc - mean_wealth_tmp^2)^(1/2);
    a_grid_size = size(agrid_fine)[1];
    z_grid_size = convert(Int64,ns_fine/size(agrid_fine)[1]);
    worker_past_dist_mat = reshape(worker_past_dist,a_grid_size,z_grid_size);
    domestic_past_dist_mat = reshape(domestic_past_dist,a_grid_size,z_grid_size);
    exporter_past_dist_mat = reshape(exporter_past_dist,a_grid_size,z_grid_size);
    dist0 = worker_past_dist_mat + domestic_past_dist_mat + exporter_past_dist_mat;
    dist1 = sum(dist0,dims = 2)/L_loc;
    cumu_wealth = cumsum(dist1,dims = 1);
    p10_index_wealth = findfirst(cumu_wealth.>0.1);
    p50_index_wealth = findfirst(cumu_wealth.>0.5);
    p90_index_wealth = findfirst(cumu_wealth.>0.9);
    p99_index_wealth = findfirst(cumu_wealth.>0.99);
    p10_wealth_tmp = sum(dist1[p10_index_wealth[1]:end].*agrid_fine[p10_index_wealth[1]:end])/mean_wealth_tmp;
    p50_wealth_tmp = sum(dist1[p50_index_wealth[1]:end].*agrid_fine[p50_index_wealth[1]:end])/mean_wealth_tmp;
    p90_wealth_tmp = sum(dist1[p90_index_wealth[1]:end].*agrid_fine[p90_index_wealth[1]:end])/mean_wealth_tmp;
    p99_wealth_tmp = sum(dist1[p99_index_wealth[1]:end].*agrid_fine[p99_index_wealth[1]:end])/mean_wealth_tmp;
    wealth_of_workers_tmp = sum(worker_past_dist.*s_fine[:,1])/mean_wealth_tmp/L_loc;
    wealth_of_domestic_tmp = sum(domestic_past_dist.*s_fine[:,1])/mean_wealth_tmp/L_loc;
    wealth_of_exporters_tmp = sum(exporter_past_dist.*s_fine[:,1])/mean_wealth_tmp/L_loc;

    income_policy = zeros(3*ns_fine);
    income_policy[1:ns_fine] = income_1_fine;
    income_policy[(ns_fine+1):(2*ns_fine)] = income_2_fine;
    income_policy[(2*ns_fine + 1):(3*ns_fine)] = income_3_fine;
    income_sort_index = sortperm(income_policy);
    income_sorted = income_policy[income_sort_index];
    current_distr_income_sorted = current_distr_store_tmp[income_sort_index]/L_loc;
    cumu_income_distr =  cumsum(current_distr_income_sorted,dims = 1);
    p10_index_income = findfirst(cumu_income_distr.>0.1);
    p50_index_income = findfirst(cumu_income_distr.>0.5);
    p90_index_income = findfirst(cumu_income_distr.>0.9);
    p99_index_income = findfirst(cumu_income_distr.>0.99);
    p10_income_tmp = sum(current_distr_income_sorted[p10_index_income[1]:end].*income_sorted[p10_index_income[1]:end])/mean_income_tmp;
    p50_income_tmp = sum(current_distr_income_sorted[p50_index_income[1]:end].*income_sorted[p50_index_income[1]:end])/mean_income_tmp;
    p90_income_tmp = sum(current_distr_income_sorted[p90_index_income[1]:end].*income_sorted[p90_index_income[1]:end])/mean_income_tmp;
    p99_income_tmp = sum(current_distr_income_sorted[p99_index_income[1]:end].*income_sorted[p99_index_income[1]:end])/mean_income_tmp;
    income_of_workers_tmp = sum(current_workers.*income_1_fine)/mean_income_tmp/L_loc;
    income_of_domestic_tmp = sum(current_domestic.*income_2_fine)/mean_income_tmp/L_loc;
    income_of_exporters_tmp = sum(current_exporter.*income_3_fine)/mean_income_tmp/L_loc;

    cons_policy = zeros(3*ns_fine);
    cons_policy[1:ns_fine] = cons_fine_avgs[1,1];
    cons_policy[(ns_fine+1):(2*ns_fine)] = cons_fine_avgs[2,1];
    cons_policy[(2*ns_fine+1):(3*ns_fine)] = cons_fine_avgs[3,1];
    cons_sort_index = sortperm(cons_policy);
    cons_sorted = cons_policy[cons_sort_index];
    current_distr_cons_sorted = distr_prev[cons_sort_index]/L_loc;
    cumu_cons_distr =  cumsum(current_distr_cons_sorted,dims = 1);
    p10_index_cons = findfirst(cumu_cons_distr.>0.1);
    p50_index_cons = findfirst(cumu_cons_distr.>0.5);
    p90_index_cons = findfirst(cumu_cons_distr.>0.9);
    p99_index_cons = findfirst(cumu_cons_distr.>0.99);
    p10_cons_tmp = sum(current_distr_cons_sorted[p10_index_cons[1]:end].*cons_sorted[p10_index_cons[1]:end])/mean_cons_tmp;
    p50_cons_tmp = sum(current_distr_cons_sorted[p50_index_cons[1]:end].*cons_sorted[p50_index_cons[1]:end])/mean_cons_tmp;
    p90_cons_tmp = sum(current_distr_cons_sorted[p90_index_cons[1]:end].*cons_sorted[p90_index_cons[1]:end])/mean_cons_tmp;
    p99_cons_tmp = sum(current_distr_cons_sorted[p99_index_cons[1]:end].*cons_sorted[p99_index_cons[1]:end])/mean_cons_tmp;
    cons_of_workers_tmp = sum(worker_past_dist.*cons_fine_avgs[1,1])/mean_cons_tmp/L_loc;
    cons_of_domestic_tmp = sum(domestic_past_dist.*cons_fine_avgs[2,1])/mean_cons_tmp/L_loc;
    cons_of_exporters_tmp = sum(exporter_past_dist.*cons_fine_avgs[3,1])/mean_cons_tmp/L_loc;
    # Zombies - experimental
    making_losses = (income_1_fine .>income_3_fine);
    exit_zombie = sum(((future_occupation_fine_local[:,3,1] .!=3.0) .*making_losses).*exporter_past_dist)/exporter_pop_tmp;
    fraction_zombie_exporter_tmp = sum(making_losses.*exporter_past_dist)/exporter_pop_tmp - exit_zombie
    return (labor_excess_demand_tmp_percent,excess_demand_final_good_tmp_percent,
    total_demand_final_good_tmp,export_price_sum_tmp, domestic_price_sum_tmp,NFA_tmp,asset_supply_tmp,asset_demand_tmp,capital_demand_tmp,
    p90_wealth_tmp,p90_cons_tmp,p90_income_tmp,fraction_zombie_exporter_tmp,domestic_prod_tmp, export_prod_tmp,export_value_tmp,
    sd_MRPK_d_tmp,sd_MRPK_x_tmp,sd_MRPK_tmp,D_d_denom_tmp, D_d_tmp, D_x_denom_tmp, D_x_tmp, D_d_denom_eff_tmp, D_d_eff_tmp, D_x_denom_eff_tmp, D_x_eff_tmp,
    domestic_pop_tmp,exporter_pop_tmp)
end
end
function Residual_transition_sequential(price_trans_actual::Array{Float64,2},capital_trans::Array{Float64,2},price_trans::Array{Float64,2}
    ,distr_store::Array{Float64,3},T::Int64,parameter_end::Parameter_type,coeff_store::Array{Float64,4},τ_trans::Array{Float64,2},
    openness_transition::Array{Int64,1})
    coeff_store_tmp = convert(SharedArray,coeff_store);
    distr_store_tmp = convert(SharedArray,distr_store);
    capital_supply_future = convert(SharedArray,capital_trans);
    #coeff_store_tmp = copy(coeff_store);
    #distr_store_tmp = copy(distr_store);
    #capital_supply_future = copy(capital_trans);
    # Extract the necessary local parameters:
    (β,α,δ,θ,α₁,α₂,σ,α₁_eff,α₂_eff,ω,L,FM,FMₓ,F,Fₓ,
        Exit,Exitₓ,iid_cost_value,iid_cost_prob,size_iid_cost_val,country_no,τ,
        a_min,a_max,fspace_a,fspace_a_fine,agrid_fine,banking_cost,
        bounds,Country_spec_p,s_cell,ns_cell,s_fine_cell,ns_fine_cell,Phi_z_cell,
        Phi_z_fine_cell,Phi_cell,Phi_aug_cell,P_kron_cell,P_kron1_cell,P_kron_fine_cell,
        exp_egrid_cell,ns_tmp,ns_tmp_fine,openness) = local_parameters(parameter_end);
    # Policy function
    (a_prime_fine,future_occupation_fine,cons_fine,rev_d_fine_store,
        rev_dx_fine_store,rev_xx_fine_store,k_x_fine,k_d_fine,
        l_x_fine,l_d_fine,Profit_fine_x,Profit_fine_d,output_d_fine_store,output_dx_fine_store,output_xx_fine_store,
        price_d_fine_store,price_dx_fine_store,price_xx_fine_store,labor_excess_demand_store_percent,
        excess_demand_final_good_store_percent,total_demand_final_good,export_price_sum_store,domestic_price_sum_store,NFA_store,
        asset_supply_store,asset_demand_store) =local_var_creator_policy(country_no,ns_tmp,ns_tmp_fine,size_iid_cost_val,T);
    #Precheck the prices
    price_check_sum = 0;
    for country = 1:country_no
        # Initialize each country, only once for each
        (exitflag_tmp,s,ns,s_fine,ns_fine,Phi_z,Phi_z_fine,Phi,Phi_aug,
        P_kron,P_kron1,P_kron_fine,z_tilde,z_tilde_fine,income_mat,income_mat_fine,
        a_prime_fine_local,future_occupation_fine_local,cons_fine_local,coeff,
        coeff_next,θ_loc,L_loc,τ_loc,banking_cost_loc,
        ones_tmp,x_tmp,V_tmp,β_loc,V_next_stacked,iterate1,Phi_prime_tmp,
        D_deriv_tmp_block,conv,Q_trans,ones_tmp_fine) = country_local_initialize_noprice(country,s_cell,
        ns_cell,s_fine_cell,ns_fine_cell,Phi_z_cell,Phi_z_fine_cell,Phi_cell,Phi_aug_cell,
        P_kron_cell,P_kron1_cell,P_kron_fine_cell,σ,size_iid_cost_val,θ, L, τ,country_no,banking_cost,δ,ω,β);
        for t in (T+1):-1:2
            # Prices and other time dependent quantities:
            τ_loc = τ_trans[country, t];
            openness = openness_transition[t];
            (price_final_prev,price_final_current,total_foreign_demand,W_loc,r_loc,avg_foreign_price,R,constant_d,
            constant_x,constant_z,constant_z_bar,constant_z_fine,constant_z_bar_fine,output_final_loc,
            residual,price_check_tmp) = price_reshaper_single(price_trans_actual[:,t],price_trans_actual[:,t-1],openness,country_no,
                bounds,country,banking_cost_loc,δ,ω,σ,τ_loc,z_tilde,z_tilde_fine);
            price_check_sum = price_check_tmp + price_check_sum;
        end
    end
    if price_check_sum>0
        println("Guess out of bounds")
        residual_store = 1000*ones((4*country_no),T);
        price_final_actual = 1000*ones(2,T);
        total_demand_final_good = 1000*ones(2,T);
    else
        #Enter the loop
        #Backward loop:
        @sync @distributed for country = 1:country_no #
            # Initialize each country, only once for each
            (exitflag_tmp,s,ns,s_fine,ns_fine,Phi_z,Phi_z_fine,Phi,Phi_aug,
            P_kron,P_kron1,P_kron_fine,z_tilde,z_tilde_fine,income_mat,income_mat_fine,
            a_prime_fine_local,future_occupation_fine_local,cons_fine_local,coeff,
            coeff_next,θ_loc,L_loc,τ_loc,banking_cost_loc,
            ones_tmp,x_tmp,V_tmp,β_loc,V_next_stacked,iterate1,Phi_prime_tmp,
            D_deriv_tmp_block,conv,Q_trans,ones_tmp_fine) = country_local_initialize_noprice(country,s_cell,
            ns_cell,s_fine_cell,ns_fine_cell,Phi_z_cell,Phi_z_fine_cell,Phi_cell,Phi_aug_cell,
            P_kron_cell,P_kron1_cell,P_kron_fine_cell,σ,size_iid_cost_val,θ, L, τ,country_no,banking_cost,δ,ω,β);
            for t in (T+1):-1:2
            # Prices and other time dependent quantities:
            τ_loc = τ_trans[country, t];
            openness = openness_transition[t];
            (price_final_prev,price_final_current,total_foreign_demand,W_loc,r_loc,avg_foreign_price,R,constant_d,
            constant_x,constant_z,constant_z_bar,constant_z_fine,constant_z_bar_fine,output_final_loc,
            residual,price_check_tmp) = price_reshaper_single(price_trans_actual[:,t],price_trans_actual[:,t-1],openness,country_no,
                bounds,country,banking_cost_loc,δ,ω,σ,τ_loc,z_tilde,z_tilde_fine);
            # Get the coefficients
            coeff_store_tmp[:,:,country,t] = Residual_transition_backward(
                country,price_final_prev,price_final_current, R,W_loc,r_loc ,s,
                constant_x,constant_d,constant_z,constant_z_bar,
                β_loc,θ_loc,τ_loc,σ,α₁_eff,α₂_eff,α₁,α₂,
                ones_tmp,ns,FM,FMₓ,V_tmp,
                F, Fₓ,size_iid_cost_val,iid_cost_value,iid_cost_prob,Exit,
                Exitₓ,a_min,a_max,α,fspace_a,
                fspace_a_fine,openness,agrid_fine,coeff_store_tmp[:,:,country,t+1],
                P_kron,Phi,Phi_z,income_mat,x_tmp,
                coeff_next);
            end
            #Forward loop
            for t in 2:(T+1)
                #println(t)
                τ_loc = τ_trans[country, t];
                openness = openness_transition[t];
                (price_final_prev,price_final_current,total_foreign_demand,W_loc,r_loc,avg_foreign_price,R,constant_d,
                constant_x,constant_z,constant_z_bar,constant_z_fine,constant_z_bar_fine,output_final_loc,
                residual,price_check_tmp) = price_reshaper_single(price_trans_actual[:,t],price_trans_actual[:,t-1],openness,country_no,
                    bounds,country,banking_cost_loc,δ,ω,σ,τ_loc,z_tilde,z_tilde_fine);
                (a_prime_fine[:,:,country,1,t-1],future_occupation_fine[:,:,country,1,t-1],cons_fine[:,:,country,1,t-1],rev_d_fine_store[:,country,t-1],
                    rev_dx_fine_store[:,country,t-1],rev_xx_fine_store[:,country,t-1],k_x_fine[:,country,t-1],k_d_fine[:,country,t-1],
                    l_x_fine[:,country,t-1],l_d_fine[:,country,t-1],Profit_fine_x[:,country,t-1],Profit_fine_d[:,country,t-1],output_d_fine_store[:,country,t-1],
                    output_dx_fine_store[:,country,t-1],output_xx_fine_store[:,country,t-1],
                    price_d_fine_store[:,country,t-1],price_dx_fine_store[:,country,t-1],price_xx_fine_store[:,country,t-1],distr_store_tmp[:,country,t],lambdda_d_fine_tmp,lambdda_x_fine_tmp) = Residual_transition_forward(country,price_final_prev,price_final_current, R,W_loc,r_loc ,s_fine,
                    constant_x,constant_d,constant_z_fine,constant_z_bar_fine,β_loc,θ_loc,τ_loc,σ,α₁_eff,α₂_eff,α₁,α₂,
                    ones_tmp_fine,ns_fine,FM,FMₓ,V_tmp,
                    F, Fₓ,size_iid_cost_val,iid_cost_value,iid_cost_prob,Exit,
                    Exitₓ,a_min,a_max,α,fspace_a,
                    fspace_a_fine,openness,agrid_fine,coeff_store_tmp[:,:,country,t+1],
                    P_kron_fine,P_kron1,Phi,Phi_z_fine,income_mat_fine,x_tmp,
                     cons_fine_local,a_prime_fine_local,future_occupation_fine_local,Q_trans,distr_store_tmp[:,country,t-1]);
            end
            #Second backward loop for the aggregates
            for t in (T+1):-1:2
                τ_loc = τ_trans[country, t];
                openness = openness_transition[t];
                (price_final_prev,price_final_current,total_foreign_demand,W_loc,r_loc,avg_foreign_price,R,constant_d,
                constant_x,constant_z,constant_z_bar,constant_z_fine,constant_z_bar_fine,output_final_loc,
                residual,price_check_tmp) = price_reshaper_single(price_trans_actual[:,t],price_trans_actual[:,t-1],openness,country_no,
                    bounds,country,banking_cost_loc,δ,ω,σ,τ_loc,z_tilde,z_tilde_fine);
                (labor_excess_demand_store_percent[country,t-1],excess_demand_final_good_store_percent[country,t-1],
                total_demand_final_good[country,t-1],
                export_price_sum_store[country,t-1], domestic_price_sum_store[country,t-1],NFA_store[country,t-1],
                asset_supply_store[country,t-1],asset_demand_store[country,t-1],capital_supply_future[country,t]
                )= Residual_transition_backward_aggregates(country,output_final_loc,price_final_prev,price_final_current, R,W_loc,r_loc ,L_loc,s_fine,
                σ,ns_fine,FM,FMₓ,banking_cost_loc,δ,F, Fₓ,size_iid_cost_val,iid_cost_value,iid_cost_prob,Exit,
                Exitₓ,capital_supply_future[country,t+1],a_prime_fine[:,:,country,1,t-1],future_occupation_fine[:,:,country,1,t-1],cons_fine[:,:,country,1,t-1],
                rev_d_fine_store[:,country,t-1],rev_dx_fine_store[:,country,t-1],rev_xx_fine_store[:,country,t-1],k_x_fine[:,country,t-1],
                k_d_fine[:,country,t-1],l_x_fine[:,country,t-1],l_d_fine[:,country,t-1],
                Profit_fine_x[:,country,t-1],Profit_fine_d[:,country,t-1],output_d_fine_store[:,country,t-1],
                output_dx_fine_store[:,country,t-1],output_xx_fine_store[:,country,t-1],
                price_d_fine_store[:,country,t-1],price_dx_fine_store[:,country,t-1],price_xx_fine_store[:,country,t-1], distr_store_tmp[:,country,t-1],capital_supply_future[country,t])
            end
        end
        import_price = zeros(country_no,T);
        price_final_actual = zeros(country_no,T);
        residual_store = zeros(4*country_no,T);
        for country = 1:country_no
            (exitflag_tmp,s,ns,s_fine,ns_fine,Phi_z,Phi_z_fine,Phi,Phi_aug,
            P_kron,P_kron1,P_kron_fine,z_tilde,z_tilde_fine,income_mat,income_mat_fine,
            a_prime_fine_local,future_occupation_fine_local,cons_fine_local,coeff,
            coeff_next,θ_loc,L_loc,τ_loc,banking_cost_loc,
            ones_tmp,x_tmp,V_tmp,β_loc,V_next_stacked,iterate1,Phi_prime_tmp,
            D_deriv_tmp_block,conv,Q_trans,ones_tmp_fine) = country_local_initialize_noprice(country,s_cell,
            ns_cell,s_fine_cell,ns_fine_cell,Phi_z_cell,Phi_z_fine_cell,Phi_cell,Phi_aug_cell,
            P_kron_cell,P_kron1_cell,P_kron_fine_cell,σ,size_iid_cost_val,θ, L, τ,country_no,banking_cost,δ,ω,β);
            for t in 1:T
                τ_loc = τ_trans[country, t+1];
                openness = openness_transition[t+1];
                residual_store[country,t] = labor_excess_demand_store_percent[country,t];
                residual_store[country_no + country,t] = excess_demand_final_good_store_percent[country,t];
                import_price[country,t] = (sum(export_price_sum_store[:,t])-  export_price_sum_store[country,t])/(country_no - 1);
                price_final_actual[country,t] = min((ω^σ * domestic_price_sum_store[country,t] + (1.0 - ω)^σ  * import_price[country,t])^(1.0
                    / (1.0 - σ)),10000);
                (price_final_prev,price_final_current,total_foreign_demand,W_loc,r_loc,avg_foreign_price,R,constant_d,
                constant_x,constant_z,constant_z_bar,constant_z_fine,constant_z_bar_fine,output_final_loc,
                residual,price_check_tmp) = price_reshaper_single(price_trans_actual[:,t+1],price_trans_actual[:,t],openness,country_no,
                    bounds,country,banking_cost_loc,δ,ω,σ,τ_loc,z_tilde,z_tilde_fine);

                residual_store[2* country_no + country,t] = (price_final_actual[country,t] - price_final_current)./(
                    price_final_actual[country,t] + price_final_current);
                if (openness == 0)
                    residual_store[3* country_no + country,t] = NFA_store[country,t]/(asset_demand_store[country,t] + (1 + banking_cost[country])*asset_supply_store[country,t]);
                end
            end
        end
        if openness == 1
            for t in 2:T
                residual_store[7,t] = sum(NFA_store[:,t])/sum((asset_demand_store[:,t]+ (ones(country_no) + banking_cost).*asset_supply_store[:,t]));
            end
        end
    end
    return residual_store,price_final_actual,total_demand_final_good
end
function Residual_transition_sequential_detailed(price_trans_actual::Array{Float64,2},capital_trans::Array{Float64,2},price_trans::Array{Float64,2}
    ,distr_store::Array{Float64,3},T::Int64,parameter_end::Parameter_type,coeff_store::Array{Float64,4},τ_trans::Array{Float64,2},
    openness_transition::Array{Int64,1})
    coeff_store_tmp = convert(SharedArray,coeff_store);
    distr_store_tmp = convert(SharedArray,distr_store);
    capital_supply_future = convert(SharedArray,capital_trans);
    #coeff_store_tmp = copy(coeff_store);
    #distr_store_tmp = copy(distr_store);
    #capital_supply_future = copy(capital_trans);
    # Extract the necessary local parameters:
    (β,α,δ,θ,α₁,α₂,σ,α₁_eff,α₂_eff,ω,L,FM,FMₓ,F,Fₓ,
        Exit,Exitₓ,iid_cost_value,iid_cost_prob,size_iid_cost_val,country_no,τ,
        a_min,a_max,fspace_a,fspace_a_fine,agrid_fine,banking_cost,
        bounds,Country_spec_p,s_cell,ns_cell,s_fine_cell,ns_fine_cell,Phi_z_cell,
        Phi_z_fine_cell,Phi_cell,Phi_aug_cell,P_kron_cell,P_kron1_cell,P_kron_fine_cell,
        exp_egrid_cell,ns_tmp,ns_tmp_fine,openness) = local_parameters(parameter_end);
    # Policy function
    (a_prime_fine,future_occupation_fine,cons_fine,rev_d_fine_store,
        rev_dx_fine_store,rev_xx_fine_store,k_x_fine,k_d_fine,
        l_x_fine,l_d_fine,Profit_fine_x,Profit_fine_d,output_d_fine_store,output_dx_fine_store,output_xx_fine_store,
        price_d_fine_store,price_dx_fine_store,price_xx_fine_store,labor_excess_demand_store_percent,
        excess_demand_final_good_store_percent,total_demand_final_good,export_price_sum_store,domestic_price_sum_store,NFA_store,
        asset_supply_store,asset_demand_store,lambdda_d_fine_store,lambdda_x_fine_store,p90_wealth_store,p90_cons_store,p90_income_store,
        fraction_zombie_exporter_store,domestic_prod_store,export_prod_store,export_value_store,sd_MRPK_d_store,sd_MRPK_x_store,sd_MRPK_store,
        D_d_denom_store,D_d_store,D_x_denom_store,D_x_store,D_d_denom_eff_store,D_d_eff_store,D_x_denom_eff_store,D_x_eff_store,domestic_pop_store,
        exporter_pop_store) =local_var_creator_policy_detailed(country_no,ns_tmp,ns_tmp_fine,size_iid_cost_val,T);
    #Precheck the prices
    price_check_sum = 0;
    for country = 1:country_no
        # Initialize each country, only once for each
        (exitflag_tmp,s,ns,s_fine,ns_fine,Phi_z,Phi_z_fine,Phi,Phi_aug,
        P_kron,P_kron1,P_kron_fine,z_tilde,z_tilde_fine,income_mat,income_mat_fine,
        a_prime_fine_local,future_occupation_fine_local,cons_fine_local,coeff,
        coeff_next,θ_loc,L_loc,τ_loc,banking_cost_loc,
        ones_tmp,x_tmp,V_tmp,β_loc,V_next_stacked,iterate1,Phi_prime_tmp,
        D_deriv_tmp_block,conv,Q_trans,ones_tmp_fine) = country_local_initialize_noprice(country,s_cell,
        ns_cell,s_fine_cell,ns_fine_cell,Phi_z_cell,Phi_z_fine_cell,Phi_cell,Phi_aug_cell,
        P_kron_cell,P_kron1_cell,P_kron_fine_cell,σ,size_iid_cost_val,θ, L, τ,country_no,banking_cost,δ,ω,β);
        for t in (T+1):-1:2
            # Prices and other time dependent quantities:
            τ_loc = τ_trans[country, t];
            openness = openness_transition[t];
            (price_final_prev,price_final_current,total_foreign_demand,W_loc,r_loc,avg_foreign_price,R,constant_d,
            constant_x,constant_z,constant_z_bar,constant_z_fine,constant_z_bar_fine,output_final_loc,
            residual,price_check_tmp) = price_reshaper_single(price_trans_actual[:,t],price_trans_actual[:,t-1],openness,country_no,
                bounds,country,banking_cost_loc,δ,ω,σ,τ_loc,z_tilde,z_tilde_fine);
            price_check_sum = price_check_tmp + price_check_sum;
        end
    end
    price_check_sum = 0;
    if price_check_sum>0
        println("Guess out of bounds")
        residual_store = 1000*ones((4*country_no),T);
        price_final_actual = 1000*ones(2,T);
        total_demand_final_good = 1000*ones(2,T);
    else
        #Enter the loop
        #Backward loop:
        @sync @distributed for country = 1:country_no #
            # Initialize each country, only once for each
            (exitflag_tmp,s,ns,s_fine,ns_fine,Phi_z,Phi_z_fine,Phi,Phi_aug,
            P_kron,P_kron1,P_kron_fine,z_tilde,z_tilde_fine,income_mat,income_mat_fine,
            a_prime_fine_local,future_occupation_fine_local,cons_fine_local,coeff,
            coeff_next,θ_loc,L_loc,τ_loc,banking_cost_loc,
            ones_tmp,x_tmp,V_tmp,β_loc,V_next_stacked,iterate1,Phi_prime_tmp,
            D_deriv_tmp_block,conv,Q_trans,ones_tmp_fine) = country_local_initialize_noprice(country,s_cell,
            ns_cell,s_fine_cell,ns_fine_cell,Phi_z_cell,Phi_z_fine_cell,Phi_cell,Phi_aug_cell,
            P_kron_cell,P_kron1_cell,P_kron_fine_cell,σ,size_iid_cost_val,θ, L, τ,country_no,banking_cost,δ,ω,β);
            for t in (T+1):-1:2
            # Prices and other time dependent quantities:
            τ_loc = τ_trans[country, t];
            openness = openness_transition[t];
            (price_final_prev,price_final_current,total_foreign_demand,W_loc,r_loc,avg_foreign_price,R,constant_d,
            constant_x,constant_z,constant_z_bar,constant_z_fine,constant_z_bar_fine,output_final_loc,
            residual,price_check_tmp) = price_reshaper_single(price_trans_actual[:,t],price_trans_actual[:,t-1],openness,country_no,
                bounds,country,banking_cost_loc,δ,ω,σ,τ_loc,z_tilde,z_tilde_fine);
            # Get the coefficients
            coeff_store_tmp[:,:,country,t] = Residual_transition_backward(
                country,price_final_prev,price_final_current, R,W_loc,r_loc ,s,
                constant_x,constant_d,constant_z,constant_z_bar,
                β_loc,θ_loc,τ_loc,σ,α₁_eff,α₂_eff,α₁,α₂,
                ones_tmp,ns,FM,FMₓ,V_tmp,
                F, Fₓ,size_iid_cost_val,iid_cost_value,iid_cost_prob,Exit,
                Exitₓ,a_min,a_max,α,fspace_a,
                fspace_a_fine,openness,agrid_fine,coeff_store_tmp[:,:,country,t+1],
                P_kron,Phi,Phi_z,income_mat,x_tmp,
                coeff_next);
            end
            #Forward loop
            for t in 2:(T+1)
                #println(t)
                τ_loc = τ_trans[country, t];
                openness = openness_transition[t];
                (price_final_prev,price_final_current,total_foreign_demand,W_loc,r_loc,avg_foreign_price,R,constant_d,
                constant_x,constant_z,constant_z_bar,constant_z_fine,constant_z_bar_fine,output_final_loc,
                residual,price_check_tmp) = price_reshaper_single(price_trans_actual[:,t],price_trans_actual[:,t-1],openness,country_no,
                    bounds,country,banking_cost_loc,δ,ω,σ,τ_loc,z_tilde,z_tilde_fine);
                (a_prime_fine[:,:,country,1,t-1],future_occupation_fine[:,:,country,1,t-1],cons_fine[:,:,country,1,t-1],rev_d_fine_store[:,country,t-1],
                    rev_dx_fine_store[:,country,t-1],rev_xx_fine_store[:,country,t-1],k_x_fine[:,country,t-1],k_d_fine[:,country,t-1],
                    l_x_fine[:,country,t-1],l_d_fine[:,country,t-1],Profit_fine_x[:,country,t-1],Profit_fine_d[:,country,t-1],output_d_fine_store[:,country,t-1],
                    output_dx_fine_store[:,country,t-1],output_xx_fine_store[:,country,t-1],
                    price_d_fine_store[:,country,t-1],price_dx_fine_store[:,country,t-1],price_xx_fine_store[:,country,t-1],distr_store_tmp[:,country,t],lambdda_d_fine_store[:,country,t-1],lambdda_x_fine_store[:,country,t-1]) = Residual_transition_forward(country,price_final_prev,price_final_current, R,W_loc,r_loc ,s_fine,
                    constant_x,constant_d,constant_z_fine,constant_z_bar_fine,β_loc,θ_loc,τ_loc,σ,α₁_eff,α₂_eff,α₁,α₂,
                    ones_tmp_fine,ns_fine,FM,FMₓ,V_tmp,
                    F, Fₓ,size_iid_cost_val,iid_cost_value,iid_cost_prob,Exit,
                    Exitₓ,a_min,a_max,α,fspace_a,
                    fspace_a_fine,openness,agrid_fine,coeff_store_tmp[:,:,country,t+1],
                    P_kron_fine,P_kron1,Phi,Phi_z_fine,income_mat_fine,x_tmp,
                     cons_fine_local,a_prime_fine_local,future_occupation_fine_local,Q_trans,distr_store_tmp[:,country,t-1]);
            end
            #Second backward loop for the aggregates - now more detailed
            for t in (T+1):-1:2
                τ_loc = τ_trans[country, t];
                openness = openness_transition[t];
                (price_final_prev,price_final_current,total_foreign_demand,W_loc,r_loc,avg_foreign_price,R,constant_d,
                constant_x,constant_z,constant_z_bar,constant_z_fine,constant_z_bar_fine,output_final_loc,
                residual,price_check_tmp) = price_reshaper_single(price_trans_actual[:,t],price_trans_actual[:,t-1],openness,country_no,
                    bounds,country,banking_cost_loc,δ,ω,σ,τ_loc,z_tilde,z_tilde_fine);
                    (labor_excess_demand_store_percent[country,t-1],excess_demand_final_good_store_percent[country,t-1],
                    total_demand_final_good[country,t-1],
                    export_price_sum_store[country,t-1], domestic_price_sum_store[country,t-1],NFA_store[country,t-1],
                    asset_supply_store[country,t-1],asset_demand_store[country,t-1],capital_supply_future[country,t],p90_wealth_store[country,t-1],
                    p90_cons_store[country,t-1],p90_income_store[country,t-1],
                    fraction_zombie_exporter_store[country,t-1],domestic_prod_store[country,t-1],export_prod_store[country,t-1],export_value_store[country,t-1],
                    sd_MRPK_d_store[country,t-1],sd_MRPK_x_store[country,t-1],sd_MRPK_store[country,t-1],
                    D_d_denom_store[country,t-1],D_d_store[country,t-1],D_x_denom_store[country,t-1],D_x_store[country,t-1],D_d_denom_eff_store[country,t-1],
                    D_d_eff_store[country,t-1],D_x_denom_eff_store[country,t-1],D_x_eff_store[country,t-1],domestic_pop_store[country,t-1],exporter_pop_store[country,t-1]
                    )= Residual_transition_backward_aggregates_detailed(country,output_final_loc,price_final_prev,price_final_current, R,W_loc,r_loc ,L_loc,s_fine,
                σ,ns_fine,FM,FMₓ,banking_cost_loc,δ,F, Fₓ,size_iid_cost_val,iid_cost_value,iid_cost_prob,Exit,
                Exitₓ,capital_supply_future[country,t+1],a_prime_fine[:,:,country,1,t-1],future_occupation_fine[:,:,country,1,t-1],cons_fine[:,:,country,1,t-1],
                rev_d_fine_store[:,country,t-1],rev_dx_fine_store[:,country,t-1],rev_xx_fine_store[:,country,t-1],k_x_fine[:,country,t-1],
                k_d_fine[:,country,t-1],l_x_fine[:,country,t-1],l_d_fine[:,country,t-1],
                Profit_fine_x[:,country,t-1],Profit_fine_d[:,country,t-1],output_d_fine_store[:,country,t-1],
                output_dx_fine_store[:,country,t-1],output_xx_fine_store[:,country,t-1],
                price_d_fine_store[:,country,t-1],price_dx_fine_store[:,country,t-1],price_xx_fine_store[:,country,t-1],
                distr_store_tmp[:,country,t-1],lambdda_d_fine_store[:,country,t-1],lambdda_x_fine_store[:,country,t-1],
                ones_tmp_fine,z_tilde_fine,agrid_fine,α₁_eff,α₂_eff,capital_supply_future[country,t])
            end
        end
        import_price = zeros(country_no,T);
        price_final_actual = zeros(country_no,T);
        residual_store = zeros(4*country_no,T);
        total_production_store = zeros(country_no,T);
        import_share_store = zeros(country_no,T);
        TOT_store = zeros(country_no,T);
        TFP_d_store = zeros(country_no,T);
        TFP_x_store = zeros(country_no,T);
        TFP_d_efficient_store = zeros(country_no,T);
        TFP_x_efficient_store = zeros(country_no,T);
        TFP_store = zeros(country_no,T);
        TFP_within_store = zeros(country_no,T);
        TFP_across_store= zeros(country_no,T);
        TFP_second_best_store= zeros(country_no,T);
        K_d_ratio_store = zeros(country_no,T);
        K_x_ratio_store = zeros(country_no,T);
        L_d_ratio_store = zeros(country_no,T);
        L_x_ratio_store = zeros(country_no,T);
        K_d_ratio_eff_store = zeros(country_no,T);
        K_x_ratio_eff_store = zeros(country_no,T);
        L_d_ratio_eff_store = zeros(country_no,T);
        L_x_ratio_eff_store = zeros(country_no,T);
        for country = 1:country_no
            (exitflag_tmp,s,ns,s_fine,ns_fine,Phi_z,Phi_z_fine,Phi,Phi_aug,
            P_kron,P_kron1,P_kron_fine,z_tilde,z_tilde_fine,income_mat,income_mat_fine,
            a_prime_fine_local,future_occupation_fine_local,cons_fine_local,coeff,
            coeff_next,θ_loc,L_loc,τ_loc,banking_cost_loc,
            ones_tmp,x_tmp,V_tmp,β_loc,V_next_stacked,iterate1,Phi_prime_tmp,
            D_deriv_tmp_block,conv,Q_trans,ones_tmp_fine) = country_local_initialize_noprice(country,s_cell,
            ns_cell,s_fine_cell,ns_fine_cell,Phi_z_cell,Phi_z_fine_cell,Phi_cell,Phi_aug_cell,
            P_kron_cell,P_kron1_cell,P_kron_fine_cell,σ,size_iid_cost_val,θ, L, τ,country_no,banking_cost,δ,ω,β);
            for t in 1:T
                τ_loc = τ_trans[country, t+1];
                openness = openness_transition[t+1];
                residual_store[country,t] = labor_excess_demand_store_percent[country,t];
                residual_store[country_no + country,t] = excess_demand_final_good_store_percent[country,t];
                import_price[country,t] = (sum(export_price_sum_store[:,t])-  export_price_sum_store[country,t])/(country_no - 1);
                price_final_actual[country,t] = min((ω^σ * domestic_price_sum_store[country,t] + (1.0 - ω)^σ  * import_price[country,t])^(1.0
                    / (1.0 - σ)),10000);
                import_curr = (sum(export_prod_store[:,t])-  export_prod_store[country,t])/(country_no - 1.0);
                total_production_store[country, t] =  (ω *domestic_prod_store[country,t] + (1.0 -ω) *import_curr)^((σ)/(σ-1));
                (price_final_prev,price_final_current,total_foreign_demand,W_loc,r_loc,avg_foreign_price,R,constant_d,
                constant_x,constant_z,constant_z_bar,constant_z_fine,constant_z_bar_fine,output_final_loc,
                residual,price_check_tmp) = price_reshaper_single(price_trans_actual[:,t+1],price_trans_actual[:,t],openness,country_no,
                    bounds,country,banking_cost_loc,δ,ω,σ,τ_loc,z_tilde,z_tilde_fine);

                residual_store[2* country_no + country,t] = (price_final_actual[country,t] - price_final_current)./(
                    price_final_actual[country,t] + price_final_current);
                # TFP calculations
                TFP_d_store[country, t] = ω * (D_d_store[country, t].^(1 - α₂_eff)./D_d_denom_store[country, t].^α₁_eff);
                TFP_x_store[country, t] = (D_x_store[country, t].^(1 - α₂_eff)./D_x_denom_store[country, t].^α₁_eff);
                TFP_d_efficient_store[country, t] = ω * (D_d_eff_store[country, t].^(1 - α₂_eff)./D_d_denom_eff_store[country, t].^α₁_eff);
                TFP_x_efficient_store[country, t] = (D_x_eff_store[country, t].^(1 - α₂_eff)./D_x_denom_eff_store[country, t].^α₁_eff);
                #total_foreign_demand = (sum(output_final) - output_final[country,1])/(country_no - 1);
                import_value_tmp = (sum(export_value_store[:,t])-  export_value_store[country,t]);
                export_value_tradeoff = import_value_tmp/export_value_store[country,t];
                import_share_store[country,t] = import_value_tmp/(country_no - 1)/price_final_actual[country,t] / total_production_store[country, t];
                TOT_store[country,t] = (ω *constant_d^(σ - 1) + (1 - ω) *constant_x^(σ - 1)*(
                 export_value_tradeoff *price_final_actual[country,t]/ avg_foreign_price)^((σ -1)/σ))/(
                constant_d^σ + (1 +τ_loc)* constant_x^σ)^((σ -1)/σ);
                TFP_x_store[country,t] = TFP_x_store[country,t] *TOT_store[country,t];
                TFP_x_efficient_store[country,t] = TFP_x_efficient_store[country,t] *TOT_store[country,t];
                D_d_K_tmp = D_d_denom_store[country,t];
                D_x_K_tmp = D_x_denom_store[country,t];
                D_d_K_eff_tmp = D_d_denom_eff_store[country,t];
                D_x_K_eff_tmp = D_x_denom_eff_store[country,t];
                K_d_ratio_store[country,t] = constant_d^σ * D_d_K_tmp/(constant_d^σ * D_d_K_tmp + (constant_d^σ
                 +(1 +τ_loc)*constant_x^σ)* D_x_K_tmp);
                K_x_ratio_store[country,t] = (constant_d^σ
                 +(1 +τ_loc)*constant_x^σ)* D_x_K_tmp/(constant_d^σ * D_d_K_tmp + (constant_d^σ
                  +(1 +τ_loc)*constant_x^σ)* D_x_K_tmp);
                L_d_ratio_store[country,t] = constant_d^σ * D_d_store[country,t]/(constant_d^σ * D_d_store[country,t] + (constant_d^σ
                    +(1 +τ_loc)*constant_x^σ)* D_x_store[country,t]);
                L_x_ratio_store[country,t] = (constant_d^σ
                    +(1 +τ_loc)*constant_x^σ)* D_x_store[country,t]/(constant_d^σ * D_d_store[country,t] + (constant_d^σ
                    +(1 +τ_loc)*constant_x^σ)* D_x_store[country,t]);
                K_d_ratio_eff_store[country,t] = constant_d^σ * D_d_K_eff_tmp/(constant_d^σ * D_d_K_eff_tmp + (constant_d^σ
                 +(1 +τ_loc)*constant_x^σ)* D_x_K_eff_tmp);
                K_x_ratio_eff_store[country,t] = (constant_d^σ
                 +(1 +τ_loc)*constant_x^σ)* D_x_K_eff_tmp/(constant_d^σ * D_d_K_eff_tmp + (constant_d^σ
                  +(1 +τ_loc)*constant_x^σ)* D_x_K_eff_tmp);
                L_d_ratio_eff_store[country,t] = constant_d^σ * D_d_eff_store[country,t]/(constant_d^σ * D_d_eff_store[country,t] + (constant_d^σ
                    +(1 +τ_loc)*constant_x^σ)* D_x_eff_store[country,t]);
                L_x_ratio_eff_store[country,t] = (constant_d^σ
                    +(1 +τ_loc)*constant_x^σ)* D_x_eff_store[country,t]/(constant_d^σ * D_d_eff_store[country,t] + (constant_d^σ
                    +(1 +τ_loc)*constant_x^σ)* D_x_eff_store[country,t]);
                TFP_store[country,t] = (TFP_d_store[country,t] * K_d_ratio_store[country,t]^α₁_eff * L_d_ratio_store[country,t]^α₂_eff
                    + TFP_x_store[country,t] * K_x_ratio_store[country,t]^α₁_eff * L_x_ratio_store[country,t]^α₂_eff)^(σ/(σ-1));
                TFP_within_store[country,t] = (TFP_d_efficient_store[country,t] * K_d_ratio_store[country,t]^α₁_eff * L_d_ratio_store[country,t]^α₂_eff
                    + TFP_x_efficient_store[country,t] * K_x_ratio_store[country,t]^α₁_eff * L_x_ratio_store[country,t]^α₂_eff)^(σ/(σ-1));
                TFP_across_store[country,t] = (TFP_d_store[country,t] * K_d_ratio_eff_store[country,t]^α₁_eff * L_d_ratio_eff_store[country,t]^α₂_eff
                    + TFP_x_store[country,t] * K_x_ratio_eff_store[country,t]^α₁_eff * L_x_ratio_eff_store[country,t]^α₂_eff)^(σ/(σ-1)); # this doesnt work
                TFP_second_best_store[country,t] =(TFP_d_efficient_store[country,t] * K_d_ratio_eff_store[country,t]^α₁_eff * L_d_ratio_eff_store[country,t]^α₂_eff
                    + TFP_x_efficient_store[country,t] * K_x_ratio_eff_store[country,t]^α₁_eff * L_x_ratio_eff_store[country,t]^α₂_eff)^(σ/(σ-1));
                if (openness == 0)
                    residual_store[3* country_no + country,t] = NFA_store[country,t]/(asset_demand_store[country,t] + (1 + banking_cost[country])*asset_supply_store[country,t]);
                end
            end
        end
        if openness == 1
            for t in 2:T
                residual_store[7,t] = sum(NFA_store[:,t])/sum((asset_demand_store[:,t]+ (ones(country_no) + banking_cost).*asset_supply_store[:,t]));
            end
        end
    end
    return (residual_store,price_final_actual,total_demand_final_good,total_production_store,import_share_store,TOT_store,TFP_d_store,
    TFP_x_store,TFP_d_efficient_store,TFP_x_efficient_store,TFP_store ,TFP_within_store,TFP_across_store,TFP_second_best_store,K_d_ratio_store,
    K_x_ratio_store,L_d_ratio_store,L_x_ratio_store,K_d_ratio_eff_store,K_x_ratio_eff_store,L_d_ratio_eff_store,L_x_ratio_eff_store,p90_wealth_store,
    p90_cons_store,p90_income_store,fraction_zombie_exporter_store,sd_MRPK_d_store,sd_MRPK_x_store,sd_MRPK_store,domestic_pop_store,
    exporter_pop_store,distr_store_tmp,coeff_store_tmp,capital_supply_future,price_trans_actual)
end
function Residual_transition_total_nonlin(vec_price_trans::Array{Float64,1},capital_trans::Array{Float64,2},price_trans::Array{Float64,2}
    ,distr_store::Array{Float64,3},T::Int64,parameter_end::Parameter_type,coeff_store::Array{Float64,4},τ_trans::Array{Float64,2},openness_transition::Array{Int64,1},country_no::Int64 = 2)
    price_trans_actual = copy(price_trans);
    price_trans_actual[:,2:(end-1)] = reshape(vec_price_trans,7,T);
    residual_store,price_final_actual,total_demand_final_good =  Residual_transition_sequential(price_trans_actual,capital_trans,
    price_trans,distr_store,T,parameter_end,coeff_store,τ_trans,openness_transition);
    if openness == 1
        vec_residual_store = reshape(residual_store,(3*country_no+1)*T,1);
    else
        vec_residual_store = reshape(residual_store,(4*country_no)*T,1);
    end
    return vec_residual_store
end
function Residual_transition_iterative(vec_price_trans::Array{Float64,1},capital_trans::Array{Float64,2},price_trans::Array{Float64,2}
    ,distr_store::Array{Float64,3},T::Int64,parameter_end::Parameter_type,coeff_store::Array{Float64,4},τ_trans::Array{Float64,2},
    openness_transition::Array{Int64,1},Jac_init::Array{Float64,2},Jac_fini::Array{Float64,2},country_no::Int64 = 2)
    price_trans_actual = copy(price_trans);
    price_trans_actual[:,2:(end-1)] = reshape(vec_price_trans,7,T);
    residual_store,price_final_actual,total_demand_final_good =  Residual_transition_sequential(price_trans_actual,capital_trans,
    price_trans,distr_store,T,parameter_end,coeff_store,τ_trans,openness_transition);
    maxval, maxval_ind = findmax(abs.(residual_store));
    println("Max residual: ",maxval," in time period: ",maxval_ind[2]," for market: ",maxval_ind[1])
    max_T = maxval_ind[2];
    price_trans_next = zeros(7,T);
    price_trans_next_slow = zeros(7,T);
    price_trans_next_maxT = copy(price_trans_actual[:,2:(end-1)]);
    price_trans_next_maxT_jac_based = copy(price_trans_actual[:,2:(end-1)]);
    price_trans_next_single_update_plus = copy(price_trans_actual[:,2:(end-1)]);
    price_trans_next_single_update_minus = copy(price_trans_actual[:,2:(end-1)]);
    res_lab = residual_store[1:2,:];
    res_tmp = res_lab.>0.2;
    res_lab[res_tmp] .= 0.2;
    res_tmp = res_lab.<-0.2;
    res_lab[res_tmp] .= -0.2;

    res_y = residual_store[3:4,:];
    res_tmp = res_y.>0.2;
    res_y[res_tmp] .= 0.2;
    res_tmp = res_y.<-0.2;
    res_y[res_tmp] .= -0.2;
    res_p = residual_store[5,:];
    res_tmp = res_p.>0.2;
    res_p[res_tmp] .= 0.2;
    res_tmp = res_p.<-0.2;
    res_p[res_tmp] .= -0.2;
    res_cap = residual_store[7:8,:];
    res_tmp = res_cap.>0.2;
    res_cap[res_tmp] .= 0.2;
    res_tmp = res_cap.<-0.2;
    res_cap[res_tmp] .= -0.2;
    within_country_time_cross_terms = zeros(4,4);
    within_country_time_cross_terms[1,1] =1;
    within_country_time_cross_terms[2,2] =1;
    within_country_time_cross_terms[3,3] =1;
    within_country_time_cross_terms[4,4] =-1;
    #within_country_time_cross_terms[1,:] = [1, 0.4 ,-0.2, 0];
    #within_country_time_cross_terms[2,:] = [0.4, 1, -0.2, 0];
    #within_country_time_cross_terms[3,:] = [-0.2, -0.2, 1, 0];
    #within_country_time_cross_terms[4,:] = [0, 0, 0,-1];
    res_p_aug = zeros(2,T);
    res_p_aug[1,:] = res_p;
    price_trans_next[1:2,:] = (ones(2,T)+within_country_time_cross_terms[1,1]*res_lab+within_country_time_cross_terms[1,2]*res_y +
    within_country_time_cross_terms[1,3]*res_p_aug + within_country_time_cross_terms[1,4]*res_cap ).*price_trans_actual[1:2,2:(end-1)];
    price_trans_next[3:4,:] =  (ones(2,T)+within_country_time_cross_terms[2,1]*res_lab+within_country_time_cross_terms[2,2]*res_y +
    within_country_time_cross_terms[2,3]*res_p_aug + within_country_time_cross_terms[2,4]*res_cap).*price_trans_actual[3:4,2:(end-1)];
    price_trans_next[5,:] =  (ones(T)+within_country_time_cross_terms[3,1]*res_lab[1,:]+within_country_time_cross_terms[3,2]*res_y[1,:] +
    within_country_time_cross_terms[3,3]*res_p + within_country_time_cross_terms[3,4]*res_cap[1,:]).*price_trans_actual[5,2:(end-1)];
    price_trans_next[6:7,:] = (ones(2,T)+within_country_time_cross_terms[4,1]*res_lab+within_country_time_cross_terms[4,2]*res_y +
    within_country_time_cross_terms[4,3]*res_p_aug + within_country_time_cross_terms[4,4]*res_cap).*price_trans_actual[6:7,2:(end-1)];
    residual_adj = zeros(7,T);
    residual_adj[1:2,:] = res_lab;
    residual_adj[3:4,:] = res_y;
    residual_adj[5,:] = res_p;
    residual_adj[6:7,:] = res_cap;
    res_tmp = (price_trans_actual[:,2:(end-1)].==0.0);
    residual_adj[res_tmp].=0;
    price_trans_next_maxT[:,1:max_T] = price_trans_next[:,1:max_T];
    jacob_improvement_init = Jac_init \ residual_adj;
    #time_weight_mat = repeat(1:T,1,7)'./T
    time_weight_mat = ones(7,T);
    if parameter_end.openness == 0
        jacob_improvement_fini = Jac_fini \ residual_adj;
        price_trans_next_jac_based =  price_trans_actual[:,2:(end-1)] - time_weight_mat.*jacob_improvement_fini - (1 .- time_weight_mat).*jacob_improvement_init;
    else
        jacob_improvement_fini = copy(jacob_improvement_init);
        jacob_improvement_fini[1:6,:] = Jac_fini \ residual_adj[1:6,:];
        price_trans_next_jac_based =  price_trans_actual[:,2:(end-1)] - time_weight_mat.*jacob_improvement_fini - (1 .- time_weight_mat).*jacob_improvement_init;
        price_trans_next_jac_based[7,2:end] .= 0;
    end
    price_trans_next_maxT_jac_based[:,1:max_T] = price_trans_next_jac_based[:,1:max_T];
    price_trans_next_single_update_plus[maxval_ind[1],maxval_ind[2]] = price_trans_next[maxval_ind[1],maxval_ind[2]];
    price_trans_next_single_update_minus[maxval_ind[1],maxval_ind[2]] = price_trans[maxval_ind[1],maxval_ind[2]+1].^2/price_trans_next[maxval_ind[1],maxval_ind[2]];
    vec_price_trans_next_maxT = reshape(price_trans_next_maxT,7*T);
    vec_price_trans_next = reshape(price_trans_next,7*T);
    vec_price_trans_next_jac_based = reshape(price_trans_next_jac_based,7*T);
    vec_price_trans_next_maxT_jac_based = reshape(price_trans_next_maxT_jac_based,7*T);
    vec_price_trans_next_single_update_plus = reshape(price_trans_next_single_update_plus,7*T);
    vec_price_trans_next_single_update_minus = reshape(price_trans_next_single_update_minus,7*T);
    return vec_price_trans_next,residual_store,vec_price_trans_next_maxT,vec_price_trans_next_jac_based,vec_price_trans_next_maxT_jac_based,vec_price_trans_next_single_update_plus,vec_price_trans_next_single_update_minus
end
#end
if load_solution == 0
    dampen_start = 0.90;
    dampen = copy(dampen_start);
    #iter_end = 10;
    #dampen_grid = range(dampen,0.995,length = iter_end);
    #for iterate1 = iter_end:-1:50
    iterate1 = 1;
    conv_res_past = 2000;
    conv_res = 1;
    vec_price_trans_smallest_res = copy(vec_price_trans);
    vec_price_trans_next_applied = copy(vec_price_trans);
    for i = 1:200
        (vec_price_trans_next,residual_store,vec_price_trans_next_maxT,vec_price_trans_next_jac_based,
        vec_price_trans_next_maxT_jac_based,vec_price_trans_next_single_update_plus,vec_price_trans_next_single_update_minus) = Residual_transition_iterative(vec_price_trans,capital_trans,price_trans,
        distr_store,T,parameter_end,coeff_store,τ_trans,openness_transition,Jac_init,Jac_fini);
        #vec_price_trans_next,residual_store,vec_price_trans_next_maxT = Residual_transition_iterative(vec_price_trans,capital_trans,price_trans,
        #distr_store,T,parameter_end,coeff_store,τ_trans,openness_transition,Jac_init,Jac_fini);
        vec_price_trans_next,residual_store,vec_price_trans_next_maxT = Residual_transition_iterative(vec_price_trans_smallest_res,capital_trans,price_trans,
        distr_store,T,parameter_end,coeff_store,τ_trans,openness_transition,Jac_init,Jac_fini);
        conv = maximum(abs.(vec_price_trans_next - vec_price_trans));
        conv_res = maximum(abs.(residual_store));
        if (conv_res<conv_res_past || dampen>0.995)
            if conv_res<conv_res_past
                println("Smaller residual found")
                conv_res_past = copy(conv_res)
                dampen = copy(dampen_start);
                vec_price_trans_smallest_res = copy(vec_price_trans);
                vec_price_trans_next_applied = copy(vec_price_trans_next);
                #vec_price_trans_next_applied = copy(vec_price_trans_next_maxT);
                #vec_price_trans_next_applied = copy(vec_price_trans_next_jac_based);
                #vec_price_trans_next_applied = copy(vec_price_trans_next_single_update_plus);
            elseif dampen>0.995
                println("Dampening too low")
                vec_price_trans_next_single_update_plus_dampened =  (1 - dampen)*copy(vec_price_trans_next_single_update_plus)+dampen*vec_price_trans_smallest_res;
                vec_price_trans_next_single_update_minus_dampened =  (1 - dampen)*copy(vec_price_trans_next_single_update_minus)+dampen*vec_price_trans_smallest_res;
                (vec_price_trans_next1,residual_store1,vec_price_trans_next_maxT1,vec_price_trans_next_jac_based,vec_price_trans_next_maxT_jac_based,
                vec_price_trans_next_single_update_plus1,vec_price_trans_next_single_update_minus1)= Residual_transition_iterative(vec_price_trans_next_single_update_plus_dampened,
                capital_trans,price_trans,distr_store,T,parameter_end,coeff_store,τ_trans,openness_transition,Jac_init,Jac_fini);
                conv_res_plus = maximum(abs.(residual_store1));
                (vec_price_trans_next2,residual_store2,vec_price_trans_next_maxT2,vec_price_trans_next_jac_based,vec_price_trans_next_maxT_jac_based,
                vec_price_trans_next_single_update_plus2,vec_price_trans_next_single_update_minus2) = Residual_transition_iterative(vec_price_trans_next_single_update_minus_dampened,
                capital_trans,price_trans, distr_store,T,parameter_end,coeff_store,τ_trans,openness_transition,Jac_init,Jac_fini);
                conv_res_minus = maximum(abs.(residual_store2));
                if (conv_res_plus< conv_res_minus && conv_res_plus <conv_res)
                    println("Positive price updating")
                    conv_res_past = copy(conv_res_plus)
                    dampen = copy(dampen_start);
                    vec_price_trans_next_applied = copy(vec_price_trans_next_single_update_plus1);
                    vec_price_trans_smallest_res = copy(vec_price_trans_next_single_update_plus_dampened);
                elseif (conv_res_minus< conv_res_plus && conv_res_minus <conv_res)
                    println("Negative price updating")
                    conv_res_past = copy(conv_res_minus)
                    dampen = copy(dampen_start);
                    vec_price_trans_next_applied = copy(vec_price_trans_next_single_update_plus2);
                    vec_price_trans_smallest_res = copy(vec_price_trans_next_single_update_minus_dampened);
                else
                    println("Move to a new guess")
                    #conv_res_past = copy(conv_res)
                    dampen = copy(dampen_start);
                    vec_price_trans_next_applied = copy(vec_price_trans_next);
                    #vec_price_trans_next_applied = copy(vec_price_trans_next_single_update_plus_dampened);
                end
            end
        else
            dampen = (1.0 + 5 * dampen)/6
            println("Current dampening: ",dampen)
        end

        #vec_price_trans = (1 - dampen_grid[iterate1])*copy(vec_price_trans_next)+dampen_grid[iterate1]*vec_price_trans;
        vec_price_trans = (1 - dampen)*copy(vec_price_trans_next_applied)+dampen*vec_price_trans_smallest_res;
        #vec_price_trans = (1 - dampen)*copy(vec_price_trans_next)+dampen*vec_price_trans_smallest_res;
        println("Iteration:",iterate1," Prices:",conv," Residuals:",conv_res)
        iterate1 = iterate1 + 1;
    end
elseif load_solution==1
    price_trans_actual = copy(price_trans);
    price_trans_actual[:,2:(end-1)] = reshape(vec_price_trans,7,T);
    if case_final == 1
        #open trade closed capital markets:
        (residual_store_open_CM_closed_trade,price_final_actual_open_CM_closed_trade,total_demand_final_good_open_CM_closed_trade,total_production_store_open_CM_closed_trade,import_share_store_open_CM_closed_trade,TOT_store_open_CM_closed_trade,TFP_d_store_open_CM_closed_trade,
        TFP_x_store_open_CM_closed_trade,TFP_d_efficient_store_open_CM_closed_trade,TFP_x_efficient_store_open_CM_closed_trade,TFP_store_open_CM_closed_trade ,TFP_within_store_open_CM_closed_trade,TFP_across_store_open_CM_closed_trade,TFP_second_best_store_open_CM_closed_trade,K_d_ratio_store_open_CM_closed_trade,
        K_x_ratio_store_open_CM_closed_trade,L_d_ratio_store_open_CM_closed_trade,L_x_ratio_store_open_CM_closed_trade,K_d_ratio_eff_store_open_CM_closed_trade,K_x_ratio_eff_store_open_CM_closed_trade,L_d_ratio_eff_store_open_CM_closed_trade,L_x_ratio_eff_store_open_CM_closed_trade,p90_wealth_store_open_CM_closed_trade,
        p90_cons_store_open_CM_closed_trade,p90_income_store_open_CM_closed_trade,fraction_zombie_exporter_store_open_CM_closed_trade,sd_MRPK_d_store_open_CM_closed_trade,sd_MRPK_x_store_open_CM_closed_trade,sd_MRPK_store_open_CM_closed_trade,domestic_pop_store_open_CM_closed_trade,
        exporter_pop_store_open_CM_closed_trade,distr_store_open_CM_closed_trade,coeff_store_open_CM_closed_trade,capital_supply_future_open_CM_closed_trade,
        price_trans_actual_open_CM_closed_trade) =  Residual_transition_sequential_detailed(price_trans_actual,capital_trans,
        price_trans,distr_store,T,parameter_end,coeff_store,τ_trans,openness_transition);
    elseif case_final == 2
        #open trade closed capital markets:
        (residual_store_closed_CM_open_trade,price_final_actual_closed_CM_open_trade,total_demand_final_good_closed_CM_open_trade,total_production_store_closed_CM_open_trade,import_share_store_closed_CM_open_trade,TOT_store_closed_CM_open_trade,TFP_d_store_closed_CM_open_trade,
        TFP_x_store_closed_CM_open_trade,TFP_d_efficient_store_closed_CM_open_trade,TFP_x_efficient_store_closed_CM_open_trade,TFP_store_closed_CM_open_trade ,TFP_within_store_closed_CM_open_trade,TFP_across_store_closed_CM_open_trade,TFP_second_best_store_closed_CM_open_trade,K_d_ratio_store_closed_CM_open_trade,
        K_x_ratio_store_closed_CM_open_trade,L_d_ratio_store_closed_CM_open_trade,L_x_ratio_store_closed_CM_open_trade,K_d_ratio_eff_store_closed_CM_open_trade,K_x_ratio_eff_store_closed_CM_open_trade,L_d_ratio_eff_store_closed_CM_open_trade,L_x_ratio_eff_store_closed_CM_open_trade,p90_wealth_store_closed_CM_open_trade,
        p90_cons_store_closed_CM_open_trade,p90_income_store_closed_CM_open_trade,fraction_zombie_exporter_store_closed_CM_open_trade,sd_MRPK_d_store_closed_CM_open_trade,sd_MRPK_x_store_closed_CM_open_trade,sd_MRPK_store_closed_CM_open_trade,domestic_pop_store_closed_CM_open_trade,
        exporter_pop_store_closed_CM_open_trade,distr_store_closed_CM_open_trade,coeff_store_closed_CM_open_trade,capital_supply_future_closed_CM_open_trade,
        price_trans_actual_closed_CM_open_trade) =  Residual_transition_sequential_detailed(price_trans_actual,capital_trans,
        price_trans,distr_store,T,parameter_end,coeff_store,τ_trans,openness_transition);
    elseif case_final == 3
        #open trade open capital markets:
        (residual_store_open_CM_open_trade,price_final_actual_open_CM_open_trade,total_demand_final_good_open_CM_open_trade,total_production_store_open_CM_open_trade,import_share_store_open_CM_open_trade,TOT_store_open_CM_open_trade,TFP_d_store_open_CM_open_trade,
        TFP_x_store_open_CM_open_trade,TFP_d_efficient_store_open_CM_open_trade,TFP_x_efficient_store_open_CM_open_trade,TFP_store_open_CM_open_trade ,TFP_within_store_open_CM_open_trade,TFP_across_store_open_CM_open_trade,TFP_second_best_store_open_CM_open_trade,K_d_ratio_store_open_CM_open_trade,
        K_x_ratio_store_open_CM_open_trade,L_d_ratio_store_open_CM_open_trade,L_x_ratio_store_open_CM_open_trade,K_d_ratio_eff_store_open_CM_open_trade,K_x_ratio_eff_store_open_CM_open_trade,L_d_ratio_eff_store_open_CM_open_trade,L_x_ratio_eff_store_open_CM_open_trade,p90_wealth_store_open_CM_open_trade,
        p90_cons_store_open_CM_open_trade,p90_income_store_open_CM_open_trade,fraction_zombie_exporter_store_open_CM_open_trade,sd_MRPK_d_store_open_CM_open_trade,sd_MRPK_x_store_open_CM_open_trade,sd_MRPK_store_open_CM_open_trade,domestic_pop_store_open_CM_open_trade,
        exporter_pop_store_open_CM_open_trade,distr_store_open_CM_open_trade,coeff_store_open_CM_open_trade,capital_supply_future_open_CM_open_trade,
        price_trans_actual_open_CM_open_trade) =  Residual_transition_sequential_detailed(price_trans_actual,capital_trans,
        price_trans,distr_store,T,parameter_end,coeff_store,τ_trans,openness_transition);
    elseif case_final == 4
        (residual_store_open_CMdelayed_open_trade,price_final_actual_open_CMdelayed_open_trade,total_demand_final_good_open_CMdelayed_open_trade,total_production_store_open_CMdelayed_open_trade,import_share_store_open_CMdelayed_open_trade,TOT_store_open_CMdelayed_open_trade,TFP_d_store_open_CMdelayed_open_trade,
        TFP_x_store_open_CMdelayed_open_trade,TFP_d_efficient_store_open_CMdelayed_open_trade,TFP_x_efficient_store_open_CMdelayed_open_trade,TFP_store_open_CMdelayed_open_trade ,TFP_within_store_open_CMdelayed_open_trade,TFP_across_store_open_CMdelayed_open_trade,TFP_second_best_store_open_CMdelayed_open_trade,K_d_ratio_store_open_CMdelayed_open_trade,
        K_x_ratio_store_open_CMdelayed_open_trade,L_d_ratio_store_open_CMdelayed_open_trade,L_x_ratio_store_open_CMdelayed_open_trade,K_d_ratio_eff_store_open_CMdelayed_open_trade,K_x_ratio_eff_store_open_CMdelayed_open_trade,L_d_ratio_eff_store_open_CMdelayed_open_trade,L_x_ratio_eff_store_open_CMdelayed_open_trade,p90_wealth_store_open_CMdelayed_open_trade,
        p90_cons_store_open_CMdelayed_open_trade,p90_income_store_open_CMdelayed_open_trade,fraction_zombie_exporter_store_open_CMdelayed_open_trade,sd_MRPK_d_store_open_CMdelayed_open_trade,sd_MRPK_x_store_open_CMdelayed_open_trade,sd_MRPK_store_open_CMdelayed_open_trade,domestic_pop_store_open_CMdelayed_open_trade,
        exporter_pop_store_open_CMdelayed_open_trade,distr_store_open_CMdelayed_open_trade,coeff_store_open_CMdelayed_open_trade,capital_supply_future_open_CMdelayed_open_trade,
        price_trans_actual_open_CMdelayed_open_trade) =  Residual_transition_sequential_detailed(price_trans_actual,capital_trans,
        price_trans,distr_store,T,parameter_end,coeff_store,τ_trans,openness_transition);
    elseif case_final == 5
        #open trade open capital markets:
        (residual_store_open_CM_closed_trade_dev,price_final_actual,total_demand_final_good,total_production_store_open_CM_closed_trade_dev,import_share_store_open_CM_closed_trade_dev,TOT_store_open_CM_closed_trade_dev,TFP_d_store_open_CM_closed_trade_dev,
        TFP_x_store_open_CM_closed_trade_dev,TFP_d_efficient_store_open_CM_closed_trade_dev,TFP_x_efficient_store_open_CM_closed_trade_dev,TFP_store_open_CM_closed_trade_dev ,TFP_within_store_open_CM_closed_trade_dev,TFP_across_store_open_CM_closed_trade_dev,TFP_second_best_store_open_CM_closed_trade_dev,K_d_ratio_store_open_CM_closed_trade_dev,
        K_x_ratio_store_open_CM_closed_trade_dev,L_d_ratio_store_open_CM_closed_trade_dev,L_x_ratio_store_open_CM_closed_trade_dev,K_d_ratio_eff_store_open_CM_closed_trade_dev,K_x_ratio_eff_store_open_CM_closed_trade_dev,L_d_ratio_eff_store_open_CM_closed_trade_dev,L_x_ratio_eff_store_open_CM_closed_trade_dev,p90_wealth_store_open_CM_closed_trade_dev,
        p90_cons_store_open_CM_closed_trade_dev,p90_income_store_open_CM_closed_trade_dev,fraction_zombie_exporter_store_open_CM_closed_trade_dev,sd_MRPK_d_store_open_CM_closed_trade_dev,sd_MRPK_x_store_open_CM_closed_trade_dev,sd_MRPK_store_open_CM_closed_trade_dev,domestic_pop_store_open_CM_closed_trade_dev,
        exporter_pop_store_open_CM_closed_trade_dev,distr_store_open_CM_closed_trade_dev,coeff_store_open_CM_closed_trade_dev) =  Residual_transition_sequential_detailed(price_trans_actual,capital_trans,
        price_trans,distr_store,T,parameter_end,coeff_store,τ_trans,openness_transition);
    elseif case_final == 6
        #open trade closed capital markets:
        (residual_store_closed_CM_open_trade_dev,price_final_actual,total_demand_final_good,total_production_store_closed_CM_open_trade_dev,import_share_store_closed_CM_open_trade_dev,TOT_store_closed_CM_open_trade_dev,TFP_d_store_closed_CM_open_trade_dev,
        TFP_x_store_closed_CM_open_trade_dev,TFP_d_efficient_store_closed_CM_open_trade_dev,TFP_x_efficient_store_closed_CM_open_trade_dev,TFP_store_closed_CM_open_trade_dev ,TFP_within_store_closed_CM_open_trade_dev,TFP_across_store_closed_CM_open_trade_dev,TFP_second_best_store_closed_CM_open_trade_dev,K_d_ratio_store_closed_CM_open_trade_dev,
        K_x_ratio_store_closed_CM_open_trade_dev,L_d_ratio_store_closed_CM_open_trade_dev,L_x_ratio_store_closed_CM_open_trade_dev,K_d_ratio_eff_store_closed_CM_open_trade_dev,K_x_ratio_eff_store_closed_CM_open_trade_dev,L_d_ratio_eff_store_closed_CM_open_trade_dev,L_x_ratio_eff_store_closed_CM_open_trade_dev,p90_wealth_store_closed_CM_open_trade_dev,
        p90_cons_store_closed_CM_open_trade_dev,p90_income_store_closed_CM_open_trade_dev,fraction_zombie_exporter_store_closed_CM_open_trade_dev,sd_MRPK_d_store_closed_CM_open_trade_dev,sd_MRPK_x_store_closed_CM_open_trade_dev,sd_MRPK_store_closed_CM_open_trade_dev,domestic_pop_store_closed_CM_open_trade_dev,
        exporter_pop_store_closed_CM_open_trade_dev,distr_store_closed_CM_open_trade_dev,coeff_store_closed_CM_open_trade_dev) =  Residual_transition_sequential_detailed(price_trans_actual,capital_trans,
        price_trans,distr_store,T,parameter_end,coeff_store,τ_trans,openness_transition);
    elseif case_final == 7
        #open trade open capital markets:
        (residual_store_open_CM_open_trade_dev,price_final_actual,total_demand_final_good,total_production_store_open_CM_open_trade_dev,import_share_store_open_CM_open_trade_dev,TOT_store_open_CM_open_trade_dev,TFP_d_store_open_CM_open_trade_dev,
        TFP_x_store_open_CM_open_trade_dev,TFP_d_efficient_store_open_CM_open_trade_dev,TFP_x_efficient_store_open_CM_open_trade_dev,TFP_store_open_CM_open_trade_dev ,TFP_within_store_open_CM_open_trade_dev,TFP_across_store_open_CM_open_trade_dev,TFP_second_best_store_open_CM_open_trade_dev,K_d_ratio_store_open_CM_open_trade_dev,
        K_x_ratio_store_open_CM_open_trade_dev,L_d_ratio_store_open_CM_open_trade_dev,L_x_ratio_store_open_CM_open_trade_dev,K_d_ratio_eff_store_open_CM_open_trade_dev,K_x_ratio_eff_store_open_CM_open_trade_dev,L_d_ratio_eff_store_open_CM_open_trade_dev,L_x_ratio_eff_store_open_CM_open_trade_dev,p90_wealth_store_open_CM_open_trade_dev,
        p90_cons_store_open_CM_open_trade_dev,p90_income_store_open_CM_open_trade_dev,fraction_zombie_exporter_store_open_CM_open_trade_dev,sd_MRPK_d_store_open_CM_open_trade_dev,sd_MRPK_x_store_open_CM_open_trade_dev,sd_MRPK_store_open_CM_open_trade_dev,domestic_pop_store_open_CM_open_trade_dev,
        exporter_pop_store_open_CM_open_trade_dev,distr_store_open_CM_open_trade_dev,coeff_store_open_CM_open_trade_dev) =  Residual_transition_sequential_detailed(price_trans_actual,capital_trans,
        price_trans,distr_store,T,parameter_end,coeff_store,τ_trans,openness_transition);
    elseif case_final == 8
        (residual_store_open_CMdelayed_open_trade_dev,price_final_actual,total_demand_final_good,total_production_store_open_CMdelayed_open_trade_dev,import_share_store_open_CMdelayed_open_trade_dev,TOT_store_open_CMdelayed_open_trade_dev,TFP_d_store_open_CMdelayed_open_trade_dev,
        TFP_x_store_open_CMdelayed_open_trade_dev,TFP_d_efficient_store_open_CMdelayed_open_trade_dev,TFP_x_efficient_store_open_CMdelayed_open_trade_dev,TFP_store_open_CMdelayed_open_trade_dev ,TFP_within_store_open_CMdelayed_open_trade_dev,TFP_across_store_open_CMdelayed_open_trade_dev,TFP_second_best_store_open_CMdelayed_open_trade_dev,K_d_ratio_store_open_CMdelayed_open_trade_dev,
        K_x_ratio_store_open_CMdelayed_open_trade_dev,L_d_ratio_store_open_CMdelayed_open_trade_dev,L_x_ratio_store_open_CMdelayed_open_trade_dev,K_d_ratio_eff_store_open_CMdelayed_open_trade_dev,K_x_ratio_eff_store_open_CMdelayed_open_trade_dev,L_d_ratio_eff_store_open_CMdelayed_open_trade_dev,L_x_ratio_eff_store_open_CMdelayed_open_trade_dev,p90_wealth_store_open_CMdelayed_open_trade_dev,
        p90_cons_store_open_CMdelayed_open_trade_dev,p90_income_store_open_CMdelayed_open_trade_dev,fraction_zombie_exporter_store_open_CMdelayed_open_trade_dev,sd_MRPK_d_store_open_CMdelayed_open_trade_dev,sd_MRPK_x_store_open_CMdelayed_open_trade_dev,sd_MRPK_store_open_CMdelayed_open_trade_dev,domestic_pop_store_open_CMdelayed_open_trade_dev,
        exporter_pop_store_open_CMdelayed_open_trade_dev,distr_store_open_CMdelayed_open_trade_dev,coeff_store_open_CMdelayed_open_trade_dev) =  Residual_transition_sequential_detailed(price_trans_actual,capital_trans,
        price_trans,distr_store,T,parameter_end,coeff_store,τ_trans,openness_transition);
    end
elseif load_solution==2
    @everywhere function f_true(vec_price_trans)
        vec_residual_store = Residual_transition_total_nonlin(vec_price_trans,capital_trans,price_trans
        ,distr_store,T,parameter_end,coeff_store,τ_trans,openness_transition);
        res = sum(vec_residual_store.^2);
        println(res, " ")
        return vec_residual_store
    end
    @everywhere function f_sum(vec_price_trans)
        vec_residual_store = Residual_transition_total_nonlin(vec_price_trans,capital_trans,price_trans
        ,distr_store,T,parameter_end,coeff_store,τ_trans,openness_transition);
        res = sum(vec_residual_store.^2);
        println(res, " ")
        return res
    end
#    Trans_simplex = Optim.AffineSimplexer(0.025,0.005)
#    optim_res = optimize(f_sum, vec_price_trans, method =NelderMead(initial_simplex = Trans_simplex),store_trace = true,show_trace = true,
#    extended_trace=true,time_limit = 10000.0);
#    prices = optim_res.minimizer
#    ls_res= LeastSquaresOptim.optimize(f_true, vec_price_trans, LevenbergMarquardt(),show_trace = true, store_trace = true,
#    x_tol = 1e-10, f_tol=  1e-10,iterations=40,time_limit = 10000.0);
    @everywhere function gradient1(f::Function,x::Array{Float64,1})
        #Forward differencing
        delta=1e-7*Matrix(I,length(x),length(x));
        J=convert(SharedArray,zeros(length(x)));
        f_x = f(x);
        @sync @distributed for i=1:length(x)
            tmp =(f(x+delta[:,i])-f_x)/1e-7;
            J[i] = tmp;
        end
        return J
    end
    @everywhere function jacob_f_sum(x)
        G = gradient1(f_sum,x);
        return G
    end
    optim_res1 = optimize(f_sum,jacob_f_sum,vec_price_trans, BFGS(alphaguess = LineSearches.InitialStatic(alpha=1e-6,scaled=true)); inplace = false,time_limit = 20000.0)

end

if save_solution==1
    # This assumes that the main script has been evaluated
    if case_final == 1
        open("vec_price_openCM_notrade.csv", "a") do io
            writedlm(io, vec_price_trans_smallest_res,',')
        end
    elseif case_final == 2
        open("vec_price_clCM.csv", "a") do io
            writedlm(io, vec_price_trans_smallest_res,',')
        end
    elseif case_final == 3
        open("vec_price_openCM.csv", "a") do io
            writedlm(io, vec_price_trans_smallest_res,',')
        end
    elseif case_final == 4
        open("vec_price_openCM_delayed.csv", "a") do io
            writedlm(io, vec_price_trans_smallest_res,',')
        end
    elseif case_final == 5
        open("vec_price_openCM_notrade_dev.csv", "a") do io
            writedlm(io, vec_price_trans_smallest_res,',')
        end
    elseif case_final == 6
        open("vec_price_clCM_dev.csv", "a") do io
            writedlm(io, vec_price_trans_smallest_res,',')
        end
    elseif case_final == 7
        open("vec_price_openCM_dev.csv", "a") do io
            writedlm(io, vec_price_trans_smallest_res,',')
        end
    elseif case_final == 8
        open("vec_price_openCM_delayed_dev.csv", "a") do io
            writedlm(io, vec_price_trans_smallest_res,',')
        end
    end
end
