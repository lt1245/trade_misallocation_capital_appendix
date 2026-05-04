include("code_new.jl")
prices_initial = [  0.24190237714043653
0.38077781106180153
0.2248600604939829
2.274816020027389
1.2945582363748034
0.05653056683161781
0.049450671591694605];

τ_closed_trade = 0.525; #Final bilateral low tariff state
Baseline_parameter.τ = τ_closed_trade * [1.0,1.0];
parameter_tmp = copy(Baseline_parameter);
parameter_tmp.σₛ = 1.00 *parameter_tmp.σₛ;
parameter_tmp.β[1] = 0.948;
parameter_tmp.θ[1] = 0.86;

prices = copy(prices_initial);
residual = Residual_stst(prices, parameter_tmp)
function f_true(prices)
    residual = Residual_stst(prices, parameter_tmp);
    return residual
end
function f_sum(prices)
    residual = Residual_stst(prices, parameter_tmp);
    return sum(residual.^2)
end

lower_guess_bound = 0.1;
upper_guess_bound = 2.0;

ls_res= LeastSquaresOptim.optimize(f_true, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
upper = upper_guess_bound * prices);
prices = ls_res.minimizer;
#        prices  = optim_res1.minimizer;
#        stst_simplex = Optim.AffineSimplexer(0.025,optim_res1.ssr);
stst_simplex = Optim.AffineSimplexer(0.025,ls_res.ssr);
#        stst_simplex = Optim.AffineSimplexer(0.025,0.005);
optim_res = optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
extended_trace=true);
prices = optim_res.minimizer;
ls_res= LeastSquaresOptim.optimize(f_true, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
upper = upper_guess_bound * prices);
prices = ls_res.minimizer;
stst_simplex = Optim.AffineSimplexer(0.025,ls_res.ssr);
#        stst_simplex = Optim.AffineSimplexer(0.025,0.005);
optim_res = optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
extended_trace=true,time_limit = 3000.0);
prices = optim_res.minimizer;


#Closed capital markets and after trade liberalization
closed_CM_open_trade_parameter = copy(parameter_tmp);
τ_open_trade = 0.354; #Final bilateral low tariff state
closed_CM_open_trade_parameter.τ = τ_open_trade * [1.0,1.0];
function f_true(prices)
    residual = Residual_stst(prices, closed_CM_open_trade_parameter);
    return residual
end
function f_sum(prices)
    residual = Residual_stst(prices, closed_CM_open_trade_parameter);
    return sum(residual.^2)
end
prices_closed_CM_open_trade = [ 0.24892225689053699
0.37975992418504045
0.23339271207738851
2.1597019583593378
1.2604151452576928
0.056426060722649676
0.05126321391584822];
lower_guess_bound = 0.1;
upper_guess_bound = 2.0;
prices = copy(prices_closed_CM_open_trade);
ls_res= LeastSquaresOptim.optimize(f_true, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
upper = upper_guess_bound * prices);
prices = ls_res.minimizer;
#        prices  = optim_res1.minimizer;
#        stst_simplex = Optim.AffineSimplexer(0.025,optim_res1.ssr);
stst_simplex = Optim.AffineSimplexer(0.025,ls_res.ssr);
#        stst_simplex = Optim.AffineSimplexer(0.025,0.005);
optim_res = optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
extended_trace=true);
prices = optim_res.minimizer;

#Open capital markets and after trade liberalization

open_CM_open_trade_parameter = copy(parameter_tmp);
open_CM_open_trade_parameter.openness = 1;

τ_open_trade = 0.354; #Final bilateral low tariff state
open_CM_open_trade_parameter.τ = τ_open_trade * [1.0,1.0];

prices_open_CM_open_trade =  [  0.2361233906432929
0.37955347959555313
0.22729335230233944
2.169395820881968
1.2286854773034048
0.05065824098026955];
function f_true(prices)
    residual = Residual_stst(prices, open_CM_open_trade_parameter);
    return residual
end
function f_sum(prices)
    residual = Residual_stst(prices, open_CM_open_trade_parameter);
    return sum(residual.^2)
end

lower_guess_bound = 0.1;
upper_guess_bound = 2.0;
prices = copy(prices_open_CM_open_trade);
residual = Residual_stst(prices, open_CM_open_trade_parameter)
ls_res= LeastSquaresOptim.optimize(f_true, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
upper = upper_guess_bound * prices);
prices = ls_res.minimizer;
#        prices  = optim_res1.minimizer;
#        stst_simplex = Optim.AffineSimplexer(0.025,optim_res1.ssr);
stst_simplex = Optim.AffineSimplexer(0.025,ls_res.ssr);
#        stst_simplex = Optim.AffineSimplexer(0.025,0.005);
optim_res = optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
extended_trace=true);
prices = optim_res.minimizer;
#Open capital markets only
#parameter_tmp = copy(Baseline_parameter);
open_CM_closed_trade_parameter = copy(parameter_tmp);
open_CM_closed_trade_parameter.openness = 1;

τ_closed_trade = 0.525; #Final bilateral low tariff state
open_CM_closed_trade_parameter.τ = τ_closed_trade * [1.0,1.0];
prices_open_CM_closed_trade =  [    0.23682352802411427,
 0.38112005081922895,
 0.23684346058879885,
 2.2709889030267885,
 1.2852541180790449,
 0.05012701009487472];
function f_true(prices)
    residual = Residual_stst(prices, open_CM_closed_trade_parameter);
    return residual
end
function f_sum(prices)
    residual = Residual_stst(prices, open_CM_closed_trade_parameter);
    return sum(residual.^2)
end

lower_guess_bound = 0.1;
upper_guess_bound = 2.0;
prices = copy(prices_open_CM_closed_trade);
residual = Residual_stst(prices, open_CM_closed_trade_parameter)
ls_res= LeastSquaresOptim.optimize(f_true, prices, LevenbergMarquardt(),show_trace = true, store_trace = true,
x_tol = 1e-9, f_tol= 1e-5,iterations=20,lower = lower_guess_bound * prices,
upper = upper_guess_bound * prices);
prices = ls_res.minimizer;
#        prices  = optim_res1.minimizer;
#        stst_simplex = Optim.AffineSimplexer(0.025,optim_res1.ssr);
stst_simplex = Optim.AffineSimplexer(0.025,ls_res.ssr);
#        stst_simplex = Optim.AffineSimplexer(0.025,0.005);
optim_res = optimize(f_sum, prices, method =NelderMead(initial_simplex = stst_simplex),store_trace = true,show_trace = true,
extended_trace=true);
prices = optim_res.minimizer;
