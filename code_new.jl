using DelimitedFiles
using Optim
using DataFrames
using CSV

using SharedArrays, Distributed
using CompEcon, QuantEcon
using LinearAlgebra, Statistics
using SparseArrays, StructArrays
using BasisMatrices, Arpack
using NLsolve, BlackBoxOptim, LeastSquaresOptim
using BenchmarkTools #,FiniteDiff
addprocs(2,exeflags="--project");
#addprocs(1);
@everywhere begin
using CompEcon,QuantEcon
using SparseArrays, StructArrays,SharedArrays
using LinearAlgebra, Statistics
using BasisMatrices, Arpack,BlackBoxOptim#, FiniteDiff
using DelimitedFiles
using Optim
using DataFrames
using CSV
end

@everywhere begin
mutable struct Parameter_type
    country_no::Int64 #number of countries
    α::Float64 # capital share
    η::Float64 # decreasing returns
    δ::Float64 # Depreciation rate
    σ::Float64 # Elasticity of substitution
    α₁::Float64 # shortcut for the C-D functions, effective capital share
    α₂::Float64 # shortcut for the C-D functions, effective labor share
    α₁_eff::Float64 # shortcut for the C-D functions, markup adjusted capital share
    α₂_eff::Float64 # shortcut for the C-D functions, markup adjusted labor share
    L::Array{Float64,1} # size of the economies
    openness::Int # 0: closed capital markets, 1: open capital markets
    β::Array{Float64,1} # Discount factor for each countries
    ω::Float64 # Home bias in final goods
    θ::Array{Float64,1} # financial development - max leverage constraint
    banking_cost::Array{Float64,1} # banking friction, wedge between the savings and lending rate -unused
    Avg_prod::Array{Float64,1} # Productivity shifter
    ρ::Array{Float64,1} # autoregressive part of the transitory productivity shock
    σₛ::Array{Float64,1} # Standard deviation of the transitory productivity shock
    κ::Float64 # Scaled entry cost
    FM::Array{Float64,1} #  Maintenance cost of the domestic sector
    FMₓ::Array{Float64,1} #  Maintenance cost of the export sector
    F::Array{Float64,1} #  Entry cost to the domestic sector
    Fₓ::Array{Float64,1} #  Entry cost to the export sector
    Exit::Array{Float64,1} #  Exit cost from the domestic sector  -unused
    Exitₓ::Array{Float64,1} #  Exit cost from the export sector  -unused
    Improv::Float64 # Productivity improvement in the domestic sector -unused
    Improvₓ::Float64 # Productivity improvement in the domestic sector -unused
    iid_cost_value::Array{Float64,1} #  grid for the iid entry cost into the production sector
    iid_cost_prob::Array{Float64,1} #  probabilities for the iid entry cost into the production sector
    n::Array{Int64,1} # joint approximation number of gridpoints for 1) assets and 2) productivity
    n_fine::Array{Int64,1} # joint simulation number of gridpoints for 1) assets and 2) productivity
    agrid::Array{Float64,1} # joint approximation number of gridpoints for 1) assets and 2) productivity
    agrid_fine::Array{Float64,1} # joint simulation number of gridpoints for 1) assets and 2) productivity
    a_min::Float64 # Lower bound for the asset grid/ borrowing constraint
    a_max::Float64 # Lower bound for the asset grid/ borrowing constraint
    spliorder::Int64 # Order of the spline approximation
    fspace_a::Dict{Symbol,Any} # function space for the approximation of the value function
    fspace_a_fine::Dict{Symbol,Any} # function space for the approximation of the distribution
    τ::Array{Float64,1} # Value added export tariff rate
# Country specific processes:
    Country_spec_p
# This stores the following objects
    #s: matrix of states, with the first column for assets, second for productivity
    #ns: total number of states, taken as a cartesian product of "s"
    #s_fine: matrix of states for the distribution, with the first column for assets, second for productivity
    #ns_fine: total number of states for the distribution, taken as a cartesian product of "s_fine"
    #Phi_z: basis matrix for the productivity
    #Phi_z_fine: basis matrix for the productivity in the distribution
    #Phi joint basis matrix
    #Phi_aug joint basis matrix with occupation choice
    #P_kron,P_kron1,P_kron_fine, are useful object for the distribution
    #exp_egrid is the grid of productivities
# Optimization bounds:
    bounds::Array{Float64,2}
end
#Redefine copy to allow for this in the parameter types.
Base.copy(x::T) where T = T([deepcopy(getfield(x, k)) for k ∈ fieldnames(T)]...)
function setup_rouwen(rho_Rouw::Number, mu_uncond::Number, sig_uncond::Number, n_R::Int)
    step_R = sig_uncond*sqrt(n_R - 1)
    z_Rouw = Array(-1:2/(n_R-1):1);
    z_Rouw = mu_uncond*ones(size(z_Rouw))+step_R.*z_Rouw;
    p=(rho_Rouw + 1)/2;
    q=p;
    P_Rouw=[ p  (1-p);
            (1-q) q];
    for i_R=2:n_R-1
        a1R=[P_Rouw zeros(i_R, 1); zeros(1, i_R+1)];
        a2R=[zeros(i_R, 1) P_Rouw; zeros(1, i_R+1)];
        a3R=[zeros(1,i_R+1); P_Rouw zeros(i_R,1)];
        a4R=[zeros(1,i_R+1); zeros(i_R,1) P_Rouw];
        P_Rouw=p*a1R+(1-p)*a2R+(1-q)*a3R+q*a4R;
        P_Rouw[2:i_R, :] = P_Rouw[2:i_R, :]/2;
    end
    P_Rouw=P_Rouw';
    for i_R = 1:n_R
        P_Rouw[:,i_R] = P_Rouw[:,i_R]/sum(P_Rouw[:,i_R]);
    end
    return (P_Rouw, z_Rouw)
end
function goldenx(f::Function,a::Array{Float64,1},b::Array{Float64,1},tol::Float64=1e-10,arg...)
    """Vectorized golden section search to maximize univariate functions simultaneously
    Returns the maximum and the maximal value of f. Closer to Matlab code
    Parameters
    ----------
    f : function
        Function that maps from R^n to R^n such that it is an augmented univariate function
    a : array_like
        The lower bound for the maximization, for each dimension
    b : array_like
        The upper bound for the maximization, for each dimension
    tol : float
        Specifies the default tolerance value (for the argument of f) - default is 10**(-10)
    Returns
    -------
    x1 : ndarray
       Returns the n dimensional solution of the maximization.
    f1 : ndarray
       Returns the n dimensional maximum values.
    Notes
    -----
    To test in Julia:
    f(x::Array{Float64,1},c::Array{Float64,1},d::Float64) = -d*(x-c).^2;
    c = Array(range(1,10,length=10));
    a = zeros(10);
    b = ones(10) * 20;
    tol = 1e-10;
    d= 2.0;
    x, val = goldenx(f,a,b,tol,c,d);
    """

    alpha1 = (3.0 - sqrt(5)) / 2.0;
    alpha2 = 1.0 - alpha1;
    d  = b - a;
    x1 = a + alpha1 * d;
    x2 = a + alpha2 * d;
    n_tmp = size(x1);
    ones_tmp = ones(n_tmp);
    s  = ones_tmp;
    f1 = f(x1,arg...);
    f2 = f(x2,arg...);
    d = alpha1 * alpha2 * d;
    conv = 2.0;
    while conv > tol
        i = f2.>f1;
        not_i = ones_tmp - i;
        x1[i] = x2[i];
        f1[i] = f2[i];
        d = alpha2 * d;
        x2 = x1 + s.*(i - not_i).*d;
        s = sign.(x2-x1);
        f2 = f(x2,arg...);
        conv = maximum(d);
    end
    return x1, f1
end
function setup_state_space(country::Int64, parameters_tmp::Parameter_type)
    # Set up the productivity process
    (P_Rouw, z_Rouw) = setup_rouwen( parameters_tmp.ρ[country],parameters_tmp.Avg_prod[country],
        parameters_tmp.σₛ[country] /sqrt(1-parameters_tmp.ρ[country]^2),parameters_tmp.n[2] );
    exp_egrid = exp.(z_Rouw);
    P = P_Rouw';
    check_nonzero = P.>1e-10;
    P = P.*check_nonzero;
    # Set up the function spaces
    fspace = fundef((:spli, parameters_tmp.agrid, 0,parameters_tmp.spliorder),(:spli, exp_egrid,0,1));# joint
    fspace_fine = fundef((:spli, parameters_tmp.agrid_fine, 0,1),(:spli, exp_egrid,0,1));# joint distribution
    fspace_x = fundef((:spli, exp_egrid,0,1));# only productivity
    # Set up the grids
    grid = funnode(fspace);
    s    = grid[1];
    ns = length(s[:,1]);
    grid_fine = funnode(fspace_fine);
    s_fine    = grid_fine[1];
    ns_fine   = length(s_fine[:,1]);
    # Set up the basis matrices
    Phi_z = funbase(fspace_x, s[:,2]);
    Phi_z_fine = funbase(fspace_x, s_fine[:,2]);
    Phi = funbase(fspace, s);
    check_nonzero = Phi.>1e-10;
    Phi = Phi.*check_nonzero;
    Phi_aug =  kron(Matrix(1.0I, 3, 3) ,Phi);
    P_kron = kron(P,Matrix(1.0I, parameters_tmp.n[1], parameters_tmp.n[1]));
    P_kron1 =  kron(Matrix(1.0I, 3, 3),P_kron);
    P_kron_fine = kron(ones(parameters_tmp.n_fine[1],1),P);
    #Eliminate small entries
    check_nonzero = P_kron.>1e-10;
    P_kron = P_kron.*check_nonzero;
    check_nonzero = P_kron1.>1e-10;
    P_kron1 = P_kron1.*check_nonzero;
    check_nonzero = P_kron_fine.>1e-10;
    P_kron_fine = P_kron_fine.*check_nonzero;
    #Force the rest to be in sparse format
    P_kron = SparseArrays.sparse(P_kron);
    P_kron1 = SparseArrays.sparse(P_kron1);
    P_kron_fine = SparseArrays.sparse(P_kron_fine);
    return (s,ns,s_fine,ns_fine,Phi_z,Phi_z_fine,Phi,Phi_aug,P_kron,P_kron1,P_kron_fine,exp_egrid)
end

function optimal_k( R::Float64,W_loc::Float64,constant_z::Array{Float64,1},
        constant_z_bar::Array{Float64,1},σ::Float64,α₁_eff::Float64,α₂_eff::Float64)
    #Unconditionally optimal k - needed for individual incomes
    k_d = (α₁_eff^(1- α₂_eff) * α₂_eff^(α₂_eff) * constant_z_bar * R^(α₂_eff - 1) * W_loc^(-α₂_eff)).^σ;
    k_x = (α₁_eff^(1- α₂_eff) * α₂_eff^(α₂_eff) * constant_z * R^(α₂_eff - 1) * W_loc^(-α₂_eff)).^σ;
    return (k_d,k_x)
end
function true_profit(price_final_prev::Float64, R::Float64,W_loc::Float64 ,s_curr::Array{Float64,2},
        constant_x::Float64,constant_d::Float64,constant_z::Array{Float64,1},constant_z_bar::Array{Float64,1},
        θ_loc::Float64,τ_loc::Float64,σ::Float64,α₁_eff::Float64,α₂_eff::Float64,α₁::Float64,α₂::Float64,
        ones_tmp::Array{Float64,1})
    # Given prices, evaluate the income of households needed for VFI
    (k_d,k_x ) = optimal_k(R,W_loc,constant_z,constant_z_bar ,σ,α₁_eff,α₂_eff);
    real_value_asset = s_curr[:,1]/((1 - θ_loc) * price_final_prev);
    which_k_d = k_d .< real_value_asset;
    which_k_x = k_x .< real_value_asset;
    k_choice_d = which_k_d .*k_d + (ones_tmp - which_k_d).*(real_value_asset);
    k_choice_x = which_k_x .*k_x + (ones_tmp - which_k_x).*(real_value_asset);
    #Implied Lagrange multiplier
    lambdda_d = ((k_choice_d./k_d).^(1/σ /(α₂_eff - 1)) - ones_tmp )*R;
    lambdda_x = ((k_choice_x./k_x).^(1/σ /(α₂_eff - 1)) - ones_tmp )*R;
    labor_d = (α₁_eff^(α₁_eff) * α₂_eff^(1 - α₁_eff) * constant_z_bar.* (R * ones_tmp + lambdda_d).^(-α₁_eff) * W_loc^(α₁_eff-1)).^σ;
    labor_x = (α₁_eff^(α₁_eff) * α₂_eff^(1 - α₁_eff) * constant_z.* (R* ones_tmp + lambdda_x).^(-α₁_eff) * W_loc^(α₁_eff-1)).^σ;
    total_produced_x = (s_curr[:,2]).* k_choice_x.^(α₁).* labor_x.^(α₂);
    output_dx =total_produced_x * constant_d^σ / (constant_d^σ + (1+τ_loc)*(constant_x)^σ);
    output_xx =total_produced_x * constant_x^σ / (constant_d^σ + (1+τ_loc)*(constant_x)^σ); # Note - iceberg cost paid here
    output_d = ( s_curr[:,2]).* k_choice_d.^(α₁).* labor_d.^(α₂);
    price_dx = constant_d .* output_dx.^(-1/σ);
    price_d = constant_d .* output_d.^(-1/σ);
    price_xx = (1+τ_loc)* constant_x .* output_xx.^(-1/σ); # final price without tariffs
    rev_d = price_d.*output_d ;
    rev_dx = price_dx.*output_dx;
    rev_xx = price_xx.*output_xx;
    Profits_d = rev_d - W_loc*labor_d - R *k_choice_d;
    Profits_x = rev_dx + rev_xx   - W_loc*labor_x - R *k_choice_x;
    return (Profits_d,Profits_x,output_dx,output_xx,output_d,labor_d,labor_x,k_choice_d,
        k_choice_x,price_dx,price_d,price_xx,lambdda_d,lambdda_x,rev_d,rev_dx,rev_xx,k_d,k_x )
end

function BellmanF_stable(x::Array{Float64,1},coeff::Array{Float64,2},Phi_z_co::SparseMatrixCSC{Float64,Int64},
        fspace_a_co::Dict{Symbol,Any},income_mat_co::Array{Float64,1},β_applied::Float64,
        future_occupation_guess::Int64,price_final_current::Float64,Return_Phi::Int64=0)
    #Current utility conditional on occupation
    cons_guess = (income_mat_co- x)/price_final_current;
    curr_util = @fastmath log.(cons_guess);
    Phi_prime_a = funbase(fspace_a_co, x);
    Phi_prime = row_kron(Phi_z_co,Phi_prime_a); # Note - check if RowKron==dprod
    #Phi_prime = row_kron(Phi_prime_a,Phi_z_co);
    V_fut_mat = Phi_prime * coeff[:,future_occupation_guess];
    util= curr_util + β_applied*V_fut_mat;
    if Return_Phi==0
        return util
    elseif Return_Phi==1
        return util , Phi_prime
    end
end

function local_parameters(parameters_tmp::Parameter_type)
    # Function that extract the structure parameters for better readability
    β = parameters_tmp.β;
    α = parameters_tmp.α;
    δ = parameters_tmp.δ;
    θ = parameters_tmp.θ;
    α₁ = parameters_tmp.α₁;
    α₂ = parameters_tmp.α₂;
    σ = parameters_tmp.σ;
    α₁_eff = parameters_tmp.α₁_eff;
    α₂_eff = parameters_tmp.α₂_eff;
    ω = parameters_tmp.ω;
    L = parameters_tmp.L;
    FM = parameters_tmp.FM;
    FMₓ = parameters_tmp.FMₓ;
    F = parameters_tmp.F;
    Fₓ = parameters_tmp.Fₓ;
    Exit = parameters_tmp.Exit;
    Exitₓ = parameters_tmp.Exitₓ;
    iid_cost_value =parameters_tmp.iid_cost_value;
    iid_cost_prob = parameters_tmp.iid_cost_prob;
    size_iid_cost_val = size(iid_cost_value,1);
    country_no = parameters_tmp.country_no;
    τ = parameters_tmp.τ;
    a_min = parameters_tmp.a_min;
    a_max = parameters_tmp.a_max;
    fspace_a = parameters_tmp.fspace_a;
    fspace_a_fine = parameters_tmp.fspace_a_fine;
    agrid_fine = parameters_tmp.agrid_fine;
    banking_cost = parameters_tmp.banking_cost;
    bounds = parameters_tmp.bounds;
    # Extract country specific parameters:
    Country_spec_p = parameters_tmp.Country_spec_p;
    s_cell = [Country_spec_p[i][1] for i=1:parameters_tmp.country_no];
    ns_cell= [Country_spec_p[i][2] for i=1:parameters_tmp.country_no];
    s_fine_cell= [Country_spec_p[i][3] for i=1:parameters_tmp.country_no];
    ns_fine_cell= [Country_spec_p[i][4] for i=1:parameters_tmp.country_no];
    Phi_z_cell= [Country_spec_p[i][5] for i=1:parameters_tmp.country_no];
    Phi_z_fine_cell= [Country_spec_p[i][6] for i=1:parameters_tmp.country_no];
    Phi_cell= [Country_spec_p[i][7] for i=1:parameters_tmp.country_no];
    Phi_aug_cell = [Country_spec_p[i][8] for i=1:parameters_tmp.country_no];
    P_kron_cell= [Country_spec_p[i][9] for i=1:parameters_tmp.country_no];
    P_kron1_cell= [Country_spec_p[i][10] for i=1:parameters_tmp.country_no];
    P_kron_fine_cell= [Country_spec_p[i][11] for i=1:parameters_tmp.country_no];
    exp_egrid_cell= [Country_spec_p[i][12] for i=1:parameters_tmp.country_no];
    #
    ns_tmp = ns_cell[1];
    ns_tmp_fine = ns_fine_cell[1];
    openness = parameters_tmp.openness;
    return (β,α,δ,θ,α₁,α₂,σ,α₁_eff,α₂_eff,ω,L,FM,FMₓ,F,Fₓ,
        Exit,Exitₓ,iid_cost_value,iid_cost_prob,size_iid_cost_val,country_no,τ,
        a_min,a_max,fspace_a,fspace_a_fine,agrid_fine,banking_cost,
        bounds,Country_spec_p,s_cell,ns_cell,s_fine_cell,ns_fine_cell,Phi_z_cell,
        Phi_z_fine_cell,Phi_cell,Phi_aug_cell,P_kron_cell,P_kron1_cell,P_kron_fine_cell,
        exp_egrid_cell,ns_tmp,ns_tmp_fine,openness)
end
function price_reshaper(prices::Array{Float64,1},openness::Int64,country_no::Int64,bounds::Array{Float64})
    if openness==0
        W      = prices[1:country_no];
        output_final = prices[(country_no + 1):(2 *country_no)];
        price_final = ones(country_no,2);
        price_final[1:(country_no - 1),1] = prices[(2*country_no + 1):(3 *country_no-1)];
        r      = prices[(3*country_no):(4 *country_no-1)];
        residual = SharedArray{Float64}(4*country_no,1);
    elseif openness ==1
        W      = prices[1:country_no];
        output_final = prices[(country_no + 1):(2 *country_no)];
        price_final = ones(country_no,2);
        price_final[1:(country_no - 1),1] = prices[(2*country_no + 1):(3 *country_no - 1)];
        r      = prices[3 *country_no] * ones(country_no);
        residual = SharedArray{Float64}(3*country_no + 1,1);
        residual[3*country_no + 1,1] = 0;
    end
    price_final[:,2] = price_final[:,1];
    # Steady state: take prices as given from "last period",
    # for steady state this is obvious
    price_check_tmp = sum((minimum(W)<bounds[1,1])  + (minimum(r)<bounds[2,1]) + (maximum(r)>bounds[2,2]
            ) + (minimum(price_final[:,1])<bounds[3,1]) + (minimum(output_final)<bounds[4,1]));
    return W,r,output_final,price_final,r,residual,price_check_tmp
end
function price_reshaper_single(prices::Array{Float64,1},prices_past::Array{Float64,1},openness::Int64,country_no::Int64,
    bounds::Array{Float64},country::Int64,banking_cost_loc::Float64,δ::Float64,ω::Float64,σ::Float64,τ_loc::Float64,
    z_tilde::Array{Float64,1},z_tilde_fine::Array{Float64,1})
    if openness==0
        W      = prices[1:country_no];
        output_final = prices[(country_no + 1):(2 *country_no)];
        price_final = ones(country_no,2);
        price_final[1:(country_no - 1),1] = prices[(2*country_no + 1):(3 *country_no-1)];
        price_final[1:(country_no - 1),2] = prices_past[(2*country_no + 1):(3 *country_no - 1)];
        r      = prices[(3*country_no):(4 *country_no-1)];
    elseif openness ==1
        W      = prices[1:country_no];
        output_final = prices[(country_no + 1):(2 *country_no)];
        price_final = ones(country_no,2);
        price_final[1:(country_no - 1),1] = prices[(2*country_no + 1):(3 *country_no - 1)];
        price_final[1:(country_no - 1),2] = prices_past[(2*country_no + 1):(3 *country_no - 1)];
        r      = prices[3 *country_no] * ones(country_no);
    end
    residual = Array{Float64}(undef,4,1);
    r_loc = r[country,1];
    W_loc = W[country,1];
    price_final_current = price_final[country,1];
    price_final_prev = price_final[country,2];
    output_final_loc = output_final[country,1]
    R = maximum([price_final_prev * (1.0 + r_loc + r_loc*banking_cost_loc + banking_cost_loc
        - price_final_current/price_final_prev * (1.0 - δ)),0]);
    avg_foreign_price = (sum(price_final[:,2] ) - price_final[country,2])/(country_no - 1);
    total_foreign_demand = (sum(output_final) - output_final_loc);
    constant_d = ω * price_final_current * output_final_loc^(1/σ);
    constant_x = (1.0 - ω)/ (1.0 + τ_loc) * avg_foreign_price * total_foreign_demand^(1.0/σ);
    constant_z = z_tilde* (constant_d^σ + (1+τ_loc)*(constant_x)^σ)^(1/σ);
    constant_z_bar = z_tilde* constant_d;
    constant_z_fine = z_tilde_fine* (constant_d^σ + (1+τ_loc)*(constant_x)^σ)^(1/σ);
    constant_z_bar_fine = z_tilde_fine * constant_d;
    price_check_tmp = sum((W_loc<bounds[1,1])+(W_loc>bounds[1,2])  + (r_loc<bounds[2,1]) + (r_loc>bounds[2,2]
            ) + (price_final_current<bounds[3,1]) + (price_final_current>bounds[3,2])  + (output_final_loc<bounds[4,1]) + (output_final_loc>bounds[4,2]));
    return (price_final_prev,price_final_current,total_foreign_demand,W_loc,r_loc,avg_foreign_price,R,constant_d,
    constant_x,constant_z,constant_z_bar,constant_z_fine,constant_z_bar_fine,output_final_loc,
    residual,price_check_tmp)
end
function local_var_creator_minimal(country_no::Int64)
    # Initialization of local variables - expect overhead:
    export_prod = SharedArray{Float64}(country_no,1);
    export_price_sum = SharedArray{Float64}(country_no,1);
    import_price = Array{Float64}(undef,country_no,1);
    domestic_price_sum = SharedArray{Float64}(country_no,1);
    price_final_actual = Array{Float64}(undef,country_no,1);
    NFA = SharedArray{Float64}(country_no,1);
    asset_demand = SharedArray{Float64}(country_no,1);
    asset_supply = SharedArray{Float64}(country_no,1);
    exitflag = SharedArray{Int64}(country_no,1);
    return (export_prod,export_price_sum,import_price,domestic_price_sum,
    price_final_actual,NFA,asset_supply,asset_demand,exitflag)
end
function local_var_creator_detailed(country_no::Int64,ns_tmp::Int64,ns_tmp_fine::Int64,size_iid_cost_val::Int64)
    a_prime_fine = SharedArray{Float64}(ns_tmp_fine,3,country_no,size_iid_cost_val);
    future_occupation_fine = SharedArray{Float64}(ns_tmp_fine,3,country_no,size_iid_cost_val);
    cons_fine = SharedArray{Float64}(ns_tmp_fine,3,country_no,size_iid_cost_val);
    rev_d_fine_store = SharedArray{Float64}(ns_tmp_fine,country_no);
    rev_dx_fine_store = SharedArray{Float64}(ns_tmp_fine,country_no);
    rev_xx_fine_store = SharedArray{Float64}(ns_tmp_fine,country_no);
    coeff_final = SharedArray{Float64}(ns_tmp,3,country_no);
    k_x_fine = SharedArray{Float64}(ns_tmp_fine,country_no);
    k_d_fine = SharedArray{Float64}(ns_tmp_fine,country_no);
    k_opt_d_fine= SharedArray{Float64}(ns_tmp_fine,country_no);
    k_opt_x_fine= SharedArray{Float64}(ns_tmp_fine,country_no);
    l_x_fine = SharedArray{Float64}(ns_tmp_fine,country_no);
    l_d_fine = SharedArray{Float64}(ns_tmp_fine,country_no);
    current_distr_store = SharedArray{Float64}(ns_tmp_fine*3,country_no);
    distr_current = SharedArray{Float64}(ns_tmp_fine*3,country_no);
    exitflag = SharedArray{Float64}(country_no,1); # Set to 1 if it goes wrong within VFI, 2 if later
    domestic_prod = SharedArray{Float64}(country_no,1);
    exporter_pop = SharedArray{Float64}(country_no,1);
    domestic_pop = SharedArray{Float64}(country_no,1);
    worker_pop = SharedArray{Float64}(country_no,1);
    entry_share_to_domestic = SharedArray{Float64}(country_no,1);
    exit_share_to_domestic = SharedArray{Float64}(country_no,1);
    entry_share_to_exporter = SharedArray{Float64}(country_no,1);
    exit_share_to_worker = SharedArray{Float64}(country_no,1);
    total_entry_cost = SharedArray{Float64}(country_no,1);
    total_incumbents_cost = SharedArray{Float64}(country_no,1);
    total_exit_cost= SharedArray{Float64}(country_no,1);
    labor_excess_demand = SharedArray{Float64}(country_no,1);
    labor_demand = SharedArray{Float64}(country_no,1);
    total_consumption = SharedArray{Float64}(country_no,1);
    capital_demand = SharedArray{Float64}(country_no,1);
    capital_next = SharedArray{Float64}(country_no,1);
    total_demand_final_good = SharedArray{Float64}(country_no,1);
    excess_demand_final_good = SharedArray{Float64}(country_no,1);
    GDP = SharedArray{Float64}(country_no,1);
    nomGDP = SharedArray{Float64}(country_no,1);
    PPI = SharedArray{Float64}(country_no,1);
    TFP = SharedArray{Float64}(country_no,1);
    TFP_within = SharedArray{Float64}(country_no,1);
    TFP_across = SharedArray{Float64}(country_no,1);
    TFP_second_best = SharedArray{Float64}(country_no,1);
    TOT = SharedArray{Float64}(country_no,1);
    total_production = SharedArray{Float64}(country_no,1);
    domestic_price_sum = SharedArray{Float64}(country_no,1);
    export_value= SharedArray{Float64}(country_no,1);
    investment = SharedArray{Float64}(country_no,1);
    export_value= SharedArray{Float64}(country_no,1);
    L_x = SharedArray{Float64}(country_no,1);
    L_d = SharedArray{Float64}(country_no,1);
    K_x = SharedArray{Float64}(country_no,1);
    K_d = SharedArray{Float64}(country_no,1);
    x_k =SharedArray{Float64}(country_no,1); # capital allocated to exports
    x_l =SharedArray{Float64}(country_no,1); # labor allocated to exports
    D_d = SharedArray{Float64}(country_no,1); #exporting tfp construction
    D_x = SharedArray{Float64}(country_no,1); #domestic tfp construction
    D_d_denom = SharedArray{Float64}(country_no,1); #exporting tfp construction
    D_x_denom = SharedArray{Float64}(country_no,1); #domestic tfp construction
    D_d_eff = SharedArray{Float64}(country_no,1); #exporting tfp construction
    D_x_eff = SharedArray{Float64}(country_no,1); #domestic tfp construction
    D_d_denom_eff = SharedArray{Float64}(country_no,1); #exporting tfp construction
    D_x_denom_eff = SharedArray{Float64}(country_no,1); #domestic tfp construction
    K_d_ratio = SharedArray{Float64}(country_no,1);
    K_x_ratio = SharedArray{Float64}(country_no,1);
    L_d_ratio = SharedArray{Float64}(country_no,1);
    L_x_ratio = SharedArray{Float64}(country_no,1);
    K_x_ratio_eff = SharedArray{Float64}(country_no,1);
    K_d_ratio_eff = SharedArray{Float64}(country_no,1);
    L_d_ratio_eff = SharedArray{Float64}(country_no,1);
    L_x_ratio_eff = SharedArray{Float64}(country_no,1);
    nomGDP_d = SharedArray{Float64}(country_no,1);
    nomGDP_xd = SharedArray{Float64}(country_no,1);
    nomGDP_xx = SharedArray{Float64}(country_no,1);
    RGDP_d = SharedArray{Float64}(country_no,1);
    RGDP_xd = SharedArray{Float64}(country_no,1);
    RGDP_xx = SharedArray{Float64}(country_no,1);
    RGDP = SharedArray{Float64}(country_no,1);
    mean_MRPK_d = SharedArray{Float64}(country_no,1);
    mean_MRPK_x = SharedArray{Float64}(country_no,1);
    sd_MRPK_d = SharedArray{Float64}(country_no,1);
    sd_MRPK_x = SharedArray{Float64}(country_no,1);
    mean_MRPK = SharedArray{Float64}(country_no,1);
    sd_MRPK = SharedArray{Float64}(country_no,1);
    mean_MRPL_d = SharedArray{Float64}(country_no,1);
    mean_MRPL_x = SharedArray{Float64}(country_no,1);
    sd_MRPL_d = SharedArray{Float64}(country_no,1);
    sd_MRPL_x = SharedArray{Float64}(country_no,1);
    mean_MRPL = SharedArray{Float64}(country_no,1);
    sd_MRPL = SharedArray{Float64}(country_no,1);
    cov_TFPR_z_d = SharedArray{Float64}(country_no,1);
    cov_TFPR_z_x = SharedArray{Float64}(country_no,1);
    cov_TFPR_z = SharedArray{Float64}(country_no,1);
    corr_TFPR_z_d = SharedArray{Float64}(country_no,1);
    corr_TFPR_z_x = SharedArray{Float64}(country_no,1);
    corr_TFPR_z = SharedArray{Float64}(country_no,1);
    mean_logProd_d = SharedArray{Float64}(country_no,1);
    mean_logProd_x = SharedArray{Float64}(country_no,1);
    mean_logProd = SharedArray{Float64}(country_no,1);
    sd_logProd_d = SharedArray{Float64}(country_no,1);
    sd_logProd_x = SharedArray{Float64}(country_no,1);
    sd_logProd = SharedArray{Float64}(country_no,1);
    import_share = SharedArray{Float64}(country_no,1);
    worker_bond_holding = SharedArray{Float64}(country_no,1);
    domestic_bond_holding= SharedArray{Float64}(country_no,1);
    exporter_bond_holding= SharedArray{Float64}(country_no,1);
    domestic_firm_debt= SharedArray{Float64}(country_no,1);
    exporter_firm_debt= SharedArray{Float64}(country_no,1);
    #Profits:
    avg_Pi_x = SharedArray{Float64}(country_no,1);
    avg_Pi_d = SharedArray{Float64}(country_no,1);
    Profit_fine_x = SharedArray{Float64}(ns_tmp_fine,country_no);
    Profit_fine_d = SharedArray{Float64}(ns_tmp_fine,country_no);
    #Inequality:
    mean_cons = SharedArray{Float64}(country_no,1);
    sd_cons = SharedArray{Float64}(country_no,1);
    mean_income =  SharedArray{Float64}(country_no,1);
    sd_income =  SharedArray{Float64}(country_no,1);
    mean_wealth = SharedArray{Float64}(country_no,1);
    sd_wealth = SharedArray{Float64}(country_no,1);
    wealth_of_workers = SharedArray{Float64}(country_no,1);
    wealth_of_domestic = SharedArray{Float64}(country_no,1);
    wealth_of_exporters = SharedArray{Float64}(country_no,1);
    income_of_workers = SharedArray{Float64}(country_no,1);
    income_of_domestic = SharedArray{Float64}(country_no,1);
    income_of_exporters = SharedArray{Float64}(country_no,1);
    cons_of_workers = SharedArray{Float64}(country_no,1);
    cons_of_domestic = SharedArray{Float64}(country_no,1);
    cons_of_exporters = SharedArray{Float64}(country_no,1);
    GINI_wealth= SharedArray{Float64}(country_no,1);
    GINI_income= SharedArray{Float64}(country_no,1);
    GINI_cons = SharedArray{Float64}(country_no,1);
    p10_wealth = SharedArray{Float64}(country_no,1);
    p50_wealth = SharedArray{Float64}(country_no,1);
    p90_wealth = SharedArray{Float64}(country_no,1);
    p99_wealth = SharedArray{Float64}(country_no,1);
    p10_income = SharedArray{Float64}(country_no,1);
    p50_income = SharedArray{Float64}(country_no,1);
    p90_income = SharedArray{Float64}(country_no,1);
    p99_income = SharedArray{Float64}(country_no,1);
    p10_cons = SharedArray{Float64}(country_no,1);
    p50_cons = SharedArray{Float64}(country_no,1);
    p90_cons = SharedArray{Float64}(country_no,1);
    p99_cons = SharedArray{Float64}(country_no,1);
    # Banking real cost of monitoring
    Aggr_bank_cost = SharedArray{Float64}(country_no,1);
    mean_leverage_d= SharedArray{Float64}(country_no,1);
    mean_leverage_x= SharedArray{Float64}(country_no,1);
    mean_leverage= SharedArray{Float64}(country_no,1);
    sd_leverage_d= SharedArray{Float64}(country_no,1);
    sd_leverage_x= SharedArray{Float64}(country_no,1);
    sd_leverage= SharedArray{Float64}(country_no,1);
    corr_MRPK_lev_d= SharedArray{Float64}(country_no,1);
    corr_MRPK_lev_x= SharedArray{Float64}(country_no,1);
    corr_MRPK_lev= SharedArray{Float64}(country_no,1);
    exit_domestic_to_work_sum= SharedArray{Float64}(country_no,1);
    exit_exporting_to_work_sum= SharedArray{Float64}(country_no,1);
    mean_growth_rev= SharedArray{Float64}(country_no,1);
    sd_growth_rev= SharedArray{Float64}(country_no,1);
    mean_growth_k= SharedArray{Float64}(country_no,1);
    sd_growth_k= SharedArray{Float64}(country_no,1);
    autocorr_rev= SharedArray{Float64}(country_no,1);
    fraction_zombie_exporter= SharedArray{Float64}(country_no,1);
    return (a_prime_fine,future_occupation_fine,cons_fine,rev_d_fine_store,
        rev_dx_fine_store,rev_xx_fine_store,coeff_final,k_x_fine,k_d_fine,
        l_x_fine,l_d_fine,current_distr_store,distr_current,exitflag,domestic_prod,exporter_pop,
        domestic_pop,worker_pop,entry_share_to_domestic,exit_share_to_domestic,
        entry_share_to_exporter,exit_share_to_worker,total_entry_cost,
        total_incumbents_cost,total_exit_cost,labor_excess_demand,labor_demand,
        total_consumption,capital_demand,capital_next,
        total_demand_final_good,excess_demand_final_good,GDP,nomGDP,PPI,TFP,
        TFP_within,TFP_across,TFP_second_best,TOT,total_production,
        domestic_price_sum,export_value,investment,
        export_value,L_x,L_d,K_x,K_d,x_k,x_l,D_d,D_x,D_d_denom,D_x_denom,
        D_d_eff,D_x_eff,D_d_denom_eff,D_x_denom_eff,K_d_ratio,K_x_ratio,
        L_d_ratio,L_x_ratio,K_x_ratio_eff,K_d_ratio_eff,L_d_ratio_eff,
        L_x_ratio_eff,nomGDP_d,nomGDP_xd,nomGDP_xx,RGDP_d,RGDP_xd,RGDP_xx,RGDP,
        mean_MRPK_d,mean_MRPK_x,sd_MRPK_d,sd_MRPK_x,mean_MRPK,sd_MRPK,
        mean_MRPL_d,mean_MRPL_x,sd_MRPL_d,sd_MRPL_x,mean_MRPL,sd_MRPL,
        cov_TFPR_z_d,cov_TFPR_z_x,cov_TFPR_z,corr_TFPR_z_d,corr_TFPR_z_x,
        corr_TFPR_z,mean_logProd_d,mean_logProd_x,mean_logProd,sd_logProd_d,
        sd_logProd_x,sd_logProd,import_share,worker_bond_holding,
        domestic_bond_holding,exporter_bond_holding,domestic_firm_debt,
        exporter_firm_debt,avg_Pi_x,avg_Pi_d,Profit_fine_x,
        Profit_fine_d,mean_cons,sd_cons,mean_income,sd_income,mean_wealth,
        sd_wealth,wealth_of_workers,wealth_of_domestic,wealth_of_exporters,
        income_of_workers,income_of_domestic,income_of_exporters,
        cons_of_workers,cons_of_domestic,cons_of_exporters,GINI_wealth,
        GINI_income,GINI_cons,p10_wealth,p50_wealth,p90_wealth,p99_wealth,p10_income,p50_income,p90_income,
        p99_income,p10_cons,p50_cons,p90_cons,p99_cons,Aggr_bank_cost,mean_leverage_d,mean_leverage_x,mean_leverage,
        sd_leverage_d,sd_leverage_x,sd_leverage,corr_MRPK_lev_d,corr_MRPK_lev_x,corr_MRPK_lev,exit_domestic_to_work_sum,
        exit_exporting_to_work_sum,mean_growth_rev,sd_growth_rev,mean_growth_k,sd_growth_k,
        autocorr_rev,fraction_zombie_exporter,k_opt_d_fine,k_opt_x_fine)
end

function country_local_initialize(country::Int64,s_cell::Array{Array{Float64,2},1},
    ns_cell::Array{Int64,1},s_fine_cell::Array{Array{Float64,2},1},
    ns_fine_cell::Array{Int64,1},
    Phi_z_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    Phi_z_fine_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    Phi_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    Phi_aug_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    P_kron_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    P_kron1_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    P_kron_fine_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    σ::Float64,size_iid_cost_val::Int64,price_final::Array{Float64,2},
    W::Array{Float64,1}, r::Array{Float64,1}, output_final::Array{Float64,1},
    θ::Array{Float64,1}, L::Array{Float64,1}, τ::Array{Float64,1},
    country_no::Int64,banking_cost::Array{Float64,1},δ::Float64,ω::Float64,
    β::Array{Float64,1})
    exitflag_tmp = 0;
    #Load country specific productivity process:
    s = s_cell[country];
    ns = ns_cell[country];
    s_fine = s_fine_cell[country];
    ns_fine = ns_fine_cell[country];
    Phi_z = Phi_z_cell[country];
    Phi_z_fine = Phi_z_fine_cell[country];
    Phi = Phi_cell[country];
    Phi_aug = Phi_aug_cell[country];
    P_kron = P_kron_cell[country];
    P_kron1 = P_kron1_cell[country];
    P_kron_fine = P_kron_fine_cell[country];
    z_tilde = ( s[:,2]).^((σ - 1)/σ);
    z_tilde_fine = (s_fine[:,2]).^((σ - 1)/σ);
    income_mat =  Array{Float64}(undef, ns,3,3,size_iid_cost_val);
    income_mat_fine = Array{Float64}(undef, ns_fine,3,3,size_iid_cost_val);
    a_prime_fine_local = Array{Float64}(undef, ns_fine,3,size_iid_cost_val);
    future_occupation_fine_local = Array{Float64}(undef, ns_fine,3,size_iid_cost_val);
    cons_fine_local = Array{Float64}(undef, ns_fine,3,size_iid_cost_val);
    #Resetting coefficients:
    coeff = zeros(ns,3);#Array{Float64}(undef, ns,3);
    coeff_next = zeros(ns,3);
    price_final_prev = price_final[country,2];
    price_final_current = price_final[country,1];
    total_foreign_demand = (sum(output_final) - output_final[country,1]);
    W_loc = W[country,1];
    θ_loc = θ[country,1];
    L_loc = L[country,1];
    τ_loc = τ[country,1];
    r_loc = r[country,1];
    β_loc = β[country,1];
    banking_cost_loc = banking_cost[country];
    avg_foreign_price = (sum(price_final[:,2] ) - price_final[country,2])/(country_no - 1);
    R = price_final_prev * (1.0 + r_loc + r_loc*banking_cost_loc + banking_cost_loc
        - price_final_current/price_final_prev * (1.0 - δ));
    constant_d = ω * price_final_current * output_final[country,1]^(1/σ);
    constant_x = (1.0 - ω)/ (1.0 + τ_loc) * avg_foreign_price * total_foreign_demand^(1.0/σ);
    constant_z = z_tilde* (constant_d^σ + (1+τ_loc)*(constant_x)^σ)^(1/σ);
    constant_z_bar = z_tilde* constant_d;
    constant_z_fine = z_tilde_fine* (constant_d^σ + (1+τ_loc)*(constant_x)^σ)^(1/σ);
    constant_z_bar_fine = z_tilde_fine * constant_d;
    # Do it once for the VFI
    ones_tmp = ones(ns);
    x_tmp = zeros(ns,3,size_iid_cost_val);
    V_tmp = zeros(ns,3,size_iid_cost_val);
    V_next_stacked = zeros(ns*3,1);
    iterate1 = 0;
    conv = 0.0
    Phi_prime_tmp = Array{SparseMatrixCSC{Float64,Int64}}(undef, 3,size_iid_cost_val);
    D_deriv_tmp_block = Array{SparseMatrixCSC{Float64,Int64}}(undef, 3,1);
    Q_trans = spzeros(ns_fine*3,ns_fine*3);
    ones_tmp_fine = ones(ns_fine);
    return (exitflag_tmp,s,ns,s_fine,ns_fine,Phi_z,Phi_z_fine,Phi,Phi_aug,
    P_kron,P_kron1,P_kron_fine,z_tilde,z_tilde_fine,income_mat,income_mat_fine,
    a_prime_fine_local,future_occupation_fine_local,cons_fine_local,coeff,
    coeff_next,price_final_prev,price_final_current,total_foreign_demand,W_loc,
    θ_loc,L_loc,τ_loc,r_loc,banking_cost_loc,avg_foreign_price,R,constant_d,
    constant_x,constant_z,constant_z_bar,constant_z_fine,constant_z_bar_fine,
    ones_tmp,x_tmp,V_tmp,β_loc,V_next_stacked,iterate1,Phi_prime_tmp,
    D_deriv_tmp_block,conv,Q_trans,ones_tmp_fine)
end
function country_local_initialize_noprice(country::Int64,s_cell::Array{Array{Float64,2},1},
    ns_cell::Array{Int64,1},s_fine_cell::Array{Array{Float64,2},1},
    ns_fine_cell::Array{Int64,1},
    Phi_z_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    Phi_z_fine_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    Phi_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    Phi_aug_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    P_kron_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    P_kron1_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    P_kron_fine_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    σ::Float64,size_iid_cost_val::Int64,
    θ::Array{Float64,1}, L::Array{Float64,1}, τ::Array{Float64,1},
    country_no::Int64,banking_cost::Array{Float64,1},δ::Float64,ω::Float64,
    β::Array{Float64,1})
    exitflag_tmp = 0;
    #Load country specific productivity process:
    s = s_cell[country];
    ns = ns_cell[country];
    s_fine = s_fine_cell[country];
    ns_fine = ns_fine_cell[country];
    Phi_z = Phi_z_cell[country];
    Phi_z_fine = Phi_z_fine_cell[country];
    Phi = Phi_cell[country];
    Phi_aug = Phi_aug_cell[country];
    P_kron = P_kron_cell[country];
    P_kron1 = P_kron1_cell[country];
    P_kron_fine = P_kron_fine_cell[country];
    z_tilde = ( s[:,2]).^((σ - 1)/σ);
    z_tilde_fine = (s_fine[:,2]).^((σ - 1)/σ);
    income_mat =  Array{Float64}(undef, ns,3,3,size_iid_cost_val);
    income_mat_fine = Array{Float64}(undef, ns_fine,3,3,size_iid_cost_val);
    a_prime_fine_local = Array{Float64}(undef, ns_fine,3,size_iid_cost_val);
    future_occupation_fine_local = Array{Float64}(undef, ns_fine,3,size_iid_cost_val);
    cons_fine_local = Array{Float64}(undef, ns_fine,3,size_iid_cost_val);
    #Resetting coefficients:
    coeff = zeros(ns,3);#Array{Float64}(undef, ns,3);
    coeff_next = zeros(ns,3);
    θ_loc = θ[country,1];
    L_loc = L[country,1];
    τ_loc = τ[country,1];
    β_loc = β[country,1];
    banking_cost_loc = banking_cost[country];
    # Do it once for the VFI
    ones_tmp = ones(ns);
    x_tmp = zeros(ns,3,size_iid_cost_val);
    V_tmp = zeros(ns,3,size_iid_cost_val);
    V_next_stacked = zeros(ns*3,1);
    iterate1 = 0;
    conv = 0.0
    Phi_prime_tmp = Array{SparseMatrixCSC{Float64,Int64}}(undef, 3,size_iid_cost_val);
    D_deriv_tmp_block = Array{SparseMatrixCSC{Float64,Int64}}(undef, 3,1);
    Q_trans = spzeros(ns_fine*3,ns_fine*3);
    ones_tmp_fine = ones(ns_fine);
    return (exitflag_tmp,s,ns,s_fine,ns_fine,Phi_z,Phi_z_fine,Phi,Phi_aug,
    P_kron,P_kron1,P_kron_fine,z_tilde,z_tilde_fine,income_mat,income_mat_fine,
    a_prime_fine_local,future_occupation_fine_local,cons_fine_local,coeff,
    coeff_next,θ_loc,L_loc,τ_loc,banking_cost_loc,
    ones_tmp,x_tmp,V_tmp,β_loc,V_next_stacked,iterate1,Phi_prime_tmp,
    D_deriv_tmp_block,conv,Q_trans,ones_tmp_fine)
end

function income_creator(W_loc::Float64,ones_tmp::Array{Float64},
    Profits_d::Array{Float64},Profits_x::Array{Float64},FM::Array{Float64},
    FMₓ::Array{Float64},country::Int64,r_loc::Float64,s::Array{Float64,2},
    iid_cost_value::Array{Float64},F::Array{Float64}, Fₓ::Array{Float64},
    size_iid_cost_val::Int64,income_mat::Array{Float64,4},ns::Int64,
    a_min::Float64,a_max::Float64)
    income_1 = W_loc* ones_tmp;
    income_2 = Profits_d - W_loc * FM[country,1] * ones_tmp;
    income_3 = Profits_x - W_loc * FMₓ[country,1]* ones_tmp;
    for j_cost = 1:size_iid_cost_val
        # Index for the iid cost
        for j = 1:3
            #Create the income matrix for all possible scenarios, transitioning from occupation j to these occupation:
            income_mat[:,1,j,j_cost] = income_1 + (1 + r_loc ) * s[:,1];
            income_mat[:,2,j,j_cost] = income_2 - W_loc*iid_cost_value[j_cost]*  F[country,1] * (j == 1
                )*ones_tmp + (1 + r_loc ) * s[:,1];
            income_mat[:,3,j,j_cost] = income_3 - W_loc *iid_cost_value[j_cost]* Fₓ[country,1] * (j != 3
                )*ones_tmp + (1 + r_loc ) * s[:,1];
        end
    end
    max_x = zeros(size(income_mat));
    for j = 1:3
        for j_cost = 1:size_iid_cost_val
            for jj =1:3
                Applied_max_tmp = Array{Float64}(undef, ns,2);
                Applied_max_tmp[:,1] = income_mat[:,jj,j,j_cost];
                Applied_max_tmp[:,2] = a_max * ones_tmp;
                max_x_tmp = minimum(Applied_max_tmp,dims=2);
                Applied_max = Array{Float64}(undef, ns,2);
                Applied_max[:,1] = max_x_tmp;
                Applied_max[:,2] = a_min * ones_tmp;
                max_x_tmp = maximum(Applied_max,dims=2)[:,1];
                max_x[:,jj,j,j_cost] = max_x_tmp;
            end
        end
    end
    return income_mat,max_x
end
function Bellman_iteration(coeff::Array{Float64,2},coeff_next::Array{Float64,2},
    ns::Int64, x_tmp::Array{Float64,3},V_tmp::Array{Float64,3},ones_tmp::Array{Float64,1},
    size_iid_cost_val::Int64,iid_cost_prob::Array{Float64},income_mat::Array{Float64,4},
    P_kron::SparseMatrixCSC{Float64,Int64},Phi::SparseMatrixCSC{Float64,Int64},
    a_min::Float64,Phi_z::SparseMatrixCSC{Float64,Int64},β_loc::Float64,
    fspace_a::Dict{Symbol,Any},price_final_current::Float64,max_x::Array{Float64,4},
    tol::Float64 = 1e-10)
    for j = 1:3
        for j_cost = 1:size_iid_cost_val
            for jj =1:3
                max_x_tmp = max_x[:,jj,j,j_cost];
                x_tmp[:,jj,j_cost],fval = goldenx(BellmanF_stable,a_min*ones_tmp,
                    max_x_tmp,tol,coeff,Phi_z,fspace_a,income_mat[:,jj,j,j_cost],β_loc,jj,
                    price_final_current);
                V_tmp[:,jj,j_cost] = BellmanF_stable(x_tmp[:,jj,j_cost],coeff,
                    Phi_z,fspace_a,income_mat[:,jj,j,j_cost],β_loc,jj,price_final_current);
            end
        end
        V_tmp[isnan.(V_tmp)] .= -10000;
        W_tmp, future_occupation = findmax(V_tmp; dims=2);
        V_next =iid_cost_prob[1]* W_tmp[:,1,1];
        if size_iid_cost_val>1
            for j_cost = 2:size_iid_cost_val
                V_next = V_next + iid_cost_prob[j_cost]* W_tmp[:,1,j_cost];
            end
        end
        E_V_next = P_kron * V_next;
        coeff_next[:,j] = Phi\E_V_next;
    end
    conv = maximum(abs.(coeff_next - coeff));
    #println("Bellman iteration convergence: ", conv)
    return (coeff_next,conv)
end
function Newton_iteration(coeff::Array{Float64,2},coeff_next::Array{Float64,2},
    ns::Int64, x_tmp::Array{Float64,3},V_tmp::Array{Float64,3},ones_tmp::Array{Float64,1},
    size_iid_cost_val::Int64,iid_cost_prob::Array{Float64},income_mat::Array{Float64,4},
    P_kron::SparseMatrixCSC{Float64,Int64},Phi::SparseMatrixCSC{Float64,Int64},
    a_min::Float64,Phi_z::SparseMatrixCSC{Float64,Int64},β_loc::Float64,
    fspace_a::Dict{Symbol,Any},price_final_current::Float64,max_x::Array{Float64,4},
    V_next_stacked::Array{Float64,2},iterate1::Int64,
    Phi_prime_tmp::Array{SparseMatrixCSC{Float64,Int64}},
    D_deriv_tmp_block::Array{SparseMatrixCSC{Float64,Int64}},
    Phi_aug::SparseMatrixCSC{Float64,Int64},
    P_kron1::SparseMatrixCSC{Float64,Int64},exitflag_tmp::Int64,
    tol::Float64 = 1e-10)
    D_deriv = spzeros(ns*3,ns*3);
    for j = 1:3
        for j_cost = 1:size_iid_cost_val
            for jj =1:3
                max_x_tmp = max_x[:,jj,j,j_cost];
                x_tmp[:,jj,j_cost],fval = goldenx(BellmanF_stable,a_min*ones_tmp,max_x_tmp,
                    tol,coeff,Phi_z,fspace_a,income_mat[:,jj,j,j_cost],β_loc,jj,price_final_current);
                (V_tmp[:,jj,j_cost], Phi_prime_tmp_loc)= BellmanF_stable(x_tmp[:,jj,j_cost],coeff,
                    Phi_z,fspace_a,income_mat[:,jj,j,j_cost],β_loc,jj,price_final_current,1);
                check_nonzero = Phi_prime_tmp_loc.>tol;
                Phi_prime_tmp[jj,j_cost] = Phi_prime_tmp_loc.*check_nonzero;
            end
        end
        V_tmp[isnan.(V_tmp)] .= -10000;
        W_tmp, future_occupation = findmax(V_tmp; dims=2);
        future_occupation_index = getindex.(future_occupation,2);
        V_next =iid_cost_prob[1]* W_tmp[:,1,1];
        for jj=1:3
            D_deriv_tmp_block[jj,1] = iid_cost_prob[1] * P_kron * row_kron(sparse(hcat(float(
                        future_occupation_index[:,1,1].==jj))),Phi_prime_tmp[jj,1]);
        end
        if size_iid_cost_val>1
            for j_cost = 2:size_iid_cost_val
                V_next = V_next + iid_cost_prob[j_cost]* W_tmp[:,1,j_cost];
                for jj=1:3
                    D_deriv_tmp_block[jj,1] = (D_deriv_tmp_block[jj,1]+
                        iid_cost_prob[j_cost] * P_kron * row_kron(sparse(hcat(float(
                        future_occupation_index[:,1,j_cost].==jj))),Phi_prime_tmp[jj,j_cost]));
                end
            end
        end
        D_deriv[((j-1)*ns + 1):(j*ns),1:ns] = D_deriv_tmp_block[1,1];
        D_deriv[((j-1)*ns + 1):(j*ns),(ns+1):(2*ns)]= D_deriv_tmp_block[2,1];
        D_deriv[((j-1)*ns + 1):(j*ns),(2*ns+1):(3*ns)] = D_deriv_tmp_block[3,1];
        V_next_stacked[((j-1)*ns + 1):(j*ns),1] = V_next;
    end
    D = Phi_aug - β_loc * D_deriv;
    E_V_next = P_kron1 * V_next_stacked;
    g = Phi_aug *[coeff[:,1];coeff[:,2];coeff[:,3]] - E_V_next;
    improvement = D\g;
    improvement_mat = zeros(size(coeff));
    improvement_mat[:,1] = improvement[((1-1)*ns + 1):(1*ns)];
    improvement_mat[:,2] = improvement[((2-1)*ns + 1):(2*ns)];
    improvement_mat[:,3] = improvement[((3-1)*ns + 1):(3*ns)];
    coeff_next = coeff - improvement_mat;
    conv = maximum(abs.(coeff_next - coeff));
    #println("Newton iteration convergence: ", conv)
    iterate1 = iterate1 + 1;
    coeff = coeff_next;
    if iterate1 ==20
        println("Converge not achieved after ", iterate1, "Newton iterations")
        if conv>10^(-4)
            # Only stop if the imprecision is very high
            exitflag_tmp = 4;
        end
        conv =0;
    end
    return (coeff,conv,iterate1,exitflag_tmp)
end
function Q_transition(coeff::Array{Float64,2},
    ns_fine::Int64,ones_tmp_fine::Array{Float64,1},
    size_iid_cost_val::Int64,iid_cost_prob::Array{Float64},income_mat_fine::Array{Float64,4},
    P_kron_fine::SparseMatrixCSC{Float64,Int64},
    a_min::Float64,Phi_z_fine::SparseMatrixCSC{Float64,Int64},β_loc::Float64,
    fspace_a::Dict{Symbol,Any},fspace_a_fine::Dict{Symbol,Any},
    price_final_current::Float64,max_x_fine::Array{Float64,4},
    P_kron1::SparseMatrixCSC{Float64,Int64},
    a_prime_fine_local::Array{Float64},
    future_occupation_fine_local::Array{Float64},
    cons_fine_local::Array{Float64},Q_trans::SparseMatrixCSC{Float64,Int64},
    tol::Float64 = 1e-10)
    #tol= 1e-8
    x_tmp = zeros(ns_fine,3,size_iid_cost_val);
    V_tmp = zeros(ns_fine,3,size_iid_cost_val);
    cons_tmp = zeros(ns_fine,3,size_iid_cost_val);
    Phi_prime_fine = Array{SparseMatrixCSC{Float64,Int64}}(undef, 3,size_iid_cost_val);
    Q_tmp_block = Array{SparseMatrixCSC{Float64,Int64}}(undef, 3,1);
    for j = 1:3
        for j_cost = 1:size_iid_cost_val
            for jj =1:3
                max_x_tmp = max_x_fine[:,jj,j,j_cost]
                xtmp_loc,V_tmp[:,jj,j_cost] = goldenx(BellmanF_stable,a_min*ones_tmp_fine,max_x_tmp,tol,coeff,
                        Phi_z_fine,fspace_a,income_mat_fine[:,jj,j,j_cost],β_loc,jj,price_final_current);
                #(V_tmp[:,jj,j_cost])= BellmanF_stable(x_tmp[:,jj,j_cost],coeff,
                #    Phi_z_fine,fspace_a,income_mat_fine[:,jj,j,j_cost],β_loc,jj,price_final_current);
                Phi_prime_fine_a_tmp = funbase(fspace_a_fine, xtmp_loc);
                #println(typeof(P_kron_fine),typeof(Phi_prime_fine_a_tmp))
                #check_nonzero = Phi_prime_fine_a_tmp.>tol;
                #Phi_prime_fine_a_tmp = check_nonzero.*Phi_prime_fine_a_tmp;
                #println(nnz(Phi_prime_fine_a_tmp))
                #Phi_prime_fine_local =  row_kron(P_kron_fine,Phi_prime_fine_a_tmp);
                #check_nonzero = Phi_prime_fine_local.>tol;
                #Phi_prime_fine[jj,j_cost] = Phi_prime_fine_local.*check_nonzero;
                Phi_prime_fine[jj,j_cost] =  row_kron(P_kron_fine,Phi_prime_fine_a_tmp);
                cons_tmp[:,jj,j_cost] = (income_mat_fine[:,jj,j,j_cost] - xtmp_loc
                    )/price_final_current;
                x_tmp[:,jj,j_cost] = xtmp_loc;
            end
        end
        V_tmp[isnan.(V_tmp)] .= -10000;
        W_tmp, future_occupation = findmax(V_tmp; dims=2);
        future_occupation_index = getindex.(future_occupation,2);
        V_next = iid_cost_prob[1]* W_tmp[:,1,1];
        cons = zeros(ns_fine,1);
        a_prime_local = zeros(ns_fine,1);
        for jj=1:3
            future_occupation_tmp_dummy = sparse(hcat(float(
                        future_occupation_index[:,1,1].==jj)));
            Q_tmp_block[jj,1] = iid_cost_prob[1]  * row_kron(future_occupation_tmp_dummy,Phi_prime_fine[jj,1]);
            cons = cons + future_occupation_tmp_dummy.*cons_tmp[:,jj,1];
            a_prime_local = a_prime_local + future_occupation_tmp_dummy.*x_tmp[:,jj,1];
        end
        cons_fine_local[:,j,1] = cons;
        a_prime_fine_local[:,j,1] = a_prime_local;
        if size_iid_cost_val>1
            for j_cost = 2:size_iid_cost_val
                cons = spzeros(ns_fine,1);
                a_prime_local = spzeros(ns_fine,1);
                V_next = V_next + iid_cost_prob[j_cost]* W_tmp[:,1,j_cost];
                for jj=1:3
                    future_occupation_tmp_dummy = sparse(hcat(float(
                        future_occupation_index[:,1,j_cost].==jj)));
                    Q_tmp_block[jj,1] = Q_tmp_block[jj,1]  + iid_cost_prob[j_cost]  * row_kron(future_occupation_tmp_dummy,
                        Phi_prime_fine[jj,j_cost]);
                    cons = cons + future_occupation_tmp_dummy.*cons_tmp[:,jj,j_cost];
                    a_prime_local = a_prime_local + future_occupation_tmp_dummy.*x_tmp[:,jj,j_cost];
                end
                cons_fine_local[:,j,j_cost] = cons;
                a_prime_fine_local[:,j,j_cost] = a_prime_local;
            end
        end
        future_occupation_fine_local[:,j,:] = future_occupation_index;
        j_prev_ns_fine = (j-1)*ns_fine+1;
        j_ns_fine = j*ns_fine;
        Q_trans[j_prev_ns_fine:j_ns_fine,1:ns_fine] = Q_tmp_block[1,1];
        Q_trans[j_prev_ns_fine:j_ns_fine,(ns_fine+1):(2*ns_fine)]= Q_tmp_block[2,1];
        Q_trans[j_prev_ns_fine:j_ns_fine,(2*ns_fine+1):(3*ns_fine)] = Q_tmp_block[3,1];
    end
    check_nonzero = Q_trans.>tol;
    Q_trans = Q_trans.*check_nonzero;
    I1,J1,V1 = findnz(Q_trans);
    I2 = [I1;3 * ns_fine];
    J2 = [J1;3 * ns_fine];
    V2 = [V1;0.0];
    Q_trans_prime = sparse(J2,I2,V2);
    return Q_trans_prime,cons_fine_local,a_prime_fine_local,future_occupation_fine_local
end
function stationary_distribution(Q_trans_prime::SparseMatrixCSC{Float64,Int64},
    L_loc::Float64,ns_fine::Int64,exitflag_tmp::Int64,tol::Float64 = 1e-10)
    try
        lambda,V,eigflag = eigs(Q_trans_prime,nev=60, v0 = ones(3*ns_fine)./(3*ns_fine), which=:LR,maxiter = 1000);
        if (real(lambda[1] - lambda[2]) < tol)
            println("Unstable solution: second dominant eigenvalue is close to 1")
        end
        stat_distr = real(V[:,1])/sum(real(V[:,1]));
        stat_distr = stat_distr .* (stat_distr.>1e-12);
        stat_distr = stat_distr * L_loc;
        return stat_distr,exitflag_tmp
    catch y
        #warn("Exception: ", y) # What to do on error.
        #if String(string(y))==".ARPACKException(1)"
        if isa(y, Arpack.ARPACKException)
            println("Other issues with the stationary distribution")
            exitflag_tmp = 5;
            return ones(3*ns_fine)./(3*ns_fine),exitflag_tmp
        end
    end
end
function predict_irreducibility(future_occupation_fine_local::Array{Float64,3},
    exitflag_tmp::Int64)
    from_worker_trans_up = maximum(future_occupation_fine_local[:,1,1]);
    from_domestic_trans_up = maximum(future_occupation_fine_local[:,2,1]);
    from_domestic_trans_down = minimum(future_occupation_fine_local[:,2,1]);
    from_exporter_trans_down = minimum(future_occupation_fine_local[:,3,1]);
    if (from_worker_trans_up<3 && from_domestic_trans_up<3)
        println("No exporters in the economy")
        exitflag_tmp = exitflag_tmp + 1;
    elseif from_worker_trans_up<2
        println("No firms in the economy")
        exitflag_tmp = exitflag_tmp + 2;
    elseif from_domestic_trans_down>1 && from_exporter_trans_down>1
        println("No workers in the economy")
            # Only stop if the imprecision is very high
        exitflag_tmp = exitflag_tmp + 3;
    else
        # Final check: unstructured irredcucible
#        Q_trans = transpose(Q_trans_prime); # tradeoff between additional
#        Q_trans1 = Q_trans./sum(Q_trans,dims =2);
#        mc = MarkovChain(Q_trans1)
#        if QuantEcon.is_irreducible(mc)==true
        exitflag_tmp = 0;
#        else
#            exitflag_tmp = exitflag_tmp + 4;
#        end
    end
    return exitflag_tmp
end
function country_residual(country::Int64,s_cell::Array{Array{Float64,2},1},
    ns_cell::Array{Int64,1},s_fine_cell::Array{Array{Float64,2},1},
    ns_fine_cell::Array{Int64,1},
    Phi_z_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    Phi_z_fine_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    Phi_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    Phi_aug_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    P_kron_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    P_kron1_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    P_kron_fine_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    σ::Float64,size_iid_cost_val::Int64,price_final::Array{Float64,2},
    W::Array{Float64,1}, r::Array{Float64,1}, output_final::Array{Float64,1},
    θ::Array{Float64,1}, L::Array{Float64,1}, τ::Array{Float64,1},
    country_no::Int64,banking_cost::Array{Float64,1},δ::Float64,ω::Float64,
    β::Array{Float64,1},
    α₁_eff::Float64,α₂_eff::Float64,α₁::Float64,α₂::Float64,FM::Array{Float64,1},FMₓ::Array{Float64,1},
    F::Array{Float64,1}, Fₓ::Array{Float64,1},
    iid_cost_value::Array{Float64,1},iid_cost_prob::Array{Float64,1},Exit::Array{Float64,1},
    Exitₓ::Array{Float64,1},a_min::Float64,a_max::Float64,α::Float64,fspace_a::Dict{Symbol,Any},
    fspace_a_fine::Dict{Symbol,Any},openness::Int64,agrid_fine::Array{Float64,1})
    # Initialize each country # Threads.@threads for country = 1:country_no
    (exitflag_tmp,s,ns,s_fine,ns_fine,Phi_z,Phi_z_fine,Phi,Phi_aug,
    P_kron,P_kron1,P_kron_fine,z_tilde,z_tilde_fine,income_mat,income_mat_fine,
    a_prime_fine_local,future_occupation_fine_local,cons_fine_local,coeff,
    coeff_next,price_final_prev,price_final_current,total_foreign_demand,W_loc,
    θ_loc,L_loc,τ_loc,r_loc,banking_cost_loc,avg_foreign_price,R,constant_d,
    constant_x,constant_z,constant_z_bar,constant_z_fine,constant_z_bar_fine,
    ones_tmp,x_tmp,V_tmp,β_loc,V_next_stacked,iterate1,Phi_prime_tmp,
    D_deriv_tmp_block,conv,Q_trans,ones_tmp_fine) = country_local_initialize(country,s_cell,ns_cell,s_fine_cell,ns_fine_cell,
    Phi_z_cell,Phi_z_fine_cell,Phi_cell,Phi_aug_cell,P_kron_cell,P_kron1_cell,
    P_kron_fine_cell,σ,size_iid_cost_val,price_final,W, r, output_final,θ,L,τ,
    country_no,banking_cost,δ,ω,β);

    (Profits_d,Profits_x,output_dx,output_xx,output_d,labor_d,labor_x,k_choice_d,k_choice_x,price_dx,price_d,
        price_xx,lambdda_d,lambdda_x,rev_d,rev_dx,rev_xx,k_d_opt,k_x_opt ) = true_profit(price_final_prev,R,W_loc,s,
        constant_x,constant_d,constant_z,constant_z_bar,θ_loc,τ_loc,σ,α₁_eff,α₂_eff,α₁,α₂,ones_tmp);
    (income_mat,max_x) = income_creator(W_loc,ones_tmp,Profits_d,Profits_x,FM,FMₓ,country,
        r_loc,s,iid_cost_value,F, Fₓ,size_iid_cost_val,income_mat,ns,a_min,a_max);
    for i = 1:3
        (coeff[:], conv) = Bellman_iteration(coeff,coeff_next,ns,x_tmp,V_tmp,ones_tmp,size_iid_cost_val,
            iid_cost_prob,income_mat,P_kron,Phi,a_min,Phi_z,β_loc,
            fspace_a,price_final_current,max_x);
    end
    while conv> 10^(-12)
        (coeff[:], conv,iterate1,exitflag_tmp) = Newton_iteration(coeff,coeff_next,ns,x_tmp,V_tmp,ones_tmp,size_iid_cost_val,
            iid_cost_prob,income_mat,P_kron,Phi,a_min,Phi_z,β_loc,
            fspace_a,price_final_current,max_x,V_next_stacked,iterate1,
            Phi_prime_tmp,D_deriv_tmp_block,Phi_aug,P_kron1,exitflag_tmp);
    end
    # Initialize return values if something goes wrong
    export_price_sum_tmp=0.0;
    domestic_price_sum_tmp=0.0;
    NFA_tmp=0.0;
    asset_supply_tmp=0.0;
    asset_demand_tmp=0.0;
    # Solving for the aggregates - stationary distribution:
    if exitflag_tmp>0
        println("Nonstationarity")
    else
        #Continue only if we have convergence in VFI
        (Profits_d_fine,Profits_x_fine,output_dx_fine,output_xx_fine,output_d_fine,labor_d_fine,labor_x_fine
            ,k_choice_d_fine,k_choice_x_fine,price_dx_fine,price_d_fine,price_xx_fine,lambdda_d_fine,
            lambdda_x_fine,rev_d_fine,rev_dx_fine,rev_xx_fine,k_opt_d_fine,k_opt_x_fine) = true_profit(price_final_prev,R,W_loc,s_fine
            ,constant_x,constant_d,constant_z_fine,constant_z_bar_fine,θ_loc,τ_loc,σ,α₁_eff,α₂_eff,α₁,α₂,ones_tmp_fine);
        (income_mat_fine,max_x_fine) = income_creator(W_loc,ones_tmp_fine,Profits_d_fine,Profits_x_fine,FM,FMₓ,country,
                r_loc,s_fine,iid_cost_value,F, Fₓ,size_iid_cost_val,income_mat_fine,ns_fine,a_min,a_max);
        (Q_trans_prime,cons_fine_local,a_prime_fine_local,
        future_occupation_fine_local) = Q_transition(coeff,
            ns_fine,ones_tmp_fine,size_iid_cost_val,iid_cost_prob,income_mat_fine,
            P_kron_fine,a_min,Phi_z_fine,β_loc,fspace_a,fspace_a_fine,
            price_final_current,max_x_fine,P_kron1,cons_fine_local,a_prime_fine_local,
            future_occupation_fine_local,Q_trans);
        exitflag_tmp = predict_irreducibility(future_occupation_fine_local,exitflag_tmp);
    end
    if exitflag_tmp==1 || exitflag_tmp==2
        labor_excess_demand_tmp_percent = -200.0 ;
        excess_demand_final_good_tmp_percent = 200.0
    elseif exitflag_tmp==3
        labor_excess_demand_tmp_percent = 200.0 ;
        excess_demand_final_good_tmp_percent = -200.0
    elseif exitflag_tmp==4
        labor_excess_demand_tmp_percent = 2000.0 ;
        excess_demand_final_good_tmp_percent = 2000.0
    else
        distr_current,exitflag_tmp = stationary_distribution(Q_trans_prime,L_loc,ns_fine,exitflag_tmp);
        if exitflag_tmp==5
            labor_excess_demand_tmp_percent = 1000.0 ;
            excess_demand_final_good_tmp_percent = 1000.0
        else
            current_distr_store_tmp = copy(distr_current);
            # Calculate aggregates
            # Distribution of current workers:
            worker_past_dist = distr_current[(ns_fine *0 + 1):(ns_fine *1)];
            domestic_past_dist = distr_current[(ns_fine *1 + 1):(ns_fine *2)];
            exporter_past_dist = distr_current[(ns_fine *2 + 1):(ns_fine *3)];

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
            investment_tmp = δ * capital_demand_tmp;
            total_demand_final_good_tmp = total_consumption_tmp  + investment_tmp;
            excess_demand_final_good_tmp = total_demand_final_good_tmp - output_final[country,1];
            labor_excess_demand_tmp_percent = labor_excess_demand_tmp./L_loc;
            excess_demand_final_good_tmp_percent = excess_demand_final_good_tmp/(output_final[country,1] + total_demand_final_good_tmp);
        end
    end
    return (labor_excess_demand_tmp_percent,excess_demand_final_good_tmp_percent,
    exitflag_tmp,export_price_sum_tmp, domestic_price_sum_tmp,NFA_tmp,asset_supply_tmp,asset_demand_tmp)
end
function country_residual_detailed(country::Int64,s_cell::Array{Array{Float64,2},1},
    ns_cell::Array{Int64,1},s_fine_cell::Array{Array{Float64,2},1},
    ns_fine_cell::Array{Int64,1},
    Phi_z_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    Phi_z_fine_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    Phi_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    Phi_aug_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    P_kron_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    P_kron1_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    P_kron_fine_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    σ::Float64,size_iid_cost_val::Int64,price_final::Array{Float64,2},
    W::Array{Float64,1}, r::Array{Float64,1}, output_final::Array{Float64,1},
    θ::Array{Float64,1}, L::Array{Float64,1}, τ::Array{Float64,1},
    country_no::Int64,banking_cost::Array{Float64,1},δ::Float64,ω::Float64,
    β::Array{Float64,1},
    α₁_eff::Float64,α₂_eff::Float64,α₁::Float64,α₂::Float64,FM::Array{Float64,1},FMₓ::Array{Float64,1},
    F::Array{Float64,1}, Fₓ::Array{Float64,1},
    iid_cost_value::Array{Float64,1},iid_cost_prob::Array{Float64,1},Exit::Array{Float64,1},
    Exitₓ::Array{Float64,1},a_min::Float64,a_max::Float64,α::Float64,fspace_a::Dict{Symbol,Any},
    fspace_a_fine::Dict{Symbol,Any},openness::Int64,agrid_fine::Array{Float64,1})
    # Initialize each country # Threads.@threads for country = 1:country_no
    (exitflag_tmp,s,ns,s_fine,ns_fine,Phi_z,Phi_z_fine,Phi,Phi_aug,
    P_kron,P_kron1,P_kron_fine,z_tilde,z_tilde_fine,income_mat,income_mat_fine,
    a_prime_fine_local,future_occupation_fine_local,cons_fine_local,coeff,
    coeff_next,price_final_prev,price_final_current,total_foreign_demand,W_loc,
    θ_loc,L_loc,τ_loc,r_loc,banking_cost_loc,avg_foreign_price,R,constant_d,
    constant_x,constant_z,constant_z_bar,constant_z_fine,constant_z_bar_fine,
    ones_tmp,x_tmp,V_tmp,β_loc,V_next_stacked,iterate1,Phi_prime_tmp,
    D_deriv_tmp_block,conv,Q_trans,ones_tmp_fine) = country_local_initialize(country,s_cell,ns_cell,s_fine_cell,ns_fine_cell,
    Phi_z_cell,Phi_z_fine_cell,Phi_cell,Phi_aug_cell,P_kron_cell,P_kron1_cell,
    P_kron_fine_cell,σ,size_iid_cost_val,price_final,W, r, output_final,θ,L,τ,
    country_no,banking_cost,δ,ω,β);

    (Profits_d,Profits_x,output_dx,output_xx,output_d,labor_d,labor_x,k_choice_d,k_choice_x,price_dx,price_d,
        price_xx,lambdda_d,lambdda_x,rev_d,rev_dx,rev_xx,k_opt_d,k_opt_x ) = true_profit(price_final_prev,R,W_loc,s,
        constant_x,constant_d,constant_z,constant_z_bar,θ_loc,τ_loc,σ,α₁_eff,α₂_eff,α₁,α₂,ones_tmp);
    (income_mat,max_x) = income_creator(W_loc,ones_tmp,Profits_d,Profits_x,FM,FMₓ,country,
        r_loc,s,iid_cost_value,F, Fₓ,size_iid_cost_val,income_mat,ns,a_min,a_max);
    for i = 1:3
        (coeff[:], conv) = Bellman_iteration(coeff,coeff_next,ns,x_tmp,V_tmp,ones_tmp,size_iid_cost_val,
            iid_cost_prob,income_mat,P_kron,Phi,a_min,Phi_z,β_loc,
            fspace_a,price_final_current,max_x);
    end
    while conv> 10^(-12)
        (coeff[:], conv,iterate1,exitflag_tmp) = Newton_iteration(coeff,coeff_next,ns,x_tmp,V_tmp,ones_tmp,size_iid_cost_val,
            iid_cost_prob,income_mat,P_kron,Phi,a_min,Phi_z,β_loc,
            fspace_a,price_final_current,max_x,V_next_stacked,iterate1,
            Phi_prime_tmp,D_deriv_tmp_block,Phi_aug,P_kron1,exitflag_tmp);
    end
    # Initialize return values if something goes wrong
    export_price_sum_tmp=0.0;domestic_price_sum_tmp=0.0;labor_excess_demand_tmp_percent=0.0;
    excess_demand_final_good_tmp_percent=0.0; NFA_tmp=0.0;asset_supply_tmp=0.0;asset_demand_tmp=0.0;
    current_distr_store_tmp=zeros(3*ns_fine);distr_current=zeros(3*ns_fine);domestic_prod_tmp=0.0;
    export_prod_tmp=0.0;export_value_tmp=0.0;exit_share_to_domestic_tmp=0.0;entry_share_to_domestic_tmp=0.0;
    exit_share_to_worker_tmp=0.0;entry_share_to_exporter_tmp=0.0;exporter_pop_tmp=0.0;domestic_pop_tmp=0.0;
    worker_pop_tmp=0.0;L_x_tmp=0.0;L_d_tmp=0.0;total_entry_cost_tmp=0.0;total_incumbents_cost_tmp=0.0;total_exit_cost_tmp=0.0;
    labor_demand_tmp=0.0;labor_excess_demand_tmp=0.0;K_x_tmp=0.0;K_d_tmp=0.0;capital_demand_tmp=0.0;worker_bond_holding_sum_tmp=0.0;
    domestic_bond_holding_sum_tmp=0.0;exporter_bond_holding_sum_tmp=0.0;domestic_firm_debt_tmp=0.0;exporter_firm_debt_tmp=0.0;
    total_consumption_tmp=0.0;Aggr_bank_cost_tmp=0.0;x_k_tmp=0.0;x_l_tmp=0.0;investment_tmp=0.0;total_demand_final_good_tmp=0.0;
    excess_demand_final_good_tmp=0.0;mean_MRPK_d_tmp=0.0;mean_MRPK_x_tmp=0.0;mean_MRPK_tmp=0.0;sd_MRPK_d_tmp=0.0; sd_MRPK_x_tmp=0.0;sd_MRPK_tmp=0.0;
    mean_logProd_d_tmp=0.0;mean_logProd_x_tmp=0.0;mean_logProd_tmp=0.0;sd_logProd_d_tmp=0.0;sd_logProd_x_tmp=0.0;sd_logProd_tmp=0.0;
    cov_TFPR_z_d_tmp=0.0;cov_TFPR_z_x_tmp=0.0;cov_TFPR_z_tmp=0.0;corr_TFPR_z_d_tmp=0.0;corr_TFPR_z_x_tmp=0.0;corr_TFPR_z_tmp=0.0;
    avg_Pi_x_tmp=0.0;avg_Pi_d_tmp=0.0;D_d_denom_tmp=0.0;D_x_denom_tmp=0.0;D_x_tmp=0.0;D_d_tmp=0.0;D_d_denom_eff_tmp=0.0;D_x_denom_eff_tmp=0.0;
    D_x_eff_tmp=0.0;D_d_eff_tmp=0.0;mean_cons_tmp=0.0;sd_cons_tmp=0.0;mean_income_tmp=0.0;sd_income_tmp=0.0;mean_wealth_tmp=0.0;sd_wealth_tmp=0.0;
    p10_wealth_tmp=0.0;p50_wealth_tmp=0.0;p90_wealth_tmp=0.0;p99_wealth_tmp=0.0;p10_income_tmp=0.0;p50_income_tmp=0.0;p90_income_tmp=0.0;p99_income_tmp=0.0;p10_cons_tmp=0.0;
    p50_cons_tmp=0.0;p90_cons_tmp=0.0;p99_cons_tmp=0.0;wealth_of_workers_tmp=0.0;wealth_of_domestic_tmp=0.0;
    wealth_of_exporters_tmp=0.0;income_of_workers_tmp=0.0;income_of_domestic_tmp=0.0;
    income_of_exporters_tmp=0.0;cons_of_workers_tmp=0.0;cons_of_domestic_tmp=0.0;
    cons_of_exporters_tmp=0.0;mean_leverage_d_tmp=0.0;mean_leverage_x_tmp=0.0;mean_leverage_tmp=0.0;
    sd_leverage_d_tmp=0.0;sd_leverage_x_tmp=0.0;sd_leverage_tmp=0.0;corr_MRPK_lev_d_tmp=0.0;corr_MRPK_lev_x_tmp=0.0;corr_MRPK_lev_tmp=0.0;
    exit_domestic_sum_tmp=0.0;exit_exporting_to_work_sum_tmp=0.0;mean_growth_rev=0.0;sd_growth_rev=0.0;mean_growth_k=0.0;sd_growth_k=0.0;autocorr_rev=0.0;fraction_zombie_exporter_tmp=0.0;
    # Solving for the aggregates - stationary distribution:
    if exitflag_tmp>0
        println("Nonstationarity")
    else
        #Continue only if we have convergence in VFI
        (Profits_d_fine,Profits_x_fine,output_dx_fine,output_xx_fine,output_d_fine,labor_d_fine,labor_x_fine
            ,k_choice_d_fine,k_choice_x_fine,price_dx_fine,price_d_fine,price_xx_fine,lambdda_d_fine,
            lambdda_x_fine,rev_d_fine,rev_dx_fine,rev_xx_fine,k_opt_d_fine,k_opt_x_fine) = true_profit(price_final_prev,R,W_loc,s_fine
            ,constant_x,constant_d,constant_z_fine,constant_z_bar_fine,θ_loc,τ_loc,σ,α₁_eff,α₂_eff,α₁,α₂,ones_tmp_fine);
        (income_mat_fine,max_x_fine) = income_creator(W_loc,ones_tmp_fine,Profits_d_fine,Profits_x_fine,FM,FMₓ,country,
                r_loc,s_fine,iid_cost_value,F, Fₓ,size_iid_cost_val,income_mat_fine,ns_fine,a_min,a_max);
        (Q_trans_prime,cons_fine_local,a_prime_fine_local,
        future_occupation_fine_local) = Q_transition(coeff,
            ns_fine,ones_tmp_fine,size_iid_cost_val,iid_cost_prob,income_mat_fine,
            P_kron_fine,a_min,Phi_z_fine,β_loc,fspace_a,fspace_a_fine,
            price_final_current,max_x_fine,P_kron1,cons_fine_local,a_prime_fine_local,
            future_occupation_fine_local,Q_trans);
        exitflag_tmp = predict_irreducibility(future_occupation_fine_local,exitflag_tmp);
    end
    if exitflag_tmp==1 || exitflag_tmp==2
        labor_excess_demand_tmp_percent = -200.0 ;
        excess_demand_final_good_tmp_percent = 200.0
    elseif exitflag_tmp==3
        labor_excess_demand_tmp_percent = 200.0 ;
        excess_demand_final_good_tmp_percent = -200.0
    elseif exitflag_tmp==4
        labor_excess_demand_tmp_percent = 2000.0 ;
        excess_demand_final_good_tmp_percent = 2000.0
    else
        distr_current,exitflag_tmp = stationary_distribution(Q_trans_prime,L_loc,ns_fine,exitflag_tmp);
        if exitflag_tmp==5
            labor_excess_demand_tmp_percent = 1000.0 ;
            excess_demand_final_good_tmp_percent = 1000.0
        else
            current_distr_store_tmp = zeros(size(distr_current));
            # Calculate aggregates
            # Distribution of current workers:
            worker_past_dist = distr_current[(ns_fine *0 + 1):(ns_fine *1)];
            domestic_past_dist = distr_current[(ns_fine *1 + 1):(ns_fine *2)];
            exporter_past_dist = distr_current[(ns_fine *2 + 1):(ns_fine *3)];
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
            investment_tmp = δ * capital_demand_tmp;
            total_demand_final_good_tmp = total_consumption_tmp  + investment_tmp;
            excess_demand_final_good_tmp = total_demand_final_good_tmp - output_final[country,1];
            labor_excess_demand_tmp_percent = labor_excess_demand_tmp./L_loc;
            excess_demand_final_good_tmp_percent = excess_demand_final_good_tmp/(output_final[country,1] + total_demand_final_good_tmp);
            # Check if the upper limit has been reached:
            #upper_lim = check_upper_limit(worker_past_dist,domestic_past_dist,exporter_past_dist,
            #ns_fine,agrid_fine);

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

        # New stuff. Leverage with TFPR, as both are measurable
            mean_leverage_d_tmp = sum((current_domestic.*((k_choice_d_fine - s_fine[:,1])./k_choice_d_fine) )/domestic_pop_tmp);
            mean_leverage_x_tmp = sum((current_exporter.*((k_choice_x_fine - s_fine[:,1])./k_choice_x_fine) )/exporter_pop_tmp);
            mean_leverage_tmp =(domestic_pop_tmp* mean_leverage_d_tmp +
             exporter_pop_tmp *mean_leverage_x_tmp)/(
             domestic_pop_tmp + exporter_pop_tmp) ;

            sd_leverage_d_tmp =    (sum(current_domestic.*((k_choice_d_fine - s_fine[:,1])./k_choice_d_fine).^2)/domestic_pop_tmp - mean_leverage_d_tmp^2)^(1/2);
            sd_leverage_x_tmp =     (sum(current_exporter.*((k_choice_x_fine - s_fine[:,1])./k_choice_x_fine).^2)/exporter_pop_tmp - mean_leverage_x_tmp^2)^(1/2);
            sd_leverage_tmp =  ((sum((current_domestic.*((k_choice_d_fine - s_fine[:,1])./k_choice_d_fine).^2)) +
            sum(current_exporter.*((k_choice_x_fine - s_fine[:,1])./k_choice_x_fine).^2))/(domestic_pop_tmp + exporter_pop_tmp) - mean_leverage_tmp^2)^(1/2);

            cov_MRPK_lev_d_tmp = @fastmath sum(current_domestic.*(((k_choice_d_fine - s_fine[:,1])./k_choice_d_fine).*log.((lambdda_d_fine +
            R*ones_tmp_fine))))/domestic_pop_tmp- mean_MRPK_d_tmp * mean_leverage_d_tmp;
            cov_MRPK_lev_x_tmp = @fastmath sum(current_exporter.*(((k_choice_x_fine - s_fine[:,1])./k_choice_x_fine).*log.((lambdda_x_fine +
            R*ones_tmp_fine))))/exporter_pop_tmp- mean_MRPK_x_tmp * mean_leverage_x_tmp;
            cov_MRPK_lev_tmp = @fastmath (sum(current_exporter.*(((k_choice_x_fine - s_fine[:,1])./k_choice_x_fine).*log.((lambdda_x_fine +
            R*ones_tmp_fine)))) + sum(current_domestic.*(((k_choice_d_fine - s_fine[:,1])./k_choice_d_fine).*log.((lambdda_d_fine +
            R*ones_tmp_fine)))))/(domestic_pop_tmp + exporter_pop_tmp) - mean_MRPK_tmp * mean_leverage_tmp;

            corr_MRPK_lev_d_tmp = cov_MRPK_lev_d_tmp / (sd_MRPK_d_tmp* sd_leverage_d_tmp );
            corr_MRPK_lev_x_tmp = cov_MRPK_lev_x_tmp / (sd_MRPK_x_tmp* sd_leverage_x_tmp );
            corr_MRPK_lev_tmp = cov_MRPK_lev_tmp / (sd_MRPK_tmp* sd_leverage_tmp );

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
            current_distr_cons_sorted = distr_current[cons_sort_index]/L_loc;
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

            #Surviving firms
            sum_prob_survive = sum(Q_trans[(ns_fine+1):end,(ns_fine+1):end],dims = 2);
            Q_trans_firm = Q_trans[(ns_fine+1):end,(ns_fine+1):end]
            surviving_index = findall(x->x>0,sum_prob_survive);
            surviving_index = first.(Tuple.(surviving_index));
            sum_prob_survive_nonzero = sum_prob_survive[surviving_index];
            Q_trans_survive_firm= Q_trans_firm[surviving_index,:]./sum_prob_survive_nonzero;
            dropzeros!(Q_trans_survive_firm);
            current_firm_distr = current_distr_store_tmp[(ns_fine *1 + 1):(ns_fine *3)] ;
            # Normalize to 1 on surviving firms
            current_firm_distr_survive = current_firm_distr[surviving_index];
            current_firm_distr_survive = current_firm_distr_survive./sum(current_firm_distr_survive);
            current_firm_distr_survive_rep = repeat(current_firm_distr_survive,1,2*ns_fine);
            distribution_of_firms = current_firm_distr_survive_rep.*Q_trans_survive_firm;
            # growth rate of revenue:
            do_do_growth = zeros(ns_fine,ns_fine);
            do_ex_growth = zeros(ns_fine,ns_fine);
            ex_do_growth = zeros(ns_fine,ns_fine);
            ex_ex_growth = zeros(ns_fine,ns_fine);

            for ii= 1:ns_fine
                do_do_growth[ii,:] = log.(rev_d_fine) - log.(rev_d_fine[ii] *ones_tmp_fine) ;
                do_ex_growth[ii,:] = log.(rev_dx_fine+rev_xx_fine) - log.(rev_d_fine[ii] *ones_tmp_fine) ;
                ex_do_growth[ii,:] = log.(rev_d_fine) - log.( (rev_dx_fine[ii]+rev_xx_fine[ii]) *ones_tmp_fine);
                ex_ex_growth[ii,:] = log.(rev_dx_fine+rev_xx_fine) -  log.((rev_dx_fine[ii]+rev_xx_fine[ii]) *ones_tmp_fine);
            end

            growth_rate_rev_mat = zeros(2*ns_fine,2*ns_fine);
            growth_rate_rev_mat[1:ns_fine,1:ns_fine] = do_do_growth;
            growth_rate_rev_mat[ns_fine+1:end,1:ns_fine] = ex_do_growth;
            growth_rate_rev_mat[1:ns_fine,ns_fine+1:end] = do_ex_growth;
            growth_rate_rev_mat[ns_fine+1:end,ns_fine+1:end] = ex_ex_growth;
            growth_rate_rev_mat_survive = growth_rate_rev_mat[surviving_index,:];
            mean_growth_rev = sum(growth_rate_rev_mat_survive.*distribution_of_firms);
            sd_growth_rev = (sum(growth_rate_rev_mat_survive.^2 .*distribution_of_firms) - mean_growth_rev^2)^(1/2);
            # growth rate of capital:
            do_do_growth = zeros(ns_fine,ns_fine);
            do_ex_growth = zeros(ns_fine,ns_fine);
            ex_do_growth = zeros(ns_fine,ns_fine);
            ex_ex_growth = zeros(ns_fine,ns_fine);

            for ii= 1:ns_fine
                do_do_growth[ii,:] = log.(k_choice_d_fine) - log.(k_choice_d_fine[ii] *ones_tmp_fine);
                do_ex_growth[ii,:] = log.(k_choice_x_fine) - log.(k_choice_d_fine[ii] *ones_tmp_fine);
                ex_do_growth[ii,:] = log.(k_choice_d_fine) - log.(k_choice_x_fine[ii] *ones_tmp_fine);
                ex_ex_growth[ii,:] = log.(k_choice_x_fine) -  log.(k_choice_x_fine[ii]*ones_tmp_fine);
            end

            growth_rate_k_mat = zeros(2*ns_fine,2*ns_fine);
            growth_rate_k_mat[1:ns_fine,1:ns_fine] = do_do_growth;
            growth_rate_k_mat[ns_fine+1:end,1:ns_fine] = ex_do_growth;
            growth_rate_k_mat[1:ns_fine,ns_fine+1:end] = do_ex_growth;
            growth_rate_k_mat[ns_fine+1:end,ns_fine+1:end] = ex_ex_growth;

            growth_rate_k_mat_survive = growth_rate_k_mat[surviving_index,:];
            mean_growth_k = sum(growth_rate_k_mat_survive.*distribution_of_firms);
            sd_growth_k = (sum(growth_rate_k_mat_survive.^2 .*distribution_of_firms) - mean_growth_k^2)^(1/2);

            # autocorr of revenue:
            do_do_mul = zeros(ns_fine,ns_fine);
            do_ex_mul = zeros(ns_fine,ns_fine);
            ex_do_mul = zeros(ns_fine,ns_fine);
            ex_ex_mul = zeros(ns_fine,ns_fine);

            for ii= 1:ns_fine
                do_do_mul[ii,:] = log.(rev_d_fine)*log(rev_d_fine[ii]) ;
                do_ex_mul[ii,:] = log.(rev_dx_fine+rev_xx_fine)*log(rev_d_fine[ii]) ;
                ex_do_mul[ii,:] = log.(rev_d_fine)*log(rev_dx_fine[ii]+rev_xx_fine[ii]);
                ex_ex_mul[ii,:] = log.(rev_dx_fine+rev_xx_fine)*log(rev_dx_fine[ii]+rev_xx_fine[ii]);
            end

            mul_rate_rev_mat = zeros(2*ns_fine,2*ns_fine);
            mul_rate_rev_mat[1:ns_fine,1:ns_fine] = do_do_mul;
            mul_rate_rev_mat[ns_fine+1:end,1:ns_fine] = ex_do_mul;
            mul_rate_rev_mat[1:ns_fine,ns_fine+1:end] = do_ex_mul;
            mul_rate_rev_mat[ns_fine+1:end,ns_fine+1:end] = ex_ex_mul;
            mul_rate_rev_mat_survive = mul_rate_rev_mat[surviving_index,:];
            mean_mul_rev = sum(mul_rate_rev_mat_survive.*distribution_of_firms);
            revenue_combined = zeros(2*ns_fine);
            revenue_combined[1:ns_fine] = log.(rev_d_fine);
            revenue_combined[ns_fine+1:end] = log.(rev_dx_fine+rev_xx_fine);
            distr_tmr_firm = transpose(sum(distribution_of_firms,dims = 1));
            mean_log_rev = sum(revenue_combined.*distr_tmr_firm);
            sd_log_rev = (sum(revenue_combined.^2 .*distr_tmr_firm) - mean_log_rev.^2)^(1/2);
            revenue_combined_survive = revenue_combined[surviving_index,:];
            mean_log_rev_survive = sum(revenue_combined_survive .*current_firm_distr_survive);
            sd_log_rev_survive = (sum(revenue_combined_survive.^2 .*current_firm_distr_survive) - mean_log_rev_survive.^2)^(1/2);
            autocorr_rev = (mean_mul_rev - mean_log_rev_survive*mean_log_rev)/sd_log_rev_survive/sd_log_rev;
            # Zombies - experimental
            making_losses = (income_1_fine .>income_3_fine);
            exit_zombie = sum(((future_occupation_fine_local[:,3,1] .!=3.0) .*making_losses).*current_exporter)/exporter_pop_tmp;
            fraction_zombie_exporter_tmp = sum(making_losses.*current_exporter)/exporter_pop_tmp - 3 * exit_zombie
        end
    end

    return (labor_excess_demand_tmp_percent,excess_demand_final_good_tmp_percent,
    exitflag_tmp,export_price_sum_tmp, domestic_price_sum_tmp,NFA_tmp,asset_supply_tmp,
    asset_demand_tmp,a_prime_fine_local,future_occupation_fine_local,cons_fine_local,
    k_choice_x_fine,labor_x_fine,k_choice_d_fine,labor_d_fine,Profits_x_fine,Profits_d_fine,
    rev_d_fine,rev_dx_fine,rev_xx_fine,current_distr_store_tmp,distr_current,coeff,domestic_prod_tmp,
    export_prod_tmp,export_value_tmp,exit_share_to_domestic_tmp,entry_share_to_domestic_tmp,
    exit_share_to_worker_tmp,entry_share_to_exporter_tmp,exporter_pop_tmp,domestic_pop_tmp,
    worker_pop_tmp,L_x_tmp,L_d_tmp,total_entry_cost_tmp,total_incumbents_cost_tmp,total_exit_cost_tmp,
    labor_demand_tmp,labor_excess_demand_tmp,K_x_tmp,K_d_tmp,capital_demand_tmp,worker_bond_holding_sum_tmp,
    domestic_bond_holding_sum_tmp,exporter_bond_holding_sum_tmp,domestic_firm_debt_tmp,exporter_firm_debt_tmp,
    total_consumption_tmp,Aggr_bank_cost_tmp,x_k_tmp,x_l_tmp,investment_tmp,total_demand_final_good_tmp,
    excess_demand_final_good_tmp,mean_MRPK_d_tmp,mean_MRPK_x_tmp,mean_MRPK_tmp,sd_MRPK_d_tmp, sd_MRPK_x_tmp,sd_MRPK_tmp,
    mean_logProd_d_tmp,mean_logProd_x_tmp,mean_logProd_tmp,sd_logProd_d_tmp,sd_logProd_x_tmp,sd_logProd_tmp,
    cov_TFPR_z_d_tmp,cov_TFPR_z_x_tmp,cov_TFPR_z_tmp,corr_TFPR_z_d_tmp,corr_TFPR_z_x_tmp,corr_TFPR_z_tmp,
    avg_Pi_x_tmp,avg_Pi_d_tmp,D_d_denom_tmp,D_x_denom_tmp,D_x_tmp,D_d_tmp,D_d_denom_eff_tmp,D_x_denom_eff_tmp,
    D_x_eff_tmp,D_d_eff_tmp,mean_cons_tmp,sd_cons_tmp,mean_income_tmp,sd_income_tmp,mean_wealth_tmp,sd_wealth_tmp,
    p10_wealth_tmp,p50_wealth_tmp,p90_wealth_tmp,p99_wealth_tmp,p10_income_tmp,p50_income_tmp,p90_income_tmp,p99_income_tmp,p10_cons_tmp,
    p50_cons_tmp,p90_cons_tmp,p99_cons_tmp,wealth_of_workers_tmp,wealth_of_domestic_tmp,
    wealth_of_exporters_tmp,income_of_workers_tmp,income_of_domestic_tmp,
    income_of_exporters_tmp,cons_of_workers_tmp,cons_of_domestic_tmp,
    cons_of_exporters_tmp,mean_leverage_d_tmp,mean_leverage_x_tmp,mean_leverage_tmp,
    sd_leverage_d_tmp,sd_leverage_x_tmp,sd_leverage_tmp,corr_MRPK_lev_d_tmp,corr_MRPK_lev_x_tmp,corr_MRPK_lev_tmp,
    exit_domestic_sum_tmp,exit_exporting_to_work_sum_tmp,mean_growth_rev,sd_growth_rev,mean_growth_k,sd_growth_k,autocorr_rev,fraction_zombie_exporter_tmp,k_opt_d_fine,k_opt_x_fine)
end

function check_upper_limit(worker_past_dist::Array{Float64,1},domestic_past_dist::Array{Float64,1},
    exporter_past_dist::Array{Float64,1},ns_fine::Int64,agrid_fine::Array{Float64,1},tol::Float64=1e-10 )
    a_grid_size = size(agrid_fine)[1];
    z_grid_size = convert(Int64,ns_fine/size(agrid_fine)[1]);
    worker_past_dist_mat = reshape(worker_past_dist,a_grid_size,z_grid_size);
    domestic_past_dist_mat = reshape(domestic_past_dist,a_grid_size,z_grid_size);
    exporter_past_dist_mat = reshape(exporter_past_dist,a_grid_size,z_grid_size);
    ii =1;
    upper_lim = 1.0;
    while (upper_lim>tol || ii ==size(agrid_fine)[1])
        upper_lim = (sum(sum(worker_past_dist_mat,dims =2)[ii:a_grid_size])+
        sum(sum(domestic_past_dist_mat,dims =2)[ii:a_grid_size]) +
        sum(sum(exporter_past_dist_mat,dims =2)[ii:a_grid_size]));
        ii = ii +1;
    end
    if ii == size(agrid_fine)[1]
        println("Upper limit reached")
    end
    return agrid_fine[ii]
end
function Residual_stst(prices::Array{Float64,1},parameters_tmp::Parameter_type)
    # Extract the necessary local parameters:
    (β,α,δ,θ,α₁,α₂,σ,α₁_eff,α₂_eff,ω,L,FM,FMₓ,F,Fₓ,
        Exit,Exitₓ,iid_cost_value,iid_cost_prob,size_iid_cost_val,country_no,τ,
        a_min,a_max,fspace_a,fspace_a_fine,agrid_fine,banking_cost,
        bounds,Country_spec_p,s_cell,ns_cell,s_fine_cell,ns_fine_cell,Phi_z_cell,
        Phi_z_fine_cell,Phi_cell,Phi_aug_cell,P_kron_cell,P_kron1_cell,P_kron_fine_cell,
        exp_egrid_cell,ns_tmp,ns_tmp_fine,openness) = local_parameters(parameters_tmp)
    # Prices:
    (W,r,output_final,price_final,r,residual,price_check_tmp) = price_reshaper(prices,openness,country_no,bounds)
    # Check prices to be reasonable:
    if price_check_tmp>0
        println("Guess out of bounds")
        residual = 1000.0 * ones(8);
        return residual
    end
    # Initialization of local variables:
    (export_prod,export_price_sum,import_price,domestic_price_sum,
    price_final_actual,NFA,asset_supply,asset_demand,exitflag) = local_var_creator_minimal(country_no);
    @sync @distributed for country = 1:country_no#@sync @distributed for country = 1:country_no#
        (residual[country,1],residual[2* country_no + country,1],exitflag[country,1],export_price_sum[country,1],
        domestic_price_sum[country,1],NFA[country,1],asset_supply[country,1],
        asset_demand[country,1])  = country_residual(country,s_cell,ns_cell,s_fine_cell,ns_fine_cell,
        Phi_z_cell,Phi_z_fine_cell,Phi_cell,Phi_aug_cell,P_kron_cell,P_kron1_cell,
        P_kron_fine_cell,σ,size_iid_cost_val,price_final,W, r, output_final,θ,L,τ,
        country_no,banking_cost,δ,ω,β,α₁_eff,α₂_eff,α₁,α₂,FM,FMₓ,F,Fₓ,iid_cost_value,
        iid_cost_prob,Exit,Exitₓ,a_min,a_max,α,fspace_a,fspace_a_fine,openness,agrid_fine);
    end
    for country=1:country_no
        if exitflag[country,1]==1 || exitflag[country,1]==2# Too large supply of assets
            if openness == 1
                residual[7,1] =  100.0;
            else
                residual[3* country_no + country,1] = 100.0;
            end
        elseif exitflag[country,1]==3# Too little supply of assets?
            if openness == 1
                residual[7,1] =  - 100.0;
            else
                residual[3* country_no + country,1] = - 100.0;
            end
        elseif exitflag[country,1]==4
            if openness == 1
                residual[7,1] =  1000.0;
            else
                residual[3* country_no + country,1] =  1000.0;
            end
        elseif exitflag[country,1]==5
            if openness == 1
                residual[7,1] = 500.0;
            else
                residual[3* country_no + country,1] =  500.0;
            end
        else
            import_price[country,1] = (sum(export_price_sum[:,1])-  export_price_sum[country,1])/(country_no - 1);
            price_final_actual[country,1] = min((ω^σ * domestic_price_sum[country,1] + (1.0 - ω)^σ  * import_price[country,1])^(1.0
                / (1.0 - σ)),10000);
            residual[country_no + country,1] = (price_final_actual[country,1] - price_final[country,1])./(
                price_final_actual[country,1] + price_final[country,1]);
            if openness == 0
                residual[3* country_no + country,1] = NFA[country,1]/(asset_demand[country,1] + (1 + banking_cost[country,1]).*asset_supply[country,1]);
            end
        end
    end
    if openness == 1
        #residual[7,1] = sum(NFA)/sum((asset_demand+ (ones(country_no) + banking_cost).*asset_supply));
        residual[7,1] = sum(NFA)/output_final[1];
    end
    replace!(residual, NaN=>100.0);
    return residual
end
function Residual_stst_detailed(prices::Array{Float64,1},parameters_tmp::Parameter_type)
    # Extract the necessary local parameters:
    (β,α,δ,θ,α₁,α₂,σ,α₁_eff,α₂_eff,ω,L,FM,FMₓ,F,Fₓ,
        Exit,Exitₓ,iid_cost_value,iid_cost_prob,size_iid_cost_val,country_no,τ,
        a_min,a_max,fspace_a,fspace_a_fine,agrid_fine,banking_cost,
        bounds,Country_spec_p,s_cell,ns_cell,s_fine_cell,ns_fine_cell,Phi_z_cell,
        Phi_z_fine_cell,Phi_cell,Phi_aug_cell,P_kron_cell,P_kron1_cell,P_kron_fine_cell,
        exp_egrid_cell,ns_tmp,ns_tmp_fine,openness) = local_parameters(parameters_tmp)
    # Prices:
    (W,r,output_final,price_final,r,residual,price_check_tmp) = price_reshaper(prices,openness,country_no,bounds)
    # Check prices to be reasonable:
    if price_check_tmp>0
        println("Guess out of bounds")
        residual = 1000.0 * ones(8);
        return residual
    end
    # Initialization of local variables:
    (export_prod,export_price_sum,import_price,domestic_price_sum,
    price_final_actual,NFA,asset_supply,asset_demand,exitflag) = local_var_creator_minimal(country_no);
    (a_prime_fine,future_occupation_fine,cons_fine,rev_d_fine_store,
    rev_dx_fine_store,rev_xx_fine_store,coeff_final,k_x_fine,k_d_fine,
    l_x_fine,l_d_fine,current_distr_store,distr_current,exitflag,domestic_prod,exporter_pop,
    domestic_pop,worker_pop,entry_share_to_domestic,exit_share_to_domestic,
    entry_share_to_exporter,exit_share_to_worker,total_entry_cost,
    total_incumbents_cost,total_exit_cost,labor_excess_demand,labor_demand,
    total_consumption,capital_demand,capital_next,
    total_demand_final_good,excess_demand_final_good,GDP,nomGDP,PPI,TFP,
    TFP_within,TFP_across,TFP_second_best,TOT,total_production,
    domestic_price_sum,export_value,investment,
    export_value,L_x,L_d,K_x,K_d,x_k,x_l,D_d,D_x,D_d_denom,D_x_denom,
    D_d_eff,D_x_eff,D_d_denom_eff,D_x_denom_eff,K_d_ratio,K_x_ratio,
    L_d_ratio,L_x_ratio,K_x_ratio_eff,K_d_ratio_eff,L_d_ratio_eff,
    L_x_ratio_eff,nomGDP_d,nomGDP_xd,nomGDP_xx,RGDP_d,RGDP_xd,RGDP_xx,RGDP,
    mean_MRPK_d,mean_MRPK_x,sd_MRPK_d,sd_MRPK_x,mean_MRPK,sd_MRPK,
    mean_MRPL_d,mean_MRPL_x,sd_MRPL_d,sd_MRPL_x,mean_MRPL,sd_MRPL,
    cov_TFPR_z_d,cov_TFPR_z_x,cov_TFPR_z,corr_TFPR_z_d,corr_TFPR_z_x,
    corr_TFPR_z,mean_logProd_d,mean_logProd_x,mean_logProd,sd_logProd_d,
    sd_logProd_x,sd_logProd,import_share,worker_bond_holding,
    domestic_bond_holding,exporter_bond_holding,domestic_firm_debt,
    exporter_firm_debt,avg_Pi_x,avg_Pi_d,Profit_fine_x,
    Profit_fine_d,mean_cons,sd_cons,mean_income,sd_income,mean_wealth,
    sd_wealth,wealth_of_workers,wealth_of_domestic,wealth_of_exporters,
    income_of_workers,income_of_domestic,income_of_exporters,
    cons_of_workers,cons_of_domestic,cons_of_exporters,GINI_wealth,
    GINI_income,GINI_cons,p10_wealth,p50_wealth,p90_wealth,p99_wealth,p10_income,p50_income,p90_income,
    p99_income,p10_cons,p50_cons,p90_cons,p99_cons,Aggr_bank_cost,mean_leverage_d,mean_leverage_x,mean_leverage,
    sd_leverage_d,sd_leverage_x,sd_leverage,corr_MRPK_lev_d,corr_MRPK_lev_x,corr_MRPK_lev,
    exit_domestic_to_work_sum,exit_exporting_to_work_sum,mean_growth_rev,sd_growth_rev,mean_growth_k,sd_growth_k,
    autocorr_rev,fraction_zombie_exporter,k_opt_d_fine,k_opt_x_fine)= local_var_creator_detailed(country_no,ns_tmp,ns_tmp_fine,size_iid_cost_val);
    @sync @distributed for country = 1:country_no#
        tmp_vector  = country_residual_detailed(country,s_cell,ns_cell,s_fine_cell,ns_fine_cell,
        Phi_z_cell,Phi_z_fine_cell,Phi_cell,Phi_aug_cell,P_kron_cell,P_kron1_cell,
        P_kron_fine_cell,σ,size_iid_cost_val,price_final,W, r, output_final,θ,L,τ,
        country_no,banking_cost,δ,ω,β,α₁_eff,α₂_eff,α₁,α₂,FM,FMₓ,F,Fₓ,iid_cost_value,
        iid_cost_prob,Exit,Exitₓ,a_min,a_max,α,fspace_a,fspace_a_fine,openness,agrid_fine);
        residual[country,1] = tmp_vector[1];
        residual[2* country_no + country,1] = tmp_vector[2];
        exitflag[country,1] = tmp_vector[3];
        export_price_sum[country,1] = tmp_vector[4];
        domestic_price_sum[country,1] = tmp_vector[5];
        NFA[country,1] = tmp_vector[6];
        asset_supply[country,1] = tmp_vector[7];
        asset_demand[country,1] = tmp_vector[8];
        a_prime_fine[:,:,country,:] = tmp_vector[9];
        future_occupation_fine[:,:,country,:] = tmp_vector[10];
        cons_fine[:,:,country]  = tmp_vector[11];
        k_x_fine[:,country] = tmp_vector[12];
        l_x_fine[:,country] = tmp_vector[13];
        k_d_fine[:,country] = tmp_vector[14];
        l_d_fine[:,country] = tmp_vector[15];
        Profit_fine_x[:,country] = tmp_vector[16];
        Profit_fine_d[:,country] = tmp_vector[17];
        rev_d_fine_store[:,country] = tmp_vector[18];
        rev_dx_fine_store[:,country] = tmp_vector[19];
        rev_xx_fine_store[:,country] = tmp_vector[20];
        current_distr_store[:,country] = tmp_vector[21];
        distr_current[:,country] = tmp_vector[22];
        coeff_final[:,:,country] = tmp_vector[23];
        domestic_prod[country,1] = tmp_vector[24];
        export_prod[country,1] = tmp_vector[25];
        export_value[country,1] = tmp_vector[26];
        exit_share_to_domestic[country,1] = tmp_vector[27];
        entry_share_to_domestic[country,1] = tmp_vector[28];
        exit_share_to_worker[country,1] = tmp_vector[29];
        entry_share_to_exporter[country,1] = tmp_vector[30];
        exporter_pop[country,1] = tmp_vector[31];
        domestic_pop[country,1] = tmp_vector[32];
        worker_pop[country,1] = tmp_vector[33];
        L_x[country,1] = tmp_vector[34];
        L_d[country,1] = tmp_vector[35];
        total_entry_cost[country,1] = tmp_vector[36];
        total_incumbents_cost[country,1] = tmp_vector[37];
        total_exit_cost[country,1] = tmp_vector[38];
        labor_demand[country,1] = tmp_vector[39];
        labor_excess_demand[country,1] = tmp_vector[40];
        K_x[country,1] = tmp_vector[41];
        K_d[country,1] = tmp_vector[42];
        capital_demand[country,1] = tmp_vector[43];
        worker_bond_holding[country,1] = tmp_vector[44];
        domestic_bond_holding[country,1] = tmp_vector[45];
        exporter_bond_holding[country,1] = tmp_vector[46];
        domestic_firm_debt[country,1] = tmp_vector[47];
        exporter_firm_debt[country,1] = tmp_vector[48];
        total_consumption[country,1] = tmp_vector[49];
        Aggr_bank_cost[country,1] = tmp_vector[50];
        x_k[country,1] = tmp_vector[51];
        x_l[country,1] = tmp_vector[52];
        investment[country,1] = tmp_vector[53];
        total_demand_final_good[country,1] = tmp_vector[54];
        excess_demand_final_good[country,1] = tmp_vector[55];
        mean_MRPK_d[country,1] = tmp_vector[56];
        mean_MRPK_x[country,1 ] = tmp_vector[57];
        mean_MRPK[country,1] = tmp_vector[58];
        sd_MRPK_d[country,1] = tmp_vector[59];
        sd_MRPK_x[country,1] = tmp_vector[60];
        sd_MRPK[country,1] = tmp_vector[61];
        mean_logProd_d[country,1] = tmp_vector[62];
        mean_logProd_x[country,1] = tmp_vector[63];
        mean_logProd[country,1] = tmp_vector[64];
        sd_logProd_d[country,1] = tmp_vector[65];
        sd_logProd_x[country,1] = tmp_vector[66];
        sd_logProd[country,1] = tmp_vector[67];
        cov_TFPR_z_d[country,1] = tmp_vector[68];
        cov_TFPR_z_x[country,1] = tmp_vector[69];
        cov_TFPR_z[country,1] = tmp_vector[70];
        corr_TFPR_z_d[country,1] = tmp_vector[71];
        corr_TFPR_z_x[country,1] = tmp_vector[72];
        corr_TFPR_z[country,1] = tmp_vector[73];
        avg_Pi_x[country,1] = tmp_vector[74];
        avg_Pi_d[country,1] = tmp_vector[75];
        D_d_denom[country,1] = tmp_vector[76];
        D_x_denom[country,1] = tmp_vector[77];
        D_x[country,1] = tmp_vector[78];
        D_d[country,1] = tmp_vector[79];
        D_d_denom_eff[country,1] = tmp_vector[80];
        D_x_denom_eff[country,1] = tmp_vector[81];
        D_x_eff[country,1] = tmp_vector[82];
        D_d_eff[country,1] = tmp_vector[83];
        mean_cons[country,1] = tmp_vector[84];
        sd_cons[country,1] = tmp_vector[85];
        mean_income[country,1] = tmp_vector[86];
        sd_income[country,1] = tmp_vector[87];
        mean_wealth[country,1] = tmp_vector[88];
        sd_wealth[country,1] = tmp_vector[89];
        p10_wealth[country,1] = tmp_vector[90];
        p50_wealth[country,1] = tmp_vector[91];
        p90_wealth[country,1] = tmp_vector[92];
        p99_wealth[country,1] = tmp_vector[93];
        p10_income[country,1] = tmp_vector[94];
        p50_income[country,1] = tmp_vector[95];
        p90_income[country,1] = tmp_vector[96];
        p99_income[country,1] = tmp_vector[97];
        p10_cons[country,1] = tmp_vector[98];
        p50_cons[country,1] = tmp_vector[99];
        p90_cons[country,1] = tmp_vector[100];
        p99_cons[country,1] = tmp_vector[101];
        wealth_of_workers[country,1] = tmp_vector[102];
        wealth_of_domestic[country,1] = tmp_vector[103];
        wealth_of_exporters[country,1] = tmp_vector[104];
        income_of_workers[country,1] = tmp_vector[105];
        income_of_domestic[country,1] = tmp_vector[106];
        income_of_exporters[country,1] = tmp_vector[107];
        cons_of_workers[country,1] = tmp_vector[108];
        cons_of_domestic[country,1] = tmp_vector[109];
        cons_of_exporters[country,1] = tmp_vector[110];
        mean_leverage_d[country,1] = tmp_vector[111];
        mean_leverage_x[country,1] = tmp_vector[112];
        mean_leverage[country,1] = tmp_vector[113];
        sd_leverage_d[country,1] = tmp_vector[114];
        sd_leverage_x[country,1] = tmp_vector[115];
        sd_leverage[country,1] = tmp_vector[116];
        corr_MRPK_lev_d[country,1] = tmp_vector[117];
        corr_MRPK_lev_x[country,1] = tmp_vector[118];
        corr_MRPK_lev[country,1] = tmp_vector[119];
        exit_domestic_to_work_sum[country,1] = tmp_vector[120];
        exit_exporting_to_work_sum[country,1] = tmp_vector[121];
        mean_growth_rev[country,1] = tmp_vector[122];
        sd_growth_rev[country,1] = tmp_vector[123];
        mean_growth_k[country,1] = tmp_vector[124];
        sd_growth_k[country,1] = tmp_vector[125];
        autocorr_rev[country,1] = tmp_vector[126];
        fraction_zombie_exporter[country,1] = tmp_vector[127];
        k_opt_d_fine[:,country] = tmp_vector[128];
        k_opt_x_fine[:,country] = tmp_vector[129];
    end
    for country=1:country_no
        if exitflag[country,1]==1 || exitflag[country,1]==2# Too large supply of assets
            if openness == 1
                residual[7,1] =  100.0;
            else
                residual[3* country_no + country,1] = 100.0;
            end
        elseif exitflag[country,1]==3# Too little supply of assets?
            if openness == 1
                residual[7,1] = - 100.0;
            else
                residual[3* country_no + country,1] = - 100.0;
            end
        elseif exitflag[country,1]==4
            if openness == 1
                residual[7,1] =  1000.0;
            else
                residual[3* country_no + country,1] =  1000.0;
            end
        elseif exitflag[country,1]==5
            if openness == 1
                residual[7,1] = 500.0;
            else
                residual[3* country_no + country,1] =  500.0;
            end
        else
            import_curr = (sum(export_prod[:,1])-  export_prod[country,1])/(country_no - 1.0);
            total_production[country,1] = (ω *domestic_prod[country,1] + (1.0 -ω) *import_curr)^((σ)/(σ-1));
            import_price[country,1] = (sum(export_price_sum[:,1])-  export_price_sum[country,1])/(country_no - 1);
            price_final_actual[country,1] = min((ω^σ * domestic_price_sum[country,1] + (1.0 - ω)^σ  * import_price[country,1])^(1.0
                / (1.0 - σ)),10000);
            residual[country_no + country,1] = (price_final_actual[country,1] - price_final[country,1])./(
                price_final_actual[country,1] + price_final[country,1]);
            if openness == 0
                residual[3* country_no + country,1] = NFA[country,1]/(asset_demand[country,1] + (1 + banking_cost[country,1]).*asset_supply[country,1]);
            end
        end
    end
    if openness == 1
        #residual[7,1] = sum(NFA)/sum((asset_demand+ (ones(country_no) + banking_cost).*asset_supply));
        residual[7,1] = sum(NFA)/output_final[1];
    end
    TFP_d = ω * (D_d.^(1 - α₂_eff)./D_d_denom.^α₁_eff);
    TFP_x = (D_x.^(1 - α₂_eff)./D_x_denom.^α₁_eff);

    #TFP_d = D_d./D_d_denom;
    #TFP_x = D_x./D_x_denom;
    TFP_d_efficient = ω * (D_d_eff.^(1 - α₂_eff)./D_d_denom_eff.^α₁_eff);
    TFP_x_efficient = (D_x_eff.^(1 - α₂_eff)./D_x_denom_eff.^α₁_eff);
    #TFP_d_efficient = D_d_eff./D_d_denom_eff;
    #TFP_x_efficient = D_x_eff./D_x_denom_eff;
    # Other definitions for TFP and some other useful measures
    for country = 1:country_no
        total_foreign_demand = (sum(total_production) - total_production[country,1])/(country_no - 1);
        #total_foreign_demand = (sum(output_final) - output_final[country,1])/(country_no - 1);
        avg_foreign_price = (sum(price_final_actual) - price_final_actual[country,1])/(country_no - 1);
        nomGDP[country,1] = price_final_actual[country,1] * total_production[country,1];
        GDP[country,1] = total_production[country,1];
        import_value_tmp = (sum(export_value[:,1])-  export_value[country,1]);
        export_value_tradeoff = import_value_tmp/export_value[country,1];
        import_share[country,1] = import_value_tmp/(country_no - 1)/nomGDP[country,1];
        #constant_d_tmp =   ω * price_final_actual[country,1] * output_final[country,1]^(1/σ);
        constant_d_tmp =   ω * price_final_actual[country,1] * total_production[country,1]^(1/σ);
        constant_x_tmp = (1 - ω)/ (1 + τ[country,1]) * avg_foreign_price * total_foreign_demand^(1/σ);
        TOT[country,1] =  (ω *constant_d_tmp^(σ - 1) + (1 - ω) *constant_x_tmp^(σ - 1)*(
         export_value_tradeoff *price_final_actual[country,1]/ avg_foreign_price)^((σ -1)/σ))/(
        constant_d_tmp^σ + (1 +τ[country])* constant_x_tmp^σ
        )^((σ -1)/σ);
        TFP_x[country,1] = TFP_x[country,1] *TOT[country,1];
        TFP_x_efficient[country,1] = TFP_x_efficient[country,1] * TOT[country,1];
        D_d_K_tmp = D_d_denom[country,1];
        D_x_K_tmp = D_x_denom[country,1];
        D_d_K_eff_tmp = D_d_denom_eff[country,1];
        D_x_K_eff_tmp = D_x_denom_eff[country,1];
        K_d_ratio[country,1] = constant_d_tmp^σ * D_d_K_tmp/(constant_d_tmp^σ * D_d_K_tmp + (constant_d_tmp^σ
         +(1 +τ[country])*constant_x_tmp^σ)* D_x_K_tmp);
        K_x_ratio[country,1] = (constant_d_tmp^σ
         +(1 +τ[country])*constant_x_tmp^σ)* D_x_K_tmp/(constant_d_tmp^σ * D_d_K_tmp + (constant_d_tmp^σ
          +(1 +τ[country])*constant_x_tmp^σ)* D_x_K_tmp);
        L_d_ratio[country,1] = constant_d_tmp^σ * D_d[country,1]/(constant_d_tmp^σ * D_d[country,1] + (constant_d_tmp^σ
            +(1 +τ[country])*constant_x_tmp^σ)* D_x[country,1]);
        L_x_ratio[country,1] = (constant_d_tmp^σ
            +(1 +τ[country])*constant_x_tmp^σ)* D_x[country,1]/(constant_d_tmp^σ * D_d[country,1] + (constant_d_tmp^σ
            +(1 +τ[country])*constant_x_tmp^σ)* D_x[country,1]);

        K_d_ratio_eff[country,1] = constant_d_tmp^σ * D_d_K_eff_tmp/(constant_d_tmp^σ * D_d_K_eff_tmp + (constant_d_tmp^σ
         +(1 +τ[country])*constant_x_tmp^σ)* D_x_K_eff_tmp);
        K_x_ratio_eff[country,1] = (constant_d_tmp^σ
         +(1 +τ[country])*constant_x_tmp^σ)* D_x_K_eff_tmp/(constant_d_tmp^σ * D_d_K_eff_tmp + (constant_d_tmp^σ
          +(1 +τ[country])*constant_x_tmp^σ)* D_x_K_eff_tmp);
        L_d_ratio_eff[country,1] = constant_d_tmp^σ * D_d_eff[country,1]/(constant_d_tmp^σ * D_d_eff[country,1] + (constant_d_tmp^σ
            +(1 +τ[country])*constant_x_tmp^σ)* D_x_eff[country,1]);
        L_x_ratio_eff[country,1] = (constant_d_tmp^σ
            +(1 +τ[country])*constant_x_tmp^σ)* D_x_eff[country,1]/(constant_d_tmp^σ * D_d_eff[country,1] + (constant_d_tmp^σ
            +(1 +τ[country])*constant_x_tmp^σ)* D_x_eff[country,1]);
        TFP[country,1] = (TFP_d[country,1] * K_d_ratio[country,1]^α₁_eff * L_d_ratio[country,1]^α₂_eff
            + TFP_x[country,1] * K_x_ratio[country,1]^α₁_eff * L_x_ratio[country,1]^α₂_eff)^(σ/(σ-1));
        TFP_within[country,1] = (TFP_d_efficient[country,1] * K_d_ratio[country,1]^α₁_eff * L_d_ratio[country,1]^α₂_eff
            + TFP_x_efficient[country,1] * K_x_ratio[country,1]^α₁_eff * L_x_ratio[country,1]^α₂_eff)^(σ/(σ-1));
        TFP_across[country,1] = (TFP_d[country,1] * K_d_ratio_eff[country,1]^α₁_eff * L_d_ratio_eff[country,1]^α₂_eff
            + TFP_x[country,1] * K_x_ratio_eff[country,1]^α₁_eff * L_x_ratio_eff[country,1]^α₂_eff)^(σ/(σ-1)); # this doesnt work
        TFP_second_best[country,1] =(TFP_d_efficient[country,1] * K_d_ratio_eff[country,1]^α₁_eff * L_d_ratio_eff[country,1]^α₂_eff
            + TFP_x_efficient[country,1] * K_x_ratio_eff[country,1]^α₁_eff * L_x_ratio_eff[country,1]^α₂_eff)^(σ/(σ-1));
    end
    Misallocation_within_d =ones(country_no) -  TFP_d./TFP_d_efficient;
    Misallocation_within_x =ones(country_no) - TFP_x./TFP_x_efficient;
    replace!(residual, NaN=>100.0);
    return (residual,a_prime_fine,future_occupation_fine,cons_fine,rev_d_fine_store,
    rev_dx_fine_store,rev_xx_fine_store,coeff_final,k_x_fine,k_d_fine,
    l_x_fine,l_d_fine,current_distr_store,distr_current,exitflag,domestic_prod,exporter_pop,
    domestic_pop,worker_pop,entry_share_to_domestic,exit_share_to_domestic,
    entry_share_to_exporter,exit_share_to_worker,total_entry_cost,
    total_incumbents_cost,total_exit_cost,labor_excess_demand,labor_demand,
    total_consumption,capital_demand,capital_next,
    total_demand_final_good,excess_demand_final_good,GDP,nomGDP,PPI,TFP,
    TFP_within,TFP_across,TFP_second_best,TOT,TFP_d,
    TFP_d_efficient,TFP_x,TFP_x_efficient,total_production,
    domestic_price_sum,export_value,investment,
    export_value,L_x,L_d,K_x,K_d,x_k,x_l,D_d,D_x,D_d_denom,D_x_denom,
    D_d_eff,D_x_eff,D_d_denom_eff,D_x_denom_eff,K_d_ratio,K_x_ratio,
    L_d_ratio,L_x_ratio,K_x_ratio_eff,K_d_ratio_eff,L_d_ratio_eff,
    L_x_ratio_eff,nomGDP_d,nomGDP_xd,nomGDP_xx,RGDP_d,RGDP_xd,RGDP_xx,RGDP,
    mean_MRPK_d,mean_MRPK_x,sd_MRPK_d,sd_MRPK_x,mean_MRPK,sd_MRPK,
    mean_MRPL_d,mean_MRPL_x,sd_MRPL_d,sd_MRPL_x,mean_MRPL,sd_MRPL,
    cov_TFPR_z_d,cov_TFPR_z_x,cov_TFPR_z,corr_TFPR_z_d,corr_TFPR_z_x,
    corr_TFPR_z,mean_logProd_d,mean_logProd_x,mean_logProd,sd_logProd_d,
    sd_logProd_x,sd_logProd,import_share,worker_bond_holding,
    domestic_bond_holding,exporter_bond_holding,domestic_firm_debt,
    exporter_firm_debt,avg_Pi_x,avg_Pi_d,Profit_fine_x,
    Profit_fine_d,mean_cons,sd_cons,mean_income,sd_income,mean_wealth,
    sd_wealth,wealth_of_workers,wealth_of_domestic,wealth_of_exporters,
    income_of_workers,income_of_domestic,income_of_exporters,
    cons_of_workers,cons_of_domestic,cons_of_exporters,GINI_wealth,
    GINI_income,GINI_cons,p10_wealth,p50_wealth,p90_wealth,p99_wealth,p10_income,p50_income,p90_income,
    p99_income,p10_cons,p50_cons,p90_cons,p99_cons,Aggr_bank_cost,mean_leverage_d,mean_leverage_x,mean_leverage,
    sd_leverage_d,sd_leverage_x,sd_leverage,corr_MRPK_lev_d,corr_MRPK_lev_x,corr_MRPK_lev,exit_domestic_to_work_sum,
    exit_exporting_to_work_sum,mean_growth_rev,sd_growth_rev,mean_growth_k,sd_growth_k,autocorr_rev,fraction_zombie_exporter,
    Misallocation_within_d,Misallocation_within_x,k_opt_d_fine,k_opt_x_fine)
end
function true_profit_frictionless(price_final_prev::Float64, R::Float64,W_loc::Float64 ,s_curr::Array{Float64,2},
        constant_x::Float64,constant_d::Float64,constant_z::Array{Float64,1},constant_z_bar::Array{Float64,1},
        θ_loc::Float64,τ_loc::Float64,σ::Float64,α₁_eff::Float64,α₂_eff::Float64,α₁::Float64,α₂::Float64,
        ones_tmp::Array{Float64,1})
    # Given prices, evaluate the income of households needed for VFI
    (k_d,k_x ) = optimal_k(R,W_loc,constant_z,constant_z_bar ,σ,α₁_eff,α₂_eff);

    #Implied Lagrange multiplier

    labor_d = (α₁_eff^(α₁_eff) * α₂_eff^(1 - α₁_eff) * constant_z_bar.* (R * ones_tmp ).^(-α₁_eff) * W_loc^(α₁_eff-1)).^σ;
    labor_x = (α₁_eff^(α₁_eff) * α₂_eff^(1 - α₁_eff) * constant_z.* (R* ones_tmp ).^(-α₁_eff) * W_loc^(α₁_eff-1)).^σ;
    total_produced_x = (s_curr[:,2]).* k_x.^(α₁).* labor_x.^(α₂);
    output_dx =total_produced_x * constant_d^σ / (constant_d^σ + (1+τ_loc)*(constant_x)^σ);
    output_xx =total_produced_x * constant_x^σ / (constant_d^σ + (1+τ_loc)*(constant_x)^σ); # Note - iceberg cost paid here
    output_d = ( s_curr[:,2]).* k_d.^(α₁).* labor_d.^(α₂);
    price_dx = constant_d .* output_dx.^(-1/σ);
    price_d = constant_d .* output_d.^(-1/σ);
    price_xx = (1+τ_loc)* constant_x .* output_xx.^(-1/σ); # final price without tariffs
    rev_d = price_d.*output_d ;
    rev_dx = price_dx.*output_dx;
    rev_xx = price_xx.*output_xx;
    Profits_d = rev_d - W_loc*labor_d - R *k_d;
    Profits_x = rev_dx + rev_xx   - W_loc*labor_x - R *k_x;
    return (Profits_d,Profits_x,output_dx,output_xx,output_d,labor_d,labor_x,k_d,
        k_x,price_dx,price_d,price_xx,rev_d,rev_dx,rev_xx)
end

function country_residual_frictionless(country::Int64,s_cell::Array{Array{Float64,2},1},
    ns_cell::Array{Int64,1},s_fine_cell::Array{Array{Float64,2},1},
    ns_fine_cell::Array{Int64,1},
    Phi_z_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    Phi_z_fine_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    Phi_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    Phi_aug_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    P_kron_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    P_kron1_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    P_kron_fine_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    σ::Float64,size_iid_cost_val::Int64,price_final::Array{Float64,2},
    W::Array{Float64,1}, r::Array{Float64,1}, output_final::Array{Float64,1},
    θ::Array{Float64,1}, L::Array{Float64,1}, τ::Array{Float64,1},
    country_no::Int64,banking_cost::Array{Float64,1},δ::Float64,ω::Float64,
    β::Array{Float64,1},
    α₁_eff::Float64,α₂_eff::Float64,α₁::Float64,α₂::Float64,FM::Array{Float64,1},FMₓ::Array{Float64,1},
    F::Array{Float64,1}, Fₓ::Array{Float64,1},
    iid_cost_value::Array{Float64,1},iid_cost_prob::Array{Float64,1},Exit::Array{Float64,1},
    Exitₓ::Array{Float64,1},a_min::Float64,a_max::Float64,α::Float64,fspace_a::Dict{Symbol,Any},
    fspace_a_fine::Dict{Symbol,Any},openness::Int64,agrid_fine::Array{Float64,1})
    # Initialize each country # Threads.@threads for country = 1:country_no
    (exitflag_tmp,s,ns,s_fine,ns_fine,Phi_z,Phi_z_fine,Phi,Phi_aug,
    P_kron,P_kron1,P_kron_fine,z_tilde,z_tilde_fine,income_mat,income_mat_fine,
    a_prime_fine_local,future_occupation_fine_local,cons_fine_local,coeff,
    coeff_next,price_final_prev,price_final_current,total_foreign_demand,W_loc,
    θ_loc,L_loc,τ_loc,r_loc,banking_cost_loc,avg_foreign_price,R,constant_d,
    constant_x,constant_z,constant_z_bar,constant_z_fine,constant_z_bar_fine,
    ones_tmp,x_tmp,V_tmp,β_loc,V_next_stacked,iterate1,Phi_prime_tmp,
    D_deriv_tmp_block,conv,Q_trans,ones_tmp_fine) = country_local_initialize(country,s_cell,ns_cell,s_fine_cell,ns_fine_cell,
    Phi_z_cell,Phi_z_fine_cell,Phi_cell,Phi_aug_cell,P_kron_cell,P_kron1_cell,
    P_kron_fine_cell,σ,size_iid_cost_val,price_final,W, r, output_final,θ,L,τ,
    country_no,banking_cost,δ,ω,β);

    (Profits_d,Profits_x,output_dx,output_xx,output_d,labor_d,labor_x,k_choice_d,k_choice_x,price_dx,price_d,price_xx
    ,rev_d,rev_dx,rev_xx ) = true_profit_frictionless(price_final_prev,R,W_loc,s,
        constant_x,constant_d,constant_z,constant_z_bar,θ_loc,τ_loc,σ,α₁_eff,α₂_eff,α₁,α₂,ones_tmp);
    (income_mat,max_x) = income_creator(W_loc,ones_tmp,Profits_d,Profits_x,FM,FMₓ,country,
        r_loc,s,iid_cost_value,F, Fₓ,size_iid_cost_val,income_mat,ns,a_min,a_max);
    for i = 1:3
        (coeff[:], conv) = Bellman_iteration(coeff,coeff_next,ns,x_tmp,V_tmp,ones_tmp,size_iid_cost_val,
            iid_cost_prob,income_mat,P_kron,Phi,a_min,Phi_z,β_loc,
            fspace_a,price_final_current,max_x);
    end
    while conv> 10^(-12)
        (coeff[:], conv,iterate1,exitflag_tmp) = Newton_iteration(coeff,coeff_next,ns,x_tmp,V_tmp,ones_tmp,size_iid_cost_val,
            iid_cost_prob,income_mat,P_kron,Phi,a_min,Phi_z,β_loc,
            fspace_a,price_final_current,max_x,V_next_stacked,iterate1,
            Phi_prime_tmp,D_deriv_tmp_block,Phi_aug,P_kron1,exitflag_tmp);
    end
    # Initialize return values if something goes wrong
    export_price_sum_tmp=0.0;
    domestic_price_sum_tmp=0.0;
    NFA_tmp=0.0;
    asset_supply_tmp=0.0;
    asset_demand_tmp=0.0;
    # Solving for the aggregates - stationary distribution:
    if exitflag_tmp>0
        println("Nonstationarity")
    else
        #Continue only if we have convergence in VFI
        (Profits_d_fine,Profits_x_fine,output_dx_fine,output_xx_fine,output_d_fine,labor_d_fine,labor_x_fine,
        k_choice_d_fine,k_choice_x_fine, price_dx_fine,price_d_fine,price_xx_fine,rev_d_fine,rev_dx_fine,rev_xx_fine
        ) = true_profit_frictionless(price_final_prev,R,W_loc,s_fine
            ,constant_x,constant_d,constant_z_fine,constant_z_bar_fine,θ_loc,τ_loc,σ,α₁_eff,α₂_eff,α₁,α₂,ones_tmp_fine);
        (income_mat_fine,max_x_fine) = income_creator(W_loc,ones_tmp_fine,Profits_d_fine,Profits_x_fine,FM,FMₓ,country,
                r_loc,s_fine,iid_cost_value,F, Fₓ,size_iid_cost_val,income_mat_fine,ns_fine,a_min,a_max);
        (Q_trans_prime,cons_fine_local,a_prime_fine_local,
        future_occupation_fine_local) = Q_transition(coeff,
            ns_fine,ones_tmp_fine,size_iid_cost_val,iid_cost_prob,income_mat_fine,
            P_kron_fine,a_min,Phi_z_fine,β_loc,fspace_a,fspace_a_fine,
            price_final_current,max_x_fine,P_kron1,cons_fine_local,a_prime_fine_local,
            future_occupation_fine_local,Q_trans);
        exitflag_tmp = predict_irreducibility(future_occupation_fine_local,exitflag_tmp);
    end
    if exitflag_tmp==1 || exitflag_tmp==2
        labor_excess_demand_tmp_percent = -200.0 ;
        excess_demand_final_good_tmp_percent = 200.0
    elseif exitflag_tmp==3
        labor_excess_demand_tmp_percent = 200.0 ;
        excess_demand_final_good_tmp_percent = -200.0
    elseif exitflag_tmp==4
        labor_excess_demand_tmp_percent = 2000.0 ;
        excess_demand_final_good_tmp_percent = 2000.0
    else
        distr_current,exitflag_tmp = stationary_distribution(Q_trans_prime,L_loc,ns_fine,exitflag_tmp);
        if exitflag_tmp==5
            labor_excess_demand_tmp_percent = 1000.0 ;
            excess_demand_final_good_tmp_percent = 1000.0
        else
            current_distr_store_tmp = copy(distr_current);
            # Calculate aggregates
            # Distribution of current workers:
            worker_past_dist = distr_current[(ns_fine *0 + 1):(ns_fine *1)];
            domestic_past_dist = distr_current[(ns_fine *1 + 1):(ns_fine *2)];
            exporter_past_dist = distr_current[(ns_fine *2 + 1):(ns_fine *3)];

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
            investment_tmp = δ * capital_demand_tmp;
            total_demand_final_good_tmp = total_consumption_tmp  + investment_tmp;
            excess_demand_final_good_tmp = total_demand_final_good_tmp - output_final[country,1];
            #Save coefficients and policy function for the transition dynamics
            labor_excess_demand_tmp_percent = labor_excess_demand_tmp./L_loc;
            excess_demand_final_good_tmp_percent = excess_demand_final_good_tmp/(output_final[country,1] + total_demand_final_good_tmp);
        end
    end
    return (labor_excess_demand_tmp_percent,excess_demand_final_good_tmp_percent,
    exitflag_tmp,export_price_sum_tmp, domestic_price_sum_tmp,NFA_tmp,asset_supply_tmp,asset_demand_tmp)
end
function Residual_stst_frictionless(prices::Array{Float64,1},parameters_tmp::Parameter_type)
    # Extract the necessary local parameters:
    (β,α,δ,θ,α₁,α₂,σ,α₁_eff,α₂_eff,ω,L,FM,FMₓ,F,Fₓ,
        Exit,Exitₓ,iid_cost_value,iid_cost_prob,size_iid_cost_val,country_no,τ,
        a_min,a_max,fspace_a,fspace_a_fine,agrid_fine,banking_cost,
        bounds,Country_spec_p,s_cell,ns_cell,s_fine_cell,ns_fine_cell,Phi_z_cell,
        Phi_z_fine_cell,Phi_cell,Phi_aug_cell,P_kron_cell,P_kron1_cell,P_kron_fine_cell,
        exp_egrid_cell,ns_tmp,ns_tmp_fine,openness) = local_parameters(parameters_tmp)
    # Prices:
    (W,r,output_final,price_final,r,residual,price_check_tmp) = price_reshaper(prices,openness,country_no,bounds)
    # Check prices to be reasonable:
    if price_check_tmp>0
        println("Guess out of bounds")
        residual = 1000.0 * ones(8);
        return residual
    end
    # Initialization of local variables:
    (export_prod,export_price_sum,import_price,domestic_price_sum,
    price_final_actual,NFA,asset_supply,asset_demand,exitflag) = local_var_creator_minimal(country_no);
    @sync @distributed for country = 1:country_no#for country = 1:country_no#
        (residual[country,1],residual[2* country_no + country,1],exitflag[country,1],export_price_sum[country,1],
        domestic_price_sum[country,1],NFA[country,1],asset_supply[country,1],
        asset_demand[country,1])  = country_residual_frictionless(country,s_cell,ns_cell,s_fine_cell,ns_fine_cell,
        Phi_z_cell,Phi_z_fine_cell,Phi_cell,Phi_aug_cell,P_kron_cell,P_kron1_cell,
        P_kron_fine_cell,σ,size_iid_cost_val,price_final,W, r, output_final,θ,L,τ,
        country_no,banking_cost,δ,ω,β,α₁_eff,α₂_eff,α₁,α₂,FM,FMₓ,F,Fₓ,iid_cost_value,
        iid_cost_prob,Exit,Exitₓ,a_min,a_max,α,fspace_a,fspace_a_fine,openness,agrid_fine);
    end
    for country=1:country_no
        if exitflag[country,1]==1 || exitflag[country,1]==2# Too large supply of assets
            if openness == 1
                residual[7,1] =  100.0;
            else
                residual[3* country_no + country,1] = 100.0;
            end
        elseif exitflag[country,1]==3# Too little supply of assets?
            if openness == 1
                residual[7,1] = - 100.0;
            else
                residual[3* country_no + country,1] = - 100.0;
            end
        elseif exitflag[country,1]==4
            if openness == 1
                residual[7,1] =  1000.0;
            else
                residual[3* country_no + country,1] =  1000.0;
            end
        elseif exitflag[country,1]==5
            if openness == 1
                residual[7,1] = 500.0;
            else
                residual[3* country_no + country,1] =  500.0;
            end
        else
            import_curr = (sum(export_prod[:,1])-  export_prod[country,1])/(country_no - 1.0);
            import_price[country,1] = (sum(export_price_sum[:,1])-  export_price_sum[country,1])/(country_no - 1);
            price_final_actual[country,1] = min((ω^σ * domestic_price_sum[country,1] + (1.0 - ω)^σ  * import_price[country,1])^(1.0
                / (1.0 - σ)),10000);
            residual[country_no + country,1] = (price_final_actual[country,1] - price_final[country,1])./(
                price_final_actual[country,1] + price_final[country,1]);
            if openness == 0
                residual[3* country_no + country,1] = NFA[country,1]/(asset_demand[country,1] + (1 + banking_cost[country,1]).*asset_supply[country,1]);
            end
        end
        if openness == 1
            residual[7,1] = sum(NFA)/sum((asset_demand+ (ones(country_no) + banking_cost).*asset_supply));
        end
    end
    replace!(residual, NaN=>100.0);
    return residual
end
function country_residual_detailed_frictionless(country::Int64,s_cell::Array{Array{Float64,2},1},
    ns_cell::Array{Int64,1},s_fine_cell::Array{Array{Float64,2},1},
    ns_fine_cell::Array{Int64,1},
    Phi_z_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    Phi_z_fine_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    Phi_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    Phi_aug_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    P_kron_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    P_kron1_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    P_kron_fine_cell::Array{SparseMatrixCSC{Float64,Int64},1},
    σ::Float64,size_iid_cost_val::Int64,price_final::Array{Float64,2},
    W::Array{Float64,1}, r::Array{Float64,1}, output_final::Array{Float64,1},
    θ::Array{Float64,1}, L::Array{Float64,1}, τ::Array{Float64,1},
    country_no::Int64,banking_cost::Array{Float64,1},δ::Float64,ω::Float64,
    β::Array{Float64,1},
    α₁_eff::Float64,α₂_eff::Float64,α₁::Float64,α₂::Float64,FM::Array{Float64,1},FMₓ::Array{Float64,1},
    F::Array{Float64,1}, Fₓ::Array{Float64,1},
    iid_cost_value::Array{Float64,1},iid_cost_prob::Array{Float64,1},Exit::Array{Float64,1},
    Exitₓ::Array{Float64,1},a_min::Float64,a_max::Float64,α::Float64,fspace_a::Dict{Symbol,Any},
    fspace_a_fine::Dict{Symbol,Any},openness::Int64,agrid_fine::Array{Float64,1})
    # Initialize each country # Threads.@threads for country = 1:country_no
    (exitflag_tmp,s,ns,s_fine,ns_fine,Phi_z,Phi_z_fine,Phi,Phi_aug,
    P_kron,P_kron1,P_kron_fine,z_tilde,z_tilde_fine,income_mat,income_mat_fine,
    a_prime_fine_local,future_occupation_fine_local,cons_fine_local,coeff,
    coeff_next,price_final_prev,price_final_current,total_foreign_demand,W_loc,
    θ_loc,L_loc,τ_loc,r_loc,banking_cost_loc,avg_foreign_price,R,constant_d,
    constant_x,constant_z,constant_z_bar,constant_z_fine,constant_z_bar_fine,
    ones_tmp,x_tmp,V_tmp,β_loc,V_next_stacked,iterate1,Phi_prime_tmp,
    D_deriv_tmp_block,conv,Q_trans,ones_tmp_fine) = country_local_initialize(country,s_cell,ns_cell,s_fine_cell,ns_fine_cell,
    Phi_z_cell,Phi_z_fine_cell,Phi_cell,Phi_aug_cell,P_kron_cell,P_kron1_cell,
    P_kron_fine_cell,σ,size_iid_cost_val,price_final,W, r, output_final,θ,L,τ,
    country_no,banking_cost,δ,ω,β);
    (Profits_d,Profits_x,output_dx,output_xx,output_d,labor_d,labor_x,k_choice_d,k_choice_x,price_dx,price_d,price_xx
    ,rev_d,rev_dx,rev_xx ) = true_profit_frictionless(price_final_prev,R,W_loc,s,
        constant_x,constant_d,constant_z,constant_z_bar,θ_loc,τ_loc,σ,α₁_eff,α₂_eff,α₁,α₂,ones_tmp);
    (income_mat,max_x) = income_creator(W_loc,ones_tmp,Profits_d,Profits_x,FM,FMₓ,country,
        r_loc,s,iid_cost_value,F, Fₓ,size_iid_cost_val,income_mat,ns,a_min,a_max);
    for i = 1:3
        (coeff[:], conv) = Bellman_iteration(coeff,coeff_next,ns,x_tmp,V_tmp,ones_tmp,size_iid_cost_val,
            iid_cost_prob,income_mat,P_kron,Phi,a_min,Phi_z,β_loc,
            fspace_a,price_final_current,max_x);
    end
    while conv> 10^(-12)
        (coeff[:], conv,iterate1,exitflag_tmp) = Newton_iteration(coeff,coeff_next,ns,x_tmp,V_tmp,ones_tmp,size_iid_cost_val,
            iid_cost_prob,income_mat,P_kron,Phi,a_min,Phi_z,β_loc,
            fspace_a,price_final_current,max_x,V_next_stacked,iterate1,
            Phi_prime_tmp,D_deriv_tmp_block,Phi_aug,P_kron1,exitflag_tmp);
    end
    # Initialize return values if something goes wrong
    export_price_sum_tmp=0.0;domestic_price_sum_tmp=0.0;labor_excess_demand_tmp_percent=0.0;
    excess_demand_final_good_tmp_percent=0.0; NFA_tmp=0.0;asset_supply_tmp=0.0;asset_demand_tmp=0.0;
    current_distr_store_tmp=zeros(3*ns_fine);distr_current=zeros(3*ns_fine);domestic_prod_tmp=0.0;
    export_prod_tmp=0.0;export_value_tmp=0.0;exit_share_to_domestic_tmp=0.0;entry_share_to_domestic_tmp=0.0;
    exit_share_to_worker_tmp=0.0;entry_share_to_exporter_tmp=0.0;exporter_pop_tmp=0.0;domestic_pop_tmp=0.0;
    worker_pop_tmp=0.0;L_x_tmp=0.0;L_d_tmp=0.0;total_entry_cost_tmp=0.0;total_incumbents_cost_tmp=0.0;total_exit_cost_tmp=0.0;
    labor_demand_tmp=0.0;labor_excess_demand_tmp=0.0;K_x_tmp=0.0;K_d_tmp=0.0;capital_demand_tmp=0.0;worker_bond_holding_sum_tmp=0.0;
    domestic_bond_holding_sum_tmp=0.0;exporter_bond_holding_sum_tmp=0.0;domestic_firm_debt_tmp=0.0;exporter_firm_debt_tmp=0.0;
    total_consumption_tmp=0.0;Aggr_bank_cost_tmp=0.0;x_k_tmp=0.0;x_l_tmp=0.0;investment_tmp=0.0;total_demand_final_good_tmp=0.0;
    excess_demand_final_good_tmp=0.0;mean_MRPK_d_tmp=0.0;mean_MRPK_x_tmp=0.0;mean_MRPK_tmp=0.0;sd_MRPK_d_tmp=0.0; sd_MRPK_x_tmp=0.0;sd_MRPK_tmp=0.0;
    mean_logProd_d_tmp=0.0;mean_logProd_x_tmp=0.0;mean_logProd_tmp=0.0;sd_logProd_d_tmp=0.0;sd_logProd_x_tmp=0.0;sd_logProd_tmp=0.0;
    cov_TFPR_z_d_tmp=0.0;cov_TFPR_z_x_tmp=0.0;cov_TFPR_z_tmp=0.0;corr_TFPR_z_d_tmp=0.0;corr_TFPR_z_x_tmp=0.0;corr_TFPR_z_tmp=0.0;
    avg_Pi_x_tmp=0.0;avg_Pi_d_tmp=0.0;D_d_denom_tmp=0.0;D_x_denom_tmp=0.0;D_x_tmp=0.0;D_d_tmp=0.0;D_d_denom_eff_tmp=0.0;D_x_denom_eff_tmp=0.0;
    D_x_eff_tmp=0.0;D_d_eff_tmp=0.0;mean_cons_tmp=0.0;sd_cons_tmp=0.0;mean_income_tmp=0.0;sd_income_tmp=0.0;mean_wealth_tmp=0.0;sd_wealth_tmp=0.0;
    p10_wealth_tmp=0.0;p50_wealth_tmp=0.0;p90_wealth_tmp=0.0;p99_wealth_tmp=0.0;p10_income_tmp=0.0;p50_income_tmp=0.0;p90_income_tmp=0.0;p99_income_tmp=0.0;p10_cons_tmp=0.0;
    p50_cons_tmp=0.0;p90_cons_tmp=0.0;p99_cons_tmp=0.0;wealth_of_workers_tmp=0.0;wealth_of_domestic_tmp=0.0;
    wealth_of_exporters_tmp=0.0;income_of_workers_tmp=0.0;income_of_domestic_tmp=0.0;
    income_of_exporters_tmp=0.0;cons_of_workers_tmp=0.0;cons_of_domestic_tmp=0.0;
    cons_of_exporters_tmp=0.0;mean_leverage_d_tmp=0.0;mean_leverage_x_tmp=0.0;mean_leverage_tmp=0.0;
    sd_leverage_d_tmp=0.0;sd_leverage_x_tmp=0.0;sd_leverage_tmp=0.0;corr_MRPK_lev_d_tmp=0.0;corr_MRPK_lev_x_tmp=0.0;corr_MRPK_lev_tmp=0.0;
    exit_domestic_sum_tmp=0.0;exit_exporting_to_work_sum_tmp=0.0;mean_growth_rev=0.0;sd_growth_rev=0.0;mean_growth_k=0.0;sd_growth_k=0.0;autocorr_rev=0.0;fraction_zombie_exporter_tmp=0.0;
    # Solving for the aggregates - stationary distribution:
    if exitflag_tmp>0
        println("Nonstationarity")
    else
        #Continue only if we have convergence in VFI
        (Profits_d_fine,Profits_x_fine,output_dx_fine,output_xx_fine,output_d_fine,labor_d_fine,labor_x_fine
            ,k_choice_d_fine,k_choice_x_fine,price_dx_fine,price_d_fine,price_xx_fine,rev_d_fine,rev_dx_fine,rev_xx_fine) = true_profit_frictionless(price_final_prev,R,W_loc,s_fine
            ,constant_x,constant_d,constant_z_fine,constant_z_bar_fine,θ_loc,τ_loc,σ,α₁_eff,α₂_eff,α₁,α₂,ones_tmp_fine);
        (income_mat_fine,max_x_fine) = income_creator(W_loc,ones_tmp_fine,Profits_d_fine,Profits_x_fine,FM,FMₓ,country,
                r_loc,s_fine,iid_cost_value,F, Fₓ,size_iid_cost_val,income_mat_fine,ns_fine,a_min,a_max);
        (Q_trans_prime,cons_fine_local,a_prime_fine_local,
        future_occupation_fine_local) = Q_transition(coeff,
            ns_fine,ones_tmp_fine,size_iid_cost_val,iid_cost_prob,income_mat_fine,
            P_kron_fine,a_min,Phi_z_fine,β_loc,fspace_a,fspace_a_fine,
            price_final_current,max_x_fine,P_kron1,cons_fine_local,a_prime_fine_local,
            future_occupation_fine_local,Q_trans);
        exitflag_tmp = predict_irreducibility(future_occupation_fine_local,exitflag_tmp);
    end
    if exitflag_tmp==1 || exitflag_tmp==2
        labor_excess_demand_tmp_percent = -200.0 ;
        excess_demand_final_good_tmp_percent = 200.0
    elseif exitflag_tmp==3
        labor_excess_demand_tmp_percent = 200.0 ;
        excess_demand_final_good_tmp_percent = -200.0
    elseif exitflag_tmp==4
        labor_excess_demand_tmp_percent = 2000.0 ;
        excess_demand_final_good_tmp_percent = 2000.0
    else
        distr_current,exitflag_tmp = stationary_distribution(Q_trans_prime,L_loc,ns_fine,exitflag_tmp);
        if exitflag_tmp==5
            labor_excess_demand_tmp_percent = 1000.0 ;
            excess_demand_final_good_tmp_percent = 1000.0
        else
            current_distr_store_tmp = zeros(size(distr_current));
            # Calculate aggregates
            # Distribution of current workers:
            worker_past_dist = distr_current[(ns_fine *0 + 1):(ns_fine *1)];
            domestic_past_dist = distr_current[(ns_fine *1 + 1):(ns_fine *2)];
            exporter_past_dist = distr_current[(ns_fine *2 + 1):(ns_fine *3)];
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
            investment_tmp = δ * capital_demand_tmp;
            total_demand_final_good_tmp = total_consumption_tmp  + investment_tmp;
            excess_demand_final_good_tmp = total_demand_final_good_tmp - output_final[country,1];
            labor_excess_demand_tmp_percent = labor_excess_demand_tmp./L_loc;
            excess_demand_final_good_tmp_percent = excess_demand_final_good_tmp/(output_final[country,1] + total_demand_final_good_tmp);
            # Check if the upper limit has been reached:
            #upper_lim = check_upper_limit(worker_past_dist,domestic_past_dist,exporter_past_dist,
            #ns_fine,agrid_fine);

            # Do additional calculations/statistics
            mean_MRPK_d_tmp=@fastmath sum(current_domestic.*log.((R*ones_tmp_fine))./domestic_pop_tmp);
            mean_MRPK_x_tmp=@fastmath sum(current_exporter.*log.((R*ones_tmp_fine))./exporter_pop_tmp);
            mean_MRPK_tmp = (domestic_pop_tmp* mean_MRPK_d_tmp +
             exporter_pop_tmp *mean_MRPK_x_tmp)/(
             domestic_pop_tmp + exporter_pop_tmp) ;
            sd_MRPK_d_tmp =@fastmath  (sum(current_domestic.*(log.(( R*ones_tmp_fine))).^2/domestic_pop_tmp) -
             mean_MRPK_d_tmp^2)^(1/2);
            sd_MRPK_x_tmp =@fastmath  (sum(current_exporter.*(log.((R*ones_tmp_fine))).^2/sum(current_exporter) )- mean_MRPK_x_tmp.^2)^(1/2);
            sd_MRPK_tmp =@fastmath ((sum(current_domestic.*(log.((R*ones_tmp_fine))).^2 )+ sum(current_exporter.*(log.(( R*ones_tmp_fine))).^2))/
             (domestic_pop_tmp + exporter_pop_tmp)
             - mean_MRPK_tmp^2)^(1/2);
            mean_logProd_d_tmp = @fastmath sum((current_domestic.*log.(z_tilde_fine))/domestic_pop_tmp);
            mean_logProd_x_tmp = @fastmath sum((current_exporter.*log.(z_tilde_fine))/exporter_pop_tmp);
            mean_logProd_tmp =@fastmath sum((current_exporter + current_domestic).*log.(z_tilde_fine))/(exporter_pop_tmp + domestic_pop_tmp);

            sd_logProd_d_tmp =@fastmath (sum((current_domestic.*(log.(z_tilde_fine)).^2)/domestic_pop_tmp) - mean_logProd_d_tmp^2)^(1/2);
            sd_logProd_x_tmp =@fastmath (sum((current_exporter.*(log.(z_tilde_fine)).^2)/exporter_pop_tmp) - mean_logProd_x_tmp^2)^(1/2);
            sd_logProd_tmp =@fastmath (sum((current_exporter + current_domestic).*(log.(z_tilde_fine)).^2)/(exporter_pop_tmp + domestic_pop_tmp)
            - mean_logProd_tmp^2)^(1/2);

            cov_TFPR_z_d_tmp = @fastmath sum(current_domestic.*(α₁_eff* log.(z_tilde_fine).*log.(
            R*ones_tmp_fine)))/domestic_pop_tmp- α₁_eff*mean_MRPK_d_tmp * mean_logProd_d_tmp;
            cov_TFPR_z_x_tmp = @fastmath sum(current_exporter.*(α₁_eff* log.(z_tilde_fine).*log.(
            R*ones_tmp_fine)))/exporter_pop_tmp - α₁_eff*mean_MRPK_x_tmp * mean_logProd_x_tmp;
            cov_TFPR_z_tmp =@fastmath ((sum( current_exporter.*(log.( R*ones_tmp_fine) .* α₁_eff.* log.(z_tilde_fine) )
            ) + sum( current_domestic.*(log.(( R*ones_tmp_fine)).* α₁_eff.* log.(z_tilde_fine) )))
            )/ (domestic_pop_tmp + exporter_pop_tmp) - α₁_eff*mean_MRPK_tmp * mean_logProd_tmp;
            corr_TFPR_z_d_tmp = cov_TFPR_z_d_tmp / (α₁_eff* sd_MRPK_d_tmp* sd_logProd_d_tmp );
            corr_TFPR_z_x_tmp = cov_TFPR_z_x_tmp / (α₁_eff* sd_MRPK_x_tmp* sd_logProd_x_tmp );
            corr_TFPR_z_tmp = cov_TFPR_z_tmp / (α₁_eff* sd_MRPK_tmp* sd_logProd_tmp );

        # New stuff. Leverage with TFPR, as both are measurable
            mean_leverage_d_tmp = sum((current_domestic.*((k_choice_d_fine - s_fine[:,1])./k_choice_d_fine) )/domestic_pop_tmp);
            mean_leverage_x_tmp = sum((current_exporter.*((k_choice_x_fine - s_fine[:,1])./k_choice_x_fine) )/exporter_pop_tmp);
            mean_leverage_tmp =(domestic_pop_tmp* mean_leverage_d_tmp +
             exporter_pop_tmp *mean_leverage_x_tmp)/(
             domestic_pop_tmp + exporter_pop_tmp) ;

            sd_leverage_d_tmp =    (sum(current_domestic.*((k_choice_d_fine - s_fine[:,1])./k_choice_d_fine).^2)/domestic_pop_tmp - mean_leverage_d_tmp^2)^(1/2);
            sd_leverage_x_tmp =     (sum(current_exporter.*((k_choice_x_fine - s_fine[:,1])./k_choice_x_fine).^2)/exporter_pop_tmp - mean_leverage_x_tmp^2)^(1/2);
            sd_leverage_tmp =  ((sum((current_domestic.*((k_choice_d_fine - s_fine[:,1])./k_choice_d_fine).^2)) +
            sum(current_exporter.*((k_choice_x_fine - s_fine[:,1])./k_choice_x_fine).^2))/(domestic_pop_tmp + exporter_pop_tmp) - mean_leverage_tmp^2)^(1/2);

            cov_MRPK_lev_d_tmp = @fastmath sum(current_domestic.*(((k_choice_d_fine - s_fine[:,1])./k_choice_d_fine).*log.((R*ones_tmp_fine))))/domestic_pop_tmp- mean_MRPK_d_tmp * mean_leverage_d_tmp;
            cov_MRPK_lev_x_tmp = @fastmath sum(current_exporter.*(((k_choice_x_fine - s_fine[:,1])./k_choice_x_fine).*log.((R*ones_tmp_fine))))/exporter_pop_tmp- mean_MRPK_x_tmp * mean_leverage_x_tmp;
            cov_MRPK_lev_tmp = @fastmath (sum(current_exporter.*(((k_choice_x_fine - s_fine[:,1])./k_choice_x_fine).*log.((R*ones_tmp_fine)))) + sum(current_domestic.*(((k_choice_d_fine - s_fine[:,1]
            )./k_choice_d_fine).*log.((R*ones_tmp_fine)))))/(domestic_pop_tmp + exporter_pop_tmp) - mean_MRPK_tmp * mean_leverage_tmp;

            corr_MRPK_lev_d_tmp = cov_MRPK_lev_d_tmp / (sd_MRPK_d_tmp* sd_leverage_d_tmp );
            corr_MRPK_lev_x_tmp = cov_MRPK_lev_x_tmp / (sd_MRPK_x_tmp* sd_leverage_x_tmp );
            corr_MRPK_lev_tmp = cov_MRPK_lev_tmp / (sd_MRPK_tmp* sd_leverage_tmp );

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
            D_d_denom_tmp =  sum(current_domestic.*(z_tilde_fine.^σ .*(R*ones_tmp_fine).^((α₂_eff-1)*σ)));
            D_d_tmp =  sum(current_domestic.*(z_tilde_fine.^σ .*(R*ones_tmp_fine).^(-α₁_eff*σ)));
            #D_d_denom_tmp = D_d_K_tmp.^α₁_eff*D_d_tmp^α₂_eff;
            D_x_denom_tmp =  sum(current_exporter.*(z_tilde_fine.^σ .*(R*ones_tmp_fine).^((α₂_eff-1)*σ)));
            D_x_tmp =  sum(current_exporter.*(z_tilde_fine.^σ .*(R*ones_tmp_fine).^(-α₁_eff*σ)));
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
            current_distr_cons_sorted = distr_current[cons_sort_index]/L_loc;
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

            #Surviving firms
            sum_prob_survive = sum(Q_trans[(ns_fine+1):end,(ns_fine+1):end],dims = 2);
            Q_trans_firm = Q_trans[(ns_fine+1):end,(ns_fine+1):end]
            surviving_index = findall(x->x>0,sum_prob_survive);
            surviving_index = first.(Tuple.(surviving_index));
            sum_prob_survive_nonzero = sum_prob_survive[surviving_index];
            Q_trans_survive_firm= Q_trans_firm[surviving_index,:]./sum_prob_survive_nonzero;
            dropzeros!(Q_trans_survive_firm);
            current_firm_distr = current_distr_store_tmp[(ns_fine *1 + 1):(ns_fine *3)] ;
            # Normalize to 1 on surviving firms
            current_firm_distr_survive = current_firm_distr[surviving_index];
            current_firm_distr_survive = current_firm_distr_survive./sum(current_firm_distr_survive);
            current_firm_distr_survive_rep = repeat(current_firm_distr_survive,1,2*ns_fine);
            distribution_of_firms = current_firm_distr_survive_rep.*Q_trans_survive_firm;
            # growth rate of revenue:
            do_do_growth = zeros(ns_fine,ns_fine);
            do_ex_growth = zeros(ns_fine,ns_fine);
            ex_do_growth = zeros(ns_fine,ns_fine);
            ex_ex_growth = zeros(ns_fine,ns_fine);

            for ii= 1:ns_fine
                do_do_growth[ii,:] = log.(rev_d_fine) - log.(rev_d_fine[ii] *ones_tmp_fine) ;
                do_ex_growth[ii,:] = log.(rev_dx_fine+rev_xx_fine) - log.(rev_d_fine[ii] *ones_tmp_fine) ;
                ex_do_growth[ii,:] = log.(rev_d_fine) - log.( (rev_dx_fine[ii]+rev_xx_fine[ii]) *ones_tmp_fine);
                ex_ex_growth[ii,:] = log.(rev_dx_fine+rev_xx_fine) -  log.((rev_dx_fine[ii]+rev_xx_fine[ii]) *ones_tmp_fine);
            end

            growth_rate_rev_mat = zeros(2*ns_fine,2*ns_fine);
            growth_rate_rev_mat[1:ns_fine,1:ns_fine] = do_do_growth;
            growth_rate_rev_mat[ns_fine+1:end,1:ns_fine] = ex_do_growth;
            growth_rate_rev_mat[1:ns_fine,ns_fine+1:end] = do_ex_growth;
            growth_rate_rev_mat[ns_fine+1:end,ns_fine+1:end] = ex_ex_growth;
            growth_rate_rev_mat_survive = growth_rate_rev_mat[surviving_index,:];
            mean_growth_rev = sum(growth_rate_rev_mat_survive.*distribution_of_firms);
            sd_growth_rev = (sum(growth_rate_rev_mat_survive.^2 .*distribution_of_firms) - mean_growth_rev^2)^(1/2);
            # growth rate of capital:
            do_do_growth = zeros(ns_fine,ns_fine);
            do_ex_growth = zeros(ns_fine,ns_fine);
            ex_do_growth = zeros(ns_fine,ns_fine);
            ex_ex_growth = zeros(ns_fine,ns_fine);

            for ii= 1:ns_fine
                do_do_growth[ii,:] = log.(k_choice_d_fine) - log.(k_choice_d_fine[ii] *ones_tmp_fine);
                do_ex_growth[ii,:] = log.(k_choice_x_fine) - log.(k_choice_d_fine[ii] *ones_tmp_fine);
                ex_do_growth[ii,:] = log.(k_choice_d_fine) - log.(k_choice_x_fine[ii] *ones_tmp_fine);
                ex_ex_growth[ii,:] = log.(k_choice_x_fine) -  log.(k_choice_x_fine[ii]*ones_tmp_fine);
            end

            growth_rate_k_mat = zeros(2*ns_fine,2*ns_fine);
            growth_rate_k_mat[1:ns_fine,1:ns_fine] = do_do_growth;
            growth_rate_k_mat[ns_fine+1:end,1:ns_fine] = ex_do_growth;
            growth_rate_k_mat[1:ns_fine,ns_fine+1:end] = do_ex_growth;
            growth_rate_k_mat[ns_fine+1:end,ns_fine+1:end] = ex_ex_growth;

            growth_rate_k_mat_survive = growth_rate_k_mat[surviving_index,:];
            mean_growth_k = sum(growth_rate_k_mat_survive.*distribution_of_firms);
            sd_growth_k = (sum(growth_rate_k_mat_survive.^2 .*distribution_of_firms) - mean_growth_k^2)^(1/2);

            # autocorr of revenue:
            do_do_mul = zeros(ns_fine,ns_fine);
            do_ex_mul = zeros(ns_fine,ns_fine);
            ex_do_mul = zeros(ns_fine,ns_fine);
            ex_ex_mul = zeros(ns_fine,ns_fine);

            for ii= 1:ns_fine
                do_do_mul[ii,:] = log.(rev_d_fine)*log(rev_d_fine[ii]) ;
                do_ex_mul[ii,:] = log.(rev_dx_fine+rev_xx_fine)*log(rev_d_fine[ii]) ;
                ex_do_mul[ii,:] = log.(rev_d_fine)*log(rev_dx_fine[ii]+rev_xx_fine[ii]);
                ex_ex_mul[ii,:] = log.(rev_dx_fine+rev_xx_fine)*log(rev_dx_fine[ii]+rev_xx_fine[ii]);
            end

            mul_rate_rev_mat = zeros(2*ns_fine,2*ns_fine);
            mul_rate_rev_mat[1:ns_fine,1:ns_fine] = do_do_mul;
            mul_rate_rev_mat[ns_fine+1:end,1:ns_fine] = ex_do_mul;
            mul_rate_rev_mat[1:ns_fine,ns_fine+1:end] = do_ex_mul;
            mul_rate_rev_mat[ns_fine+1:end,ns_fine+1:end] = ex_ex_mul;
            mul_rate_rev_mat_survive = mul_rate_rev_mat[surviving_index,:];
            mean_mul_rev = sum(mul_rate_rev_mat_survive.*distribution_of_firms);
            revenue_combined = zeros(2*ns_fine);
            revenue_combined[1:ns_fine] = log.(rev_d_fine);
            revenue_combined[ns_fine+1:end] = log.(rev_dx_fine+rev_xx_fine);
            distr_tmr_firm = transpose(sum(distribution_of_firms,dims = 1));
            mean_log_rev = sum(revenue_combined.*distr_tmr_firm);
            sd_log_rev = (sum(revenue_combined.^2 .*distr_tmr_firm) - mean_log_rev.^2)^(1/2);
            revenue_combined_survive = revenue_combined[surviving_index,:];
            mean_log_rev_survive = sum(revenue_combined_survive .*current_firm_distr_survive);
            sd_log_rev_survive = (sum(revenue_combined_survive.^2 .*current_firm_distr_survive) - mean_log_rev_survive.^2)^(1/2);
            autocorr_rev = (mean_mul_rev - mean_log_rev_survive*mean_log_rev)/sd_log_rev_survive/sd_log_rev;
            # Zombies - experimental
            making_losses = (income_1_fine .>income_3_fine);
            exit_zombie = sum(((future_occupation_fine_local[:,3,1] .!=3.0) .*making_losses).*current_exporter)/exporter_pop_tmp;
            fraction_zombie_exporter_tmp = sum(making_losses.*current_exporter)/exporter_pop_tmp - 3 * exit_zombie
        end
    end

    return (labor_excess_demand_tmp_percent,excess_demand_final_good_tmp_percent,
    exitflag_tmp,export_price_sum_tmp, domestic_price_sum_tmp,NFA_tmp,asset_supply_tmp,
    asset_demand_tmp,a_prime_fine_local,future_occupation_fine_local,cons_fine_local,
    k_choice_x_fine,labor_x_fine,k_choice_d_fine,labor_d_fine,Profits_x_fine,Profits_d_fine,
    rev_d_fine,rev_dx_fine,rev_xx_fine,current_distr_store_tmp,distr_current,coeff,domestic_prod_tmp,
    export_prod_tmp,export_value_tmp,exit_share_to_domestic_tmp,entry_share_to_domestic_tmp,
    exit_share_to_worker_tmp,entry_share_to_exporter_tmp,exporter_pop_tmp,domestic_pop_tmp,
    worker_pop_tmp,L_x_tmp,L_d_tmp,total_entry_cost_tmp,total_incumbents_cost_tmp,total_exit_cost_tmp,
    labor_demand_tmp,labor_excess_demand_tmp,K_x_tmp,K_d_tmp,capital_demand_tmp,worker_bond_holding_sum_tmp,
    domestic_bond_holding_sum_tmp,exporter_bond_holding_sum_tmp,domestic_firm_debt_tmp,exporter_firm_debt_tmp,
    total_consumption_tmp,Aggr_bank_cost_tmp,x_k_tmp,x_l_tmp,investment_tmp,total_demand_final_good_tmp,
    excess_demand_final_good_tmp,mean_MRPK_d_tmp,mean_MRPK_x_tmp,mean_MRPK_tmp,sd_MRPK_d_tmp, sd_MRPK_x_tmp,sd_MRPK_tmp,
    mean_logProd_d_tmp,mean_logProd_x_tmp,mean_logProd_tmp,sd_logProd_d_tmp,sd_logProd_x_tmp,sd_logProd_tmp,
    cov_TFPR_z_d_tmp,cov_TFPR_z_x_tmp,cov_TFPR_z_tmp,corr_TFPR_z_d_tmp,corr_TFPR_z_x_tmp,corr_TFPR_z_tmp,
    avg_Pi_x_tmp,avg_Pi_d_tmp,D_d_denom_tmp,D_x_denom_tmp,D_x_tmp,D_d_tmp,D_d_denom_eff_tmp,D_x_denom_eff_tmp,
    D_x_eff_tmp,D_d_eff_tmp,mean_cons_tmp,sd_cons_tmp,mean_income_tmp,sd_income_tmp,mean_wealth_tmp,sd_wealth_tmp,
    p10_wealth_tmp,p50_wealth_tmp,p90_wealth_tmp,p99_wealth_tmp,p10_income_tmp,p50_income_tmp,p90_income_tmp,p99_income_tmp,p10_cons_tmp,
    p50_cons_tmp,p90_cons_tmp,p99_cons_tmp,wealth_of_workers_tmp,wealth_of_domestic_tmp,
    wealth_of_exporters_tmp,income_of_workers_tmp,income_of_domestic_tmp,
    income_of_exporters_tmp,cons_of_workers_tmp,cons_of_domestic_tmp,
    cons_of_exporters_tmp,mean_leverage_d_tmp,mean_leverage_x_tmp,mean_leverage_tmp,
    sd_leverage_d_tmp,sd_leverage_x_tmp,sd_leverage_tmp,corr_MRPK_lev_d_tmp,corr_MRPK_lev_x_tmp,corr_MRPK_lev_tmp,
    exit_domestic_sum_tmp,exit_exporting_to_work_sum_tmp,mean_growth_rev,sd_growth_rev,mean_growth_k,sd_growth_k,autocorr_rev,fraction_zombie_exporter_tmp)
end
function Residual_stst_detailed_frictionless(prices::Array{Float64,1},parameters_tmp::Parameter_type)
    # Extract the necessary local parameters:
    (β,α,δ,θ,α₁,α₂,σ,α₁_eff,α₂_eff,ω,L,FM,FMₓ,F,Fₓ,
        Exit,Exitₓ,iid_cost_value,iid_cost_prob,size_iid_cost_val,country_no,τ,
        a_min,a_max,fspace_a,fspace_a_fine,agrid_fine,banking_cost,
        bounds,Country_spec_p,s_cell,ns_cell,s_fine_cell,ns_fine_cell,Phi_z_cell,
        Phi_z_fine_cell,Phi_cell,Phi_aug_cell,P_kron_cell,P_kron1_cell,P_kron_fine_cell,
        exp_egrid_cell,ns_tmp,ns_tmp_fine,openness) = local_parameters(parameters_tmp)
    # Prices:
    (W,r,output_final,price_final,r,residual,price_check_tmp) = price_reshaper(prices,openness,country_no,bounds)
    # Check prices to be reasonable:
    if price_check_tmp>0
        println("Guess out of bounds")
        residual = 1000.0 * ones(8);
        return residual
    end
    # Initialization of local variables:
    (export_prod,export_price_sum,import_price,domestic_price_sum,
    price_final_actual,NFA,asset_supply,asset_demand,exitflag) = local_var_creator_minimal(country_no);
    (a_prime_fine,future_occupation_fine,cons_fine,rev_d_fine_store,
    rev_dx_fine_store,rev_xx_fine_store,coeff_final,k_x_fine,k_d_fine,
    l_x_fine,l_d_fine,current_distr_store,distr_current,exitflag,domestic_prod,exporter_pop,
    domestic_pop,worker_pop,entry_share_to_domestic,exit_share_to_domestic,
    entry_share_to_exporter,exit_share_to_worker,total_entry_cost,
    total_incumbents_cost,total_exit_cost,labor_excess_demand,labor_demand,
    total_consumption,capital_demand,capital_next,
    total_demand_final_good,excess_demand_final_good,GDP,nomGDP,PPI,TFP,
    TFP_within,TFP_across,TFP_second_best,TOT,total_production,
    domestic_price_sum,export_value,investment,
    export_value,L_x,L_d,K_x,K_d,x_k,x_l,D_d,D_x,D_d_denom,D_x_denom,
    D_d_eff,D_x_eff,D_d_denom_eff,D_x_denom_eff,K_d_ratio,K_x_ratio,
    L_d_ratio,L_x_ratio,K_x_ratio_eff,K_d_ratio_eff,L_d_ratio_eff,
    L_x_ratio_eff,nomGDP_d,nomGDP_xd,nomGDP_xx,RGDP_d,RGDP_xd,RGDP_xx,RGDP,
    mean_MRPK_d,mean_MRPK_x,sd_MRPK_d,sd_MRPK_x,mean_MRPK,sd_MRPK,
    mean_MRPL_d,mean_MRPL_x,sd_MRPL_d,sd_MRPL_x,mean_MRPL,sd_MRPL,
    cov_TFPR_z_d,cov_TFPR_z_x,cov_TFPR_z,corr_TFPR_z_d,corr_TFPR_z_x,
    corr_TFPR_z,mean_logProd_d,mean_logProd_x,mean_logProd,sd_logProd_d,
    sd_logProd_x,sd_logProd,import_share,worker_bond_holding,
    domestic_bond_holding,exporter_bond_holding,domestic_firm_debt,
    exporter_firm_debt,avg_Pi_x,avg_Pi_d,Profit_fine_x,
    Profit_fine_d,mean_cons,sd_cons,mean_income,sd_income,mean_wealth,
    sd_wealth,wealth_of_workers,wealth_of_domestic,wealth_of_exporters,
    income_of_workers,income_of_domestic,income_of_exporters,
    cons_of_workers,cons_of_domestic,cons_of_exporters,GINI_wealth,
    GINI_income,GINI_cons,p10_wealth,p50_wealth,p90_wealth,p99_wealth,p10_income,p50_income,p90_income,
    p99_income,p10_cons,p50_cons,p90_cons,p99_cons,Aggr_bank_cost,mean_leverage_d,mean_leverage_x,mean_leverage,
    sd_leverage_d,sd_leverage_x,sd_leverage,corr_MRPK_lev_d,corr_MRPK_lev_x,corr_MRPK_lev,
    exit_domestic_to_work_sum,exit_exporting_to_work_sum,mean_growth_rev,sd_growth_rev,mean_growth_k,sd_growth_k,
    autocorr_rev,fraction_zombie_exporter,k_opt_d_fine,k_opt_x_fine)= local_var_creator_detailed(country_no,ns_tmp,ns_tmp_fine,size_iid_cost_val);
    @sync @distributed for country = 1:country_no#
        tmp_vector  = country_residual_detailed_frictionless(country,s_cell,ns_cell,s_fine_cell,ns_fine_cell,
        Phi_z_cell,Phi_z_fine_cell,Phi_cell,Phi_aug_cell,P_kron_cell,P_kron1_cell,
        P_kron_fine_cell,σ,size_iid_cost_val,price_final,W, r, output_final,θ,L,τ,
        country_no,banking_cost,δ,ω,β,α₁_eff,α₂_eff,α₁,α₂,FM,FMₓ,F,Fₓ,iid_cost_value,
        iid_cost_prob,Exit,Exitₓ,a_min,a_max,α,fspace_a,fspace_a_fine,openness,agrid_fine);
        residual[country,1] = tmp_vector[1];
        residual[2* country_no + country,1] = tmp_vector[2];
        exitflag[country,1] = tmp_vector[3];
        export_price_sum[country,1] = tmp_vector[4];
        domestic_price_sum[country,1] = tmp_vector[5];
        NFA[country,1] = tmp_vector[6];
        asset_supply[country,1] = tmp_vector[7];
        asset_demand[country,1] = tmp_vector[8];
        a_prime_fine[:,:,country,:] = tmp_vector[9];
        future_occupation_fine[:,:,country,:] = tmp_vector[10];
        cons_fine[:,:,country]  = tmp_vector[11];
        k_x_fine[:,country] = tmp_vector[12];
        l_x_fine[:,country] = tmp_vector[13];
        k_d_fine[:,country] = tmp_vector[14];
        l_d_fine[:,country] = tmp_vector[15];
        Profit_fine_x[:,country] = tmp_vector[16];
        Profit_fine_d[:,country] = tmp_vector[17];
        rev_d_fine_store[:,country] = tmp_vector[18];
        rev_dx_fine_store[:,country] = tmp_vector[19];
        rev_xx_fine_store[:,country] = tmp_vector[20];
        current_distr_store[:,country] = tmp_vector[21];
        distr_current[:,country] = tmp_vector[22];
        coeff_final[:,:,country] = tmp_vector[23];
        domestic_prod[country,1] = tmp_vector[24];
        export_prod[country,1] = tmp_vector[25];
        export_value[country,1] = tmp_vector[26];
        exit_share_to_domestic[country,1] = tmp_vector[27];
        entry_share_to_domestic[country,1] = tmp_vector[28];
        exit_share_to_worker[country,1] = tmp_vector[29];
        entry_share_to_exporter[country,1] = tmp_vector[30];
        exporter_pop[country,1] = tmp_vector[31];
        domestic_pop[country,1] = tmp_vector[32];
        worker_pop[country,1] = tmp_vector[33];
        L_x[country,1] = tmp_vector[34];
        L_d[country,1] = tmp_vector[35];
        total_entry_cost[country,1] = tmp_vector[36];
        total_incumbents_cost[country,1] = tmp_vector[37];
        total_exit_cost[country,1] = tmp_vector[38];
        labor_demand[country,1] = tmp_vector[39];
        labor_excess_demand[country,1] = tmp_vector[40];
        K_x[country,1] = tmp_vector[41];
        K_d[country,1] = tmp_vector[42];
        capital_demand[country,1] = tmp_vector[43];
        worker_bond_holding[country,1] = tmp_vector[44];
        domestic_bond_holding[country,1] = tmp_vector[45];
        exporter_bond_holding[country,1] = tmp_vector[46];
        domestic_firm_debt[country,1] = tmp_vector[47];
        exporter_firm_debt[country,1] = tmp_vector[48];
        total_consumption[country,1] = tmp_vector[49];
        Aggr_bank_cost[country,1] = tmp_vector[50];
        x_k[country,1] = tmp_vector[51];
        x_l[country,1] = tmp_vector[52];
        investment[country,1] = tmp_vector[53];
        total_demand_final_good[country,1] = tmp_vector[54];
        excess_demand_final_good[country,1] = tmp_vector[55];
        mean_MRPK_d[country,1] = tmp_vector[56];
        mean_MRPK_x[country,1 ] = tmp_vector[57];
        mean_MRPK[country,1] = tmp_vector[58];
        sd_MRPK_d[country,1] = tmp_vector[59];
        sd_MRPK_x[country,1] = tmp_vector[60];
        sd_MRPK[country,1] = tmp_vector[61];
        mean_logProd_d[country,1] = tmp_vector[62];
        mean_logProd_x[country,1] = tmp_vector[63];
        mean_logProd[country,1] = tmp_vector[64];
        sd_logProd_d[country,1] = tmp_vector[65];
        sd_logProd_x[country,1] = tmp_vector[66];
        sd_logProd[country,1] = tmp_vector[67];
        cov_TFPR_z_d[country,1] = tmp_vector[68];
        cov_TFPR_z_x[country,1] = tmp_vector[69];
        cov_TFPR_z[country,1] = tmp_vector[70];
        corr_TFPR_z_d[country,1] = tmp_vector[71];
        corr_TFPR_z_x[country,1] = tmp_vector[72];
        corr_TFPR_z[country,1] = tmp_vector[73];
        avg_Pi_x[country,1] = tmp_vector[74];
        avg_Pi_d[country,1] = tmp_vector[75];
        D_d_denom[country,1] = tmp_vector[76];
        D_x_denom[country,1] = tmp_vector[77];
        D_x[country,1] = tmp_vector[78];
        D_d[country,1] = tmp_vector[79];
        D_d_denom_eff[country,1] = tmp_vector[80];
        D_x_denom_eff[country,1] = tmp_vector[81];
        D_x_eff[country,1] = tmp_vector[82];
        D_d_eff[country,1] = tmp_vector[83];
        mean_cons[country,1] = tmp_vector[84];
        sd_cons[country,1] = tmp_vector[85];
        mean_income[country,1] = tmp_vector[86];
        sd_income[country,1] = tmp_vector[87];
        mean_wealth[country,1] = tmp_vector[88];
        sd_wealth[country,1] = tmp_vector[89];
        p10_wealth[country,1] = tmp_vector[90];
        p50_wealth[country,1] = tmp_vector[91];
        p90_wealth[country,1] = tmp_vector[92];
        p99_wealth[country,1] = tmp_vector[93];
        p10_income[country,1] = tmp_vector[94];
        p50_income[country,1] = tmp_vector[95];
        p90_income[country,1] = tmp_vector[96];
        p99_income[country,1] = tmp_vector[97];
        p10_cons[country,1] = tmp_vector[98];
        p50_cons[country,1] = tmp_vector[99];
        p90_cons[country,1] = tmp_vector[100];
        p99_cons[country,1] = tmp_vector[101];
        wealth_of_workers[country,1] = tmp_vector[102];
        wealth_of_domestic[country,1] = tmp_vector[103];
        wealth_of_exporters[country,1] = tmp_vector[104];
        income_of_workers[country,1] = tmp_vector[105];
        income_of_domestic[country,1] = tmp_vector[106];
        income_of_exporters[country,1] = tmp_vector[107];
        cons_of_workers[country,1] = tmp_vector[108];
        cons_of_domestic[country,1] = tmp_vector[109];
        cons_of_exporters[country,1] = tmp_vector[110];
        mean_leverage_d[country,1] = tmp_vector[111];
        mean_leverage_x[country,1] = tmp_vector[112];
        mean_leverage[country,1] = tmp_vector[113];
        sd_leverage_d[country,1] = tmp_vector[114];
        sd_leverage_x[country,1] = tmp_vector[115];
        sd_leverage[country,1] = tmp_vector[116];
        corr_MRPK_lev_d[country,1] = tmp_vector[117];
        corr_MRPK_lev_x[country,1] = tmp_vector[118];
        corr_MRPK_lev[country,1] = tmp_vector[119];
        exit_domestic_to_work_sum[country,1] = tmp_vector[120];
        exit_exporting_to_work_sum[country,1] = tmp_vector[121];
        mean_growth_rev[country,1] = tmp_vector[122];
        sd_growth_rev[country,1] = tmp_vector[123];
        mean_growth_k[country,1] = tmp_vector[124];
        sd_growth_k[country,1] = tmp_vector[125];
        autocorr_rev[country,1] = tmp_vector[126];
        fraction_zombie_exporter[country,1] = tmp_vector[127];
    end
    for country=1:country_no
        if exitflag[country,1]==1 || exitflag[country,1]==2# Too large supply of assets
            if openness == 1
                residual[7,1] =  100.0;
            else
                residual[3* country_no + country,1] = 100.0;
            end
        elseif exitflag[country,1]==3# Too little supply of assets?
            if openness == 1
                residual[7,1] = - 100.0;
            else
                residual[3* country_no + country,1] = - 100.0;
            end
        elseif exitflag[country,1]==4
            if openness == 1
                residual[7,1] =  1000.0;
            else
                residual[3* country_no + country,1] =  1000.0;
            end
        elseif exitflag[country,1]==5
            if openness == 1
                residual[7,1] = 500.0;
            else
                residual[3* country_no + country,1] =  500.0;
            end
        else
            import_curr = (sum(export_prod[:,1])-  export_prod[country,1])/(country_no - 1.0);
            total_production[country,1] = (ω *domestic_prod[country,1] + (1.0 -ω) *import_curr)^((σ)/(σ-1));
            import_price[country,1] = (sum(export_price_sum[:,1])-  export_price_sum[country,1])/(country_no - 1);
            price_final_actual[country,1] = min((ω^σ * domestic_price_sum[country,1] + (1.0 - ω)^σ  * import_price[country,1])^(1.0
                / (1.0 - σ)),10000);
            residual[country_no + country,1] = (price_final_actual[country,1] - price_final[country,1])./(
                price_final_actual[country,1] + price_final[country,1]);
            if openness == 0
                residual[3* country_no + country,1] = NFA[country,1]/(asset_demand[country,1] + (1 + banking_cost[country,1]).*asset_supply[country,1]);
            end
        end
        if openness == 1
            residual[7,1] = sum(NFA)/sum((asset_demand+ (ones(country_no) + banking_cost).*asset_supply));
        end
    end
    TFP_d = ω * (D_d.^(1 - α₂_eff)./D_d_denom.^α₁_eff);
    TFP_x = (D_x.^(1 - α₂_eff)./D_x_denom.^α₁_eff);

    #TFP_d = D_d./D_d_denom;
    #TFP_x = D_x./D_x_denom;
    TFP_d_efficient = ω * (D_d_eff.^(1 - α₂_eff)./D_d_denom_eff.^α₁_eff);
    TFP_x_efficient = (D_x_eff.^(1 - α₂_eff)./D_x_denom_eff.^α₁_eff);
    #TFP_d_efficient = D_d_eff./D_d_denom_eff;
    #TFP_x_efficient = D_x_eff./D_x_denom_eff;
    # Other definitions for TFP and some other useful measures
    for country = 1:country_no
        total_foreign_demand = (sum(total_production) - total_production[country,1])/(country_no - 1);
        #total_foreign_demand = (sum(output_final) - output_final[country,1])/(country_no - 1);
        avg_foreign_price = (sum(price_final_actual) - price_final_actual[country,1])/(country_no - 1);
        nomGDP[country,1] = price_final_actual[country,1] * total_production[country,1];
        GDP[country,1] = total_production[country,1];
        import_value_tmp = (sum(export_value[:,1])-  export_value[country,1]);
        export_value_tradeoff = import_value_tmp/export_value[country,1];
        import_share[country,1] = import_value_tmp/(country_no - 1)/nomGDP[country,1];
        #constant_d_tmp =   ω * price_final_actual[country,1] * output_final[country,1]^(1/σ);
        constant_d_tmp =   ω * price_final_actual[country,1] * total_production[country,1]^(1/σ);
        constant_x_tmp = (1 - ω)/ (1 + τ[country,1]) * avg_foreign_price * total_foreign_demand^(1/σ);
        TOT[country,1] =  (ω *constant_d_tmp^(σ - 1) + (1 - ω) *constant_x_tmp^(σ - 1)*(
         export_value_tradeoff *price_final_actual[country,1]/ avg_foreign_price)^((σ -1)/σ))/(
        constant_d_tmp^σ + (1 +τ[country])* constant_x_tmp^σ
        )^((σ -1)/σ);
        TFP_x[country,1] = TFP_x[country,1] *TOT[country,1];
        TFP_x_efficient[country,1] = TFP_x_efficient[country,1] * TOT[country,1];
        D_d_K_tmp = D_d_denom[country,1];
        D_x_K_tmp = D_x_denom[country,1];
        D_d_K_eff_tmp = D_d_denom_eff[country,1];
        D_x_K_eff_tmp = D_x_denom_eff[country,1];
        K_d_ratio[country,1] = constant_d_tmp^σ * D_d_K_tmp/(constant_d_tmp^σ * D_d_K_tmp + (constant_d_tmp^σ
         +(1 +τ[country])*constant_x_tmp^σ)* D_x_K_tmp);
        K_x_ratio[country,1] = (constant_d_tmp^σ
         +(1 +τ[country])*constant_x_tmp^σ)* D_x_K_tmp/(constant_d_tmp^σ * D_d_K_tmp + (constant_d_tmp^σ
          +(1 +τ[country])*constant_x_tmp^σ)* D_x_K_tmp);
        L_d_ratio[country,1] = constant_d_tmp^σ * D_d[country,1]/(constant_d_tmp^σ * D_d[country,1] + (constant_d_tmp^σ
            +(1 +τ[country])*constant_x_tmp^σ)* D_x[country,1]);
        L_x_ratio[country,1] = (constant_d_tmp^σ
            +(1 +τ[country])*constant_x_tmp^σ)* D_x[country,1]/(constant_d_tmp^σ * D_d[country,1] + (constant_d_tmp^σ
            +(1 +τ[country])*constant_x_tmp^σ)* D_x[country,1]);

        K_d_ratio_eff[country,1] = constant_d_tmp^σ * D_d_K_eff_tmp/(constant_d_tmp^σ * D_d_K_eff_tmp + (constant_d_tmp^σ
         +(1 +τ[country])*constant_x_tmp^σ)* D_x_K_eff_tmp);
        K_x_ratio_eff[country,1] = (constant_d_tmp^σ
         +(1 +τ[country])*constant_x_tmp^σ)* D_x_K_eff_tmp/(constant_d_tmp^σ * D_d_K_eff_tmp + (constant_d_tmp^σ
          +(1 +τ[country])*constant_x_tmp^σ)* D_x_K_eff_tmp);
        L_d_ratio_eff[country,1] = constant_d_tmp^σ * D_d_eff[country,1]/(constant_d_tmp^σ * D_d_eff[country,1] + (constant_d_tmp^σ
            +(1 +τ[country])*constant_x_tmp^σ)* D_x_eff[country,1]);
        L_x_ratio_eff[country,1] = (constant_d_tmp^σ
            +(1 +τ[country])*constant_x_tmp^σ)* D_x_eff[country,1]/(constant_d_tmp^σ * D_d_eff[country,1] + (constant_d_tmp^σ
            +(1 +τ[country])*constant_x_tmp^σ)* D_x_eff[country,1]);
        TFP[country,1] = (TFP_d[country,1] * K_d_ratio[country,1]^α₁_eff * L_d_ratio[country,1]^α₂_eff
            + TFP_x[country,1] * K_x_ratio[country,1]^α₁_eff * L_x_ratio[country,1]^α₂_eff)^(σ/(σ-1));
        TFP_within[country,1] = (TFP_d_efficient[country,1] * K_d_ratio[country,1]^α₁_eff * L_d_ratio[country,1]^α₂_eff
            + TFP_x_efficient[country,1] * K_x_ratio[country,1]^α₁_eff * L_x_ratio[country,1]^α₂_eff)^(σ/(σ-1));
        TFP_across[country,1] = (TFP_d[country,1] * K_d_ratio_eff[country,1]^α₁_eff * L_d_ratio_eff[country,1]^α₂_eff
            + TFP_x[country,1] * K_x_ratio_eff[country,1]^α₁_eff * L_x_ratio_eff[country,1]^α₂_eff)^(σ/(σ-1)); # this doesnt work
        TFP_second_best[country,1] =(TFP_d_efficient[country,1] * K_d_ratio_eff[country,1]^α₁_eff * L_d_ratio_eff[country,1]^α₂_eff
            + TFP_x_efficient[country,1] * K_x_ratio_eff[country,1]^α₁_eff * L_x_ratio_eff[country,1]^α₂_eff)^(σ/(σ-1));
    end
    Misallocation_within_d =ones(country_no) -  TFP_d./TFP_d_efficient;
    Misallocation_within_x =ones(country_no) - TFP_x./TFP_x_efficient;
    replace!(residual, NaN=>100.0);
    return (residual,a_prime_fine,future_occupation_fine,cons_fine,rev_d_fine_store,
    rev_dx_fine_store,rev_xx_fine_store,coeff_final,k_x_fine,k_d_fine,
    l_x_fine,l_d_fine,current_distr_store,distr_current,exitflag,domestic_prod,exporter_pop,
    domestic_pop,worker_pop,entry_share_to_domestic,exit_share_to_domestic,
    entry_share_to_exporter,exit_share_to_worker,total_entry_cost,
    total_incumbents_cost,total_exit_cost,labor_excess_demand,labor_demand,
    total_consumption,capital_demand,capital_next,
    total_demand_final_good,excess_demand_final_good,GDP,nomGDP,PPI,TFP,
    TFP_within,TFP_across,TFP_second_best,TOT,TFP_d,
    TFP_d_efficient,TFP_x,TFP_x_efficient,total_production,
    domestic_price_sum,export_value,investment,
    export_value,L_x,L_d,K_x,K_d,x_k,x_l,D_d,D_x,D_d_denom,D_x_denom,
    D_d_eff,D_x_eff,D_d_denom_eff,D_x_denom_eff,K_d_ratio,K_x_ratio,
    L_d_ratio,L_x_ratio,K_x_ratio_eff,K_d_ratio_eff,L_d_ratio_eff,
    L_x_ratio_eff,nomGDP_d,nomGDP_xd,nomGDP_xx,RGDP_d,RGDP_xd,RGDP_xx,RGDP,
    mean_MRPK_d,mean_MRPK_x,sd_MRPK_d,sd_MRPK_x,mean_MRPK,sd_MRPK,
    mean_MRPL_d,mean_MRPL_x,sd_MRPL_d,sd_MRPL_x,mean_MRPL,sd_MRPL,
    cov_TFPR_z_d,cov_TFPR_z_x,cov_TFPR_z,corr_TFPR_z_d,corr_TFPR_z_x,
    corr_TFPR_z,mean_logProd_d,mean_logProd_x,mean_logProd,sd_logProd_d,
    sd_logProd_x,sd_logProd,import_share,worker_bond_holding,
    domestic_bond_holding,exporter_bond_holding,domestic_firm_debt,
    exporter_firm_debt,avg_Pi_x,avg_Pi_d,Profit_fine_x,
    Profit_fine_d,mean_cons,sd_cons,mean_income,sd_income,mean_wealth,
    sd_wealth,wealth_of_workers,wealth_of_domestic,wealth_of_exporters,
    income_of_workers,income_of_domestic,income_of_exporters,
    cons_of_workers,cons_of_domestic,cons_of_exporters,GINI_wealth,
    GINI_income,GINI_cons,p10_wealth,p50_wealth,p90_wealth,p99_wealth,p10_income,p50_income,p90_income,
    p99_income,p10_cons,p50_cons,p90_cons,p99_cons,Aggr_bank_cost,mean_leverage_d,mean_leverage_x,mean_leverage,
    sd_leverage_d,sd_leverage_x,sd_leverage,corr_MRPK_lev_d,corr_MRPK_lev_x,corr_MRPK_lev,exit_domestic_to_work_sum,
    exit_exporting_to_work_sum,mean_growth_rev,sd_growth_rev,mean_growth_k,sd_growth_k,autocorr_rev,fraction_zombie_exporter,
    Misallocation_within_d,Misallocation_within_x)
end

Simple_function_space = fundefn(:spli, 10, 0.1, 1.0,1)
Baseline_parameter = Parameter_type(2,1/3,1,0.06,4,0,0,0,0,[1,4],0,[0.84,0.93],
0.5,[0.6,0.86],[0,0],[0,0],[0.9,0.9],[0.06,0.06],1,
    [0,0],[0,0],[0,0],[0,0],[0,0],[0,0],1.0,2.0,[0.99,1.01],[0.5,0.5],
[50,36],[200,36],[0,0],[0,0], 0.0005,10.0,1#[50,20],[200,20],[0,0],[0,0], 0.0005,10.0,1#
    ,Simple_function_space,Simple_function_space,
    [0.00,0.00],0,zeros(4,2));
Baseline_parameter.α₁=Baseline_parameter.α * Baseline_parameter.η;
Baseline_parameter.α₂=(1 - Baseline_parameter.α) * Baseline_parameter.η;
Baseline_parameter.α₁_eff = Baseline_parameter.α₁ * (Baseline_parameter.σ - 1)/Baseline_parameter.σ;
Baseline_parameter.α₂_eff = Baseline_parameter.α₂ * (Baseline_parameter.σ - 1)/Baseline_parameter.σ;

# Unused but implemented:
Baseline_parameter.FM = [0.0,0.0]*Baseline_parameter.κ;
Baseline_parameter.FMₓ = [0.0,0.0]*Baseline_parameter.κ;
Baseline_parameter.Exit = [0.0,0.0]*Baseline_parameter.κ;
Baseline_parameter.Exitₓ = [0.0,0.0]*Baseline_parameter.κ;
Baseline_parameter.β = [0.85,0.948];
Baseline_parameter.θ = [0.66,0.86];
Baseline_parameter.F = [0.0,0.0];
Baseline_parameter.Fₓ = [4.5977,0.75];
# For the fixed cost, ratio of outputs
Baseline_parameter.τ = [0.354,0.354];
Baseline_parameter.ρ = [0.92,0.92];
Baseline_parameter.σₛ = [0.045,0.045];
Baseline_parameter.δ = 0.06;
curve = 1.0;
curve_fine  = 0.5;
Baseline_parameter.agrid = range(Baseline_parameter.a_min^curve,Baseline_parameter.a_max^curve,length=
    Baseline_parameter.n[1]).^(1/curve);

#Baseline_parameter.agrid_fine = range((Baseline_parameter.a_min*2 )^curve_fine,(Baseline_parameter.a_max-1.0)^curve_fine,length=
#    (Baseline_parameter.n_fine[1]-Baseline_parameter.n[1])).^(1/curve_fine);
#Baseline_parameter.agrid_fine = sort(union(Baseline_parameter.agrid,Baseline_parameter.agrid_fine));
function a_grid_fine_gen(agrid_fine_tmp,n_tmp,n_fine_tmp,agrid,a_min,a_max,curve_fine)
    size_agrid_fine = size(agrid_fine_tmp)[1]
    ii = 1;
    while n_fine_tmp> size_agrid_fine
        agrid_fine_tmp1 = range(a_min^curve_fine,a_max^curve_fine,length=
            (n_fine_tmp-n_tmp+ii)).^(1.0/curve_fine);
        agrid_fine_tmp = sort(union(agrid,agrid_fine_tmp1));
        size_agrid_fine = size(agrid_fine_tmp)[1];
        ii = ii+1;
    end
    return agrid_fine_tmp
end
function a_grid_fine_gen_midpoints(agrid,a_min,a_max,multiplier,n)
    agrid_tmp = zeros(multiplier*n[1]);
    agrid_tmp[1] =a_min;
    tmp_col_it_lo = 2;
    tmp_col_it_up = multiplier+1;
    for ii = 1:(n[1]-2)
        agrid_tmp[tmp_col_it_lo:tmp_col_it_up] =  range(agrid[ii],agrid[ii+1],
        length = tmp_col_it_up - tmp_col_it_lo+2)[2:end];
        tmp_col_it_lo = tmp_col_it_lo + multiplier;
        tmp_col_it_up = tmp_col_it_up + multiplier;
    end
    top_end_node_no = multiplier*n[1];
    agrid_tmp[tmp_col_it_lo:top_end_node_no]= (range(agrid[end-1],agrid[end],
    length = top_end_node_no-tmp_col_it_lo+2)[2:end]);
    return agrid_tmp
end
#Baseline_parameter.agrid_fine = a_grid_fine_gen(Baseline_parameter.agrid_fine, Baseline_parameter.n[1],
#Baseline_parameter.n_fine[1], Baseline_parameter.agrid,Baseline_parameter.a_min,Baseline_parameter.a_max,curve_fine)
Baseline_parameter.agrid_fine = a_grid_fine_gen_midpoints(Baseline_parameter.agrid,
Baseline_parameter.a_min,Baseline_parameter.a_max,4,Baseline_parameter.n);
Baseline_parameter.n_fine[1] = size(Baseline_parameter.agrid_fine)[1];

Baseline_parameter.fspace_a = fundef((:spli, Baseline_parameter.agrid, 0,Baseline_parameter.spliorder
        ));# define function space for the approximation of the solution of vfi.
Baseline_parameter.fspace_a_fine = fundef((:spli, Baseline_parameter.agrid_fine, 0,1
        ));# define function space for the approximation of the stationary distribution.

A_bar_foreign_improvement = 0.0; # Productivity improvements abroad mimicking that
    #Western Europe opened up to trade to ROW.

# Optimization bounds, where the excess demand functions would struggle:
#Wages:
Baseline_parameter.bounds[1,1] = 0.05
Baseline_parameter.bounds[1,2] = 0.5
#Interest rate
Baseline_parameter.bounds[2,1] = 0.001
Baseline_parameter.bounds[2,2] = 1/minimum(Baseline_parameter.β) - 1+ 0.05
#Output
Baseline_parameter.bounds[3,1] = 0.01
Baseline_parameter.bounds[3,2] = 5.0
# Price of the final good
Baseline_parameter.bounds[4,1] = 0.01
Baseline_parameter.bounds[4,2] = 2.0

prices = [0.166324882228137,  0.380045802996483,  0.144597061099265 ,  2.203464218966640, 1.554473456396248,
0.097211630924466,0.068443870634479];

Country_spec_p = StructArray([setup_state_space(i, Baseline_parameter) for i=1:Baseline_parameter.country_no]);
Baseline_parameter.Country_spec_p = Country_spec_p;
Baseline_parameter.iid_cost_value = [1.00];
Baseline_parameter.iid_cost_prob = [1.0];

#@btime residual = Residual_stst(prices, Baseline_parameter);

Default_parameter = copy(Baseline_parameter);
Default_parameter.agrid = range(Default_parameter.a_min,Default_parameter.a_max,length=
    Default_parameter.n[1]);
Default_parameter.agrid_fine = range((Default_parameter.a_min*2 )^curve_fine,(Default_parameter.a_max-1.0)^curve_fine,length=
    (Default_parameter.n_fine[1]-Default_parameter.n[1])).^(1/curve_fine);
Default_parameter.agrid_fine = sort(union(Default_parameter.agrid,Default_parameter.agrid_fine));
end
