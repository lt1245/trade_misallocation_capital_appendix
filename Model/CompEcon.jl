# ==============================================================================
#
#   CompEcon.jl — Self-contained manual installation of the CompEcon package
#
#   Original package authored by:
#       Spencer Lyon 
#       Thomas Sargent & John Stachurski  (QuantEcon project)
#       and contributors at EconForge
#
#   The original source is archived at:
#       https://github.com/EconForge/CompEcon.jl
#   and was built on top of BasisMatrices.jl:
#       https://github.com/QuantEcon/BasisMatrices.jl
#
#   We are deeply grateful to the original authors for making these tools
#   freely available to the research community. This consolidated file exists
#   solely because the package can no longer be installed through the standard
#   Julia package registry. No changes to the mathematical content have been
#   made; the sub-files (core.jl, cheb.jl, spli.jl, lin.jl, compat.jl) are
#   reproduced here verbatim and in full, separated by clearly marked sections,
#   in order to respect the structure and authorship of the original work.
#
# ==============================================================================

# ------------------------------------------------------------------------------
#   Dependencies
#   (Reexport is intentionally NOT used; BasisMatrices names are available
#    in the including scope via the `using` statement below.)
# ------------------------------------------------------------------------------

using BasisMatrices

using QuantEcon: gridmake, gridmake!, ckron, fix, fix!, qnwlege, qnwcheb,
                 qnwsimp, qnwtrap, qnwbeta, qnwgamma, qnwequi, qnwnorm,
                 qnwunif, qnwlogn, quadrect, do_quad

import BasisMatrices: funeval, funfitxy, funfitf

# old API exports
export fundef, fundefn, funnode, funbase, funbasex, funeval, funbconv,
       chebdef, chebnode, chebbase, chebbasex, chebdop,
       splidef, splinode, splibase, splibasex, splidop,
       lindef, linnode, linbase, lindop

# quad name exports
export qnwlege, qnwcheb, qnwsimp, qnwtrap, qnwbeta, qnwgamma, qnwequi, qnwnorm,
       qnwunif, qnwlogn, quadrect, do_quad


# ==============================================================================
#   BEGIN: core.jl  (originally src/core.jl)
#   — Dispatcher functions and generic translated functions from the Matlab API
# ==============================================================================

# -------------------- #
# Dispatcher functions #
# -------------------- #

const BASE_TYPES = [:spli, :cheb, :lin]
const ABSR_MAP = Dict(
    :none => Direct(),
    :direct => Direct(),
    :tensor => Tensor(),
    :expanded => Expanded(),
)
get_bformat(b::T) where T<:BasisMatrix{Direct} = :direct
get_bformat(b::T) where T<:BasisMatrix{Expanded} = :expanded
get_bformat(b::T) where T<:BasisMatrix{Tensor} = :tensor

function to_dict(bm::BasisMatrix)
    B = Dict{Symbol, Any}()
    B[:order] = bm.order
    B[:format] = get_bformat(bm)
    B[:vals] = bm.vals
    B
end

function bm_from_dict(B::Dict)
    arr_type = eltype(B[:vals])
    bm = BasisMatrix{typeof(ABSR_MAP[B[:format]]),arr_type}(B[:order], B[:vals])
    bm
end

base_exists(s::Symbol) = s in BASE_TYPES

basedef(s::Symbol, args...) =
    s == :spli ? splidef(args...) :
    s == :cheb ? chebdef(args...) :
    s == :lin  ? lindef(args...)  :
    error("somehow you snuck through here you 👺")

basenode(s::Symbol, args...) =
    s == :spli ? splinode(args...) :
    s == :cheb ? chebnode(args...) :
    s == :lin  ? linnode(args...)  :
    error("somehow you snuck through here you 👺")

BasisMatrices.evalbase(s::Symbol, args...) =
    s == :spli ? splibase(args...) :
    s == :cheb ? chebbase(args...) :
    s == :lin  ? linbase(args...)  :
    error("somehow you snuck through here you 👺")

# Helper function
function squeeze_trail(x::AbstractArray)
    sz = size(x)
    squeezers = Int[]
    n = length(sz)
    for i=n:-1:1
        if sz[i] == 1
            push!(squeezers, i)
        else
            break
        end
    end
    squeeze(x, tuple(squeezers...))
end


# ---------------------------- #
# Generic translated functions #
# ---------------------------- #

# from fundef.m  -- DONE
function fundef(foo...)
    d = length(foo)  # 89
    n = zeros(Int, d)  # 93
    b = zeros(d)  # 94
    a = zeros(d)  # 95
    p = Array{Any}(undef, d)  # 96
    _params = Array{BasisMatrices.BasisParams}(undef, d)

    basetype = Array{Symbol}(undef, d)
    for j=1:d
        basetype[j] = foo[j][1]  # 99
        !(base_exists(basetype[j])) && error("Unknown basis $(foo[j][1])")
        n[j], a[j], b[j], p[j], _params[j] = basedef(basetype[j], foo[j][2:end]...)  # 124
    end

    # package output. Lines 143-150
    g = Dict{Symbol, Any}()
    g[:d] = d
    g[:n] = n
    g[:a] = a
    g[:b] = b
    g[:basetype] = basetype
    g[:params] = p
    g[:_basis_params] = _params
    g[:_basis] = Basis(_params...)
    g
end

# fundefn.m
function fundefn(basistype::Symbol, n, a, b, order=3)
    d = length(n)
    length(a) != d && error("a must be same dimension as n")
    length(b) != d && error("b must be same dimension as n")
    any(a .> b) && error("left endpoints must be less than right endpoints")
    any(n .< 2) && error("n(i) must be greater than 1")

    params = Array{Any}(undef, 1, d)
    if basistype == :cheb
        for i=1:d params[i] = Any[:cheb, n[i], a[i], b[i]] end
    elseif basistype == :spli
        for i=1:d params[i] = Any[:spli, [a[i], b[i]], n[i]-order+1, order] end
    elseif basistype == :lin
        for i=1:d params[i] = Any[:lin, [a[i], b[i]], n[i]] end
    end

    fundef(params...)
end

# funnode.m -- DONE
funnode(basis::Dict) = nodes(basis[:_basis])

# funbase.m -- DONE
function funbase(basis::Dict, x=funnode(basis)[1], order=fill(0, 1, basis[:d]))
    BasisMatrix(basis[:_basis], Expanded(), x, order).vals[1]
end

# funbasex.m -- DONE
function funbasex(basis::Dict{Symbol}, x=funnode(basis)[1], order=0,
                  bformat::Symbol=:none)
    to_dict(BasisMatrix(basis[:_basis], ABSR_MAP[bformat], x, order))
end

# funfitf.m -- DONE
funfitf(basis, f::Function, args...) = funfitf(basis[:_basis], f, args...)

# funfitxy.m -- DONE
function funfitxy(basis, x, y)
    c, bm = funfitxy(basis[:_basis], x, y)
    return c, to_dict(bm)
end

# funeval.m -- DONE
function funeval(c, basis::Dict, B, _order=0)
    isempty(c) && error("missing basis coefficients")
    order = BasisMatrices._check_order(basis[:d], _order)

    if isa(B, Dict)  # B is a basis structure
        bm = bm_from_dict(B)
        y = funeval(c, bm, order)
        return y, B
    else
        bm = BasisMatrix(basis[:_basis], B, order)
        y = funeval(c, bm, bm.order)
        return y, to_dict(bm)
    end
end

# fund.m
function fund(c, basis, x, hess_opt)
    # TODO: come back when I need this. I think I should probably do something
    #       like what optim does and write functions `f`, `fg!` and `fgh!` to
    #       replicate this for the type of basis instead of this function
    nothing
end

# funbconv.m  -- DONE
function funbconv(b::Dict, order=fill(0, 1, size(b[:order], 2)),
                  format::Symbol=:expanded)
    bm = bm_from_dict(b)
    new_bm = convert(typeof(ABSR_MAP[format]), bm, order)
    to_dict(new_bm)
end

# ==============================================================================
#   END: core.jl
# ==============================================================================


# ==============================================================================
#   BEGIN: cheb.jl  (originally src/cheb.jl)
#   — Chebyshev basis original API
# ==============================================================================

# ------------ #
# Original API #
# ------------ #

# chebdef.m -- DONE
function chebdef(n::Int, a::Real, b::Real)
    p = ChebParams(n, a, b)
    return (p.n, p.a, p.b, Any[p.n, p.a, p.b], p)
end


# chebnode.m -- DONE
chebnode(n::Int, a::Real, b::Real, nodetype=0) =
    nodes(ChebParams(n, a, b), nodetype)

# chebdop.m -- DONE
function chebdop(n, a, b, order=1)
    D, p = derivative_op(ChebParams(n, a, b), order)
    D, p.n-order, p.a, p.b, Any[p.n-order, p.a, p.b]
end


# chebbase.m -- DONE
function chebbase(n, a, b, x=chebnode(n, a, b, 1), order=0, nodetype=1)
    B = evalbase(ChebParams(n, a, b), x, order, nodetype)
    B, x
end

# chebbasex.m -- DONE
chebbasex(n, a, b, x) = evalbasex(ChebParams(n, a, b), x)

# ==============================================================================
#   END: cheb.jl
# ==============================================================================


# ==============================================================================
#   BEGIN: spli.jl  (originally src/spli.jl)
#   — Spline basis original API
# ==============================================================================

# ------------ #
# Original API #
# ------------ #

# splidef.m -- DONE
function splidef(breaks, evennum=0, k::Int=3)
    p = SplineParams(breaks, evennum, k)
    return (
        length(p.breaks),
        minimum(p.breaks),
        maximum(p.breaks),
        Any[p.breaks, 0, k],
        p
    )
end

# splinode.m  -- DONE
splinode(breaks::AbstractVector, evennum::Int, k::Int=3) =
    nodes(SplineParams(breaks, evennum, k))

function splidop(breaks, evennum=0, k=3, order=1)
    D, p = derivative_op(SplineParams(breaks, evennum, k), order)
    n = length(breaks) + k - 1

    D, n-order, breaks[1], breaks[end], Any[breaks, evennum, k-order]
end

# splibas.m -- DONE
function splibase(breaks::AbstractVector, evennum, k=3, x=splinode(breaks, evennum, k),
                  order=0)
    B = evalbase(SplineParams(breaks, evennum, k), x, order)
    B, x
end

# ==============================================================================
#   END: spli.jl
# ==============================================================================


# ==============================================================================
#   BEGIN: lin.jl  (originally src/lin.jl)
#   — Linear basis original API
# ==============================================================================

# ------------ #
# Original API #
# ------------ #

# lindef.m -- DONE
function lindef(breaks::AbstractVector, evennum::Int=0)
    p = LinParams(breaks, evennum)
    return (
        length(p.breaks),
        minimum(p.breaks),
        maximum(p.breaks),
        Any[p.breaks, evennum],
        p
    )
end

# linnode.m -- DONE
linnode(breaks, evennum) = nodes(LinParams(breaks, evennum))

# lindop.m
function lindop(breaks, evennum, order)
    D, params = derivative_op(LinParams(breaks, evennum), order)
    n, a, b = length(params.breaks), params.breaks[1], params.breaks[end]
    D, n, a, b, params
end

# linbase.m -- DONE
function linbase(breaks, evennum=0, x=breaks, order=0)
    B = evalbase(LinParams(breaks, evennum), x, order)
    B, x
end

# ==============================================================================
#   END: lin.jl
# ==============================================================================


# ==============================================================================
#   BEGIN: compat.jl  (originally src/compat.jl)
#   — Compatibility shims to maintain the Matlab API
# ==============================================================================

# Stuff We need to maintain compatibility with the matlab api

## BasisFamily stuff

# needed for the `convert(::Basis, ::Dict)` method below
Base.convert(::Type{BasisFamily}, s::Symbol) =
    s == :spli ? Spline() :
    s == :cheb ? Cheb() :
    s == :lin  ? Lin() :
    error("Unknown basis type")

old_name(::Lin) = :lin
old_name(::Cheb) = :cheb
old_name(::Spline) = :spli

## BasisParams stuff

# needed for the `convert(::Basis, ::Dict)` method below
Base.convert(::Type{BasisParams}, p::Vector{Any}) =
    length(p) == 3 && isa(p[1], Number) ? ChebParams(p...) :
    length(p) == 3 ? SplineParams(p...) :
    length(p) == 2 ? LinParams(p...) :
    error("Unknown parameter type")

old_name(::LinParams) = :lin
old_name(::ChebParams) = :cheb
old_name(::SplineParams) = :spli

old_params(p::LinParams) = Any[p.breaks, p.evennum]
old_params(p::ChebParams) = Any[p.n, p.a, p.b]
old_params(p::SplineParams) = Any[p.breaks, p.evennum, p.k]

## Basis stuff

# convert old API to new API
function Base.convert(::Type{Basis}, b::Dict{Symbol, Any})
    btype = map(x->convert(BasisFamily, x), b[:basetype])
    param = BasisParams[convert(BasisParams, x) for x in b[:params]]
    Basis(param...)
end

# convert new API to old API
function revert(b::Basis)
    B = Dict{Symbol, Any}()
    B[:d] = ndims(b)
    the_nodes = nodes.(b.params)
    B[:n] = [length(vec) for vec in the_nodes]
    B[:a] = [minimum(vec) for vec in the_nodes]
    B[:b] = [maximum(vec) for vec in the_nodes]
    B[:basetype] = Symbol[old_name(bt) for bt in b.params]
    B[:params] = Any[old_params(p) for p in b.params]
    B
end

# add method to funbasex that creates a BasisStructure
function funbasex(basis::Basis, x=nodes(basis)[1], order=0,
                  bformat::BasisMatrices.ABSR=Direct())
    BasisMatrix(Nothing, basis, bformat, x, order)
end

funbase(basis::Basis, x=nodes(basis)[1], order=fill(0, 1, ndims(basis))) =
    funbasex(basis, x, order, Expanded()).vals[1]

# ==============================================================================
#   END: compat.jl
# ==============================================================================
