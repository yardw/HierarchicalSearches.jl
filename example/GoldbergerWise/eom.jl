module GoldbergerWiseEoM
# export paramsearch, M_IR, γ²₀, k, u, ϕP, yₘ, solveODE, getφ, ϕ0
using DifferentialEquations

module Consts 
    export yₘ, u, ϕP, ϕT, k, M_IR, γ²₀
    # parameters
    const yₘ = 1e0π #overall normalization s.t. y_m * Lambda = pi
    const u  = 1.0e-1 #  u = log(ϕT / ϕP)/yₘ, this parameter should be fine-tuned to satisfy k*ym ~ O(50), but this causes instability by inrtoducing such a big hierarchy in numerical computation. Therefore here k*ym is set at O(10)
    const ϕP = 1.e-1 # The scalar field value at Plank-brane
    const ϕT = exp(-u * yₘ)*ϕP
    const k = 37u #pp13 below eq(6.6)
    const M_IR = exp(-k*yₘ) #IR brane scale;(with M_Pl=1)
    const γ²₀ = 4k+2u
    # l² = kappa^2 * phiP^2 / 2 reflects the strength of backreaction
    # γ² initially is at large gamma limit
end

module AffliatedFunctions
    using .Consts: yₘ, u, ϕP, ϕT
    export ϕ0, A, A′, A′′, V, V′, V′′, λP′, λT′, λP′′, λT′′, W, W′
    #static profile
    @inline ϕ0(  y::Number) = ϕP * (ϕT/ϕP)^(y/yₘ) #ϕ0' = -u ϕ0
    @inline A(   y::Number, l²::Number, k::Number, γ²::Number) = k * y + l²/6 * (ϕT/ϕP)^(2y/yₘ)
    @inline A′(  y::Number, l²::Number, k::Number, γ²::Number) = k     + l²/6 * (ϕT/ϕP)^(2y/yₘ)* (-2u)
    @inline A′′( y::Number, l²::Number, k::Number, γ²::Number) =         l²/6 * (ϕT/ϕP)^(2y/yₘ)* 4u^2
    @inline V′(  y::Number, l²::Number, k::Number, γ²::Number) = u*(4k + u)*ϕ0(y) - 2/3*u^2 * 2l²/ϕP^2*ϕ0(y)^3 #κ^2 = 2l²/ϕP^2
    @inline V′′( y::Number, l²::Number, k::Number, γ²::Number) = u*(4k + u -2u*2l²*(ϕT/ϕP)^(2y/yₘ))
    @inline λP′( φ::Number, l²::Number, k::Number, γ²::Number) = -2u * (ϕP + φ)
    # @inline λP′( φ::Number, l²::Number, k::Number, γ²::Number) = -2u * ϕP  + 2γ² * φ # no varphi needed, ϕ≠ϕₚ
    # @inline λP′( φ::Number, l²::Number, k::Number, γ²::Number) = -2u * ϕP * φ + 2γ² * φ
    @inline λT′( φ::Number, l²::Number, k::Number, γ²::Number) =  2u * (ϕT + φ)
    # @inline λT′( φ::Number, l²::Number, k::Number, γ²::Number) =  2u * ϕT  + 2γ² * φ # no varphi needed, ϕ≠ϕₜ
    # @inline λT′( φ::Number, l²::Number, k::Number, γ²::Number) =  2u * ϕT * φ + 2γ² * φ
    @inline λP′′(φ::Number, l²::Number, k::Number, γ²::Number) = 2γ²
    @inline λT′′(φ::Number, l²::Number, k::Number, γ²::Number) = 2γ²

    # superpotential
    W(ϕ, κ, k, u) = 6k/κ^2  - u*ϕ^2
    W′(ϕ, κ, k, u) =        - 2u*ϕ
    V(W, ϕ, κ, k, u) = 1/8 * (W′(ϕ, κ, k, u))^2 - κ^2 / 6 * W(ϕ, κ, k, u)^2
end # module AffliatedVariables

module EoMs
    using .Consts: u
    using .AffliatedFunctions: A, A′, A′′
    export P, Q
    # F'' + P F' + Q F = 0
    """
        P
    P for the EoM: F'' + P F' + Q F = 0
    """
    P(y, l², k, γ²₀) = 2u - 2A′(y, l², k, γ²₀)
    """
        Q
    Q for the EoM: F'' + P F' + Q F = 0
    """
    Q(y, l², k, γ²₀, m²) = m² * exp(2A(y, l², k, γ²₀)) - 4A′′(y, l², k, γ²₀) - 4u*A′(y, l², k, γ²₀)
end

module BCs
    using .Consts: u, ϕP, ϕT, yₘ
    using .AffliatedFunctions
    export dφP, dφT, dFP, dFT
    #BCs
    # non-perturbative gamma:
    """
        dφP
    given the value of F and φ at the Plank brane, return φ′ at Plank brane
    """
    @inline dφP(φ, F, l², k, γ²) = 0.5λP′′(φ, l², k, γ²) * φ - 2u * ϕP * F  #(3.14)
    """
        dφT
    given the value of F and φ at the TeV brane, return φ′ at TeV brane
    """
    @inline dφT(φ, F, l², k, γ²) =-0.5λT′′(φ, l², k, γ²) * φ - 2u * ϕT * F  #(3.14)
    @inline dFP(φ, F, l², k, γ²) = 2A′(0 , l², k, γ²)*F - 2u*l²*ϕ0(0)/ϕP^2/3 * φ # same as the bulk equation of φ (3.12)
    @inline dFT(φ, F, l², k, γ²) = 2A′(yₘ, l², k, γ²)*F - 2u*l²*ϕ0(yₘ)/ϕP^2/3 * φ

end # module BCs

# EoM of φ, obtained from the zero solution of the Einstein eq #(3.12)
function getφ(solF::ODESolution, params) 
    _, l², γ²= params
    F(y::Number)  = solF(y)[2]
    F′(y::Number) = solF(y)[1]
    φ(y::Number)  = -3ϕP^2/(2u*l²*ϕ0(y)) * (F′(y) - 2A′(y, l², k, γ²)*F(y))  #(3.12)
    return φ
end
function getφ(m2, l2, g2; FP=1., φP=1.)
    params = (m2, l2, g2)
    Fsol = eom.solveODE(FP, φP, params)
    eom.getφ(Fsol, params)
end

# EoM of F(perturbation on redshift factor A in metric)
function radionSpectrum_secondOrder!(ddf, df, f, params, y)
    F  = f[1]
    F′ = df[1]
    mF2, l², γ²= params
    dF′ = 2A′(y, l², k, γ²)*F′ + 4A′′(y, l², k, γ²)*F - 2u*F′ + 4u*A′(y, l², k, γ²)*F - mF2 * exp(2A(y, l², k, γ²))*F #(3.17)
    ddf[1]=dF′
end
# EoM of F(perturbation up to first order of l^2 and gamma^2)
function radionSpectrum_pert_secondOrder!(ddf, df, f, params, y)
    mF2, l², γ²= params
    M⁰₀, M¹₀, M⁰₁, M¹₁ = mF2
    F⁰₀, F¹₀, F⁰₁, F¹₁ = f
    F⁰₀′, F¹₀′, F⁰₁′, F¹₁′ = df

    dF⁰₀′ = -2(u-k)*F⁰₀′ + 4k*u*F⁰₀ - M⁰₀*exp(2k*y)*F⁰₀
    dF¹₀′ = -2(u-k)*F¹₀′ + 4k*u*F¹₀ - M¹₀*exp(2k*y)*F⁰₀ - 2u/3*exp(-2u*y)*(F⁰₀′-2u*F⁰₀)
    dF⁰₁′ = -2(u-k)*F⁰₁′ + 4k*u*F⁰₁ - M⁰₁*exp(2k*y)*F⁰₀
    dF¹₁′ = -2(u-k)*F¹₁′ + 4k*u*F¹₁ - M¹₁*exp(2k*y)*F⁰₀ - M¹₀*exp(2k*y)*F⁰₁ - 2u/3*exp(-2u*y)*(F⁰₁′-2u*F⁰₁)

    ddf[1], ddf[2], ddf[3], ddf[4] = dF⁰₀′, dF¹₀′, dF⁰₁′, dF¹₁′
end



function solveODE(FP, φP, params)
    _, l², γ²= params
    yspan = (0.0,yₘ)
    F′P= dFP(φP, FP, l², k, γ²)
    prob = SecondOrderODEProblem(radionSpectrum_secondOrder!,[F′P], [FP],yspan, params)
    # return solve(prob, Tsit5())
    return solve(prob, ImplicitEuler())
    # return solve(prob, Rosenbrock23())
end
function solveODE(settings::Settings, params)
    model_type = settings.model_type
    eom = nothing
    if model_type == Unperturbed
        eom = radionSpectrum_secondOrder!
    elseif model_type == Perturbed00to11Order
        eom = radionSpectrum_pert_secondOrder!
    else
        error("Unknown model type.")
    end
    error("Not implemented yet: solveODE for Settings inputs.")
end

function calculateΔφT(Fsol, params)
    _, l², γ²= params
    _, FT = Fsol(yₘ)
    φ = getφ(Fsol, params)
    φT = φ(yₘ)
    ys = range(yₘ*(1-1e-6) , yₘ, 2)
    φ′T = (diff(φ.(ys))/diff(ys))[1]
    # return (dφT(φT, FT, l², k, γ²) + φ′T)/sqrt((λT′′(0, l², k, γ²)/2)^2+ λT′(φT, l², k, γ² )^2+1)
    return dφT(φT, FT, l², k, γ²) + φ′T
    # return (dφT(φT, FT, l², k, γ²) + φ′T)/sqrt((λT′′(0, l², k, γ²)/2)^2+ λT′(0, l², k, γ² )^2+1)
end

function errBCwithφ(FP, params; φP = 1.)
    Fsol = solveODE(FP, φP, params)
    return calculateΔφT(Fsol, params)
end

function paramsearch(;l2=nothing, g2=nothing, FP = 1., φP = 1.)
    if !isnothing(g2) && isnothing(l2)
        function paramsearch_l2_m2(l2, m2)
            params = (m2, l2, g2)
            return errBCwithφ(FP, params, φP = φP )
        end
        return paramsearch_l2_m2
    elseif isnothing(g2) && !isnothing(l2)
        function paramserch_g2_m2(g2, m2)
            params = (m2, l2, g2)
            return errBCwithφ(FP, params, φP = φP)
        end
        return paramserch_g2_m2
    end
end
"""
    getseq(bounds, islogscaled, n = 100)
    return a sequence of n numbers between bounds, either linearly or logarithmically spaced
# Arguments
- `bounds::Tuple{Number, Number}`: the lower and upper bounds of the sequence
- `islogscaled::Bool`: whether the sequence is logarithmically or linearly spaced
- `n::Int`: the number of elements in the sequence
# Examples
```julia
julia> getseq((1, 10), true, 5)
5-element Array{Float64,1}:
  1.0
  2.51188643150958
  6.309573444801933
 15.848931924611133
 39.810717055349734
 ```
"""
function getseq(bounds::Tuple{Number,Number}, islogscaled::Bool, n = 100)
    if islogscaled
        exp10.(range(log10.(bounds)..., n))
    else
        range(bounds..., n)
    end
end
"""
    get_g2_m2_perturb_analytic_sol(s::eom.Settings; g2_logscaled = true)
    return a sequence of analytic solutions of γ² and m² to the 0-1 order perturbed Goldberger-Wise model.
# Arguments
- `s::eom.Settings`: the settings of the model
- `g2_logscaled::Bool`: whether the sequence of γ² is logarithmically or linearly spaced
# Examples
```julia
julia> get_g2_m2_perturb_analytic_sol(settings)
2-element Tuple{Array{Float64,1},Array{Float64,1}}:
 [0.001, 10.0]
 [1.0e-6, 0.1]
 ```
"""
function get_g2_m2_perturb_analytic_sol(s::Settings; g2_logscaled = true)
    l2 = s.l2; g2 = s.g2; u = s.u, k = s.k, yₘ = s.yₘ
    @assert l2 isa Number && g2 isa Tuple{Number, Number} "l² should be specified as a number, and γ² should be a range serving as the lower and upper bounds of γ²."
    g2s = getseq(g2, g2_logscaled)
    m2s = 4l2*(2k+u)*u^2/(3k)*(1-exp(2k*yₘ))/(1-exp((4k+2u)*yₘ)) * (1 .- (4k+2u) ./g2s)
    return g2s, m2s
end
end # module GoldbergerWiseEoM