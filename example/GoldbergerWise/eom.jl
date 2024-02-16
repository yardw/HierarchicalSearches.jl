module GoldbergerWiseEoM
# export paramsearch, M_IR, γ²₀, k, u, ϕP, yₘ, solveODE, getφ, ϕ0
module Consts 
    export yₘ, u, ϕP, ϕT, k, M_IR, γ²₀
    # parameters
    ### [m]=-1
    const yₘ = 1e0π #overall normalization s.t. y_m * M_Pl = pi
    
    ### [m]=0
    # l² = kappa^2 * phiP^2 / 2 reflects the strength of backreaction
    
    ### [m]=1
    const M_Pl = 1.0 #Plank mass
    const u  = 1.0e-1 #  u = log(ϕT / ϕP)/yₘ, this parameter should be fine-tuned to satisfy k*ym ~ O(50), but this causes instability by inrtoducing such a big hierarchy in numerical computation. Therefore here k*ym is set at O(10)
    const k = 37u #pp13 below eq(6.6)
    const M_IR = exp(-k*yₘ) * M_Pl #IR brane scale;(with M_Pl=1)
    const γ²₀ = 4k+2u
    # γ² initially is at large gamma limit
    
    ### [m]=3/2
    const ϕP = 1.e-1 # The scalar field value at Plank-brane
    const ϕT = exp(-u * yₘ)*ϕP
    

end

module AffliatedFunctions
    using ..Consts: yₘ, u, ϕP, ϕT
    export ϕ0, A, A′, A′′, V, V′, V′′, λP′, λT′, λP′′, λT′′, W, W′, α
    #static profile
    @inline ϕ0(  y) = ϕP * (ϕT/ϕP)^(y/yₘ) #ϕ0' = -u ϕ0
    @inline A(   y, l², k, γ²) = k * y + l²/6 * (ϕT/ϕP)^(2y/yₘ)
    @inline A′(  y, l², k, γ²) = k     + l²/6 * (ϕT/ϕP)^(2y/yₘ)* (-2u)
    @inline A′′( y, l², k, γ²) =         l²/6 * (ϕT/ϕP)^(2y/yₘ)* 4u^2
    @inline V′(  y, l², k, γ²) = u*(4k + u)*ϕ0(y) - 2/3*u^2 * 2l²/ϕP^2*ϕ0(y)^3 #κ^2 = 2l²/ϕP^2
    @inline V′′( y, l², k, γ²) = u*(4k + u -2u*2l²*(ϕT/ϕP)^(2y/yₘ))
    @inline λP′( φ, l², k, γ²) = -2u * (ϕP + φ)
    # @inline λP′( φ, l², k, γ²) = -2u * ϕP  + 2γ² * φ # no varphi needed, ϕ≠ϕₚ
    # @inline λP′( φ, l², k, γ²) = -2u * ϕP * φ + 2γ² * φ
    @inline λT′( φ, l², k, γ²) =  2u * (ϕT + φ)
    # @inline λT′( φ, l², k, γ²) =  2u * ϕT  + 2γ² * φ # no varphi needed, ϕ≠ϕₜ
    # @inline λT′( φ, l², k, γ²) =  2u * ϕT * φ + 2γ² * φ
    @inline λP′′(φ, l², k, γ²) = 2γ²
    @inline λT′′(φ, l², k, γ²) = 2γ²

    # superpotential
    W(ϕ, κ, k, u) = 6k/κ^2  - u*ϕ^2
    W′(ϕ, κ, k, u) =        - 2u*ϕ
    V(W, ϕ, κ, k, u) = 1/8 * (W′(ϕ, κ, k, u))^2 - κ^2 / 6 * W(ϕ, κ, k, u)^2

    # for convenience
    α(l²) = 3/2/(l²*u)
end # module AffliatedVariables

module EoMs
    using ..Consts: u, ϕP, yₘ
    using ..AffliatedFunctions: A, A′, A′′, α, ϕ0
    export P, Q, φtest, φ′test#, RR, SS, getφ′T, getφ′P
    # F'' + P F' + Q F = 0
    """
        P
    P for the EoM: F'' + P F' + Q F = 0
    """
    P(y, l², k, γ², m²) = 2u - 2A′(y, l², k, γ²)
    P(y, params) = P(y, params...)
    """
        Q
    Q for the EoM: F'' + P F' + Q F = 0
    """
    Q(y, l², k, γ², m²) = m² * exp(2A(y, l², k, γ²)) - 4A′′(y, l², k, γ²) - 4u*A′(y, l², k, γ²)
    Q(y, params) = Q(y, params...)
    
    # # \varphi 'P = \phi_P * (  (R1-P)F' + (S1-Q)F  )
    # R1(y, l², k, γ², m²) = u - 2A′(y, l², k, γ²)
    # R1(y,params) = R1(y, params...)
    # S1(y, l², k, γ², m²) = exp(-2u*y) / α(l²) - 2k - 2A(y, l², k, γ²) * u
    # S1(y,params) = S1(y, params...)
    # # \varphi 'T =  RR * F'T + SS * FT
    # """
    #     RR
    # RR for the EoM: varphi'T =  RR * F'T + SS * FT
    # """
    # RR(params, y) = (R1(y,params) - P(y,params))*ϕP
    # """
    #     SS
    # SS for the EoM: varphi'T =  RR * F'T + SS * FT
    # """
    # SS(params, y) = (S1(y,params) - Q(y,params))*ϕP
    # """
    #     getφ′T
    # given the value of F and F′ at the TeV brane(obtained by solving ODEs), return φ′ at TeV brane wrt the bulk equation of φ
    # """
    # getφ′T(FT, F′T, params) = RR(params, yₘ)*F′T + SS(params, yₘ)*FT
    # """
    # getφ′P
    # given the value of F and F′ at the Plank brane(obtained by solving ODEs), return φ′ at Plank brane wrt the bulk equation of φ
    # """
    # getφ′P(FP, F′P, params) = RR(params, 0)*F′P + SS(params, 0)*FP

    include("eom_addon.jl")
end # module EoMs

module BCs
    using ..Consts: u, ϕP, ϕT, yₘ
    using ..AffliatedFunctions
    using ..EoMs: φtest, φ′test #getφ′T, getφ′P, 
    export dφT, Δφ′Ttest, Δφ′Ptest#, Δφ′T, Δφ′P #, dFT, dφP, dFP,
    #BCs
    # non-perturbative gamma:
    """
        dφP
    given the value of F and F′ at the Plank brane, return φ′ wrt to the BC at Plank brane
    """
    function dφP(F′, F, l², k, γ², m²)
        φ = φP(F′, F, l², k, γ²)
        return 0.5λP′′(φ, l², k, γ²) * φ - 2u * ϕP * F  #(3.14)
    end
    dφP(F′, F, params) = dφP(F′, F, params...)
    φP(F′P, FP, l², k, γ²)  = -3ϕP^2/(2u*l²*ϕ0(0)) * (F′P - 2A′(0, l², k, γ²)*FP)  #(3.12)
    """
        Δφ′P
    given the value of F and F′ at the Plank brane, return the error of φ′to satisfy the BC at Plank brane
    """
    # function Δφ′P(F′, F, params)
    #     φ′P = dφP(F′, F, params...)
    #     return φ′P - getφ′P(F, F′, params)
    # end
    """
    dφT
    given the value of F and F′ at the TeV brane, return φ′ satisfying the BC at TeV brane
    """
    function dφT(F′, F, l², k, γ², m²) 
        φ = φT(F′, F, l², k, γ²)
        return -0.5λT′′(φ, l², k, γ²) * φ - 2u * ϕT * F  #(3.14)
    end
    dφT(F′, F, params) = dφT(F′, F, params...)
    φT(F′T, FT, l², k, γ²)  = -3ϕP^2/(2u*l²*ϕ0(yₘ)) * (F′T - 2A′(yₘ, l², k, γ²)*FT)  #(3.12)

    """
        Δφ′T
    given the value of F and F′ at the TeV brane, return the error of φ′to satisfy the BC at TeV brane
    """
    # function Δφ′T(F′, F, params)
    #     φ′T = dφT(F′, F, params...)
    #     return φ′T - getφ′T(F, F′, params)
    # end

    function Δφ′Ttest(F′, F, params)
        φ′T = dφT(F′, F, params...)
        return φ′T - φ′test(F′, F, params)(yₘ)        
    end
    function Δφ′Ptest(F′, F, params)
        φ′P = dφP(F′, F, params...)
        return φ′P - φ′test(F′, F, params)(0)        
    end



    # dφT(φ, F, l², k, γ²) =-0.5λT′′(φ, l², k, γ²) * φ - 2u * ϕT * F  #(3.14)

    # @inline dFP(φ, F, l², k, γ²) = 2A′(0 , l², k, γ²)*F - 2u*l²*ϕ0(0)/ϕP^2/3 * φ # same as the bulk equation of φ (3.12)
    # @inline dFT(φ, F, l², k, γ²) = 2A′(yₘ, l², k, γ²)*F - 2u*l²*ϕ0(yₘ)/ϕP^2/3 * φ


    # FP' + a F = 0
    # a(l², k, γ², m²) = ( 
    #     1/(α(l²)*γ²)
    #     *
    #     (-2k - exp(2A(0, l², k, γ²)) * m² + 1/α(l²) + u - 2A(0, l², k, γ²)*u + 4u*A(0, l², k, γ²) + 4A′′(0, l², k, γ²) )
    #     - 
    #     2A(0, l², k, γ²)
    #     )/(
    #     1 - 1/(α(l²)*γ²)*u
    #     )
    """
        dFP
    given the value of F at the Plank brane, return F′ satisfying the BC at Plank brane
    """
    # dFP(FP, l², k, γ², m²) = a(l², k, γ², m²)*FP
    # dFP(FP, params) = dFP(FP, params...)


end # module BCs

end # module GoldbergerWiseEoM
