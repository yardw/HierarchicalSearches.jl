using DifferentialEquations: ODESolution, ODEProblem, SecondOrderODEProblem, solve, ImplicitEuler
using Plots
include("convolve.jl")
include("storage.jl")
include("eom.jl")
using .GoldbergerWiseEoM.EoMs
import .GoldbergerWiseEoM.Consts: yₘ, k, M_IR, γ²₀, u
using .GoldbergerWiseEoM.AffliatedFunctions: ϕ0, A, A′
using .GoldbergerWiseEoM.BCs

function eom!(ddf, df, f, params, y)
    F  = f[1]
    F′ = df[1]
    dF′ = -P(y,params)*F′ -Q(y,params)*F#(3.17)
    ddf[1]=dF′
end

"""
    solveODE_LR(F0, params, y0)
# Arguments
- `F0`: initial value of F′ and f
- `params`: parameters of the model(l², k, γ², m²)
- `y0`: the initial value of y
"""
function solveODE_LR(F0, params, y0)
    F′, F = F0
    yspanL = (y0, 0)
    yspanR = (y0, yₘ)

    probL = SecondOrderODEProblem(eom!,[F′], [F],yspanL, params)
    probR = SecondOrderODEProblem(eom!,[F′], [F],yspanR, params)
    
    solL = solve(probL, ImplicitEuler())
    solR = solve(probR, ImplicitEuler())
    sol(y) = y < y0 ? solL(y) : solR(y)

    return sol
end
function errBCwithφtest(params; F0′, F0=1., y0 )
    Fsol = solveODE_LR([F0′, F0], params, y0)
    FT′,FT = Fsol(yₘ)
    FP′,FP = Fsol(0)
    errR = Δφ′Ttest(FT′, FT, params)
    errL = Δφ′Ptest(FP′, FP, params)
    return errL, errR
end
function searchYMtest(l2, g2, y0, k=k)
    (Fp, m2) -> errBCwithφtest((l2, k, g2*γ²₀, m2*M_IR^2); F0′=Fp, F0=1., y0=y0)
end
function m2s_analytical(l2, g2)
    return m2s = 4l2*(2k+u)*u^2/(3k)*(1-exp(2k*yₘ))/(1-exp((4k+2u)*yₘ)) * (1 .- 1 /g2) / M_IR^2
end
function search(l2=1e-3, g2=1e2, y0=0.99yₘ)
    # both Fp and m2 leaner scaled
    prob = searchYMtest(l2, g2, y0, k)
    n = 20
    Fpscale = 1e-6
    m2scale = 3e-5
    Fp = 7.3999643 .+ Fpscale*range(-1, stop=1, length=n)
    # Fp = 7.4 - 3.55e-5 .+ 1e-1range(-1, stop=1, length=n)
    m2 = m2scale*range(-1e-1, stop=1, length=2n)
    X = Fp*ones(2n)';
    Y = ones(n)*m2';
    Z = prob.(X,Y);
    # extract the hypersurface constrained by BCs
    Z1 = map(first, Z) .|> sign
    Z2 = map(last, Z)  .|> sign
    z1s = pick_xy((Z1 |> with_sobel_edge_detect) .> 0, Fp, m2)
    slope1, bias1 = [z1s[2] ones(length(z1s[2]))]\z1s[1]
    z2s = pick_xy((Z2 |> with_sobel_edge_detect) .> 0, Fp, m2)
    slope2, bias2 = [z2s[2] ones(length(z2s[2]))]\z2s[1]
    Fp0, m20 = [1 -slope1; 1 -slope2]\[bias1; bias2]
    # refine the search
    m2 = m20 .+ 0.1m2scale*range(-1, stop=1, length=2n)
    Fp = Fp0 .+ 0.1Fpscale*range(-1, stop=1, length=n)
    X = Fp*ones(2n)';
    Y = ones(n)*m2';
    Z = prob.(X,Y);
    # extract the hypersurface constrained by BCs
    Z1 = map(first, Z) .|> sign
    Z2 = map(last, Z)  .|> sign
    Z3 = Z1 .+ 2Z2 .+ 3
    z1s = pick_xy((Z1 |> with_sobel_edge_detect) .> 0, Fp, m2)
    slope1, bias1 = [z1s[2] ones(length(z1s[2]))]\z1s[1]
    z2s = pick_xy((Z2 |> with_sobel_edge_detect) .> 0, Fp, m2)
    slope2, bias2 = [z2s[2] ones(length(z2s[2]))]\z2s[1]
    Fp0, m20 = [1 -slope1; 1 -slope2]\[bias1; bias2]
    # println("Fp0=$Fp0, m20=$m20")
    let 
        heatmap(m2, Fp, Z3)
        plot!(x->slope1*x+bias1, m2[1], m2[end], ylims=(Fp[1], Fp[end]), label="left boundary", lw=3)
        plot!(x->slope2*x+bias2, m2[1], m2[end], ylims=(Fp[1], Fp[end]), label="right boundary", lw=3)
        scatter!([m20], [Fp0], label="critical point at $(m20)")
        xlabel!("m²/M_IR²")
        ylabel!("F′")
        title!("l²=$(l2), γ²=$(g2)γ²₀")
    end
    # return m20
end

using LaTeXStrings
g2s = exp10.(range(-1e-1, stop=3, length=30));
m2s_num = search.(1e-3, g2s, 0.98yₘ);
m2s_num2 = search.(1e0, g2s, 0.98yₘ);
m2s_ana = m2s_analytical.(1e-3, g2s);
begin
    scatter(g2s, m2s_num, label="numerical", markersize=4, size=(800, 600), dpi=300)
    plot!(g2s, m2s_ana, label="analytical", xaxis=:log, lw=2, color=:black)
    xlabel!(L"\gamma^2/\gamma^2_0")
    ylabel!(L"m^2/M_{IR}^2")
    title!("Radion Mass with " * L"l^2=1\times10^{-3}" * "\n where " * L"M_{IR} = e^{-k y_m} M_{Pl},\quad \gamma^2_0 = 4k + 2u")
    gui()
    savefig("example/GoldbergerWise/figs/m2_g2.png")
end

# solve for Ax = b => x = A\b
let A = rand(2,2)
    A \ (A*[3,2])
end

# sol(t) = du, u
let eom!(ddu, du, u, p, x) = ddu[1] = -u[1], du0 = [0.], u0 = [1.], tspan = (0., 2π), tspan_rev = (2π, 0.), p = nothing
    prob = SecondOrderODEProblem(eom!, du0, u0, tspan, p)
    sol = solve(prob)
    plot(sol)
    prob_rev = SecondOrderODEProblem(eom!, du0, u0, tspan_rev, p)
    sol_rev = solve(prob_rev)
    plot!(sol_rev)
    prob_mid = SecondOrderODEProblem(eom!, du0, -u0, (π, 2π), p)
    sol_mid = solve(prob_mid)
    plot!(sol_mid)
end