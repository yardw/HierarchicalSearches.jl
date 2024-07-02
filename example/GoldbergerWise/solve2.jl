using DifferentialEquations: ODESolution, ODEProblem, SecondOrderODEProblem, solve, ImplicitEuler
using Plots
using LaTeXStrings
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
function extract_sol(props)
    l2, g2, m2, dm2, F′ = props
    params = (l2, k, g2*γ²₀, m2*M_IR^2)
    yspan = (yₘ, 0)
    # return solve(SecondOrderODEProblem(eom!,[F′], [1.],yspan, params), ImplicitEuler())
    return solveODE_LR([F′, 1.], params, 0.99yₘ)
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
include("search_alg.jl")
function search(l2=1e-3, g2=1e2; y0=0.99yₘ, fig=false, niters=5, scale0 = 10, m2scale = 1e-4, Fprange = (2., 13.), tol=1e-6, nresol=10)
    # both Fp and m2 leaner scaled
    prob = searchYMtest(l2, g2, y0, k)
    n = nresol
    m20 = 0
    Fp1 = bisection(Fp->prob(Fp, m20)|>first, Fprange...)# Fp0 = 7.3999643
    Fp2 = bisection(Fp->prob(Fp, m20)|>last, Fprange...)# Fp0 = 7.3999643
    Fp0 = 0.5(Fp1+Fp2)
    Fpscale = 0.8abs(Fp1-Fp2)
    # extract the hypersurface constrained by BCs
    for i in 1:niters
        Fp, m2, err = sample(prob, Fp0, Fpscale, m20, m2scale, n)
        Fp0_old, m20_old = Fp0, m20
        Fp0, m20, hypersurfaces = fit(Fp, m2, err)
        # Fp0, m20, slopes, biass = fit(Fp, m2, err)
        if isnan(m20) || isnan(Fp0)  
            scale=1/scale0
            Fp0, m20 = Fp0_old, m20_old 
        else
            scale=scale0
        end
        m2scale, Fpscale = rescale(m2scale, Fpscale, m20, Fp0, hypersurfaces, scale)
        # m2scale, Fpscale = rescale(m2scale, Fpscale, m20, Fp0, slopes, biass, scale)
        m2scale, Fpscale = max(m2scale, abs(m20-m20_old)), max(Fpscale, abs(Fp0-Fp0_old)) # avoid to jump too far from the previous region
        tol > m2scale && break
    end
    Fp, m2, err = sample(prob, Fp0, Fpscale, m20, m2scale, n)
    Fp0, m20, hypersurfaces = fit(Fp, m2, err)
    
    if fig
        Z1 = map(first, err) .|> sign
        Z2 = map(last, err)  .|> sign
        Z3 = Z1 .+ 2Z2 .+ 3
        f = heatmap(m2, Fp, Z3, clims=(0,6))
        scatter!([m20], [Fp0], label="critical point at $(m20)")
        xlabel!("m²/M_IR²")
        ylabel!("F′")
        title!("l²=$(l2), γ²=$(g2)γ²₀")
        return m20, m2scale/n, f
    end
    return m20, m2scale/n, Fp0
end

include("storage.jl")
sto = storage("example/GoldbergerWise/data", "r")
storage("example/GoldbergerWise/data")
sto["l2g2m2"] = data
sto["l2g2m2"]
sto["l2g2m2dm"] = data
data = sto["l2g2m2dm"]

data = []
# l2 > 1 is not physical since l = κϕₚ/√2 and κ ≈ 1/Mₚ; shader it
l2, g2 = 1e-1, 1e1
m20, dm2, f = search(; fig = true); f
push!(data, (l2, g2, m20, dm2));
dm2

using ProgressMeter
newdata_fixtol_withFp = []
newdata_fixtol_withFp = let l2=1e-1, g2=1e3
        m20, dm2, Fp = search(l2, g2; m2scale=1e-2*l2)
        newdata_fixtol_withFp = [[l2, g2, m20, dm2, Fp]]
end
newdata_fixtol_withFp
push!(newdata_fixtol_withFp, newdata_fixtol_withFp[1] .+ [0, 0, newdata_fixtol_withFp[1][4], 0, 0])
push!(newdata_fixtol_withFp, newdata_fixtol_withFp[1] .+ [0, 0, -newdata_fixtol_withFp[1][4], 0, 0])


newdata_fixtol_withFp = let data=[], l2=1e-1
    # @showprogress for l2 in exp10.(range(-2, 0, 20)), g2 in [1.5, 3, 1e3]
    @showprogress for g2 in [0.5, 1.5, 2, 1e3]
        m20, dm2, Fp = search(l2, g2; m2scale=1e-2*l2)
        push!(data, (l2, g2, m20, dm2, Fp))
    end
    data
end
newdata_fixtol_withFp = let data=[], g2=1e3
    # @showprogress for l2 in exp10.(range(-2, 0, 20)), g2 in [1.5, 3, 1e3]
    @showprogress for l2 in [1e-3, 1e-2, 1e-1, 2e-1, 3e-1, 4e-1]
        m20, dm2, Fp = search(l2, g2; m2scale=1e-2*l2)
        push!(data, (l2, g2, m20, dm2, Fp))
    end
    data
end
f= let
    plot()
    F0(y) = extract_sol(newdata_fixtol_withFp[1])(y)[2]
    F_upper(y) = extract_sol(newdata_fixtol_withFp[2])(y)[2] 
    F_lower(y) = extract_sol(newdata_fixtol_withFp[3])(y)[2]
    x = range(1e-2,0.4yₘ,150)
    plot!(x, F_upper.(x), fillrange=F_lower.(x), label=L"l^2="*"$(newdata_fixtol_withFp[1][1])"*L"\,, \gamma^2="*"$(newdata_fixtol_withFp[1][2])"*L" \gamma^2_0,\,"*L" \Delta m^2 ="*" $((newdata_fixtol_withFp[1][4]))"*L"M^2_{\rm IR}", fillalpha=0.5)
    plot!(x, F0.(x), label=nothing)
    plot!(xlabel=L"y/R", ylabel=L"F(y)", yaxis=:log, size=(600, 450), dpi=300, legend=:topleft)
    savefig("example/GoldbergerWise/figs/F_profile.png")
    # plot!(xlabel="y/R", ylabel="F", yaxis=:log, xaxis=:log)
end
f= let
    plot()
    for i in eachindex(newdata_fixtol_withFp)
        sol = extract_sol(newdata_fixtol_withFp[i])
        F(y) = sol(y)[2]
        l2, g2, m2, dm2, Fp = newdata_fixtol_withFp[i]
        plot!(F, range(1e-2,0.9yₘ,150), label="l²=$l2, γ²=$g2 γ²₀")
    end
    plot!(xlabel="y/R", ylabel="F", yaxis=:log, legend=:topleft)
    # plot!(xlabel="y/R", ylabel="F", yaxis=:log, xaxis=:log)
end
display(f)
let
    plot()
    for i in eachindex(newdata_fixtol_withFp)
        sol = extract_sol(newdata_fixtol_withFp[i])
        F(y) = sol(y)[2]
        l2, g2, m2, dm2, Fp = newdata_fixtol_withFp[i]
        plot!(F, range(1e-2,0.3yₘ,150), label="l²=$l2, γ²=$g2 γ²₀")
    end
    plot!(xlabel="y/R", ylabel="F", yaxis=:log)
    # plot!(xlabel="y/R", ylabel="F", yaxis=:log, xaxis=:log)
end
plot(y->extract_sol(newdata_fixtol_withFp[1])(y)[2], range(0,yₘ, length=100))
plot!(extract_sol(newdata_fixtol_withFp[3]), idxs=(0,2))
plot!(extract_sol(newdata_fixtol_withFp[4]), idxs=(0,2))
newdata_fixtol_withFp



l2s0 = filter(x->x[2] ≈ 1e3, newdata_fixtol).|>x->x[1]
let data=newdata_fixtol
    cmap = :Dark2_3#Lapaz10,Paired_3
    f = plot(palette=cmap, xlabel=L"\gamma^2/\gamma^2_0", ylabel=L"\frac{m^2}{l^2}/M_{IR}^2", size=(600, 450), dpi=300, legend=:bottomright)
    xlim = (4e-1, 1e3)
    ylim = (-5e-3, 1.5e-2)
    plot!([1e-1xlim[1],1],[ylim[1],ylim[1]], fillrange=[ylim[2],ylim[2]], c=:gray, fillalpha=0.3, label="tachyonic")
    plot!(x->m2s_analytical(1, x), exp10.(range(-0.3, stop=3, length=40)), label="perturbative", lw=2, color=:black, ls=:dash)
    for i in eachindex(l2s0)
        l2 = l2s0[i]
        dat = filter(x->x[1] ≈ l2, data)
        g2s = dat.|>x->x[2]
        m2s = dat.|>x->x[3]
        dms = dat.|>x->x[4]

        plot!(g2s, m2s./(l2), label=L"l^2="*"$(round(l2, sigdigits=5)) numerical", lw=2, color=i, xaxis=:log, xlim = xlim, ylim = ylim, ribbon=dms./(l2), fillalpha=0.3)
    end
    # plot!()
    savefig("example/GoldbergerWise/figs/m2_g2_2.png")
end
let l2 = l2s0[1]
    dat = filter(x->x[1] ≈ l2, data)
    g2s = dat.|>x->x[2]
    m2s = dat.|>x->x[3]
    dms = dat.|>x->x[4]
    xlim = (1e-4, 1e3)
    ylim = 2e-5.*(-1, 1) .+ extrema(m2s)
    f = plot([1e-1xlim[1],1],[ylim[1],ylim[1]], fillrange=[ylim[2],ylim[2]], 
        label="tachyonic", fillalpha=0.3, color=:gray, lw=0, 
        xlim = xlim, ylim = ylim, legend=:topleft,
        xlabel=L"\gamma^2/\gamma^2_0", ylabel=L"m^2/M_{IR}^2",
        size=(600, 450), dpi=300)
    plot!(f, g2s, m2s, label=L"l^2="*"$(round(l2, digits=3))", xaxis=:log, lw=2, color=:brown)
    plot!(x->m2s_analytical(l2, x), exp10.(range(-4, stop=3, length=40)), label=nothing, lw=2, color=:brown, ls=:dash)
    plot!([xlim[1], 1],[0, 0], lw=1, color=:black, label=nothing, linestyle=:dash)

    plot!(f, inset=(1, bbox(0.545, 0.15, 0.46, 0.75)))

    g2s = g2s[17:30]
    m2s = m2s[17:30]
    xlim = extrema(g2s) .* (0.9, 1.1) 
    ylim = extrema(m2s) .* (0.9, 1.1)
    plot!(f[2], [1e-1xlim[1],1],[ylim[1],ylim[1]], fillrange=[ylim[2],ylim[2]], label=nothing, fillalpha=0.3, color=:gray, lw=0, xlim = xlim, ylim = ylim)
    plot!(f[2], g2s, m2s, ribbon=2 .*dms, xaxis=:log, lw=2, label=nothing, fillalpha=0.5, color=:brown, alpha=0.4)
    plot!(f[2], x->m2s_analytical(l2, x), exp10.(range(-1, stop=3, length=20)), label=nothing, lw=2, color=:brown, ls=:dash)
    plot!(f[2], [xlim[1], 1],[0, 0], lw=1, color=:black, label=nothing, linestyle=:dash)
    savefig("example/GoldbergerWise/figs/m2_g2_1.png")
end
let l2 = l2s0[20]
    dat = filter(x->x[1] ≈ l2, data)
    g2s = dat.|>x->x[2]
    m2s = dat.|>x->x[3]
    dms = dat.|>x->x[4]
    xlim = (1e-4, 1e3)
    ylim = 2e-5.*(-1, 1) .+ extrema(m2s)
    f = plot([1e-1xlim[1],1],[ylim[1],ylim[1]], fillrange=[ylim[2],ylim[2]], 
        label="tachyonic", fillalpha=0.3, color=:gray, lw=0, 
        xlim = xlim, ylim = ylim, legend=:topleft,
        xlabel=L"\gamma^2/\gamma^2_0", ylabel=L"m^2/M_{IR}^2",
        size=(600, 450), dpi=300)
    plot!(f, g2s, m2s, label=L"l^2="*"$(round(l2, digits=3))", xaxis=:log, lw=2, color=:brown)
    plot!(x->m2s_analytical(l2, x), exp10.(range(-4, stop=3, length=40)), label=nothing, lw=2, color=:brown, ls=:dash)
    plot!([xlim[1], 1],[0, 0], lw=1, color=:black, label=nothing, linestyle=:dash)

    plot!(f, inset=(1, bbox(0.545, 0.15, 0.46, 0.75)))

    g2s = g2s[17:30]
    m2s = m2s[17:30]
    xlim = extrema(g2s) .* (0.9, 1.1) 
    ylim = extrema(m2s) .* (0.9, 1.1)
    plot!(f[2], [1e-1xlim[1],1],[ylim[1],ylim[1]], fillrange=[ylim[2],ylim[2]], label=nothing, fillalpha=0.3, color=:gray, lw=0, xlim = xlim, ylim = ylim)
    plot!(f[2], g2s, m2s, ribbon=2 .*dms, xaxis=:log, lw=2, label=nothing, fillalpha=0.5, color=:brown, alpha=0.5)
    plot!(f[2], x->m2s_analytical(l2, x), exp10.(range(-1, stop=3, length=20)), label=nothing, lw=2, color=:brown, ls=:dash)
    plot!(f[2], [xlim[1], 1],[0, 0], lw=1, color=:black, label=nothing, linestyle=:dash)
    savefig("example/GoldbergerWise/figs/m2_g2_2.png")
end

g2s0 = filter(x->x[1] ≈ 1e-2, newdata_fixtol).|>x->x[2]
let data=newdata_fixtol
    cmap = :Dark2_3#Lapaz10,Paired_3
    f = plot(palette=cmap, xlabel=L"l^2", ylabel=L"m^2/M_{IR}^2", size=(600, 450), dpi=300, legend=:bottomright)
    for i in eachindex(g2s0)
        g2 = g2s0[i]
        dat = filter(x->x[2] ≈ g2, data)
        l2s = dat.|>x->x[1]
        m2s = dat.|>x->x[3]
        dms = dat.|>x->x[4]
        xlim = (1e-2, 1)
        ylim = (7e-5, 0.01)
        plot!(l2s, m2s, label=L"\gamma^2/\gamma^2_0="*"$(round(g2*(1+1e-5), sigdigits=5)) numerical", lw=2, color=i, xaxis=:log, xlim = xlim, yaxis=:log, ylim = ylim, ribbon=1.5e-1 .*m2s, fillalpha=0.3)
        plot!(x->m2s_analytical(x, g2), exp10.(range(-3, stop=xlim[2], length=40)), label=L"\gamma^2/\gamma^2_0="*"$(round(g2*(1+1e-5), sigdigits=5)) perturbative", lw=2, color=i, ls=:dash)
    end
    savefig("example/GoldbergerWise/figs/m2_l2_2.png")
end
let 
    g2 = g2s0[3]; c = :brown
    dat = filter(x->x[2] ≈ g2, data)
    l2s = dat.|>x->x[1]
    m2s = dat.|>x->x[3]
    dms = dat.|>x->x[4]
    xlim = (1e-3, 5)
    # ylim = 2e-5.*(-1, 1) .+ extrema(m2s)
    ylim = (-5e-3, 0.015)

    f = plot([1, xlim[2]], [ylim[1], ylim[1]], 
    fillrange=[ylim[2], ylim[2]], 
    color=:gray, label=L"\langle \phi \rangle > M_{Pl}", fillalpha=0.3,
    xlabel=L"l^2", ylabel=L"m^2/M_{IR}^2", 
    legend=:bottomleft, 
    size=(600, 450), dpi=300)
    plot!(l2s, m2s, label=L"\gamma^2/\gamma^2_0="*"$(round(g2, digits=5))", 
    lw=2, color=c, xaxis=:log,
    xlim = xlim, ylim = ylim)
    plot!(x->m2s_analytical(x, g2), exp10.(range(-3, stop=xlim[2], length=40)), label=nothing, lw=2, color=c, ls=:dash)

    plot!(f, inset=(1, bbox(0.15, 0.05, 0.45, 0.5)))
    l2s = l2s[1:10]
    m2s = m2s[1:10]
    plot!(f[2], l2s, m2s, xaxis=:log, lw=2, label=nothing, fillalpha=0.5, color=c, alpha=0.5)
    plot!(f[2], x->m2s_analytical(x, g2), exp10.(range(log10(l2s[1]), stop=log10(l2s[end]), length=40)), label=nothing, lw=2, color=c, ls=:dash)

    g2 = g2s0[18]; c=:chocolate
    dat = filter(x->x[2] ≈ g2, data)
    l2s = dat.|>x->x[1]
    m2s = dat.|>x->x[3]
    dms = dat.|>x->x[4]

    plot!(l2s, m2s, label=L"\gamma^2/\gamma^2_0="*"$(round(g2, digits=5))", 
    lw=2, color=c, xaxis=:log)
    plot!(x->m2s_analytical(x, g2), exp10.(range(-3, stop=xlim[2], length=40)), label=nothing, lw=2, color=c, ls=:dash)

    l2s = l2s[1:10]
    m2s = m2s[1:10]
    plot!(f[2], l2s, m2s, xaxis=:log, lw=2, label=nothing, fillalpha=0.5, color=c, alpha=0.5)
    plot!(f[2], x->m2s_analytical(x, g2), exp10.(range(log10(l2s[1]), stop=log10(l2s[end]), length=40)), label=nothing, lw=2, color=c, ls=:dash)

    g2 = g2s0[17]; c=:brown2
    dat = filter(x->x[2] ≈ g2, data)
    l2s = dat.|>x->x[1]
    m2s = dat.|>x->x[3]
    dms = dat.|>x->x[4]

    plot!(l2s, m2s, label=L"\gamma^2/\gamma^2_0="*"$(round(g2, digits=5))", 
    lw=2, color=c, xaxis=:log)
    plot!(x->m2s_analytical(x, g2), exp10.(range(-3, stop=xlim[2], length=40)), label=nothing, lw=2, color=c, ls=:dash)

    l2s = l2s[1:10]
    m2s = m2s[1:10]
    plot!(f[2], l2s, m2s, xaxis=:log, lw=2, label=nothing, fillalpha=0.5, color=c, alpha=0.5)
    plot!(f[2], x->m2s_analytical(x, g2), exp10.(range(log10(l2s[1]), stop=log10(l2s[end]), length=40)), label=nothing, lw=2, color=c, ls=:dash)
    
    savefig("example/GoldbergerWise/figs/m2_l2_1.png")
end

solveODE_LR()
data[1:3]
filter(x->x[2] ≈ 1e2, data)
filter(x->x[1] ≈ 1e-1, data)
data = sort(data; lt = (x,y)->x[1]<=y[1]&&x[2]<y[2])
l2view = Dict("$l2"=>filter(x->x[1] ≈ l2, data) for l2 in exp10.(range(0,-3, 4)))
sort!(l2view["1.0"], by = x->x[2])
let ylim = (-3e-1, 1e-1), xlim = (1e-7, 1e2)
    f = plot([1e-1xlim[1],1],[ylim[1],ylim[1]], fillrange=[ylim[2],ylim[2]], label="tachyonic", fillalpha=0.3, color=:gray, lw=0)
    
    for val in values(l2view)
        g2s = val.|>x->x[2]
        m2s = val.|>x->x[3]
        l2 = val[1][1]
        plot!(g2s, m2s, label="l²=$l2", lw=2)
        # scatter!(g2s, m2s, label="l²=$l2", markersize=4)
    end
    plot!(xaxis=:log, ylim=ylim, xlim=xlim)
    plot!(f, inset=(1, bbox(0.1, 0.1, 0.4, 0.4)))
    val = l2view["0.001"]
    g2s = val.|>x->x[2]
    m2s = val.|>x->x[3]
    plot!(f[2], g2s, m2s, lw=2, xaxis=:log)
end

m20, f = search(1e0, 1e3, 0.9yₘ, true, 4); f
m20, f = search(5e0, 1e3, 0.9yₘ, true, 4); f
m20, f = search(1e1, 1e3, 0.9yₘ, true, 4); f
m20, f = search(5e1, 1e3, 0.9yₘ, true, 4); f
m20, f = search(1e2, 1e3, 0.9yₘ, true, 2, 4, 1e-8); f



using LaTeXStrings
let l2 = 1e-3, g2s(n) = exp10.(range(-1e-1, stop=3, length=n))
    g2s_num = g2s(10);
    m2s_num  = search.(l2, g2s_num, 0.98yₘ);
    g2s_ana = g2s(50);
    m2s_ana = m2s_analytical.(l2, g2s_ana);
    sto["m2g2"][l2] = g2s_num, m2s_num, g2s_ana, m2s_ana
    begin
        scatter(g2s_num, m2s_num, label="numerical", markersize=4, size=(800, 600), dpi=300)
        plot!(g2s_ana, m2s_ana, label="analytical", xaxis=:log, lw=2, color=:black)
        xlabel!(L"\gamma^2/\gamma^2_0")
        ylabel!(L"m^2/M_{IR}^2")
        # title!("Radion Mass with " * L"l^2=1\times10^{-3}" * "\n where " * L"M_{IR} = e^{-k y_m} M_{Pl},\quad \gamma^2_0 = 4k + 2u")
        savefig("example/GoldbergerWise/figs/m2_g2.png")
    end
end
let g2 = 1e3, l2s(n) = exp10.(range(-3, stop=log10(1), length=n))
    l2s_num = l2s(10);
    m2s_num  = search.(l2s_num, g2, 0.98yₘ);
    l2s_ana = l2s(50);
    m2s_ana = m2s_analytical.(l2s_ana, g2);
    sto["m2l2/$g2"] = l2s_num, m2s_num, l2s_ana, m2s_ana
    begin
        scatter(l2s_num, m2s_num, label="numerical", markersize=4, size=(800, 600), dpi=300)
        plot!(l2s_ana, m2s_ana, label="analytical", xaxis=:log, yaxis=:log, lw=2, color=:black)
        xlabel!(L"l^2")
        ylabel!(L"m^2/M_{IR}^2")
        # title!("Radion Mass with " * L"\gamma^2=1\times10^{3}\gamma^2_0" * "\n where " * L"M_{IR} = e^{-k y_m} M_{Pl},\quad \gamma^2_0 = 4k + 2u")
        savefig("example/GoldbergerWise/figs/m2_l2.png")
    end
end
let g2 = 1e1, l2s(n) = exp10.(range(-3, stop=log10(1), length=n))
    l2s_num = l2s(10);
    m2s_num  = search.(l2s_num, g2, 0.98yₘ);
    l2s_ana = l2s(50);
    m2s_ana = m2s_analytical.(l2s_ana, g2);
    sto["m2l2/$g2"] = l2s_num, m2s_num, l2s_ana, m2s_ana
    begin
        scatter(l2s_num, m2s_num, label="numerical", markersize=4, size=(800, 600), dpi=300)
        plot!(l2s_ana, m2s_ana, label="analytical", xaxis=:log, yaxis=:log, lw=2, color=:black)
        xlabel!(L"l^2")
        ylabel!(L"m^2/M_{IR}^2")
        title!("Radion Mass with " * L"\gamma^2=1\times10^{3}\gamma^2_0" * "\n where " * L"M_{IR} = e^{-k y_m} M_{Pl},\quad \gamma^2_0 = 4k + 2u")
        savefig("example/GoldbergerWise/figs/m2_l2.png")
    end
end
sto = storage("data")
@pack! sto = [1,2]
sto

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