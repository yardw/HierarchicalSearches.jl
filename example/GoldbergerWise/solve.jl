using DifferentialEquations: ODESolution, ODEProblem, SecondOrderODEProblem, solve, ImplicitEuler
include("eom.jl")
using .GoldbergerWiseEoM.EoMs
import .GoldbergerWiseEoM.Consts: yₘ, k, M_IR, γ²₀, u
using .GoldbergerWiseEoM.AffliatedFunctions: ϕ0, A, A′
using .GoldbergerWiseEoM.BCs
# params = l², k, γ², m²
# k = 37u = k 

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

function errBCwithφ(params; F0′, F0=1., y0 = yₘ/2)
    Fsol = solveODE_LR([F0′, F0], params, y0)
    FT′,FT = Fsol(yₘ)
    FP′,FP = Fsol(0)
    errR = Δφ′T(FT, FT′, params)
    errL = Δφ′P(FP, FP′, params)
    return errL, errR
end
function errBCwithφtest(params; F0′, F0=1., y0 = yₘ/2)
    Fsol = solveODE_LR([F0′, F0], params, y0)
    FT′,FT = Fsol(yₘ)
    FP′,FP = Fsol(0)
    errR = Δφ′Ttest(FT, FT′, params)
    errL = Δφ′Ptest(FP, FP′, params)
    return errL, errR
end


function isparallel(e1, e2)
    @show e1, e2
    vec1 = e1[2] .- e1[1]
    vec2 = e2[2] .- e2[1]
    return vec1[1]*vec2[2] ≈ vec1[2]*vec2[1]
end
function maxdist(vs)
    es = ((v1, v2) for (i,v1) in enumerate(vs), (j,v2) in enumerate(vs) if i<j)
    emins = Dict()
    for e in es
        vec = (e[1] .- e[2])
        dis = abs2.(vec) |> sum
        dir = iszero(vec[1]) ? (0,1.) : (1,vec[2]/vec[1])
        if haskey(emins, dir) 
            if dis < emins[dir][1]
                emins[dir] = (dis, e)
            end
        else
            emins[dir] = (dis, e)
        end
    end
    edis, emax = emins[findmax(x->values(x)[1], emins)[2]]
    return emax, edis
end
function bisection2d(f,(x1,x2),(y1,y2), tol=1e-8)
    level = Dict((i,j)=>(NaN,NaN) for i in (true, false), j in (true, false))
    flag  = Dict((i,j)=>false for i in (true, false), j in (true, false))
    function updatelevel(v, level=level, flag=flag)
        key = Tuple(f(v[1],v[2]) .> 0)
        level[key] = v
        flag[key] = true
    end
    x0 = (x1+x2)/2
    y0 = (y1+y2)/2
    for v in ((x0,y1),(x0,y2),(x1,y0),(x2,y0))
        print(sign.(f(v[1],v[2])))
        @show v
        updatelevel(v)
    end
    valid = sum(values(flag)) > 2
    @assert valid "two or more initial points are of the same sign"
    loopbound = 10
    vertices = (level[key] for key in keys(level) if flag[key])
    emax, edis = maxdist(vertices)
    while edis > tol^2
        if loopbound == 0
            @warn "loopbound reached"
            break
        end
        loopbound -= 1
          
        # get new verteces by bisection
        vertices = (level[key] for key in keys(level) if flag[key])
        emax, edis = maxdist(vertices)
        updatelevel(0.5 .* (emax[1] .+ emax[2]))
        @show values(level)
    end
    @show collect(vertices)
    return v = reduce(.+, vertices)./sum(values(flag))
end
function bisection2d(f, vertices, tol=1e-8)
    level = Dict((i,j)=>(NaN,NaN) for i in (true, false), j in (true, false))
    flag  = Dict((i,j)=>false for i in (true, false), j in (true, false))
    function updatelevel(v, level=level, flag=flag)
        key = Tuple(f(v[1],v[2]) .> 0)
        level[key] = v
        flag[key] = true
    end

    for v in vertices
        updatelevel(v)
    end
    valid = sum(values(flag)) > 3
    @assert valid "two or more initial points are of the same sign"
    loopbound = 100
    vertices = (level[key] for key in keys(level) if flag[key])
    emax, edis = maxdist(vertices)
    while edis > tol^2
        if loopbound == 0
            @warn "loopbound reached"
            break
        end
        loopbound -= 1
          
        # get new verteces by bisection
        vertices = (level[key] for key in keys(level) if flag[key])
        emax, edis = maxdist(vertices)
        updatelevel(0.5 .* (emax[1] .+ emax[2]))
        # @show values(level)
    end
    @show collect(vertices)
    return v = reduce(.+, vertices)./sum(values(flag))
end


include("convolve.jl")

"""
    pick_xy(Z, x, y)

pick all data points of true in Z, return the coordinates with the axes x and y
"""
function pick_xy(Z, x, y)
    idxs = findall(x-> x>0, Z)
    idx1s = map(x->x[1], idxs)
    idx2s = map(x->x[2], idxs)
    return x[idx1s], y[idx2s]
end
# f(x,y) = x+y,y-x
# f(0,1) .> 0
# f(1,0) .> 0
# f(0,-1) .> 0
# f(-1,0) .> 0
# bisection2d(f, ((1,0),(-1,0),(0,1),(0,-1) ))

function searchYM(l2, g2, k=k)
    (Fp, m2) -> errBCwithφ((l2, k, g2*γ²₀, m2*M_IR^2); F0′=Fp, F0=1., y0 = 0.99*yₘ)
end
function searchYMtest(l2, g2, k=k)
    (Fp, m2) -> errBCwithφtest((l2, k, g2*γ²₀, m2*M_IR^2); F0′=Fp, F0=1., y0 = 0.99*yₘ)
end

function search(l2=1e-3, g2=1e8)
    # both Fp and m2 leaner scaled
    prob = searchYMtest(l2, g2)
    n = 20
    Fp = 7.4 - 3.55e-5 .+ 1e-6range(-1, stop=1, length=n)
    m2 = range(-5e-5, stop=1e-4, length=2n)
    X = Fp*ones(2n)';
    Y = ones(n)*m2';
    Z = prob.(X,Y);
    # extract the hypersurface constrained by BCs
    Z1 = map(first, Z) .|> sign
    Z2 = map(last, Z)  .|> sign
    Z3 = Z1+Z2 |> with_sobel_edge_detect
    z1s = pick_xy((Z1 |> with_sobel_edge_detect) .> 0, Fp, m2)
    slope1, bias1 = [z1s[2] ones(length(z1s[2]))]\z1s[1]
    z2s = pick_xy((Z2 |> with_sobel_edge_detect) .> 0, Fp, m2)
    slope2, bias2 = [z2s[2] ones(length(z2s[2]))]\z2s[1]
    Fp0, m20 = [1 -slope1; 1 -slope2]\[bias1; bias2]
    # let 
    #     heatmap(m2, Fp, Z3)
    #     scatter!([m20], [Fp0], label="critical point")
    #     xlabel!("m²/M_IR²")
    #     ylabel!("F′(0.99yₘ)")
    #     title!("l²=1e-3, γ²=1e8γ²₀")
    # end
    return m20
end

search(1e-3, 1e6)
m2s_analytical(1e-3, 1e6)
search(1e-3, 1e7)
m2s_analytical(1e-3, 1e7)
search(1e-3, 1e8)
m2s_analytical(1e-3, 1e8)

m2s_numerical(g2) = search(1e-3, g2)

plot([m2s_numerical, m2s_analytical], exp10.(range(5.5,10,10)), label=["numerical" "analytical"], xlabel="γ²/γ²₀", ylabel="m²/M_IR²", title="l²=1e-3", xaxis=:log10)

function m2s_analytical(l2, g2)
    return m2s = 4l2*(2k+u)*u^2/(3k)*(1-exp(2k*yₘ))/(1-exp((4k+2u)*yₘ)) * (1 .- 1 /g2) / M_IR^2
end
m2s_analytical(g2) = m2s_analytical(1e-3, g2)

let l2 = 1e-3, g2 = 1e6
    # both Fp and m2 leaner scaled
    prob = searchYMtest(l2, g2)
    n = 20
    Fp = 7.4 - 3.55e-5 .+ 1e-3range(-1, stop=1, length=n)
    m2 = range(-5e-5, stop=1e-4, length=2n)
    X = Fp*ones(2n)';
    Y = ones(n)*m2';
    Z = prob.(X,Y);
    # Z0 = map(Z) do (z1, z2)
    #     2sign(z1) + sign(z2)+3
    # end;
    # extract the hypersurface constrained by BCs
    Z1 = map(first, Z) .|> sign
    Z2 = map(last, Z)  .|> sign
    Z3 = Z1+Z2 |> with_sobel_edge_detect
    # z1s = pick_xy((Z1 |> with_sobel_edge_detect) .> 0, Fp, m2)
    # slope1, bias1 = [z1s[2] ones(length(z1s[2]))]\z1s[1]
    # z2s = pick_xy((Z2 |> with_sobel_edge_detect) .> 0, Fp, m2)
    # slope2, bias2 = [z2s[2] ones(length(z2s[2]))]\z2s[1]
    # Fp0, m20 = [1 -slope1; 1 -slope2]\[bias1; bias2]
    let 
        heatmap(m2, Fp, Z3)
        # scatter!([m20], [Fp0], label="critical point")
        xlabel!("m²/M_IR²")
        ylabel!("F′(0.99yₘ)")
        title!("l²=$(l2), γ²=$(g2)γ²₀")
    end
end










prob = searchYMtest(1e-3, 1e2)
n = 40
# lgFp = range(-1, stop=3, length=n)
Fp = range(-1e2, stop=1e2, length=n)
lgm2 = range(0, stop=4, length=n)
X = Fp*ones(n)';
Y = ones(n)*exp10.(lgm2)';
Z = prob.(X,Y) 
findfirst(x-> (x.>0) == (false,false), Z)
findfirst(x-> (x.>0) == (false,true), Z)
findfirst(x-> (x.>0) == (true,false), Z)
findfirst(x-> (x.>0) == (true,true), Z)
Z0 = map(Z) do (z1, z2)
    (2sign(z1) + sign(z2)+3)/2
end;
using  Plots
heatmap(lgm2, Fp, Z0)

# l2 = 1e-3, g2 = 1e2, m2 = 10^2.5 M_IR^2
Fp = range(10, stop=25, length=n)
lgm2 = range(2, stop=3, length=n)
X = Fp*ones(n)';
Y = ones(n)*exp10.(lgm2)';
Z = prob.(X,Y);
Z0 = map(Z) do (z1, z2)
    2sign(z1) + sign(z2)+3
end;
heatmap(lgm2, Fp, Z0)

prob = searchYMtest(1e0, 1e2)
Fp = range(10, stop=25, length=n)
lgm2 = range(2, stop=3, length=n)
X = Fp*ones(n)';
Y = ones(n)*exp10.(lgm2)';
Z = prob.(X,Y);
Z0 = map(Z) do (z1, z2)
    2sign(z1) + sign(z2)+3
end;
heatmap(lgm2, Fp, Z0)

# Fp linear, m2 log scaled
prob = searchYMtest(1e-3, 1e8)
n = 40
Fp = 7.4 - 3.55e-5 .+ 1e-6range(-1, stop=1, length=n)
lgm2 = range(-6, stop=-4, length=n)
X = Fp*ones(n)';
Y = ones(n)*exp10.(lgm2)';
Z = prob.(X,Y);
Z0 = map(Z) do (z1, z2)
    2sign(z1) + sign(z2)+3
end;
let 
    heatmap(lgm2, Fp, Z0, clims=(0,6))
    xlabel!("lg m²/M_IR²")
    ylabel!("F′(0.99yₘ)")
    title!("l²=1e-3, γ²=1e8γ²₀")
end






heatmap(m2, Fp, Z3)
idxs = findall(x-> x>0.7maximum(Z3), Z3)
idx1s = map(x->x[1], idxs)
idx2s = map(x->x[2], idxs)
Fp = range(Fp[minimum(idx1s)], stop=Fp[maximum(idx1s)], length=n)
m2 = range(m2[minimum(idx2s)], stop=m2[maximum(idx2s)], length=n)


prob = searchYM(1e-3, 1e2)
n = 40
lgFp = range(-1, stop=5, length=n)
lgm2 = range(-0, stop=4, length=n)
X = exp10.(lgFp)*ones(n)';
Y = ones(n)*exp10.(lgm2)';
Z = prob.(X,Y) 
findfirst(x-> (x.>0) == (false,false), Z)
findfirst(x-> (x.>0) == (false,true), Z)
findfirst(x-> (x.>0) == (true,false), Z)
findfirst(x-> (x.>0) == (true,true), Z)
Z0 = map(Z) do (z1, z2)
    2sign(z1) + sign(z2) + 3
end;
heatmap(lgm2, lgFp, Z0)

idx = CartesianIndex(33, 16)
X[idx], Y[idx], Z[idx].>0
idx = CartesianIndex(1, 16)
X[idx], Y[idx], Z[idx].>0



let F0=[1e1, 1], y0 = 0.8yₘ, params = (1e-3,k,1e1,1e-4M_IR^2)
    sol = solveODE_LR(F0, params, y0)
    plot(y->sol(y)[2] - exp(2k*(y-y0)), 0, yₘ, label="F")
end

# let F0=[1e-1, 1], y0 = 0.8yₘ, params = (1,1,1,1)
#     sol = solveODE_LR(F0, params, y0)
#     plot(y->sol(y)[1], 0, yₘ, label="F'")
# end

let F0=[1e-1, 1], y0 = 0.8yₘ, params = (1,1,1,1)
    errBCwithφ(params; F0′=F0[1], F0=F0[2], y0=y0)
end
let F0=[1e-1, 1], y0 = 0.8yₘ, params = (1,1,1,1)
    errBCwithφtest(params; F0′=F0[1], F0=F0[2], y0=y0)
end
let F0=[1e-1, 1], y0 = 0.8yₘ, params = (1e-1,k,1e1γ²₀,1e-3M_IR^2)
    errBCwithφ(params; F0′=F0[1], F0=F0[2], y0=y0)
end
let F0=[1e-1, 1], y0 = 0.8yₘ, params = (1e-1,k,1e1γ²₀,1e-3M_IR^2)
    errBCwithφtest(params; F0′=F0[1], F0=F0[2], y0=y0)
end
let F0=[1e1, 1], y0 = 0.8yₘ, params = (1e-1,k,1e1γ²₀,1e-3M_IR^2)
    errBCwithφ(params; F0′=F0[1], F0=F0[2], y0=y0)
end
let F0=[1e1, 1], y0 = 0.8yₘ, params = (1e-1,k,1e1γ²₀,1e-3M_IR^2)
    errBCwithφtest(params; F0′=F0[1], F0=F0[2], y0=y0)
end
let F0=[1e2, 1], y0 = 0.8yₘ, params = (1e-1,k,1e1γ²₀,1e-3M_IR^2)
    errBCwithφ(params; F0′=F0[1], F0=F0[2], y0=y0)
end
let F0=[1e2, 1], y0 = 0.8yₘ, params = (1e-1,k,1e1γ²₀,1e-3M_IR^2)
    errBCwithφtest(params; F0′=F0[1], F0=F0[2], y0=y0)
end
exit()