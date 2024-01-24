using DifferentialEquations: ODESolution, ODEProblem, SecondOrderODEProblem, solve, ImplicitEuler
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

function errBCwithφ(params; F0′, F0=1., y0 = yₘ/2)
    Fsol = solveODE_LR([F0′, F0], params, y0)
    FT′,FT = Fsol(yₘ)
    FP′,FP = Fsol(0)
    errR = Δφ′T(FT, FT′, params)
    errL = Δφ′P(FP, FP′, params)
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

    for v in ((x1,y1),(x1,y2),(x2,y2),(x2,y1))
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

# f(x,y) = x+y,y-x
# f(0,1) .> 0
# f(1,0) .> 0
# f(0,-1) .> 0
# f(-1,0) .> 0
# bisection2d(f, ((1,0),(-1,0),(0,1),(0,-1) ))

function searchYM(l2, g2, k=k)
    (Fp, m2) -> errBCwithφ((l2, k, g2, m2*M_IR^2); F0′=Fp, F0=1., y0 = 0.9*yₘ)
end

prob = searchYM(1e-3, 1e1)

n = 40
X = exp10.(range(-10, stop=2, length=n))*ones(n)';
Y = ones(n)*exp10.(range(-7, stop=1, length=n))';
Z = prob.(X,Y) 
findfirst(x-> (x.>0) == (false,false), Z)
findfirst(x-> (x.>0) == (false,true), Z)
findfirst(x-> (x.>0) == (true,false), Z)
findfirst(x-> (x.>0) == (true,true), Z)
Z0 = map(Z) do (z1, z2)
    (sign(z1) - sign(z2))/2
end
contourf(Z0)
idx = CartesianIndex(33, 16)
X[idx], Y[idx], Z[idx].>0
idx = CartesianIndex(1, 16)
X[idx], Y[idx], Z[idx].>0

using Plots



let F0=[1e1, 1], y0 = 0.8yₘ, params = (1e-3,k,1e1,1e-4M_IR^2)
    sol = solveODE_LR(F0, params, y0)
    plot(y->sol(y)[2] - exp(2k*(y-y0)), 0, yₘ, label="F")
end

# let F0=[1e-1, 1], y0 = 0.8yₘ, params = (1,1,1,1)
#     sol = solveODE_LR(F0, params, y0)
#     plot(y->sol(y)[1], 0, yₘ, label="F'")
# end

# let F0=[1e-1, 1], y0 = 0.8yₘ, params = (1,1,1,1)
#     errBCwithφ(params; F0′=F0[1], F0=F0[2], y0=y0)
# end
