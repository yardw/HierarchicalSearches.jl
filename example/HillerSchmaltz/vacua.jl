using IntervalArithmetic, IntervalRootFinding
roots(x -> x^2 - 2x, 0..10)
roots(x -> (x^2 - 2)^2 * (x^2 - 3), 0..10)
using StaticArrays
let g((x1,x2,x3)) = @SVector[
                          x1^2 + x2^2 + x3^2 - 1,
                          x1^2 + x3^2 - 0.25,
                          x1^2 + x2^2 - 4x3
                            ],
    x = 0..10, X = x × x × x
    @show roots(g, X)
end

let g(x) = n=length(x);@SVector[-x[i]^2 - 1 for i in 1:n],
    x0 = 0..10, X = x0 × x0 × x0 × x0
    @show roots(g, X)
end
let g0(n::Int) = x->SA[-x[i]^2 + 1 for i in 1:n]
    g = g0(4)
    x0 = 0..10, X = x0 × x0 × x0 × x0
    @show roots(g, X)
end

using Plots
let x=-1..1, X=x×x, g(x) = SA[x[1]^2 - 1, x[2]^2 - 0.25]
    rts = roots(g, X)    
    meanpoints = mid.(interval.(rts))
    xs = first.(meanpoints)
    ys = last.(meanpoints)
    scatter(xs, ys, label="roots")
end
let m=2 , n=2
    using Plots; gr()

    xp = [0, 0, 0, 0, 1, 1, 1, 1]
    yp = [0, 1, 0, 1, 0, 0, 1, 1]
    zp = [0, 0, 1, 1, 1, 0, 0, 1]
    connections = [(1,2,3), (4,2,3), (4,7,8), (7,5,6), (2,4,7), (1,6,2), (2,7,6), (7,8,5), (4,8,5), (4,5,3), (1,6,3), (6,3,5)][1:m]
    
    xe = [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0][1:n]
    ye = [0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1][1:n]
    ze = [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1][1:n]
    
    plot(xe,ye,ze; lc=:black, lw=0.5, lims=(-0.25,1.25))
    # scatter!((0.5,0.5,0.5); c=:red, ms=6, msw=0.1)
    mesh3d!(xp,yp,zp; connections, proj_type=:persp, fc=:lime, lc=:lime, fa=0.3, lw=0)
end

let x=1..2

    plot(x::Interval) = 
end