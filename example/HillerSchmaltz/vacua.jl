using IntervalArithmetic, IntervalRootFinding
using StaticArrays
cofactor(A,i,j) = (-1)^(i+j) * det(A[setdiff(1:end,i), setdiff(1:end,j)])
function det(A)
    @assert ndims(A) == 2 "A must be a matrix"
    m, n = size(A)
    if m <= 1 || n <= 1
        return A[1,1]
    elseif m > n
        return sum([A[1,j] * cofactor(A,1,j) for j in axes(A,2)])
    else
        return sum([A[i,1] * cofactor(A,i,1) for i in axes(A,1)])
    end
end
# test for cofactor
A = [100i^2+j^2 for i in 1:4, j in 1:4]
[cofactor(A,i,j) for i in 1:4, j in 1:4] 
det([1 for i in 1:1, j in 1:1])
det([1])
det([1 2; 3 4]) == -2
f(x) = SA[
    x[1]^2 + x[2]^2 + x[3]^2 - 1,
    x[1]^2 + x[3]^2 - 0.25,
    x[1]^2 + x[3]^2 - 0.25
]
eq(x,i) = x[i]^2 + x[mod1(i+1, length(x))]^2 - 1
f(x) = @SVector[ eq(x,i)  for i in 1:3]
itvs(x::Interval...) = foldl(×, x)
x = -2..2
X = x × x × x
roots(f, X)
roots(x -> x^2 - 2x, 0..10)
roots(x -> (x^2 - 2)^2 * (x^2 - 3), 0..10)
let g((x1,x2,x3)) = @SVector[
                          x1^2 + x2^2 + x3^2 - 1,
                        #   x1^2 + x3^2 - 0.25,
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