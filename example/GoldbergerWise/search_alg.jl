include("convolve.jl")
function rescale(xscale, yscale, x0, y0, hypersurfaces, scale = 20)
    k1, b1, kx1pb_style, k2, b2, kx2pb_style = hypersurfaces
    xrange = [x0-xscale/scale, x0+xscale/scale]
    yrange = [y0-yscale/scale, y0+yscale/scale]
    if kx1pb_style
        append!(yrange, k1 .* xrange .+ b1)
    else
        append!(xrange, k1 .* yrange .+ b1)
    end
    if kx2pb_style
        append!(yrange, k2 .* xrange .+ b2)
    else
        append!(xrange, k2 .* yrange .+ b2)
    end
    # for (slope, bias) in zip(slopes, biass)
    #     append!(yrange, slope .* xrange .+ bias)
    # end
    # yrange = [minimum(yrange), maximum(yrange)]
    # xrange = extrema(xrange)
    # yrange = extrema(yrange)
    if scale > 1
        xscale = 1.1maximum(abs.(x0 .- xrange))
        yscale = 1.1maximum(abs.(y0 .- yrange))
    else
        xscale = 1.1minimum(abs.(x0 .- xrange))
        yscale = 1.1minimum(abs.(y0 .- yrange))
    end
    return xscale, yscale
end
function sample(prob, x0, xsc, y0, ysc, n)
    x = x0 .+ xsc.*range(-1, 1, n)
    y = y0 .+ ysc.*range(-1, 1, n)
    X = x*ones(n)'
    Y = ones(n)*y'
    Z = prob.(X, Y)
    return x, y, Z
end
function fit(x, y, Z, isdebug=false)
    Z1 = map(first, Z) .|> sign
    Z2 = map(last, Z)  .|> sign
    # return Z1 .+ 2Z2 .+ 3
    y1, x1 = pick_xy((Z1 |> with_sobel_edge_detect) .> 0, x, y)
    y2, x2 = pick_xy((Z2 |> with_sobel_edge_detect) .> 0, x, y)
    @assert length(y1) > 1 && length(y2) > 1 "Unable to find the hypersurface within the given parameter region.\nx: $(x|>extrema)\ny: $(y|>extrema)" 
    #solve for a*x + b*y = 1
    # a1, b1 = [x1 y1]\ones(length(x1))
    # a2, b2 = [x2 y2]\ones(length(x2))
    # x0, y0 = [a1 b1; a2 b2]\[1; 1]
    #solve for x = k*y + b
    if abs(-(extrema(x1)...)) < 1e-4abs(-(extrema(y1)...))   
        kx1pb_style = false
        k1, b1 = [y1 ones(length(y1))]\x1
    else
        kx1pb_style = true
        k1, b1 = [x1 ones(length(x1))]\y1
    end
    if abs(-(extrema(x2)...)) < 1e-4abs(-(extrema(y2)...))   
        kx2pb_style = false
        k2, b2 = [y2 ones(length(y2))]\x2
    else
        kx2pb_style = true
        k2, b2 = [x2 ones(length(x2))]\y2
    end
    hypersurfaces = k1, b1, kx1pb_style, k2, b2, kx2pb_style
    if kx1pb_style == kx2pb_style && abs(k2-k1) <= 1e-4abs(b2 - b1) 
        isdebug && println("Two hypersurfaces are almost parallel.")
        return NaN, NaN, hypersurfaces
    end
    if kx1pb_style
        if kx2pb_style
            y0, x0 = [1 -k1; 1 -k2]\[b1; b2]
        else
            y0, x0 = [1 -k1; -k2 1]\[b1; b2]
        end
    else
        if kx2pb_style
            y0, x0 = [-k1 1; 1 -k2]\[b1; b2]
        else
            y0, x0 = [-k1 1; -k2 1]\[b1; b2]
        end
    end
    return y0, x0, hypersurfaces
end
function bisection(f, a, b, tol=1e-9, max_iter=100)
    fa = f(a)
    fb = f(b)
    
    @assert fa * fb < 0 "The function values at the interval endpoints must have opposite signs."
    
    c = (a + b) / 2
    fc = f(c)
    
    iter = 0
    while abs(fc) > tol && iter < max_iter
        if fa * fc < 0
            b = c
            fb = fc
        else
            a = c
            fa = fc
        end
        
        c = (a + b) / 2
        fc = f(c)
        
        iter += 1
    end
    
    return c
end