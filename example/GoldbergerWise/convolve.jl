# https://www.math.pku.edu.cn/teachers/lidf/docs/Julia/html/_book/exa-science.html#exasci-convimg2d
"""
    extend(M, i, j)

extend matrix M with its boundary values
"""
function extend(M::AbstractMatrix, i, j)
    i0 = firstindex(M, 1)
    i1 = lastindex(M,1)
    j0 = firstindex(M, 2)
    j1 = lastindex(M,2)

    ii = min(max(i, i0), i1)
    jj = min(max(j, j0), j1)
    
    return M[ii,jj]
end

"""
    convolve(M, K)

convolve matrix M with kernel K
"""
function convolve(M::AbstractMatrix, K::AbstractMatrix)
    l = (size(K,1)-1) ÷ 2
    [ sum([extend(M, s-i, t-j) for i=(-l):l, j=(-l):l] .* K) 
        for s in axes(M,1), t in axes(M,2)]
end

"""
    with_sobel_edge_detect(image)

edge detection with sobel operator
"""
function with_sobel_edge_detect(image)
    Kx = [1 0 -1; 2 0 -2; 1 0 -1]
    Ky = Kx'
    Gx = convolve(image, Kx)
    Gy = convolve(image, Ky)
    G = sqrt.(Gx .^2 + Gy .^2)
    
    return G
end

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