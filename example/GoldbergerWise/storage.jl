# store data
using JLD2
using UnPack
function storage(filename, op="a+")
    @assert filename isa String "filename is not a string"
    if !endswith(filename, ".jld2")
        filename = filename * ".jld2"
    end
    jldopen(filename, op)
end