mutable struct FixedStack{T}
    data::Array{T,1}
    next_index::UInt32
    length::UInt32

    FixedStack{T}(size::Int) where {T} = new{T}(Array{T,1}(undef, UInt32(size)), 1, 0)
end

function push!(stack::FixedStack{T}, item::T) where {T}
    stack.length < length(stack.data) && stack.length += 1
    stack.data[stack.next_index] = item
    stack.next_index = mod1(stack.next_index + 1, length(stack.data))
    return stack
end
function push!(stack::FixedStack{T}, item::T; trushbin::T0) where {T, T0<:AbstractArray{T,1}}
    stack.length < length(stack.data) && stack.length += 1
    push!(trushbin, stack.data[stack.next_index])
    stack.data[stack.next_index] = item
    stack.next_index = mod1(stack.next_index + 1, length(stack.data))
    return stack
end

function pop!(stack::FixedStack{T}) where {T}
    if stack.length <= 0
        @warn "trying to pop from an empty stack"
        return nothing
    end
    stack.length -= 1
    stack.next_index = mod1(stack.next_index - 1, length(stack.data))
    return stack.data[stack.next_index]
end