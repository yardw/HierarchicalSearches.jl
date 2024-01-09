module ModuliSpaces

abstract type AbstractClass end
import Base.:∈
Base.:∈(y, Y0::T) where {T<:AbstractClass}  = println("$(Y0) is an AbstractClass")

function priori_dfs(f::Function, X::T1, Y0, x0::T2) where {T1<:AbstractArray, T2}
    # create a closure for the function priori_dfs!
    get_priority(x) = (f(x) ∈ Y0)|>sum
    function init_stack_from(x0)
        list = [x0]
        sizehint!(list, sum(size(X)))
        return list
    end
    function init_path()
        list = T2[]
        sizehint!(list, sum(size(X)))
        return list
    end
    N = ndims(X)
    adjacents = (CartesianIndex{N}([n==n0 ? sgn : 0 for n in 1:N]...) for sgn in (1,-1) for n0 in 1:N)
    adjacent(x0) = (x0 + x for x in adjacents)
    isvisited = falses(size(X)...)
    paths = []

    function priori_dfs(x0; goal_priority=0)
        isvisited[x0] = true # mark x0 as visited
        path = init_path() # initialize a new path starting with x0
        
        # start searching
        stack = init_stack_from(x0) # initialize a stack with x0
        while !isempty(stack) # while the stack is not empty
            x = pop!(stack) # pop the last element from the stack, and name it as x
            push!(path, x) # push x to the path
            isendpoint = true # indicate that current position is an endpoint of path
    
            # push all valid adjacent vertices of x into the stack, and check if x is an endpoint
            for y in adjacent(x) # for all adjacent vertices of x
                # check if y is reachable
                if !checkbounds(Bool, isvisited, y) continue end # if y is out of bound, skip it
                if isvisited[y] continue end # if y has been visited, skip it
                isvisited[y] = true # mark y as visited
                isendpoint = false # since y is reachable, the current position is not an endpoint
    
                # check if y is a goal
                priority = get_priority(y) # compute the priority of y
                if priority < goal_priority continue end # if the priority of y is smaller than the global priority, skip it 
                if priority > goal_priority # if the priority of y is larger than the global priority
                    push!(paths, (path=path, goal_priority=goal_priority)) # archive the current path
                    goal_priority = priority # update the global priority
                    stack = init_stack_from(y) # initialize a new stack starting with y
                    path  = init_path() # initialize a new path starting with y
                    break # stop searching for adjacent vertices of x   
                end
    
                # since y have to be the goal, we just push it to the current stack
                push!(stack, y) # push y to the stack
            end
    
            if isendpoint # if the current position is an endpoint
                push!(paths, (path=path, goal_priority=goal_priority)) # archive the current path
                path = init_path() # initialize a new path starting with x
            end
        end
        
        return paths, isvisited
    end
    return priori_dfs
end

end