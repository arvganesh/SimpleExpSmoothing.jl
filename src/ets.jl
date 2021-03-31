#=
ETS model class
- fit, predict
- Constructor:
=#

abstract type ETSModel end

mutable struct ExponentialSmoothing <: ETSModel
    y::Vector{T} where {T <: AbstractFloat}
    h::Integer
    alpha::Union{Nothing, Float64}
    beta::Union{Nothing, Float64}
    init_level::Union{Nothing, Float64}
    function ExponentialSmoothing(y = []; h = 10, alpha = nothing, beta = nothing, init_level = nothing) # Constuctor 
        # Input Cleaning
        if length(y) <= 0
            error("The input array is empty.")
        end
        if h <= 0
            error("The number of prediction steps must be positive")
        end

        cleanParams_(alpha, "alpha")
        cleanParams_(beta, "beta")

        return new(y, h, alpha, beta, init_level)
    end 
end

function fit(model::ETSModel) # Set alpha + init_level
    if (isnothing(model.alpha))
        @warn "Since no value was entered for alpha, it will be chosen"
    end
    if (isnothing(model.init_level))
        @warn "Since no value was entered for init_level, it will be chosen"
    end
    if (isnothing(model.beta))
        @warn "Since no value was entered for beta, it will be chosen"
    end

    # Replace with heuristic approach to set them
    model.alpha = isnothing(model.alpha) ? 0.25 : model.alpha
    model.init_level = isnothing(model.init_level) ? model.y[1] : model.init_level
    return
end

function cleanParams_(param::Union{Nothing, Float64}, name::String, lb::Integer = 0, ub::Integer = 1)
    if !isnothing(param) && (param < lb || param > ub)
        error(name * " needs to be in the range ($lb, $ub)")
    end
end


function predict(model::ETSModel) # Return vector of fitted values of length h
    h = model.h
    alpha = model.alpha
    init_level = model.init_level
    y = model.y
    data_size = length(y)

    lvls = zeros(data_size + 1) # Stores fitted values
    lvls[1] = init_level # Set l_0

    for i in 2:(length(lvls)) # Compute l_t's
        lvls[i] = alpha * y[i-1] + (1 - alpha) * lvls[i-1];
    end

    for i in 1:(h-1)
        append!(lvls, lvls[data_size + 1])
    end
    return lvls
end