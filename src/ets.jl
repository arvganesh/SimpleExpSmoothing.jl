abstract type ETSModel end
import Optim: optimize, LBFGS, Options, minimizer, Fminbox
using LinearAlgebra
import LineSearches

"""
ExponentialSmoothing

Description:
    - Contains parameters and data to apply Simple Exponential Smoothing (SES)
    - Subtype of ETSModel, for future proofing

Constuctor:
    - Vector{AbstractFloat} y: contains observations for forecasting
    - Integer h: Stores the number of steps after time T to forecast (T = latest time in the observed data)
        - Default = 10, must be > 0
    - Union{Nothing, Float64} alpha: smoothing parameter [1]
        - 0 < alpha < 1
        - If no value is specified, an optimal value will be chosen using Brent's method.
    - Union{Nothing, Float64} init_level: initial level value [1]
        -Inf < init_level < Inf
        - If no value is specified, an optimal value will be chosen uisng a heuristic method.

References: 
    [1] Hyndman, R.J., & Athanasopoulos, G. (2019) *Forecasting:
        principles and practice*, 3rd edition, OTexts: Melbourne,
        Australia. OTexts.com/fpp3. Accessed on April 19th 2020.

EXAMPLE USAGE:
observations = [1.0, 2.0, 3.0]

mdl = ExponentialSmoothing(observations, h = 5)
mdl = ExponentialSmoothing(observations, alpha = 0.4)
mdl = ExponentialSmoothing(observations, alpha = 0.25, init_level = 500.0)
mdl = ExponentialSmoothing(observations, h = 15, alpha = 0.3, init_level = 750.0)

"""
mutable struct ExponentialSmoothing <: ETSModel
    y::Vector{AbstractFloat}
    h::Integer
    alpha::Union{Nothing, Float64}
    init_level::Union{Nothing, Float64}
    function ExponentialSmoothing(y = []; h = 10, alpha = nothing, init_level = nothing) # Constuctor 
        # Input Cleaning
        if length(y) <= 0
            error("The input array is empty.")
        end
        if h <= 0
            error("The number of prediction steps must be positive")
        end

        cleanParams_(alpha, "alpha")

        return new(y, h, alpha, init_level)
    end 
end
"""
predict_()

Description:
    - Helper function
    - computes a forecast with Simple Exponential Smoothing (SES) [1]

Parameters
    - model: ExponentialSmoothing struct (defined above)
    - Note: when predict_ is called, all parameters have been calibrated.

Returns:
    A tuple containing the following:
    - Vector{Float64} forecast: contains fitted values
    - SSE::Float64 SSE: the Sum of the Squares (error term)
        - used for parameter optimzation

References: 
    [1] Hyndman, R.J., & Athanasopoulos, G. (2019) *Forecasting:
        principles and practice*, 3rd edition, OTexts: Melbourne,
        Australia. OTexts.com/fpp3. Accessed on April 19th 2020.
"""
function predict_(model::ETSModel) # Return vector of fitted values of length h
    h = model.h
    alpha = model.alpha
    init_level = model.init_level
    y = model.y
    data_size = length(y)

    lvls = zeros(data_size + 1) # Stores smoothed values
    forecast = zeros(data_size + h) # Stores fitted values
    lvls[1] = init_level # Set l_0
    forecast[1] = init_level # Set yhat_0
    SSE = 0.0 # Compute Sum of Squared Errors

    for i in 2:(length(lvls)) # Compute l_t's
        lvls[i] = alpha * y[i-1] + (1 - alpha) * lvls[i-1];
        forecast[i] = lvls[i]
        SSE += (forecast[i] - y[i-1])^2
    end

    for i in 1:(h-1)
        forecast[length(lvls) + i] = lvls[end]
    end

    return forecast, SSE
end

"""
SSE_()

- model: ExponentialSmoothing struct (defined above)
- alpha: value of alpha use for computation

Note: when SSE_ is used, alpha ≂̸ nothing, so model.alpha is always defined.

"""
function SSE_(model::ETSModel; alpha = nothing, verbose = 0) # Helper function for optimization
    # If parameter values have been specified, keep them constant.
    model.alpha = alpha
    res = predict_(model)[2]
    if verbose > 0
        println("Error: $res, Alpha: $(model.alpha), L0: $(model.init_level)")
    end
    return res
end

"""
compute_init_heuristic_(vector{T})

Description:
    Computes the initial level (l_0) according to the heuristic approach described below 
    init_level = intercept of the linear trendline for the first 10 values (Non-seasonal)

Parameters:
    - Vector{AbstractFloat} y: observations; inputted by user

Returns:
    - Float64: intercept of the linear trendline for the first 10 values (Non-seasonal)   
"""
function compute_init_heuristic_(y::Vector{AbstractFloat}) 
    nx = collect(1.0:min(10.0, length(y))) # Computed for time values from 1.0 to 10.0 (or smaller depending on the size of y)
    ny = [y[i] for i in 1:min(10, length(y))] # Get y values for trendline
    n = length(ny) # Compute number of data points for trendline

    slope = (n * dot(nx, ny) - sum(nx) * sum(ny)) / (n * sum(nx.^2) - sum(nx)^2) # Compute trendline slope
    intercept = (sum(ny) - slope * sum(nx)) / n # Compute trendline intercept

    return intercept
end

"""
cleanParams_()

Description:
    Throws an error if a given parameters is out of the specified bounds.

Parameters:
    - param: value of parameter to be checked
    - name: name of parameter to be checked
    - lb: lower bound of parameter, default = 0.0
    - ub: upper bound of parameter, defualt = 1.0

Note: if param is not of type Nothing, then param must be in (0.0, 1.0)

Returns:
    - nothing

"""
function cleanParams_(param::Union{Nothing, Float64}, name::String, lb::Float64 = 0.0, ub::Float64 = 1.0)
    if !isnothing(param) && (param < lb || param > ub)
        error(name * " needs to be in the range ($lb, $ub)")
    end
end

"""
fit!(model::ETSModel)

Description:
    fits an ETSModel object
    - computes initial values of parameters: 'alpha' and 'init_value'
    - modifies the 'model' directly.
    - after fit! is called, predict(model) is ready to be called

Parameters:
    - model: ExponentialSmoothing struct
    - v: verbosity, default = 0. if v > 0, verbose will be used.

Computation of model parameters:
'alpha': computed using Brent's method
'init_value': computed using a heuristic method

Returns:
    - nothing
"""
function fit!(model::ETSModel; v=0) # Set alpha + init_level
    if (isnothing(model.alpha))
        @warn "Since no value was entered for 'alpha', it will be chosen"
    end
    if (isnothing(model.init_level))
        model.init_level = compute_init_heuristic_(model.y)
        @warn "Since no value was entered for 'init_level', it will be chosen"
    end 

    # Optimize predict(model) - directly updates model
    if isnothing(model.alpha) # if atleast one value needs to be optimized
        opt = LBFGS()
        res = optimize(x -> SSE_(model, alpha=x, verbose=v), 0.0, 1.0)
    end

    return
end

"""
predict()

Description:
    - Exposed function
    - computes a forecast with Simple Exponential Smoothing (SES) [1]
    - calls predict_ a helper function, which returns a vector of fitted values and SSE error term.

Parameters
    - model: ExponentialSmoothing struct (defined above)
    - Note: when predict is called, all parameters have been calibrated.

Returns:
    A tuple containing the following:
    - Vector{Float64} forecast: contains fitted values

References: 
    [1] Hyndman, R.J., & Athanasopoulos, G. (2019) *Forecasting:
        principles and practice*, 3rd edition, OTexts: Melbourne,
        Australia. OTexts.com/fpp3. Accessed on April 19th 2020.
"""
function predict(model::ETSModel)
    return predict_(model)[1]
end