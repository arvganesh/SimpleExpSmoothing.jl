using Plots
"""
plot_ts()

Description:
    Plots a time series given observations and fitted values

Parameters:
    y: observed values
    yhat: fitted values, from SES prediction

or

Parameters:
    model: given an ETSModel object, it computes fitted values, and plots observed and fitted values

Returns:
    Single plot of the observed and predicted values.

"""
function plot_ts(y::Vector{T} where {T <: AbstractFloat}, yhat::Vector{T} where {T <: AbstractFloat})
    data_size = size(yhat)
    p = plot(y, title = "Time Series Example", lw = 2, legend=:bottomright, label = "Observed")
    plot!(p, yhat, lw = 2, label = "Forecasted")
    xlabel!("Unit Time")
end

function plot_ts(model::ETSModel) # Takes in ETSModel object and plots the time series prediction
    fit!(model)
    yhat = predict(model)
    y = model.y
    plot_ts(y, yhat)
end