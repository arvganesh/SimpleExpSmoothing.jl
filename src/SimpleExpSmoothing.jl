module SimpleExpSmoothing

# Export from ets.jl
export ExponentialSmoothing, fit!, predict, plot_ts

# Include Files
include("ets.jl")
include("plot_ts.jl")

end # module
