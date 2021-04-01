# SimpleExpSmoothing.jl
Implementing Simple Exponential Smoothing (SES) in Julia.

| Linux/OSX | Windows | Coverage |
| :----: | :----: | :----: |
| [![Build Status](https://travis-ci.com/arvganesh/SimpleExpSmoothing.jl.svg?branch=master)](https://travis-ci.com/arvganesh/SimpleExpSmoothing.jl) | [![Build status](https://ci.appveyor.com/api/projects/status/lmbbqp2tf46ccvyd?svg=true)](https://ci.appveyor.com/project/arvganesh/simpleexpsmoothing-jl) | [![codecov](https://codecov.io/gh/arvganesh/SimpleExpSmoothing.jl/branch/master/graph/badge.svg?token=V7ZS8LCMKU)](https://codecov.io/gh/arvganesh/SimpleExpSmoothing.jl) |

## Installation
First, `git clone` into the package directory.

Then, in Julia, do:
```[activate .``` and ```instantiate```

## Basic Tour
```julia
using SimpleExpSmoothing

y = [1.0, 2.0, 3.0, 4.0, 5.0] 
mdl = ExponentialSmoothing(y) 

fit!(mdl) # Fit the model. Find optimal parameters for SES.
yhat = predict(mdl) # Compute the forecast based on determined parameters
```
At this point, `yhat` will look something like this:
```julia
0.0
0.9999999776901562
1.9999999776901558
2.999999977690156
3.999999977690156
4.999999977690155
4.999999977690155
4.999999977690155
4.999999977690155
4.999999977690155
```
To visualize the forecast, do this:
```julia
plot_ts(y, yhat) # Plots y and yhat against time.
```
Alternatively, to save some space, do this:
```julia
plot_ts(mdl) # Takes the ExponentialSmoothing object, makes predictions, and plots them.
```
Both methods of plotting will yield something like this:

<img src="https://user-images.githubusercontent.com/21336191/113225417-32f83400-9253-11eb-94e0-a54e5fb334b4.png" width="600" alt="example_plot">

## Customizing the Model

While the example above showed the `ExponentialSmoothing` object being used with only 1 parameter, it is possible to use more: 
```julia
# Inputted values will be used if they are specified. Otherwise, they will be computed.
mdl = ExponentialSmoothing(observations, h = 3) 
mdl = ExponentialSmoothing(observations, alpha = 0.4)
mdl = ExponentialSmoothing(observations, alpha = 0.25, init_level = 500.0)
mdl = ExponentialSmoothing(observations, h = 15, alpha = 0.3, init_level = 750.0)

Key:
h: Number of time steps to forecast (h > 0, default: 5)
alpha: Smoothing parameter used in SES. (0 < alpha < 1, default: nothing)
init_level: initial level value. (default: nothing)
```
If `alpha` or `init_level` is `nothing`, an estimated optimal value will be computed.
### Example
```julia
y = rand(1.:100., 25)
mdl = ExponentialSmoothing(y, h = 3, init_level = 2) # Since alpha is not specified, it will be computed.

plot_ts(mdl)
```
The output of the code above may look like this:

<img src="https://user-images.githubusercontent.com/21336191/113226601-35a85880-9256-11eb-9009-46977e0c6245.png" width="600" alt="example_plot">

