# SimpleExpSmoothing.jl
Implementing [Simple Exponential Smoothing](https://otexts.com/fpp2/ses.html) (SES) in Julia.

## Table of Contents
- [Installation](#installation)
- [Basic Tour](#basic-tour)
- [Customization](#customization)

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

y = [445.36, 453.2, 454.41, 422.38, 456.04, 440.39, 425.19, 486.21, 500.43, 521.28, 
     508.95, 488.89, 509.87, 456.72, 473.82, 525.95, 549.83, 542.32] # Data to forecast on. Cite: [1]
mdl = ExponentialSmoothing(y) # Initialize SES model

fit!(mdl) # Fit the model. Find optimal parameters for SES.
yhat = predict(mdl) # Compute the forecast based on determined parameters
```
At this point, `yhat` will look something like this:
```julia
446.5735923427146
445.561818502492
451.9297824011145
453.9975437047159
427.6379479108166
451.31678329980605
442.2071069188588
428.01991818111446
. . .
544.3769207888416
542.6620626999179
542.6620626999179 
543.542.6620626999179
542.6620626999179
542.6620626999179
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

<img src="https://user-images.githubusercontent.com/21336191/113236346-4bc01400-926a-11eb-8863-8024ecfed6e8.png" width="600" alt="example_plot">

## Customization

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

<img src="https://user-images.githubusercontent.com/21336191/113234520-0817db00-9267-11eb-9ebe-369ae38a3a36.png" width="600" alt="example_plot">

## References
[[1](https://otexts.com/fpp2/)] Hyndman, R.J., & Athanasopoulos, G. (2018) Forecasting: principles and practice, 2nd edition, OTexts: Melbourne, Australia. OTexts.com/fpp2. Accessed on 03/31/2021.
