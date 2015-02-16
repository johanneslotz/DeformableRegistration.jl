#using ImageRegistration.ImageProcessing
using ImageRegistration.Regularizer
using ImageRegistration.Transformation
using Base.Test

# setup logging
Logging.configure(level = Logging.INFO) # set to DEBUG to see error tables
Logging.info("Regularizer: Testing regularizers for basic plausibility...")

## check if constant deformation yields zero
omega = [0.0,1.0,0.0,1.0];
m = [5; 7];
h = [1.0,1.0];
u = ones(size(getStaggeredGrid(omega,m)))
B = createElasticOperatorStaggered(h,m)
S, dS, d2S = regularizer(u, B)
@test_approx_eq(S[1], 0)
@test sum(dS) == 0.0
Logging.info("Regularizer: Elastic Staggered: Constant deformation yields zero ✔")

## check if random deformation yields nonzero result
omega = [0.0,1.0,0.0,1.0];
m = [5; 7];
h = [1.0,1.0];
u = rand(size(getStaggeredGrid(omega,m)))
B = createElasticOperatorStaggered(h,m)
S, dS, d2S = regularizer(u, B)
@test S[1] > 0
@test (dS'*dS)[1] > 0
Logging.info("Regularizer: Elastic Staggered: Random deformation yields nonzero ✔")

## check if constant deformation yields zero
omega = [0.0,1.0,0.0,1.0];
m = [5; 7];
h = [1.0,1.0];
u = ones(size(getCellCenteredGrid(omega,m)))
B = createElasticOperatorCentered(h,m)
S, dS, d2S = regularizer(u, B)
@test_approx_eq(S[1], 0)
@test sum(dS) == 0.0
Logging.info("Regularizer: Elastic Centered: Constant deformation yields zero ✔")

## check if random deformation yields nonzero result
omega = [0.0,1.0,0.0,1.0];
m = [5; 7];
h = [1.0,1.0];
u = rand(size(getCellCenteredGrid(omega,m)))
B = createElasticOperatorCentered(h,m)
S, dS, d2S = regularizer(u, B)
@test S[1] > 0
@test (dS'*dS)[1] > 0
Logging.info("Regularizer: Elastic Centered: Random deformation yields nonzero ✔")

## check if constant deformation yields zero
omega = [0.0,1.0,0.0,1.0];
m = [5; 7];
h = [1.0,1.0];
u = ones(size(getCellCenteredGrid(omega,m)))
B = createDiffusiveOperatorCentered(h,m)
S, dS, d2S = regularizer(u, B)
@test_approx_eq(S[1], 0)
@test sum(dS) == 0.0
Logging.info("Regularizer: Diffusive Centered: Constant deformation yields zero ✔")

## check if random deformation yields nonzero result
omega = [0.0,1.0,0.0,1.0];
m = [5; 7];
h = [1.0,1.0];
u = rand(size(getCellCenteredGrid(omega,m)))
B = createDiffusiveOperatorCentered(h,m)
S, dS, d2S = regularizer(u, B)
@test S[1] > 0
@test (dS'*dS)[1] > 0
Logging.info("Regularizer: Diffusive Centered: Random deformation yields nonzero ✔")

# speed test
#omega = [0.0,1.0,0.0,1.0];
#m = [512; 512];
#h = [1.0,1.0];
#u = rand(size(getStaggeredGrid(omega,m)))
#@time B = createElasticOperatorStaggered(h,m)
#@time S, dS, d2S = regularizer(u, B)

# visualization
#using PyPlot
#spy(B)
