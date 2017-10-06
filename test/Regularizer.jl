#using DeformableRegistration.ImageProcessing
using DeformableRegistration: Regularizer, Transformation
 using Base.Test


# @testset "check if constant deformation yields zero" begin
#     omega = [0.0,1.0,0.0,1.0];
#     m = [5; 7];
#     h = [1.0,1.0];
#     u = ones(size(getStaggeredGrid(omega,m)))
#     B = createElasticOperatorStaggered(h,m)
#     S, dS, d2S = regularizer(u, B)
#     @test_approx_eq(S[1], 0)
#     @test sum(dS) == 0.0
# end
#
# @testset "check if random deformation yields nonzero result" begin
#     omega = [0.0,1.0,0.0,1.0];
#     m = [5; 7];
#     h = [1.0,1.0];
#     u = rand(size(getStaggeredGrid(omega,m)))
#     B = createElasticOperatorStaggered(h,m)
#     S, dS, d2S = regularizer(u, B)
#     @test S[1] > 0
#     @test (dS'*dS)[1] > 0
# end

@testset "Regularizer" begin
@testset "check if constant deformation yields zero" begin
    m = [5; 7];
    h = [1.0,1.0];
    u = ones(2*m[1],m[2])[:]
    B = createElasticOperatorCentered(h,m)
    S, dS, d2S = regularizer(u, B)
    @test  S[1] ≈ 0
    @test sum(dS) == 0.0
end

@testset "check if random deformation yields nonzero result" begin
    omega = [0.0,1.0,0.0,1.0];
    m = [5; 7];
    h = [1.0,1.0];
    u = rand(2*m[1],m[2])[:]
    B = createElasticOperatorCentered(h,m)
    S, dS, d2S = regularizer(u, B)
    @test S[1] > 0
    @test (dS'*dS)[1] > 0
end

@testset "check if constant deformation yields zero" begin
    omega = [0.0,1.0,0.0,1.0];
    m = [5; 7];
    h = [1.0,1.0];
    u = ones(2*m[1],m[2])[:]
    B = createDiffusiveOperatorCentered(h,m)
    S, dS, d2S = regularizer(u, B)
    @test S[1] ≈ 0
    @test sum(dS) == 0.0
end

@testset "check if constant deformation yields zero" begin
    omega = [0.0,1.0,0.0,1.0];
    m = [5; 7];
    h = [1.0,1.0];
    u = ones(m[1],2*m[2])
    #u[3,3] = 1
    u = u[:]
    B = createCurvatureOperatorCentered(h,m)
    S, dS, d2S = regularizer(u, B)
    @test S[1] ≈ 0
    @test sum(dS) == 0.0

    #using PyPlot
    #imshow(reshape(dS[1:35],(5,7)))
    #plot(dS[1:25])
end

@testset "check if random deformation yields nonzero result" begin
    omega = [0.0,1.0,0.0,1.0];
    m = [5; 7];
    h = [1.0,1.0];
    u = rand(2*m[1],m[2])[:]
    B = createCurvatureOperatorCentered(h,m)
    S, dS, d2S = regularizer(u, B)
    @test S[1] > 0
    @test (dS'*dS)[1] > 0
end

@testset "check if random deformation yields nonzero result" begin
    omega = [0.0,1.0,0.0,1.0];
    m = [5; 7];
    h = [1.0,1.0];
    u = rand(2*m[1],m[2])[:]
    B = createDiffusiveOperatorCentered(h,m)
    S, dS, d2S = regularizer(u, B)
    @test S[1] > 0
    @test (dS'*dS)[1] > 0
end
end #testset regularizer
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
