using DeformableRegistration.Visualization
using DeformableRegistration.Examples
using DeformableRegistration.ImageProcessing
using DeformableRegistration.Transformation
using DeformableRegistration.Interpolation
using DeformableRegistration.Distance
using DeformableRegistration.Types
using DeformableRegistration.Regularizer
using DeformableRegistration.Optimization: optimizeGaussNewtonAugmentedLagrangian
using DeformableRegistration:regOptions
using Test
using LinearAlgebra
using SparseArrays
using Interpolations

include("/Users/jo/uhw/experiments/experimental-julia/twoResolutionRegistration.jl")

function constructTestImages()
    dataT = zeros(120,240);
 	dataT[26:65,26:65] .= 1
 	dataT[56:95,96:125] .= 1
     # data[41:100, 181:220] = 1
    dataT[56:95, 176:215] .= 1 # 75,195
    temImg = createImage(dataT) 
    dataR = zeros(120,240);

 	m1 = [60, 60]
 	m2 = [55, 125]
 	m3 = [60, 180]

 	r = 15
 	r2 = 20

 	for i=1:size(dataR,1)
 		for j=1:size(dataR,2)
 			if norm([i,j] - m1) < r  || norm([i,j] - m2) < r2 || norm([i,j] - m3) < r
  				dataR[i,j] = 1
 			end
 		end
 	end

     for i=178:7:216
         dataT[:,i:i+4] .= 0
     end

     for i=(m3[2]-15+3):5:(m3[2]+15)
         dataR[:,i:i+2] .= 0
     end

     refImg = createImage(dataR)
     options = regOptions()
     return refImg, temImg, options
end


@testset "registerNonParametricOnGrid" begin

# "GLOBAL MULTILEVEL"
options.levels = [3,2]
options.maxIterGaussNewton = 1000
options.regularizerWeight = 10
options.interpolateToReferenceImage = true

reg = createCurvatureOperatorCentered
R, T, options = constructTestImages()
interpolationScheme = BSpline(Cubic(Line()))
#Bf=Bf'*Bf;
targetGrid =  getCellCenteredGrid(restrictResolutionToLevel(R,2));


fineDisplacement0 = registerNonParametricOnGrid(R, T, targetGrid, options,
    regularizerOperator=reg, interpolationScheme=interpolationScheme)

ssdGlobal = ssdDistance(R, T, (fineDisplacement0+getCellCenteredGrid(R)).data, doDerivative=true)[1]
#regularizerGlobal = fineDisplacement0.data[:]'*Bf*fineDisplacement0.data[:]


@info fineDisplacement0.dimensions
@info "ssd ", ssdGlobal
#@info "reg ", α() *regularizerGlobal

visualizeResults(refImgFine, temImgFine, displacement= fineDisplacement0, showDeformationImage=false, suptitle="",
    filename = "", plotInitialState=true)

end
