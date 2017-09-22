module Visualization

using PyPlot
using Images

using DeformableRegistration.Transformation
using DeformableRegistration.Interpolation

export showImage, showGrid, visualizeResults

# plotting properties
PyPlot.rc("legend",fontsize=20)
PyPlot.rc("xtick", labelsize=10, color="black", direction="out")
PyPlot.rc("ytick", labelsize=10, color="black", direction="out")

function showImage(image::ImageMeta;name="")
    showImage(image.data,spatialDomain=image["spatialdomain"],name=name)
end

function showImage(imageData::Array{Float64,2};spatialDomain=[0.0,size(imageData)[1],0.0,size(imageData)[2]],name="")
    imshow(imageData,
           interpolation="none",
           cmap="gray",
           extent=[spatialDomain[3:4],
                   spatialDomain[2],spatialDomain[1]])
    xlabel("x",fontsize=10)
    ylabel("y",fontsize=10)
    title(name,fontsize=12)
end

function showGrid(spatialDomain,gridSize,centeredGrid;
                  deformationField=zeros(2*prod(gridSize)),
                  affineParameters=[1,0,0,0,1,0.0],
                  showPoints=false, showIndices=false,
                  numberOfGridLinesX = 20, numberOfGridLinesY = 20, gridColor="red")

    # divide deformation field into x and y components
    deformationFieldX = reshape(deformationField[1:prod(gridSize)]    ,gridSize[1],gridSize[2])
    deformationFieldY = reshape(deformationField[prod(gridSize)+1:end],gridSize[1],gridSize[2])

    # small gridSize -> appropiate number of grid lines
    if(prod(gridSize) <= 64)
        numberOfGridLinesX = gridSize[2] + 1
        numberOfGridLinesY = gridSize[1] + 1
    end

    # define where grid lines should be drawn
    gridLineX = linspace(spatialDomain[3],spatialDomain[4],numberOfGridLinesX)
    gridLineY = linspace(spatialDomain[1],spatialDomain[2],numberOfGridLinesY)

    # draw each line of the grid separately
    for i=1:numberOfGridLinesX
        # define N points for a line
        N = numberOfGridLinesY*10
        line = zeros(2*N)
        line[1:N] = gridLineX[i]*ones(N)
        line[N+1:end] = linspace(spatialDomain[1],spatialDomain[2],N)
        # affine transformation
        transformedLine = transformGridAffine(line'[:],affineParameters)
        # deform line
        transformedLine[1:N] = transformedLine[1:N] + interpolateDeformationField(deformationFieldX,spatialDomain,line)
        transformedLine[N+1:end] = transformedLine[N+1:end] + interpolateDeformationField(deformationFieldY,spatialDomain,line)
        # plot line
        plot(transformedLine[1:N],transformedLine[N+1:end],gridColor)
    end
    for i=1:numberOfGridLinesY
        # define N points for a line
        N = numberOfGridLinesX*10
        line = zeros(2*N)
        line[1:N] = linspace(spatialDomain[3],spatialDomain[4],N)
        line[N+1:end] = gridLineY[i]*ones(N)
        # affine transformation
        transformedLine = transformGridAffine(line'[:],affineParameters)
        # deform line
        transformedLine[1:N] = transformedLine[1:N] + interpolateDeformationField(deformationFieldX,spatialDomain,line)
        transformedLine[N+1:end] = transformedLine[N+1:end] + interpolateDeformationField(deformationFieldY,spatialDomain,line)
        # plot line
        plot(transformedLine[1:N],transformedLine[N+1:end],gridColor)
    end

    # show grid points
    if(showPoints)
        transformedGridPoints = transformGridAffine(centeredGrid,affineParameters) + deformationField
        plot(transformedGridPoints[1:prod(gridSize)],transformedGridPoints[prod(gridSize)+1:end], "ko",markersize=2)
    end

    # show indices
    if (showIndices & showPoints)
        for i=1:prod(gridSize)
            text(transformedGridPoints[i]+0.01,transformedGridPoints[prod(gridSize)+i]-0.02,i,fontsize=8)
        end
    end

    # set axis
    axis([spatialDomain[3],spatialDomain[4],spatialDomain[2],spatialDomain[1]])

end

function visualizeResults(referenceImage,templateImage;
                          deformationField=zeros(2*prod(size(referenceImage))),
                          affineParameters=[1,0,0,0,1,0.0],
                          numberOfGridLinesX = 10, numberOfGridLinesY = 10)

    subplot(2,3,1)
    showImage(referenceImage)
    title("Reference")
    subplot(2,3,2)
    showImage(templateImage)
    title("Template")
    subplot(2,3,3)
    showImage(abs(referenceImage.data-templateImage.data),spatialDomain=referenceImage["spatialdomain"])
    title("Template-Reference")
    subplot(2,3,4)
    showImage(templateImage)
    centeredGrid = getCellCenteredGrid(referenceImage)
    showGrid(referenceImage["spatialdomain"],
             size(referenceImage),centeredGrid,
             deformationField=deformationField,
             affineParameters=affineParameters,
             numberOfGridLinesX=numberOfGridLinesX,numberOfGridLinesY=numberOfGridLinesY)
    title("y (def. Grid)")
    subplot(2,3,5)
    transformedGrid = transformGridAffine(centeredGrid,affineParameters) + deformationField
    transformedTemplate = interpolateImage(templateImage,transformedGrid)
    showImage(transformedTemplate,spatialDomain=referenceImage["spatialdomain"])
    title("Template[y]")
    subplot(2,3,6)
    differenceImage = abs(referenceImage.data-transformedTemplate)
    showImage(differenceImage,spatialDomain=referenceImage["spatialdomain"])
    xlabel("max: $(maximum(differenceImage[:]))")
    title("Template[y]-Reference")

end

end
