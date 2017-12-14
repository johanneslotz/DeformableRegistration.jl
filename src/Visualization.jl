module Visualization

using PyPlot
using Images

using DeformableRegistration: Transformation, Interpolation, ImageProcessing

export showImage, plotGrid, visualizeResults, clf, figure, suptitle

# plotting properties
PyPlot.rc("legend",fontsize=20)
PyPlot.rc("xtick", labelsize=10, color="black", direction="out")
PyPlot.rc("ytick", labelsize=10, color="black", direction="out")

function showImage(image::regImage;name="")
    extent = [image.shift[2],
    image.shift[2] + image.voxelsize[2]*size(image.data,2),
    image.shift[1],
    image.shift[1] + image.voxelsize[1]*size(image.data,1)]
    imshow(image.data,
           interpolation="none",
           cmap="gray",
           extent=extent)
    xlabel("x",fontsize=10)
    ylabel("y",fontsize=10)
    title(name,fontsize=12)
end

function plotGrid(transformedGrid::Array{Float64,1}, shape::Tuple{Int64,Int64}; numberOfGridLines=20)
    tX = reshape(transformedGrid[1:Int(length(transformedGrid)/2)],shape)
    tY = reshape(transformedGrid[Int(length(transformedGrid)/2+1):end],shape)

    for i = convert(Array{Int}, round.(linspace(1,size(tX,2), numberOfGridLines)))
        plot(tY[:,i], tX[:,i], "k", linewidth=1.0)
    end
    for i = convert(Array{Int}, round.(linspace(1,size(tX,1), numberOfGridLines)))
        plot(tY[i,:], tX[i,:],"k", linewidth=1.0)
    end
    axis("equal")
    return tX, tY
end

function plotGrid(g::scaledArray; numberOfGridLines=20)
    return plotGrid(g.data, g.dimensions, numberOfGridLines=numberOfGridLines)
end

function setFigSize()
    f = PyPlot.gcf()
    f["set_size_inches"](12,7)
    f["set_dpi"](150)
end

function visualizeResults(referenceImage,templateImage;
                          displacement=scaledArray(
                            zeros(2*prod(size(referenceImage.data))),
                            size(referenceImage.data),
                            referenceImage.voxelsize,
                            referenceImage.shift),
                          affineParameters=[1,0,0,0,1,0.0],
                          numberOfGridLines = 20,
                          showDeformationImage = false,
                          suptitle = "", cmap = "gray", filename = "")


    subplot(2,3,1)
    showImage(referenceImage)
    title("Reference")

    subplot(2,3,2)
    showImage(templateImage)
    title("Template")

    subplot(2,3,3)
    showImage(createImage(abs.(Array(referenceImage.data)-Array(templateImage.data))))
    title("Template-Reference")

    subplot(2,3,4)
    centeredGrid = getCellCenteredGrid(referenceImage)
    transformedGrid = transformGridAffine(centeredGrid,affineParameters) + displacement
    plotGrid(transformedGrid, numberOfGridLines=numberOfGridLines)
    title("y (def. Grid)")

    subplot(2,3,5)
    transformedTemplate = interpolateImage(templateImage,transformedGrid)
    showImage(createImage(transformedTemplate))
    title("Template[y]")

    subplot(2,3,6)
    differenceImage = abs.(Array(referenceImage.data)-Array(transformedTemplate))
    showImage(createImage(differenceImage))
    xlabel("max: $(maximum(differenceImage[:]))")
    title("Template[y]-Reference")
    PyPlot.suptitle(suptitle)
    setFigSize()
    if filename != ""
        savefig(filename[1:end-4]*"-A-"*filename[end-4:end])
    end

    if showDeformationImage && (displacement.data != zeros(2*prod(size(referenceImage.data))))
        shape = displacement.dimensions
        tX = reshape(displacement.data[1:Int(length(displacement.data)/2)],shape)
        tY = reshape(displacement.data[Int(length(displacement.data)/2+1):end],shape)

        figure()
        subplot(1,2,1)
        imshow(tX, cmap=cmap)
        colorbar()
        subplot(1,2,2)
        imshow(tY, cmap=cmap)
        colorbar()
        PyPlot.suptitle(suptitle)
        setFigSize()
        if filename != ""
            savefig(filename[1:end-4]*"-B-"*filename[end-4:end])
        end
    end

end

end
