function constructPatches(nPatches::Tuple{Int64, Int64}, refImg::regImage, temImg::regImage, initialDisplacement::scaledArray)
	@assert nPatches[1] == 1 "patches only supported in one (the second) dimension"
	sizePxRef = Int(size(refImg.data,2) / nPatches[2])
	sizePxDispl = Int(initialDisplacement.dimensions[2] / nPatches[2])
	iDdataX = reshape(initialDisplacement.data[1:Int(end/2)],initialDisplacement.dimensions)
	iDdataY = reshape(initialDisplacement.data[1+Int(end/2):end],initialDisplacement.dimensions)
	refImgs = []
	temImgs = []
	iDs = []
	for i = 1:nPatches[2]
		newShift = refImg.shift + [0, Int((i-1) * sizePxRef * refImg.voxelsize[2])]
		newRefImg = regImage(refImg.data[:, (i-1) * sizePxRef+1 : i * sizePxRef], refImg.voxelsize, newShift)
		newTemImg = regImage(temImg.data[:, (i-1) * sizePxRef+1 : i * sizePxRef], temImg.voxelsize, newShift)

        newShiftD = initialDisplacement.shift + [0,  Int(round((i-1) * sizePxDispl * initialDisplacement.voxelsize[2]))]
		newDisplacement = scaledArray(vcat(iDdataX[:, (i-1) * sizePxDispl+1 : i * sizePxDispl][:], iDdataY[:, (i-1) * sizePxDispl+1 : i * sizePxDispl][:]),
		 								(initialDisplacement.dimensions[1],sizePxDispl), initialDisplacement.voxelsize, newShiftD)


		refImgs = vcat(refImgs, newRefImg)
		temImgs = vcat(temImgs, newTemImg)
		iDs = vcat(iDs, newDisplacement)
	end

	return(refImgs, temImgs, iDs)
end


function constructTestImages()
    data = zeros(120,240);
	data[41:80,21:60] = 1
    data[41:80, 181:220] = 1
    temImg = createImage(data) 
    dataR = zeros(120,240);
    dataR[41:80,41:80] = 1
    dataR[41:80,  161:200] = 1
    refImg = createImage(dataR)
    options = regOptions()
    options.matrixFree = true;
    options.interpolateToReferenceImage = true
    options.regularizerWeight = 100
    options.stopping["tolQ"] = 1e-4
    options.maxIterCG = 4000
    return refImg, temImg, options
end

function combineDisplacementsNaive(displacements::Array{scaledArray},nPatches)
    @assert nPatches[1]==1 "two-dimensional patching not supported"
    finalsize = Tuple([nPatches[i]*displacements[1].dimensions[i] for i=1:2])
    finalDisplacementX = hcat([reshape(d.data[1:Int(end/2)],d.dimensions) for d in displacements]...)
    finalDisplacementY = hcat([reshape(d.data[1+Int(end/2):end],d.dimensions) for d in displacements]...)
    finalDisplacement = scaledArray(vcat(finalDisplacementX[:], finalDisplacementY[:]), finalsize, displacements[1].voxelsize, displacements[1].shift)
    return finalDisplacement
end