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

function constructTestImagesAndPatches(;patchLevel=1, imageLevel=4)
	refImg, temImg, options = constructTestImages()
	refPatches=[]
	temPatches=[]
	refPatch = createImage(refImg.data[1:120,1:120])
	temPatch = createImage(temImg.data[1:120,1:120])
	refPatch = restrictResolutionToLevel(refPatch, patchLevel)
	temPatch = restrictResolutionToLevel(temPatch, patchLevel)
	push!(refPatches, refPatch)
	push!(temPatches, temPatch)
	refPatch = createImage(refImg.data[1:120,121:240], shift=[0,120])
	temPatch = createImage(temImg.data[1:120,121:240], shift=[0,120])
	refPatch = restrictResolutionToLevel(refPatch, patchLevel)
	temPatch = restrictResolutionToLevel(temPatch, patchLevel)
	push!(refPatches, refPatch)
	push!(temPatches, temPatch)
	refImg = restrictResolutionToLevel(refImg, imageLevel)
	temImg = restrictResolutionToLevel(temImg, imageLevel)
	return refPatches, temPatches, refImg, temImg, options, patchLevel, imageLevel
end


function constructTestImages()
    dataT = zeros(120,240);
	dataT[41:100,31:90] = 1
    # data[41:100, 181:220] = 1
    dataT[21:80, 151:210] = 1
    temImg = createImage(dataT) 
    dataR = zeros(120,240);
	
	m1 = [50, 60]
	m2 = [70, 180]
	r = 20

	for i=1:size(dataR,1)
		for j=1:size(dataR,2)
			if norm([i,j] - m1) < r || norm([i,j] - m2) < r
 				dataR[i,j] = 1
			end
		end
	end
	
	# dataR[31:50,39:78] = 1
    # dataR[51:70,43:82] = 1
    # dataR[31:40,29:78] = 1
    # dataR[51:70,  163:202] = 1
    # dataR[71:90,  159:198] = 1
    # dataR[51:60,  159:210] = 1
    refImg = createImage(dataR)
    options = getSomeStandardOptions()
    return refImg, temImg, options
end

function constructTestImages2()
    dataT = zeros(120,240);
	dataT[26:65,26:65] = 1
	dataT[56:95,86:125] = 1
    # data[41:100, 181:220] = 1
    dataT[56:95, 176:215] = 1
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
	
    refImg = createImage(dataR)
    options = getSomeStandardOptions()
    return refImg, temImg, options
end

function getSomeStandardOptions()
	options = regOptions()
    options.matrixFree = true;
    options.interpolateToReferenceImage = true
    options.regularizerWeight = 100
    options.stopping["tolQ"] = 1e-4
    options.maxIterCG = 4000
	return options
end

function stripDisplacement(d::scaledArray,idx::Tuple{UnitRange{Int64},UnitRange{Int64}})
	dx = reshape(d.data[1:Int(end/2)],d.dimensions)[idx[1],idx[2]]
	dy = reshape(d.data[1+Int(end/2):end],d.dimensions)[idx[1],idx[2]]
	dims = (length(idx[1]),length(idx[2]))
	shift = [d.shift[1], d.shift[2]+d.voxelsize[2]*(idx[2][1]-1)]
	return scaledArray(vcat(dx[:], dy[:]), dims, d.voxelsize, shift)
end

function combineDisplacementsNaive(displacements::Array{scaledArray},nPatches)
    @assert nPatches[1]==1 "two-dimensional patching not supported"
	finalsizeX = displacements[1].dimensions[1]
	finalsizeY = sum([d.dimensions[2] for d in displacements])
    finalsize = (finalsizeX, finalsizeY)
    finalDisplacementX = hcat([reshape(d.data[1:Int(end/2)],d.dimensions) for d in displacements]...)
    finalDisplacementY = hcat([reshape(d.data[1+Int(end/2):end],d.dimensions) for d in displacements]...)
    finalDisplacement = scaledArray(vcat(finalDisplacementX[:], finalDisplacementY[:]), finalsize, displacements[1].voxelsize, displacements[1].shift)
    return finalDisplacement
end
