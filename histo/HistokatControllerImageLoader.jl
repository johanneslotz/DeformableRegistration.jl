module HistokatControllerImageLoader
export getImageObject, getTile
using PyCall

using DeformableRegistration.ImageProcessing: createImage
@pyimport histokat
# if HistokatController does not load:
# cp -r HistokatController /usr/local/anaconda3/lib/python3.6/site-packages/histokat/backend/lib
# cd /usr/local/anaconda3/lib/python3.6/site-packages/histokat/backend/lib
# install_name_tool -add_rpath "@loader_path/HistokatController/lib" libHistokatController.dylib

    tile_size = [1000, 1000]

    function getImageObject(filename::String)
        return histokat.controller[:open_session](filename)[:get_image]()
    end


    function getVoxelSize(image; level = 0)
        x = image[:get_voxel_size](level)[:x]/1000000
        y = image[:get_voxel_size](level)[:y]/1000000
        return (x,y)
    end
    function getExtent(image; level = 0)
        x = image[:get_extent](level)[:x]/1000000
        y = image[:get_extent](level)[:y]/1000000
        return (x,y)
    end

    function getTile(image, level, x, y, z)
        # get correct level and position
        max_level = image[:get_num_levels]() - 1
        level = max_level - level
        x = tile_size[1] * x
        y = tile_size[2] * y
        region = [[x, y, z], [x + tile_size[2] - 1, y + tile_size[2] - 1, z]]
        # load tile as numpy array
        tile = image[:open_region](level, region)
        # normalize and convert to uint8 image
        tile = tile .* (255.0/image[:get_value_range]()[:max])
        #tile = Image.fromarray(tile.astype(np.uint8))
        tile = sum(Array(tile)[1,:,:,:],3)[:,:,1]/3
        vx = getVoxelSize(image,level=level)
        tileIm = createImage(tile, voxelsize = [vx[1],vx[2]], shift = [tile_size[1] * x * vx[1], tile_size[2] * y * vx[2]])
        return tileIm
    end


end
