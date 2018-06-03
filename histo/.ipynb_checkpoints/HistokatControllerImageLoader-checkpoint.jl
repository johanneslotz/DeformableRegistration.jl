module HistokatControllerImageLoader
export getImageObject, getTile, removeBlack!
using PyCall
@pyimport histokat
# if HistokatController does not load:
# cp -r HistokatController /usr/local/anaconda3/lib/python3.6/site-packages/histokat/backend/lib
# cd /usr/local/anaconda3/lib/python3.6/site-packages/histokat/backend/lib
# install_name_tool -add_rpath "@loader_path/HistokatController/lib" libHistokatController.dylib

    tile_size = [1000, 1000]

    function getImageObject(filename::String)
        return histokat.controller[:open_session](filename)[:get_image]()
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
        
        return tile
    end

    function removeBlack!(arr::Array{Number2)
        for i = 1:length(arr)
            if arr[i]<5
                arr[i]=255
            end
        end
    end
end
