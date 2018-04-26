from histokat.backend.functions import (c_image_get_extent,
                                        c_image_get_group_id, c_image_get_hash,
                                        c_image_get_id,
                                        c_image_get_num_channels,
                                        c_image_get_num_levels,
                                        c_image_get_slice_name,
                                        c_image_get_tile_extent,
                                        c_image_get_url,
                                        c_image_get_value_range,
                                        c_image_get_voxel_size,
                                        c_image_open_region)
from .utils import Box


class Image:
    def __init__(self, image_id):
        self._image_id = image_id

    def get_url(self):
        return c_image_get_url(self._image_id)

    def get_extent(self, level=0):
        return c_image_get_extent(self._image_id, level)

    def get_num_channels(self):
        return c_image_get_num_channels(self._image_id)

    def get_num_levels(self):
        return c_image_get_num_levels(self._image_id)

    def get_value_range(self):
        return c_image_get_value_range(self._image_id)

    def get_voxel_size(self, level=0):
        return c_image_get_voxel_size(self._image_id, level)

    def get_tile_extent(self):
        return c_image_get_tile_extent(self._image_id)

    def open_region(self, level, region):
        region = Box(region)
        return c_image_open_region(self._image_id, level, region)

    def get_slice_name(self, z):
        return c_image_get_slice_name(self._image_id, z)

    def get_id(self):
        return c_image_get_id(self._image_id)

    def get_group_id(self):
        return c_image_get_group_id(self._image_id)

    def get_hash(self):
        return c_image_get_hash(self._image_id)
