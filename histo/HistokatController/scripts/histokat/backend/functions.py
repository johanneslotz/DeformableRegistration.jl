import ctypes
import shutil
import tempfile

import numpy as np

from histokat.backend import (_HISTOKAT_LIB, HistokatException,
                              InvalidIdException, NotInitializedException,
                              _cfunc_cancellation_callback,
                              _workspace_dir_generated)
from histokat.backend.py_c_converters import (convert_c_ptr_to_py_analysis_id,
                                              convert_c_ptr_to_py_class_info_list,
                                              convert_c_ptr_to_py_image_group,
                                              convert_c_ptr_to_py_value_range,
                                              convert_c_to_py_analysis_id_list_ptr,
                                              convert_c_to_py_vector_ptr,
                                              convert_py_to_c_analysis_id,
                                              convert_py_to_c_analysis_mode,
                                              convert_py_to_c_box,
                                              convert_py_to_c_vector)
from histokat.backend.types import (RETURN_CODE, c_analysis_id,
                                    c_analysis_id_list, c_class_info_list,
                                    c_image_group, c_string_list,
                                    c_value_range, c_vector)


def c_get_last_error():
    s = ctypes.c_char_p()
    _HISTOKAT_LIB.get_last_error(ctypes.byref(s))
    err = s.value.decode('utf-8')
    check(_HISTOKAT_LIB.destroy_string(ctypes.byref(s)))
    return err


def check(return_code):
    if return_code == RETURN_CODE.C_OK:
        return
    if return_code == RETURN_CODE.C_CONTROLLER_NOT_INITIALIZED:
        raise NotInitializedException()
    if return_code == RETURN_CODE.C_NONEXISTING_ID:
        raise InvalidIdException()
    if return_code == RETURN_CODE.C_CANNOT_CREATE_SESSION:
        raise HistokatException(
            'Cannot create Session: %s' % c_get_last_error())
    if return_code == RETURN_CODE.C_ERROR:
        raise HistokatException(c_get_last_error())


def c_get_controller_workspace_dir():
    char_ptr = ctypes.c_char_p()
    check(_HISTOKAT_LIB.get_controller_workspace_dir(ctypes.byref(char_ptr)))
    wrk_dir = char_ptr.value.decode('utf-8')
    check(_HISTOKAT_LIB.destroy_string(ctypes.byref(char_ptr)))
    return wrk_dir


def c_init_controller(config_dir='', workspace_dir='', plugin_path=''):
    old_workspace_dir_generated = _workspace_dir_generated.dir
    _workspace_dir_generated.update(None)

    # create temp folder if no workspace dir given
    if workspace_dir == '':
        workspace_dir = tempfile.mkdtemp()
        _workspace_dir_generated.update(workspace_dir)

    # perform initialization
    check(_HISTOKAT_LIB.init_controller(config_dir.encode('utf-8'),
                                        workspace_dir.encode('utf-8'), plugin_path.encode('utf-8'), _cfunc_cancellation_callback))

    # remove old workspace if it was generated
    if old_workspace_dir_generated is not None:
        shutil.rmtree(old_workspace_dir_generated)


def c_read_image_group(url):
    c_image_grp = ctypes.POINTER(c_image_group)()
    check(_HISTOKAT_LIB.read_image_group(
        url.encode('utf-8'), ctypes.byref(c_image_grp)))
    image_group = convert_c_ptr_to_py_image_group(c_image_grp)
    check(_HISTOKAT_LIB.destroy_image_group(ctypes.byref(c_image_grp)))
    return image_group


def c_get_analysis_ids():
    c_analysis_id_list_ptr = ctypes.POINTER(c_analysis_id_list)()
    check(_HISTOKAT_LIB.get_analysis_ids(ctypes.byref(c_analysis_id_list_ptr)))
    analysis_ids = convert_c_to_py_analysis_id_list_ptr(c_analysis_id_list_ptr)
    check(_HISTOKAT_LIB.destroy_analysis_id_list(
        ctypes.byref(c_analysis_id_list_ptr)))
    return analysis_ids


def c_is_image_openable(url):
    is_openable = ctypes.c_bool()
    check(_HISTOKAT_LIB.is_image_openable(
        url.encode('utf-8'), ctypes.byref(is_openable)))
    return bool(is_openable.value)


def c_open_session(url):
    sess_id = ctypes.c_int64()
    check(_HISTOKAT_LIB.open_session(url.encode('utf-8'), ctypes.byref(sess_id)))
    session_id = int(sess_id.value)
    if session_id == 0:
        raise HistokatException('Could not open session')
    return session_id


def c_session_exists(session_id):
    exists = ctypes.c_bool()
    check(_HISTOKAT_LIB.session_exists(
        ctypes.c_int64(session_id), ctypes.byref(exists)))
    return bool(exists.value)


def c_close_session(session_id):
    check(_HISTOKAT_LIB.close_session(ctypes.c_int64(session_id)))


def c_is_analyzable(session_id, analysis_id):
    analyzable = ctypes.c_bool()
    check(_HISTOKAT_LIB.is_analyzable(ctypes.c_int64(session_id), ctypes.byref(
        convert_py_to_c_analysis_id(analysis_id)), ctypes.byref(analyzable)))
    return bool(analyzable.value)


def c_analysis_was_performed_before(session_id, analysis_id):
    was_analyzed = ctypes.c_bool()
    check(_HISTOKAT_LIB.analysis_was_performed_before(ctypes.c_int64(session_id),
                                                      ctypes.byref(convert_py_to_c_analysis_id(analysis_id)), ctypes.byref(was_analyzed)))
    return bool(was_analyzed.value)


def c_get_default_analysis(session_id):
    c_analysis_id_ptr = ctypes.POINTER(c_analysis_id)()
    check(_HISTOKAT_LIB.get_default_analysis(
        ctypes.c_int64(session_id), ctypes.byref(c_analysis_id_ptr)))
    analysis_id = convert_c_ptr_to_py_analysis_id(c_analysis_id_ptr)
    check(_HISTOKAT_LIB.destroy_analysis_id(ctypes.byref(c_analysis_id_ptr)))
    return analysis_id


def c_get_image(session_id, viewable_representation):
    img_id = ctypes.c_int64()
    check(_HISTOKAT_LIB.get_image(ctypes.c_int64(session_id),
                                  ctypes.c_bool(viewable_representation), ctypes.byref(img_id)))
    image_id = int(img_id.value)
    if image_id == 0:
        raise HistokatException('Could not get image')
    return image_id


def c_image_exists(image_id):
    exists = ctypes.c_bool()
    check(_HISTOKAT_LIB.image_exists(
        ctypes.c_int64(image_id), ctypes.byref(exists)))
    return bool(exists.value)


def c_save_image_at_level(session_id, path, z_pos, level):
    check(_HISTOKAT_LIB.save_image_at_level(ctypes.c_int64(session_id),
                                            path.encode('utf-8'), ctypes.c_int64(z_pos), ctypes.c_int64(level)))


def c_save_visualization_at_level(session_id, path, z_pos, level, analysis_mode):
    c_mode = convert_py_to_c_analysis_mode(analysis_mode)
    check(_HISTOKAT_LIB.save_visualization_at_level(
        ctypes.c_int64(session_id),
        path.encode('utf-8'),
        ctypes.c_int64(z_pos),
        ctypes.c_int64(level),
        ctypes.byref(c_mode)))


def c_create_analysis(session_id, analysis_id):
    c_id = convert_py_to_c_analysis_id(analysis_id)
    check(_HISTOKAT_LIB.create_analysis(
        ctypes.c_int64(session_id), ctypes.byref(c_id)))


def c_close_analysis(session_id):
    check(_HISTOKAT_LIB.close_analysis(ctypes.c_int64(session_id)))


def c_is_analysis_open(session_id):
    is_open = ctypes.c_bool()
    check(_HISTOKAT_LIB.is_analysis_open(
        ctypes.c_int64(session_id), ctypes.byref(is_open)))
    return bool(is_open.value)


def c_get_visualization(session_id, analysis_mode):
    c_analysis_mode = convert_py_to_c_analysis_mode(analysis_mode)
    c_image_id = ctypes.c_int64()
    check(_HISTOKAT_LIB.get_visualization(
        ctypes.c_int64(session_id),
        ctypes.byref(c_analysis_mode),
        ctypes.byref(c_image_id)))
    image_id = int(c_image_id.value)
    return image_id


# Image #########################


def c_image_get_url(image_id):
    char_ptr = ctypes.c_char_p()
    check(_HISTOKAT_LIB.image_get_url(
        ctypes.c_int64(image_id), ctypes.byref(char_ptr)))
    url = char_ptr.value.decode('utf-8')
    check(_HISTOKAT_LIB.destroy_string(ctypes.byref(char_ptr)))
    return url


def c_image_get_extent(image_id, level):
    c_vector_ptr = ctypes.POINTER(c_vector)()
    check(_HISTOKAT_LIB.image_get_extent(ctypes.c_int64(image_id),
                                         ctypes.c_int64(level), ctypes.byref(c_vector_ptr)))
    extent = convert_c_to_py_vector_ptr(c_vector_ptr)
    check(_HISTOKAT_LIB.destroy_vector(ctypes.byref(c_vector_ptr)))
    return extent


def c_image_get_num_channels(image_id):
    int_ptr = ctypes.c_int64()
    check(_HISTOKAT_LIB.image_get_num_channels(
        ctypes.c_int64(image_id), ctypes.byref(int_ptr)))
    return int(int_ptr.value)


def c_image_get_num_levels(image_id):
    int_ptr = ctypes.c_int64()
    check(_HISTOKAT_LIB.image_get_num_levels(
        ctypes.c_int64(image_id), ctypes.byref(int_ptr)))
    return int(int_ptr.value)


def c_image_get_value_range(image_id):
    c_value_range_ptr = ctypes.POINTER(c_value_range)()
    check(_HISTOKAT_LIB.image_get_value_range(
        ctypes.c_int64(image_id), ctypes.byref(c_value_range_ptr)))
    value_range = convert_c_ptr_to_py_value_range(c_value_range_ptr)
    check(_HISTOKAT_LIB.destroy_value_range(ctypes.byref(c_value_range_ptr)))
    return value_range


def c_image_get_voxel_size(image_id, level):
    c_vector_ptr = ctypes.POINTER(c_vector)()
    check(_HISTOKAT_LIB.image_get_voxel_size(ctypes.c_int64(
        image_id), ctypes.c_int64(level), ctypes.byref(c_vector_ptr)))
    voxel_size = convert_c_to_py_vector_ptr(c_vector_ptr)
    check(_HISTOKAT_LIB.destroy_vector(ctypes.byref(c_vector_ptr)))
    return voxel_size


def c_image_get_tile_extent(image_id):
    c_vector_ptr = ctypes.POINTER(c_vector)()
    check(_HISTOKAT_LIB.image_get_tile_extent(
        ctypes.c_int64(image_id), ctypes.byref(c_vector_ptr)))
    tile_extent = convert_c_to_py_vector_ptr(c_vector_ptr)
    check(_HISTOKAT_LIB.destroy_vector(ctypes.byref(c_vector_ptr)))
    return tile_extent


def c_image_open_region(image_id, level, region):
    extent = region.get_extent()
    num_channels = c_image_get_num_channels(image_id)
    c_region = convert_py_to_c_box(region)
    buf = np.ndarray((num_channels, extent.z, extent.y, extent.x))
    check(_HISTOKAT_LIB.image_open_region(
        image_id, ctypes.c_int64(level), ctypes.byref(c_region), buf.ctypes))
    # move c-dimension to end
    return np.rollaxis(buf, 0, 4)


def c_image_get_slice_name(image_id, z):
    char_ptr = ctypes.c_char_p()
    check(_HISTOKAT_LIB.image_get_slice_name(ctypes.c_int64(
        image_id), ctypes.c_int64(z), ctypes.byref(char_ptr)))
    slice_name = char_ptr.value.decode('utf-8')
    check(_HISTOKAT_LIB.destroy_string(ctypes.byref(char_ptr)))
    return slice_name


def c_image_get_id(image_id):
    char_ptr = ctypes.c_char_p()
    check(_HISTOKAT_LIB.image_get_id(
        ctypes.c_int64(image_id), ctypes.byref(char_ptr)))
    img_id = char_ptr.value.decode('utf-8')
    check(_HISTOKAT_LIB.destroy_string(ctypes.byref(char_ptr)))
    return img_id


def c_image_get_group_id(image_id):
    char_ptr = ctypes.c_char_p()
    check(_HISTOKAT_LIB.image_get_group_id(
        ctypes.c_int64(image_id), ctypes.byref(char_ptr)))
    grp_id = char_ptr.value.decode('utf-8')
    check(_HISTOKAT_LIB.destroy_string(ctypes.byref(char_ptr)))
    return grp_id


def c_image_get_hash(image_id):
    int_ptr = ctypes.c_int64()
    check(_HISTOKAT_LIB.image_get_hash(
        ctypes.c_int64(image_id), ctypes.byref(int_ptr)))
    return int(int_ptr.value)

# Analysis #######################


def c_analysis_get_config_dir(session_id):
    char_ptr = ctypes.c_char_p()
    check(_HISTOKAT_LIB.analysis_get_config_dir(
        ctypes.c_int64(session_id), ctypes.byref(char_ptr)))
    config_dir = char_ptr.value.decode('utf-8')
    check(_HISTOKAT_LIB.destroy_string(ctypes.byref(char_ptr)))
    return config_dir


def c_analysis_get_workspace_dir(session_id):
    char_ptr = ctypes.c_char_p()
    check(_HISTOKAT_LIB.analysis_get_workspace_dir(
        ctypes.c_int64(session_id), ctypes.byref(char_ptr)))
    wrk_dir = char_ptr.value.decode('utf-8')
    check(_HISTOKAT_LIB.destroy_string(ctypes.byref(char_ptr)))
    return wrk_dir


def c_analysis_get_id(session_id):
    c_analysis_id_ptr = ctypes.POINTER(c_analysis_id)()
    check(_HISTOKAT_LIB.analysis_get_id(ctypes.c_int64(
        session_id), ctypes.byref(c_analysis_id_ptr)))
    analysis_id = convert_c_ptr_to_py_analysis_id(c_analysis_id_ptr)
    check(_HISTOKAT_LIB.destroy_analysis_id(ctypes.byref(c_analysis_id_ptr)))
    return analysis_id


def c_analysis_get_minor_modes(session_id, major_mode):
    string_list_ptr = ctypes.POINTER(c_string_list)()
    check(_HISTOKAT_LIB.analysis_get_minor_modes(
        session_id, ctypes.c_int64(major_mode), ctypes.byref(string_list_ptr)))
    modes = []
    for i in range(string_list_ptr.contents.num_strings):
        modes.append(string_list_ptr.contents.string_list[i].decode('utf-8'))
    check(_HISTOKAT_LIB.destroy_string_list(ctypes.byref(string_list_ptr)))
    return modes


def c_analysis_get_class_infos(session_id, input_mode):
    c_class_info_list_ptr = ctypes.POINTER(c_class_info_list)()
    check(_HISTOKAT_LIB.analysis_get_class_infos(
        session_id, ctypes.c_int64(input_mode), ctypes.byref(c_class_info_list_ptr)))
    classes = convert_c_ptr_to_py_class_info_list(c_class_info_list_ptr)
    check(_HISTOKAT_LIB.destroy_class_info_list(
        ctypes.byref(c_class_info_list_ptr)))
    return classes


def c_analysis_perform(session_id, analysis_mode):
    c_mode = convert_py_to_c_analysis_mode(analysis_mode)
    check(_HISTOKAT_LIB.analysis_perform(
        ctypes.c_int64(session_id), ctypes.byref(c_mode)))


def c_analysis_inspect_object(session_id, inspect_mode, position, level):
    char_ptr = ctypes.c_char_p()
    c_pos = convert_py_to_c_vector(position)
    check(_HISTOKAT_LIB.analysis_inspect_object(ctypes.c_int64(session_id), ctypes.c_int64(
        inspect_mode), ctypes.byref(c_pos), ctypes.c_int64(level), ctypes.byref(char_ptr)))
    inspection_value = char_ptr.value.decode('utf-8')
    check(_HISTOKAT_LIB.destroy_string(ctypes.byref(char_ptr)))
    return inspection_value
