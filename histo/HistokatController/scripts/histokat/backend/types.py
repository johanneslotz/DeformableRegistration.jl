import ctypes


class RETURN_CODE:
    C_OK = 0
    C_CONTROLLER_NOT_INITIALIZED = 1
    C_NONEXISTING_ID = 2
    C_CANNOT_CREATE_SESSION = 3
    C_ERROR = 4


class c_image_info(ctypes.Structure):
    _fields_ = \
        [
            ('image_id', ctypes.c_char_p),
            ('image_url', ctypes.c_char_p)
        ]


class c_image_group(ctypes.Structure):
    _fields_ = \
        [
            ('group_id', ctypes.c_char_p),
            ('group_url', ctypes.c_char_p),
            ('num_image_infos', ctypes.c_int64),
            ('image_info_list', ctypes.POINTER(c_image_info))
        ]


class c_analysis_id(ctypes.Structure):
    _fields_ = \
        [
            ('name', ctypes.c_char_p),
            ('version', ctypes.c_int64)
        ]


class c_analysis_id_list(ctypes.Structure):
    _fields_ = \
        [
            ('num_ids', ctypes.c_int64),
            ('analysis_ids', ctypes.POINTER(c_analysis_id))
        ]


class c_analysis_mode(ctypes.Structure):
    _fields_ = \
        [
            ('major', ctypes.c_int64),
            ('minor', ctypes.c_int64)
        ]


class c_vector(ctypes.Structure):
    _fields_ = \
        [
            ('x', ctypes.c_int64),
            ('y', ctypes.c_int64),
            ('z', ctypes.c_int64)
        ]


class c_box(ctypes.Structure):
    _fields_ = \
        [
            ('v1', c_vector),
            ('v2', c_vector),
        ]


class c_value_range(ctypes.Structure):
    _fields_ = \
        [
            ('min', ctypes.c_double),
            ('max', ctypes.c_double)
        ]


class c_string_list(ctypes.Structure):
    _fields_ = \
        [
            ('num_strings', ctypes.c_int64),
            ('string_list', ctypes.POINTER(ctypes.c_char_p))
        ]


class c_class_info(ctypes.Structure):
    _fields_ = \
        [
            ('name', ctypes.c_char_p),
            ('color', ctypes.c_int64),
            ('type', ctypes.c_int64)
        ]


class c_class_info_list(ctypes.Structure):
    _fields_ = \
        [
            ('num_class_infos', ctypes.c_int64),
            ('class_info_list', ctypes.POINTER(c_class_info))
        ]
