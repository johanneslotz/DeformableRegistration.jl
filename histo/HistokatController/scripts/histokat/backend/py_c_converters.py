from histokat.backend.types import c_analysis_id, c_box, c_analysis_mode, c_vector
from histokat.utils import AnalysisId, ImageInfo, ImageGroup, Vector, ValueRange, ClassInfo


def convert_py_to_c_analysis_id(analysis_id):
    c_id = c_analysis_id()
    c_id.name = analysis_id.name.encode('utf-8')
    c_id.version = analysis_id.version
    return c_id


def convert_c_ptr_to_py_analysis_id(c_analysis_id_ptr):
    return AnalysisId(
        c_analysis_id_ptr.contents.name.decode('utf-8'),
        c_analysis_id_ptr.contents.version)


def convert_c_to_py_analysis_id_list_ptr(c_analysis_id_list_ptr):
    analysis_ids = []
    for i in range(c_analysis_id_list_ptr.contents.num_ids):
        analysis_ids.append(AnalysisId(
            c_analysis_id_list_ptr.contents.analysis_ids[i].name.decode(
                'utf-8'),
            c_analysis_id_list_ptr.contents.analysis_ids[i].version))
    return analysis_ids


def convert_py_to_c_analysis_mode(analysis_mode):
    c_mode = c_analysis_mode()
    c_mode.major = analysis_mode.major
    c_mode.minor = analysis_mode.minor
    return c_mode


def convert_c_ptr_to_py_image_group(c_image_group_ptr):
    group_id = c_image_group_ptr.contents.group_id.decode('utf-8')
    group_url = c_image_group_ptr.contents.group_url.decode('utf-8')
    image_list = []
    for i in range(c_image_group_ptr.contents.num_image_infos):
        image_list.append(ImageInfo(
            c_image_group_ptr.contents.image_info_list[i].image_id.decode(
                'utf-8'),
            c_image_group_ptr.contents.image_info_list[i].image_url.decode('utf-8')))
    return ImageGroup(group_id, group_url, image_list)


def convert_py_to_c_vector(vector):
    c_vec = c_vector()
    c_vec.x = vector[0]
    c_vec.y = vector[1]
    c_vec.z = vector[2]
    return c_vec


def convert_c_to_py_vector_ptr(c_vector_ptr):
    return Vector((c_vector_ptr.contents.x, c_vector_ptr.contents.y, c_vector_ptr.contents.z))


def convert_py_to_c_box(box):
    c_bx = c_box()
    c_bx.v1.x, c_bx.v1.y, c_bx.v1.z = box.v1.x, box.v1.y, box.v1.z
    c_bx.v2.x, c_bx.v2.y, c_bx.v2.z = box.v2.x, box.v2.y, box.v2.z
    return c_bx


def convert_c_ptr_to_py_value_range(c_value_range_ptr):
    return ValueRange(c_value_range_ptr.contents.min, c_value_range_ptr.contents.max)


def convert_c_ptr_to_py_class_info_list(c_class_info_list_ptr):
    classes = []
    for i in range(c_class_info_list_ptr.contents.num_class_infos):
        classes.append(ClassInfo(
            c_class_info_list_ptr.contents.class_info_list[i].name.decode(
                'utf-8'),
            c_class_info_list_ptr.contents.class_info_list[i].color,
            c_class_info_list_ptr.contents.class_info_list[i].type
        ))
    return classes
