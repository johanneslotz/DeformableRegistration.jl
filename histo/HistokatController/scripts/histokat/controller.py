import histokat.backend.functions as c_funcs
from histokat.backend.functions import c_get_analysis_ids as get_analysis_ids
from histokat.backend.functions import \
    c_get_controller_workspace_dir as get_workspace_dir
from histokat.backend.functions import c_init_controller as init
from histokat.backend.functions import c_is_image_openable as is_image_openable
from histokat.backend.functions import c_read_image_group as read_image_group
from histokat.backend.functions import c_close_session
from histokat.session import Session
from histokat.utils import ImageNotAcceptedException


def open_session(url):
    if not is_image_openable(url):
        raise ImageNotAcceptedException()
    return Session(c_funcs.c_open_session(url))


def close_session(session):
    session.close_analysis()
    c_funcs.c_close_session(session.get_id())


__all__ = ('init', 'read_image_group',
           'get_analysis_ids', 'is_image_openable',
           'open_session', 'close_session')
