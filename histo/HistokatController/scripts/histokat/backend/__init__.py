import atexit
import ctypes
import os
import shutil
import sys

#
# Add MLAB environment variables to PATH if present
#


def add_if_exists(env, postfix):
    if env in os.environ:
        path = os.path.abspath(os.path.join(os.environ[env], postfix))
        os.environ['PATH'] = os.pathsep.join([path, os.environ['PATH']])


add_if_exists('MLAB_ROOT', 'MeVis/ThirdParty/lib')
add_if_exists('MLAB_ROOT', 'MeVisLab/IDE/bin')
add_if_exists('MLAB_FMEwork_MeVis', 'lib')
add_if_exists('MLAB_FMEwork_ThirdParty', 'lib')
add_if_exists('MLAB_FMEwork_Internal', 'lib')

#
# Load the histkat controller library into the global variable
#
lib_name = 'HistokatController'
if 'linux' in sys.platform:
    lib_name = 'lib' + lib_name + '.so'
elif 'darwin' in sys.platform:
    lib_name = 'lib' + lib_name + '.dylib'

_HISTOKAT_LIB = ctypes.cdll.LoadLibrary(lib_name)

#
# workspace dir was generated
#


class _Workspace_dir_generated:
    def __init__(self):
        self.dir = None

    def update(self, new_dir):
        self.dir = new_dir


_workspace_dir_generated = _Workspace_dir_generated()

#
# remove generated dir at exit
#


def remove_generated_dir():
    if _workspace_dir_generated.dir is not None:
        shutil.rmtree(_workspace_dir_generated.dir)


atexit.register(remove_generated_dir)

#
# Define exceptions
#


class HistokatException(Exception):
    """Raise if an exception within the histokat module occurred."""


class NotInitializedException(HistokatException):
    """Raise if histokat.init has not been called."""


class InvalidIdException(HistokatException):
    """Raise if given ID is not valid."""


#
# Init Cancellation Checker
#
_is_cancelled = False


def cancel_current_activity():
    global _is_cancelled
    _is_cancelled = True


def is_cancelled():
    global _is_cancelled
    return _is_cancelled
    _is_cancelled = False


_cfunc_cancellation_callback = ctypes.CFUNCTYPE(ctypes.c_bool)(is_cancelled)
