from histokat.image import Image
from histokat.backend.functions import \
    c_is_analyzable, \
    c_get_image, \
    c_analysis_was_performed_before, \
    c_get_default_analysis, \
    c_create_analysis, \
    c_close_analysis, \
    c_is_analysis_open, \
    c_save_image_at_level, \
    c_save_visualization_at_level, \
    c_analysis_get_config_dir, \
    c_analysis_get_workspace_dir, \
    c_analysis_get_id, \
    c_analysis_get_minor_modes, \
    c_analysis_get_class_infos, \
    c_analysis_perform, \
    c_analysis_inspect_object, \
    c_get_visualization


class Session:
    class _Analysis:
        def __init__(self, session_id):
            if session_id <= 0:
                raise Exception(
                    'Cannot create session from invalid session id')
            self._session_id = session_id

        def get_config_dir(self):
            return c_analysis_get_config_dir(self._session_id)

        def get_workspace_dir(self):
            return c_analysis_get_workspace_dir(self._session_id)

        def get_id(self):
            return c_analysis_get_id(self._session_id)

        def get_minor_modes(self, major_mode):
            return c_analysis_get_minor_modes(self._session_id, major_mode)

        def get_class_infos(self, input_mode):
            return c_analysis_get_class_infos(self._session_id, input_mode)

        def perform(self, analysis_mode):
            return c_analysis_perform(self._session_id, analysis_mode)

        def inspect_object(self, inspect_mode, position, level):
            return c_analysis_inspect_object(self._session_id, inspect_mode, position, level)

    def __init__(self, session_id):
        if session_id <= 0:
            raise Exception('Cannot create session from invalid session id')
        self._session_id = session_id
        self.analysis = None

    def get_id(self):
        return self._session_id

    def get_image(self):
        return Image(c_get_image(self._session_id, False))

    def is_analyzable(self, analysis_id):
        return c_is_analyzable(self._session_id, analysis_id)

    def analysis_was_performed_before(self, analysis_id):
        return c_analysis_was_performed_before(self._session_id, analysis_id)

    def get_default_analysis(self):
        return c_get_default_analysis(self._session_id)

    def create_analysis(self, analysis_id):
        c_create_analysis(self._session_id, analysis_id)
        self.is_analysis_open()

    def close_analysis(self):
        c_close_analysis(self._session_id)
        self.is_analysis_open()

    def is_analysis_open(self):
        is_open = c_is_analysis_open(self._session_id)
        # correct analysis member
        if is_open and not self.analysis:
            self.analysis = Session._Analysis(self._session_id)
        elif not is_open and self.analysis:
            self.analysis = None
        return is_open

    def saveImageAtLevel(self, image_path, z_pos, level):
        c_save_image_at_level(self._session_id, image_path, z_pos, level)

    def saveVisualizationAtLevel(self, image_path, z_pos, level, analysis_mode):
        c_save_visualization_at_level(
            self._session_id, image_path, z_pos, level, analysis_mode)

    def get_visualization(self, analysis_mode):
        return Image(c_get_visualization(self._session_id, analysis_mode))
