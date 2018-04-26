import collections


class ImageNotAcceptedException(Exception):
    """Raise if session is created with file not being an accepted image."""


class AnalysisId:
    def __init__(self, name=None, version=None):
        if name is None or version is None:
            self.name = '<no analysis>'
            self.version = 0
        self.name = name
        self.version = version

    def __repr__(self):
        return '{0} {1}'.format(self.name, self.version)

    def __gt__(self, other):
        return (self.name, self.version) > (other.name, other.version)

    def __eq__(self, other):
        return self.name == other.name and self.version == other.version


class AnalysisMode:
    class MajorMode:
        MODE_ANALYSE = 0
        MODE_INPUT = 1
        MODE_INSPECT = 2

    def __init__(self, major, minor=0):
        self.major = major
        self.minor = minor

    @classmethod
    def analyse(cls):
        return cls(AnalysisMode.MajorMode.MODE_ANALYSE)

    def __gt__(self, other):
        return (self.major, self.minor) > (other.major, other.minor)


class ImageInfo:
    def __init__(self, image_id, image_url):
        self.image_id = image_id
        self.image_url = image_url


class ImageGroup:
    def __init__(self, group_id, group_url, image_list):
        self.group_id = group_id
        self.group_url = group_url
        self.image_list = image_list


class Vector:
    def __init__(self, vec):
        """
        Vector is initialized by a tuple or list of 3 values or by another Vector.
        """
        if not isinstance(vec, collections.Iterable):
            raise TypeError(
                'Vector cannot be instantiated with object of type %s.' % type(vec))
        elif len(vec) < 3:
            raise TypeError(
                'Vector cannot be instantiated with less than 3 values.')
        self.x = int(vec[0])
        self.y = int(vec[1])
        self.z = int(vec[2])

    def __getitem__(self, key):
        return list([self.x, self.y, self.z]).__getitem__(key)

    def get_prod(self):
        return self.x * self.y * self.z

    def __eq__(self, other):
        return (self.x, self.y, self.z) == (other.x, other.y, other.z)

    def __repr__(self):
        return 'x: %d, y: %d, z: %d' % tuple(self)


class Box:
    def __init__(self, vectors):
        """
        Box is initialized by a tuple or list of two vectors.
        The two vectors can also be tuples or lists of 3 values.
        """
        if not isinstance(vectors, collections.Iterable):
            raise TypeError(
                'Box cannot be instantiated with object of type %s.' % type(vectors))
        elif len(vectors) < 2:
            raise TypeError(
                'Box cannot be instantiated with less than 2 vectors.')

        self.v1 = vectors[0] if isinstance(vectors[0], Vector) else Vector(vectors[0])
        self.v2 = vectors[1] if isinstance(vectors[1], Vector) else Vector(vectors[1])

    def __getitem__(self, key):
        return list([self.v1, self.v2]).__getitem__(key)

    def get_extent(self):
        return Vector((
            self.v2.x - self.v1.x + 1,
            self.v2.y - self.v1.y + 1,
            self.v2.z - self.v1.z + 1))

    def __repr__(self):
        return 'v1: %s, v2: %s' % tuple(self)


class ValueRange:
    def __init__(self, min_value, max_value):
        self.min = min_value
        self.max = max_value


class ClassInfo:
    def __init__(self, name, color, class_type):
        self.name = name
        self.color = color
        self.class_type = class_type
