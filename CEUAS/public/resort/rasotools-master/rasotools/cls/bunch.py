import numpy as np


def dim_summary(obj):
    """ Auxilliary class function for printing
    """
    if hasattr(obj, 'data_vars'):
        return "%d vars [%s]" % (len(obj.data_vars), ", ".join(["%s(%s)" % (i, j) for i, j in obj.dims.items()]))

    if hasattr(obj, 'shape'):
        i = np.shape(obj)
        if i != ():
            return i

    if hasattr(obj, 'size'):
        i = np.size(obj)
        if i == 1:
            return obj
        else:
            return i

    if isinstance(obj, (list, tuple, dict)):
        i = len(obj)
        if i == 1:
            return obj
        else:
            return i

    return obj


def formatting(obj):
    return u'<%s (%s)>' % (type(obj).__name__, dim_summary(obj))


class Bunch(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __repr__(self):
        return u"\n".join("%-10s : %s" % (i, formatting(getattr(self, i))) for i in self.__dict__.keys())

    def __iter__(self):
        # iter(self.__dict__)   # does the same
        return (x for x in self.__dict__.keys())

    def __getitem__(self, item):
        return getattr(self, item)

    def __setitem__(self, key, value):
        setattr(self, key, value)

    def __delitem__(self, key):
        if key in self.__dict__.keys():
            delattr(self, key)
