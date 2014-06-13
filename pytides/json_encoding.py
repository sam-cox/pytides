# encoding: utf-8

import json

from pytides.constituent import BaseConstituent
from pytides.tide import Tide


class PytidesJsonEncoder(json.JSONEncoder):
    def default(self, obj):
        if hasattr(obj, 'serialize'):
            return obj.serialize()
        return json.JSONEncoder.default(self, obj)


class PytidesJsonDecoder(json.JSONDecoder):
    def __init__(self, *args, **kwargs):
        json.JSONDecoder.__init__(self, object_hook=self.object_hook,
                                  *args, **kwargs)

    def object_hook(self, d):
        type_name = d.get('__type__')

        cls = {
            'pytides.Constituent': BaseConstituent,
            'pytides.Tide': Tide
        }.get(type_name)

        if cls is not None:
            return cls.deserialize(d)
        return d
