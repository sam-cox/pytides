import abc
import numpy as np


class BaseConstituent(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def speed(self, astro):
        pass

    @abc.abstractmethod
    def V(self, astro):
        pass

    @abc.abstractmethod
    def u(self, astro):
        pass

    @abc.abstractmethod
    def f(self, astro):
        pass

    def __init__(self, name):
        self.name = name

    def astro_xdo(self, a):
        return [a['T+h-s'], a['s'], a['h'], a['p'], a['N'], a['pp'], a['90']]

    def astro_speeds(self, a):
        return np.array([each.speed for each in self.astro_xdo(a)])

    def astro_values(self, a):
        return np.array([each.value for each in self.astro_xdo(a)])

    # Consider two out of phase constituents which travel at the same speed to
    # be identical
    def __eq__(self, c):
        return np.all(self.coefficients[:-1] == c.coefficients[:-1])

    def __hash__(self):
        return hash(tuple(self.coefficients[:-1]))

    def serialize(self):
        return {
            '__type__': 'pytides.Constituent',
            'name': self.name,
        }

    @staticmethod
    def deserialize(dictionary):
        from pytides.constituent import get_constituent
        return get_constituent(dictionary['name'])
