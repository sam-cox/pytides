import operator as op

from .base_constituent import BaseConstituent


class CompoundConstituent(BaseConstituent):
    def __init__(self, name, members):
        super(CompoundConstituent, self).__init__(name)
        self.members = members

        self.coefficients = reduce(
            op.add,
            [c.coefficients * n for (c, n) in members])

    def speed(self, a):
        return reduce(op.add, [n * c.speed(a) for (c, n) in self.members])

    def V(self, a):
        return reduce(op.add, [n * c.V(a) for (c, n) in self.members])

    def u(self, a):
        return reduce(op.add, [n * c.u(a) for (c, n) in self.members])

    def f(self, a):
        return reduce(op.mul, [c.f(a) ** abs(n) for (c, n) in self.members])

