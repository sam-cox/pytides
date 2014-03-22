
import string
import operator as op
import numpy as np
import nodal_corrections as nc

class BaseConstituent(object):
	xdo_int = {
		'A': 1, 'B': 2, 'C': 3, 'D': 4, 'E': 5, 'F': 6, 'G': 7, 'H': 8, 'I': 9,
		'J': 10, 'K': 11, 'L': 12, 'M': 13, 'N': 14, 'O': 15, 'P': 16, 'Q': 17,
		'R': -8, 'S': -7, 'T': -6, 'U': -5, 'V': -4, 'W': -3, 'X': -2, 'Y': -1,
		'Z': 0
	}

	int_xdo = {v:k for k, v in xdo_int.items()}

	def __init__(self, name, xdo='', coefficients=[], u=nc.u_zero, f=nc.f_unity):
		if xdo == '':
			self.coefficients = np.array(coefficients)
		else:
			self.coefficients = np.array(self.xdo_to_coefficients(xdo))
		self.name = name
		self.u = u
		self.f = f

	def xdo_to_coefficients(self, xdo):
		return [self.xdo_int[l.upper()] for l in xdo if l in string.ascii_letters]

	def coefficients_to_xdo(self, coefficients):
		return ''.join([self.int_xdo[c] for c in cooefficients])

	def V(self, astro):
		return np.dot(self.coefficients, self.astro_values(astro))

	def xdo(self):
		return self.coefficients_to_xdo(self.coefficients)

	def speed(self, a):
		return np.dot(self.coefficients, self.astro_speeds(a))

	def astro_xdo(self, a):
		return [a['T+h-s'], a['s'], a['h'], a['p'], a['N'], a['pp'], a['90']]

	def astro_speeds(self, a):
		return np.array([each.speed for each in self.astro_xdo(a)])

	def astro_values(self, a):
		return np.array([each.value for each in self.astro_xdo(a)])

	#Consider two out of phase constituents which travel at the same speed to
	#be identical
	def __eq__(self, c):
		return np.all(self.coefficients[:-1] == c.coefficients[:-1])

	def __hash__(self):
		return hash(tuple(self.coefficients[:-1]))

class CompoundConstituent(BaseConstituent):

	def __init__(self, members = [], **kwargs):
		self.members = members

		if 'u' not in kwargs:
			kwargs['u'] = self.u
		if 'f' not in kwargs:
			kwargs['f'] = self.f

		super(CompoundConstituent,self).__init__(**kwargs)

		self.coefficients = reduce(op.add,[c.coefficients * n for (c,n) in members])

	def speed(self, a):
		return reduce(op.add, [n * c.speed(a) for (c,n) in self.members])

	def V(self, a):
		return reduce(op.add, [n * c.V(a) for (c,n) in self.members])

	def u(self, a):
		return reduce(op.add, [n * c.u(a) for (c,n) in self.members])

	def f(self, a):
		return reduce(op.mul, [c.f(a) ** abs(n) for (c,n) in self.members])

###### Base Constituents
#Long Term
_Z0      = BaseConstituent(name = 'Z0',      xdo = 'Z ZZZ ZZZ', u = nc.u_zero, f = nc.f_unity)
_Sa      = BaseConstituent(name = 'Sa',      xdo = 'Z ZAZ ZZZ', u = nc.u_zero, f = nc.f_unity)
_Ssa     = BaseConstituent(name = 'Ssa',     xdo = 'Z ZBZ ZZZ', u = nc.u_zero, f = nc.f_unity)
_Mm      = BaseConstituent(name = 'Mm',      xdo = 'Z AZY ZZZ', u = nc.u_zero, f = nc.f_Mm)
_Mf      = BaseConstituent(name = 'Mf',      xdo = 'Z BZZ ZZZ', u = nc.u_Mf, f = nc.f_Mf)

#Diurnals
_Q1      = BaseConstituent(name = 'Q1',      xdo = 'A XZA ZZA', u = nc.u_O1, f = nc.f_O1)
_O1      = BaseConstituent(name = 'O1',      xdo = 'A YZZ ZZA', u = nc.u_O1, f = nc.f_O1)
_K1      = BaseConstituent(name = 'K1',      xdo = 'A AZZ ZZY', u = nc.u_K1, f = nc.f_K1)
_J1      = BaseConstituent(name = 'J1',      xdo = 'A BZY ZZY', u = nc.u_J1, f = nc.f_J1)

#M1 is a tricky business for reasons of convention, rather than theory.  The
#reasons for this are best summarised by Schureman paragraphs 126, 127 and in
#the comments found in congen_input.txt of xtides, so I won't go over all this
#again here.

_M1      = BaseConstituent(name = 'M1',      xdo = 'A ZZZ ZZA', u = nc.u_M1, f = nc.f_M1)
_P1      = BaseConstituent(name = 'P1',      xdo = 'A AXZ ZZA', u = nc.u_zero, f = nc.f_unity)
_S1      = BaseConstituent(name = 'S1',      xdo = 'A AYZ ZZZ', u = nc.u_zero, f = nc.f_unity)
_OO1     = BaseConstituent(name = 'OO1',     xdo = 'A CZZ ZZY', u = nc.u_OO1, f = nc.f_OO1)

#Semi-Diurnals
_2N2     = BaseConstituent(name = '2N2',     xdo = 'B XZB ZZZ', u = nc.u_M2, f = nc.f_M2)
_N2      = BaseConstituent(name = 'N2',      xdo = 'B YZA ZZZ', u = nc.u_M2, f = nc.f_M2)
_nu2     = BaseConstituent(name = 'nu2',     xdo = 'B YBY ZZZ', u = nc.u_M2, f = nc.f_M2)
_M2      = BaseConstituent(name = 'M2',      xdo = 'B ZZZ ZZZ', u = nc.u_M2, f = nc.f_M2)
_lambda2 = BaseConstituent(name = 'lambda2', xdo = 'B AXA ZZB', u = nc.u_M2, f = nc.f_M2)
_L2      = BaseConstituent(name = 'L2',      xdo = 'B AZY ZZB', u = nc.u_L2, f = nc.f_L2)
_T2      = BaseConstituent(name = 'T2',      xdo = 'B BWZ ZAZ', u = nc.u_zero, f = nc.f_unity)
_S2      = BaseConstituent(name = 'S2',      xdo = 'B BXZ ZZZ', u = nc.u_zero, f = nc.f_unity)
_R2      = BaseConstituent(name = 'R2',      xdo = 'B BYZ ZYB', u = nc.u_zero, f = nc.f_unity)
_K2      = BaseConstituent(name = 'K2',      xdo = 'B BZZ ZZZ', u = nc.u_K2, f = nc.f_K2)

#Third-Diurnals
_M3      = BaseConstituent(name = 'M3',      xdo = 'C ZZZ ZZZ', u = lambda a: nc.u_Modd(a,3), f = lambda a: nc.f_Modd(a,3))

###### Compound Constituents
#Long Term
_MSF     = CompoundConstituent(name = 'MSF',  members = [(_S2, 1), (_M2, -1)])

#Diurnal
_2Q1     = CompoundConstituent(name = '2Q1',  members = [(_N2, 1), (_J1, -1)])
_rho1    = CompoundConstituent(name = 'rho1', members = [(_nu2, 1), (_K1, -1)])

#Semi-Diurnal

_mu2     = CompoundConstituent(name = 'mu2',  members = [(_M2, 2), (_S2, -1)]) #2MS2
_2SM2    = CompoundConstituent(name = '2SM2', members = [(_S2, 2), (_M2, -1)])

#Third-Diurnal
_2MK3    = CompoundConstituent(name = '2MK3', members = [(_M2, 1), (_O1, 1)])
_MK3     = CompoundConstituent(name = 'MK3',  members = [(_M2, 1), (_K1, 1)])

#Quarter-Diurnal
_MN4     = CompoundConstituent(name = 'MN4',  members = [(_M2, 1), (_N2, 1)])
_M4      = CompoundConstituent(name = 'M4',   members = [(_M2, 2)])
_MS4     = CompoundConstituent(name = 'MS4',  members = [(_M2, 1), (_S2, 1)])
_S4      = CompoundConstituent(name = 'S4',   members = [(_S2, 2)])

#Sixth-Diurnal
_M6      = CompoundConstituent(name = 'M6',   members = [(_M2, 3)])
_S6      = CompoundConstituent(name = 'S6',   members = [(_S2, 3)])

#Eighth-Diurnals
_M8      = CompoundConstituent(name = 'M8',   members = [(_M2, 4)])


noaa = [
	_M2, _S2, _N2, _K1, _M4, _O1, _M6, _MK3, _S4, _MN4, _nu2, _S6, _mu2,
	_2N2, _OO1, _lambda2, _S1, _M1, _J1, _Mm, _Ssa, _Sa, _MSF, _Mf,
	_rho1, _Q1, _T2, _R2, _2Q1, _P1, _2SM2, _M3, _L2, _2MK3, _K2,
	_M8, _MS4
]

