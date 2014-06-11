from .constituent import Constituent
from .compound_constituent import CompoundConstituent

from ..nodal_corrections import (
    f_unity, f_Modd, u_zero, u_Modd,
    f_Mm, f_Mf, f_O1, f_K1, f_J1, f_M1, f_M2, f_OO1, f_L2, f_K2,
    u_Mf, u_O1, u_K1, u_J1, u_M1, u_M2, u_OO1, u_L2, u_K2)


def get_constituent(name):
    return CONSTITUENTS[name]


def get_constituent_set(name):
    return CONSTITUENT_SETS[name]


# ## ## # Constituents
# Long Term
_Z0 = Constituent(name='Z0',  xdo='Z ZZZ ZZZ', u=u_zero, f=f_unity)
_Sa = Constituent(name='Sa',  xdo='Z ZAZ ZZZ', u=u_zero, f=f_unity)
_Ssa = Constituent(name='Ssa', xdo='Z ZBZ ZZZ', u=u_zero, f=f_unity)
_Mm = Constituent(name='Mm',  xdo='Z AZY ZZZ', u=u_zero, f=f_Mm)
_Mf = Constituent(name='Mf',  xdo='Z BZZ ZZZ', u=u_Mf, f=f_Mf)

# Diurnals
_Q1 = Constituent(name='Q1',      xdo='A XZA ZZA', u=u_O1, f=f_O1)
_O1 = Constituent(name='O1',      xdo='A YZZ ZZA', u=u_O1, f=f_O1)
_K1 = Constituent(name='K1',      xdo='A AZZ ZZY', u=u_K1, f=f_K1)
_J1 = Constituent(name='J1',      xdo='A BZY ZZY', u=u_J1, f=f_J1)

# M1 is a tricky business for reasons of convention, rather than theory.  The
# reasons for this are best summarised by Schureman paragraphs 126, 127 and in
# the comments found in congen_input.txt of xtides, so I won't go over all this
# again here.

_M1 = Constituent(name='M1',   xdo='A ZZZ ZZA', u=u_M1, f=f_M1)
_P1 = Constituent(name='P1',   xdo='A AXZ ZZA', u=u_zero, f=f_unity)
_S1 = Constituent(name='S1',   xdo='A AYZ ZZZ', u=u_zero, f=f_unity)
_OO1 = Constituent(name='OO1', xdo='A CZZ ZZY', u=u_OO1, f=f_OO1)

# Semi-Diurnals
_2N2 = Constituent(name='2N2', xdo='B XZB ZZZ', u=u_M2, f=f_M2)
_N2 = Constituent(name='N2',   xdo='B YZA ZZZ', u=u_M2, f=f_M2)
_nu2 = Constituent(name='nu2', xdo='B YBY ZZZ', u=u_M2, f=f_M2)
_M2 = Constituent(name='M2',   xdo='B ZZZ ZZZ', u=u_M2, f=f_M2)
_lambda2 = Constituent(name='lambda2', xdo='B AXA ZZB',
                       u=u_M2, f=f_M2)
_L2 = Constituent(name='L2',   xdo='B AZY ZZB', u=u_L2, f=f_L2)
_T2 = Constituent(name='T2',   xdo='B BWZ ZAZ', u=u_zero, f=f_unity)
_S2 = Constituent(name='S2',   xdo='B BXZ ZZZ', u=u_zero, f=f_unity)
_R2 = Constituent(name='R2',   xdo='B BYZ ZYB', u=u_zero, f=f_unity)
_K2 = Constituent(name='K2',   xdo='B BZZ ZZZ', u=u_K2, f=f_K2)

# Third-Diurnals
_M3 = Constituent(name='M3', xdo='C ZZZ ZZZ',
                  u=lambda a: u_Modd(a, 3),
                  f=lambda a: f_Modd(a, 3))

# ## ## # Compound Constituents
# Long Term
_MSF = CompoundConstituent(name='MSF',  members=[(_S2, 1), (_M2, -1)])

# Diurnal
_2Q1 = CompoundConstituent(name='2Q1',  members=[(_N2, 1), (_J1, -1)])
_rho1 = CompoundConstituent(name='rho1', members=[(_nu2, 1), (_K1, -1)])

# Semi-Diurnal

_mu2 = CompoundConstituent(name='mu2',  members=[(_M2, 2), (_S2, -1)])  # 2MS2
_2SM2 = CompoundConstituent(name='2SM2', members=[(_S2, 2), (_M2, -1)])

# Third-Diurnal
_2MK3 = CompoundConstituent(name='2MK3', members=[(_M2, 1), (_O1, 1)])
_MK3 = CompoundConstituent(name='MK3',  members=[(_M2, 1), (_K1, 1)])

# Quarter-Diurnal
_MN4 = CompoundConstituent(name='MN4',  members=[(_M2, 1), (_N2, 1)])
_M4 = CompoundConstituent(name='M4',   members=[(_M2, 2)])
_MS4 = CompoundConstituent(name='MS4',  members=[(_M2, 1), (_S2, 1)])
_S4 = CompoundConstituent(name='S4',   members=[(_S2, 2)])

# Sixth-Diurnal
_M6 = CompoundConstituent(name='M6',   members=[(_M2, 3)])
_S6 = CompoundConstituent(name='S6',   members=[(_S2, 3)])

# Eighth-Diurnals
_M8 = CompoundConstituent(name='M8',   members=[(_M2, 4)])


noaa = [
    _M2, _S2, _N2, _K1, _M4, _O1, _M6, _MK3, _S4, _MN4, _nu2, _S6, _mu2,
    _2N2, _OO1, _lambda2, _S1, _M1, _J1, _Mm, _Ssa, _Sa, _MSF, _Mf,
    _rho1, _Q1, _T2, _R2, _2Q1, _P1, _2SM2, _M3, _L2, _2MK3, _K2,
    _M8, _MS4
]

CONSTITUENTS = {
    'Z0': _Z0,
    'Sa': _Sa,
    'Ssa': _Ssa,
    'Mm': _Mm,
    'Mf': _Mf,
    'Q1': _Q1,
    'O1': _O1,
    'K1': _K1,
    'J1': _J1,
    'M1': _M1,
    'P1': _P1,
    'S1': _S1,
    'OO1': _OO1,
    '2N2': _2N2,
    'N2': _N2,
    'nu2': _nu2,
    'M2': _M2,
    'lambda2': _lambda2,
    'L2': _L2,
    'T2': _T2,
    'S2': _S2,
    'R2': _R2,
    'K2': _K2,
    'M3': _M3,
    'MSF': _MSF,
    '2Q1': _2Q1,
    'rho1': _rho1,
    'mu2': _mu2,
    '2SM2': _2SM2,
    '2MK3': _2MK3,
    'MK3': _MK3,
    'MN4': _MN4,
    'M4': _M4,
    'MS4': _MS4,
    'S4': _S4,
    'M6': _M6,
    'S6': _S6,
    'M8': _M8,
}

CONSTITUENT_SETS = {
    'NOAA': noaa,
}
