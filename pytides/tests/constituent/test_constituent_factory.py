from nose.tools import assert_equal, assert_is_instance

from pytides.constituent import (
    get_constituent, get_constituent_set, Constituent, CompoundConstituent)


class TestGetConstituent(object):
    def test_get_constituent_gets_a_normal_constituent_type(self):
        z0 = get_constituent('Z0')
        assert_is_instance(z0, Constituent)

    def test_get_constituent_gets_a_compound_constituent_type(self):
        z0 = get_constituent('mu2')
        assert_is_instance(z0, CompoundConstituent)

    def test_get_z0_constituent_has_correct_name(self):
        z0 = get_constituent('Z0')
        assert_equal('Z0', z0.name)

    def test_get_mu2_constituent_has_correct_name(self):
        z0 = get_constituent('mu2')
        assert_equal('mu2', z0.name)


class TestGetConstituentSet(object):
    def test_get_noaa_constituents_returns_correct_number(self):
        noaa_constituents = get_constituent_set('NOAA')
        assert_equal(37, len(noaa_constituents))

    def test_get_noaa_constituents_returns_correct_names(self):
        names = set([c.name for c in get_constituent_set('NOAA')])
        assert_equal(
            set([
                'Q1', 'MS4', 'S2', 'rho1', 'S6', 'S4', 'M4', 'K1', 'M6',
                'M1', 'M3', 'M2', 'MSF', 'Sa', 'M8', 'O1', '2Q1', 'lambda2',
                '2N2', 'Ssa', 'S1', 'Mf', 'R2', 'Mm', 'MN4', 'T2', 'nu2',
                '2MK3', 'N2', 'OO1', 'P1', 'K2', 'mu2', 'J1', '2SM2', 'MK3',
                'L2']),
            names)
