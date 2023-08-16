import datetime

import pytz

from nose.tools import assert_equal

from pytides.observations import Observations
from ..utils import fixture_filename


class TestCsvObservations(object):

    @classmethod
    def setup_class(cls):
        cls.csv_obs = Observations.from_csv(fixture_filename(
            'example_observations_2.csv'))
        cls.np_datetimes, cls.np_heights = cls.csv_obs.to_numpy_arrays()

    @classmethod
    def _get_first_five_datetimes_from_csv(cls):
        it = iter(cls.csv_obs)
        return [next(it)[0] for i in range(5)]

    @classmethod
    def _get_first_five_heights_from_csv(cls):
        it = iter(cls.csv_obs)
        return [next(it)[1] for i in range(5)]

    def test_that_correct_number_of_observations_are_loaded(self):
        assert_equal(8784, len(list(self.csv_obs)))

    def test_that_first_five_datetimes_are_correct(self):
        datetimes = self._get_first_five_datetimes_from_csv()
        expected = [
            datetime.datetime(2013, 1, 1, 0, 0, 0, tzinfo=pytz.UTC),
            datetime.datetime(2013, 1, 1, 1, 0, 0, tzinfo=pytz.UTC),
            datetime.datetime(2013, 1, 1, 2, 0, 0, tzinfo=pytz.UTC),
            datetime.datetime(2013, 1, 1, 3, 0, 0, tzinfo=pytz.UTC),
            datetime.datetime(2013, 1, 1, 4, 0, 0, tzinfo=pytz.UTC)
        ]
        assert_equal(expected, datetimes)

    def test_that_first_five_heights_are_correct(self):
        heights = self._get_first_five_heights_from_csv()
        expected = [1.225, 1.275, 1.334, 1.390, 1.411]
        assert_equal(expected, heights)

    def test_that_observations_can_be_iterated(self):
        datetime, height = next(iter(self.csv_obs))

    def test_that_observations_iterators_are_independent(self):
        datetime_1, height_1 = next(iter(self.csv_obs))
        datetime_2, height_2 = next(iter(self.csv_obs))

        assert_equal(datetime_1, datetime_2)

    def test_that_numpy_datetimes_array_is_correct_size(self):
        assert_equal(8784, len(self.np_datetimes))

    def test_that_numpy_heights_array_is_correct_size(self):
        assert_equal(8784, len(self.np_heights))

    def test_that_first_five_numpy_datetimes_are_correct(self):
        expected_datetimes = self._get_first_five_datetimes_from_csv()
        numpy_datetimes = list(self.np_datetimes[0:5])

        assert_equal(expected_datetimes, numpy_datetimes)

    def test_that_first_five_numpy_heights_are_correct(self):
        expected_heights = self._get_first_five_heights_from_csv()
        numpy_heights = list(self.np_heights[0:5])

        assert_equal(expected_heights, numpy_heights)
