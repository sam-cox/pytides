from nose.tools import assert_almost_equal

from pytides.observations import Observations
from pytides.tide import Tide

from ..utils import fixture_filename


def test_example_observations_2():
    observations = Observations.from_csv(
        fixture_filename('example_observations_2.csv'))

    datetimes, heights = observations.to_numpy_arrays()

    my_tide = Tide.decompose(heights, datetimes)

    predictions = my_tide.at(datetimes)

    for i, expected in [
            (0,    1.2301709375405498),
            (5,    1.4042197729365355),
            (10,   1.25361670249086),
            (50,   1.252261186358141),
            (100,  1.2792179066476739),
            (500,  1.3432777519478505),
            (1000, 1.7457095044631059),
            (1500, 1.472996756326626),
            (2000, 1.3003693570959916),
            (2500, 1.6887094920566177),
            (7500, 1.2907397986524105),
            ]:
        yield _assert_prediction_is_correct, expected, predictions[i]


def _assert_prediction_is_correct(expected, actual):
    assert_almost_equal(expected, actual)
