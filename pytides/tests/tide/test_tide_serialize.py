# encoding: utf-8

import datetime
import json
import numpy

from nose.tools import assert_equal, assert_is_instance

from pytides.constituent import get_constituent
from pytides.tide import Tide

from pytides import PytidesJsonEncoder, PytidesJsonDecoder

from ..utils import fixture_filename


def test_dump_tide_as_json():
    data = {'foo': _make_fake_tide()}
    string = json.dumps(data, cls=PytidesJsonEncoder, indent=2)

    assert_equal({'foo': {
        '__type__': 'pytides.Tide',
        'constituents': [
            {
                'constituent': {
                    '__type__': 'pytides.Constituent',
                    'name': 'Z0'},
                'amplitude': 1.2,
                'phase': 0.0
            },
            {
                'constituent': {
                    '__type__': 'pytides.Constituent',
                    'name': 'M2'},
                'amplitude': 0.25,
                'phase': 12.3
            },
            {
                'constituent': {
                    '__type__': 'pytides.Constituent',
                    'name': 'mu2'},
                'amplitude': 0.56,
                'phase': 45.6
            }
        ]}},
        json.loads(string))


def test_load_tide_from_json():
    string = """
    {
       "foo":{
          "__type__":"pytides.Tide",
          "constituents":[
             {
                "constituent":{
                   "__type__":"pytides.Constituent",
                   "name":"Z0"
                },
                "amplitude":1.2,
                "phase":0.0
             },
             {
                "constituent":{
                   "__type__":"pytides.Constituent",
                   "name":"M2"
                },
                "amplitude":0.25,
                "phase":12.3
             },
             {
                "constituent":{
                   "__type__":"pytides.Constituent",
                   "name":"mu2"
                },
                "amplitude":0.56,
                "phase":45.6
             }
          ]
       }
    }
    """
    actual_tide = json.loads(string, cls=PytidesJsonDecoder)['foo']

    expected_tide = _make_fake_tide()

    assert_is_instance(actual_tide, type(expected_tide))

    assert_equal(
        list(expected_tide.model),
        list(actual_tide.model))


def test_that_a_deserialized_tide_model_makes_correct_predictions():
    with open(fixture_filename('example_model_2.json'), 'r') as f:
        tide = json.loads(f.read(), cls=PytidesJsonDecoder)

        for dt, expected in [
                (datetime.datetime(2013, 1, 2, 12, 0, 0), 1.355),
                (datetime.datetime(2013, 1, 20, 0, 0, 0), 1.279),
                (datetime.datetime(2013, 2, 3, 13, 0, 0), 1.445),
                (datetime.datetime(2013, 5, 2, 8, 0, 0), 1.722),
                (datetime.datetime(2013, 7, 25, 16, 0, 0), 1.536),
                ]:
            yield assert_equal, expected, round(tide.at([dt])[0], 3)


def _make_fake_tide():
    model = numpy.zeros(3, dtype=Tide.dtype)
    model[0] = (get_constituent('Z0'), 1.2, 0.0)
    model[1] = (get_constituent('M2'), 0.25, 12.3)
    model[2] = (get_constituent('mu2'), 0.56, 45.6)

    return Tide(model=model)
