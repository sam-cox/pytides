import json

from nose.tools import assert_equal, assert_is_instance, assert_is

from pytides.constituent import get_constituent

from pytides import PytidesJsonEncoder, PytidesJsonDecoder


def test_dump_simple_constituent_as_json():
    data = {'my_constituents': [get_constituent('Z0'), get_constituent('M2')]}
    string = json.dumps(data, cls=PytidesJsonEncoder, indent=2)
    assert_equal(
        {'my_constituents': [
            {
                '__type__': 'pytides.Constituent',
                'name': 'Z0',
            },
            {
                '__type__': 'pytides.Constituent',
                'name': 'M2',
            }
        ]},
        json.loads(string))


def test_load_simple_constituent_from_json():
    string = ('{"foo": {'
              '  "__type__": "pytides.Constituent",'
              '  "name": "M2"'
              '  }'
              '}')
    expected_m2 = get_constituent('M2')
    actual_m2 = json.loads(string, cls=PytidesJsonDecoder)['foo']

    assert_is_instance(actual_m2, type(expected_m2))
    assert_is(expected_m2, actual_m2)
