import abc
import csv
import datetime

from itertools import tee

import pytz
import numpy


class Observations(object):
    __metaclass__ = abc.ABCMeta

    @staticmethod
    def from_csv(filename):
        return CsvObservations(filename)

    @staticmethod
    def separate_datetimes_heights(observations):
        """
        Return a pair of generators: datetimes and heights from a single
        generator of (datetime, height) tuples.
        """
        x, y = tee(observations, 2)  # copy the generator, can't rewind it.
        return (tup[0] for tup in x), (tup[1] for tup in y)

    @abc.abstractmethod
    def to_numpy_arrays(self):
        """
        Return (datetimes, heights) where both are numpy arrays.
        """
        pass

    @abc.abstractmethod
    def __iter__(self):
        """
        Yield observations values as (datetime, height) tuples.

        This must iterate from the start of the observations every time it's
        called.
        """
        pass


class CsvObservations(Observations):
    def __init__(self, filename):
        self.filename = filename

    def __iter__(self):
        return self.observations()

    def to_numpy_arrays(self):
        datetimes, heights = Observations.separate_datetimes_heights(
            self.observations())
        return (
            numpy.array(list(datetimes)),
            numpy.array(list(heights))
        )

    def observations(self):
        with open(self.filename, 'r') as f:
            reader = csv.DictReader(
                line for line in f if not line.startswith('#'))
            for row in reader:
                yield (CsvObservations.parse_datetime(row['datetime']),
                       float(row['height']))

    @staticmethod
    def parse_datetime(string):
        return datetime.datetime.strptime(
            string, '%Y-%m-%dT%H:%M:%S+00:00').replace(tzinfo=pytz.UTC)
