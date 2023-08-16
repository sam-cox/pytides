from contextlib import contextmanager
from os.path import abspath, dirname, join as pjoin

FIXTURE_DIR = pjoin(abspath(dirname(__file__)), 'fixtures')


def fixture_filename(basename):
    return pjoin(FIXTURE_DIR, basename)


@contextmanager
def open_fixture(basename, *args):
    with open(fixture_filename(basename), *args) as f:
        yield f
