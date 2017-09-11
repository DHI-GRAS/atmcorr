import warnings
import functools
from contextlib import contextmanager

import numpy as np


@contextmanager
def ignore_nan_warnings():
    old_settings = np.seterr(invalid='ignore')
    yield
    np.seterr(**old_settings)


def wrap_ignore_nan_warnings(func):
    """Ignore NaN warnings"""

    @functools.wraps(func)
    def wrapped(*args, **kwargs):
        with ignore_nan_warnings():
            return func(*args, **kwargs)


def wrap_nowarn(func):
    @functools.wraps(func)
    def wrapped(*args, **kwargs):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return func(*args, **kwargs)
    return wrapped
