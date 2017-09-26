import numpy as np


def do_dos(data):
    """DOS correction using histogram with 2048 bins since WV2 has 11 bit radiometric resolution"""
    nonzero = data != 0
    n = np.sum(nonzero)
    tiny_fraction = n - n * 0.999999

    nbands = data.shape[0]
    dosDN = np.zeros(nbands)
    for i in range(nbands):
        hist, edges = np.histogram(
                data[i], bins=2048, range=(1, 2048), density=False)
        for k in range(1, len(hist)):
            if hist[k] - hist[i - 1] > tiny_fraction:
                dosDN[i] = k - 1
                break
    return dosDN
