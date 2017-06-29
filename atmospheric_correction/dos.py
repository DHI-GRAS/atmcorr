import numpy as np


def do_dos(img):
    """DOS correction using histogram with 2048 bins since WV2 has 11 bit radiometric resolution"""
    dosDN = []
    tempData = img.GetRasterBand(1).ReadAsArray()
    numElements = np.size(tempData[tempData != 0])

    tiny_fraction = numElements - numElements * 0.999999

    nbands = img.RasterCount
    for b in range(nbands):
        hist, edges = np.histogram(
                img.GetRasterBand(b + 1).ReadAsArray(),
                bins=2048, range=(1, 2048), density=False)
        for i in range(1, len(hist)):
            if hist[i] - hist[i - 1] > tiny_fraction:
                dosDN.append(i - 1)
                break
    return dosDN
