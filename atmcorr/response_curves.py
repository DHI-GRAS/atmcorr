import sensor_response_curves as srcurves
import sensor_response_curves.resample as srcurves_resample


def get_response_curves(sensor, band_ids):
    # Set 6S band filters
    wavelength, rcurves = srcurves.get_response_curves(
            sensor, band_ids=band_ids)
    wavelength = wavelength.astype('float') / 1e3
    # Also need to resample the band filters from 1nm to 2.5nm
    # as this is the highest spectral resolution supported by 6S
    wavelength, rcurves = srcurves_resample.resample_response_curves(
                wavelength, rcurves, resolution=0.0025)
    return dict(
        rcurves=rcurves,
        start_wv=wavelength[0],
        end_wv=wavelength[-1])
