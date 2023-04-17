from atmcorr.processing import main
import json

from bathymetry.utils import inout

dn_data_bandfirst, dn_profile = inout.read_file(
    infile='./outdir/Input/mosaic.tif',
    aoi_file='./outdir/SHPs/AOI.geojson',
    bandlast=False)

o_atmcorr = {'sixs_params': {'atm': {'AOT': 0.00021172959435072913, 'PWV': 1.0, 'ozone': 0.15},
                            'aeroProfile': 'Maritime'},
             'adjCorr': 0,
             'aotMultiplier': 1.0,
             'tileSize': 2000,
             'sensor': 'PNEO',
             'mtdFile': './metadata/IMG_01_PNEO3_MS-FS/DIM_PNEO3_202210271846151_MS-FS_ORT_bcda267352f4468bc6de24b84b0b33f3.XML',
             'band_ids': [0, 1, 2, 5],
             'interpolation': True
             }

o_atmcorr.update(
    data=dn_data_bandfirst,
    profile=dn_profile
)
print(dn_profile)
refl_data_bandfirst, profile = main(**o_atmcorr)
refl_data = inout.to_bandlast(refl_data_bandfirst)

inout.save(data=refl_data, outfile='outfile.tif', profile=profile, is_bandlast=True)
