


def read_landsat(img, sensor):

    if img.RasterCount > 1:
        # Panchromatic should only be one band but this way the isPan option can also
        # be used to processed L8 images which are stacked in one file.
        rawdata = np.zeros((img.RasterYSize, img.RasterXSize, img.RasterCount))
        for i in range(img.RasterCount):
            rawdata[:, :, i] = img.GetRasterBand(i+1).ReadAsArray()
    else:
        visnirbands = metamod.l78.get_visnirbands(sensor)
        rawdata = np.zeros((img.RasterYSize, img.RasterXSize, len(visnirbands)))

        # Raw Landsat 8/7 data has each band in a separate image.
        # Therefore first open images with all the required band data.
        imgdir = os.path.dirname(img.GetFileList()[0])
        pattern = os.path.join(imgdir, '*.TIF')
        bandfiles = glob.glob(pattern)
        for bf in sorted(bandfiles):
            try:
                bandstr = re.search('_B(\d+)\.TIF$', os.path.basename(bf)).group(1)
                band = int(bandstr)
            except AttributeError:
                continue
            if band in visnirbands:
                logger.info(band)
                bandimg = gdal.Open(os.path.join(imgdir, bf), gdal.GA_ReadOnly)
                rawdata[:, :, band-1] = bandimg.GetRasterBand(1).ReadAsArray()
                rawdata = np.int_(rawdata)
    return rawdata
