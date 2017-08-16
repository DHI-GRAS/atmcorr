import logging
import click
import yamlconfig.click_option

logger = logging.getLogger(__name__)


@click.command()
@yamlconfig.click_option.yaml_config_option(
        keys=[
            'sensor', 'dnFile', 'mtdFile', 'roiFile',
            'method', 'atm', 'aeroProfile', 'tileSizePixels',
            'outfile',
            ('isPan', False),
            ('adjCorr', False),
            ('mtdFile_tile', None),
            ('nprocs', None),
            ('band_ids', None)],
        allow_missing=False,
        parse_kwargs=dict(
            join_rootdir=True),
        help='Atmospheric correction config file')
def cli(outfile, **config):
    """Run atmospheric correction"""
    from atmospheric_correction.logs import set_cli_logger
    set_cli_logger(level='INFO')
    from atmospheric_correction.processing import main
    from gdal_utils import gdal_utils as gu
    img = main(**config)
    logger.info('Saving to %s', outfile)
    gu.dump_gtiff(img, outfile)
