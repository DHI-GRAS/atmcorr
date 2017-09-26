import logging
import click
import yamlconfig.click_option

logger = logging.getLogger(__name__)


@click.command()
@yamlconfig.click_option.yaml_config_option(
        keys=[
            'sensor', 'dnFile', 'mtdFile',
            'method', 'atm', 'aeroProfile', 'tileSizePixels',
            'outfile',
            ('adjCorr', False),
            ('mtdFile_tile', None),
            ('nprocs', None),
            ('band_ids', None),
            ('use_modis', False),
            ('earthdata_credentials', {})],
        allow_missing=False,
        parse_kwargs=dict(
            join_rootdir=True),
        help='Atmospheric correction config file')
def cli(**config):
    """Run atmospheric correction"""
    from atmcorr.logs import set_cli_logger
    set_cli_logger(level='INFO')
    from atmcorr.processing import main
    main(**config)
