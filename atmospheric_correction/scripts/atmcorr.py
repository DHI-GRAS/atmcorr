import click

import yamlconfig.click_option


@click.command()
@yamlconfig.click_option.yaml_config_option(
        keys=[
            'sensor', 'dnFile', 'mtdFile', 'roiFile',
            'method', 'atm', 'aeroProfile', 'tileSizePixels',
            ('isPan', False),
            ('adjCorr', False),
            ('mtdFile_tile', None),
            ('nprocs', None),
            ('band_ids', None)],
        allow_missing=False,
        parse_kwargs=dict(
            join_rootdir=True),
        help='Atmospheric correction config file')
def cli(**config):
    """Run atmospheric correction"""
    from atmospheric_correction.logs import set_cli_logger
    set_cli_logger(level='INFO')
    from atmospheric_correction.processing import main
    main(**config)
