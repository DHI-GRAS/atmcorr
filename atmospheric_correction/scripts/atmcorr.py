import click

import yamlconfig.click_option


@click.command()
@yamlconfig.click_option.yaml_config_option(
        keys=[
            'sensor', 'dnFile', 'mtdfile',
            'method', 'atm', 'isPan', 'adjCorr',
            'aeroProfile', 'tileSizePixels',
            'isPan', 'adjCorr',
            'tile', 'mtdfile_tile', 'roiFile',
            'nprocs', 'band_ids'],
        allow_missing=True,
        help='Atmospheric correction config file')
def cli(**config):
    """Run atmospheric correction"""
    from atmospheric_correction.processing import main
    main(**config)
