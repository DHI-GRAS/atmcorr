import click

import yamlconfig.click_option


@click.command()
@yamlconfig.click_option.yaml_config_option(
        keys=[
            'sensor', 'dnFile', 'metadataFile',
            'atmCorrMethod', 'atm', 'isPan', 'adjCorr',
            'aeroProfile', 'tileSizePixels'],
        help='Atmospheric correction config file')
def cli(**config):
    """Run atmospheric correction"""
    from atmospheric_correction.processing import main
    main(**config)
