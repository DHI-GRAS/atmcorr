import click
import yamlconfig.click_option


@click.command()
@yamlconfig.click_option.yaml_config_option(
        keys=[
            'sensor', 'dnFile', 'mtdFile', 'band_ids',
            'atm', 'aeroProfile', 'tileSize',
            'outfile',
            ('adjCorr', False),
            ('mtdFile_tile', None),
            ('use_modis', False),
            ('earthdata_credentials', None)],
        allow_missing=False,
        parse_kwargs=dict(
            join_rootdir=True),
        help='Atmospheric correction config file')
def cli(dnFile, outfile, **config):
    """Run atmospheric correction"""
    from atmcorr.logs import set_cli_logger
    set_cli_logger(level='INFO')
    from atmcorr.processing import main
    import rasterio
    click.echo('Reading input data ...')
    with rasterio.open(dnFile) as src:
        data = src.read()
        profile = src.profile
    click.echo('Correcting ...')
    outdata, outprofile = main(data=data, profile=profile, **config)
    click.echo('Writing output data ...')
    with rasterio.open(outfile, 'w', **outprofile) as dst:
        dst.write(outdata)
    click.echo('All done.')
