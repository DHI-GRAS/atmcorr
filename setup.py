from setuptools import setup, find_packages

setup(
    name='atmcorr',
    version='0.13',
    description='Atmospheric Correction using 6S',
    author='Jonas Solvsteen',
    author_email='josl@dhi-gras.com',
    packages=find_packages(),
    include_package_data=True,
    entry_points="""
    [console_scripts]
    atmcorr=atmcorr.scripts.atmcorr:cli
    """,
    install_requires=[
        'dg_calibration',
        'sensor_response_curves>=0.4',
        'satmeta'],
    dependency_links=[
        'https://github.com/DHI-GRAS/dg-calibration/archive/master.zip',
        'https://github.com/DHI-GRAS/sensor_response_curves/archive/v0.4.tar.gz#egg=sensor_response_curves-0.4',
        'https://github.com/DHI-GRAS/satmeta/archive/v0.11.tar.gz#egg=satmeta-0.11'])
