from setuptools import setup, find_packages

setup(
    name='atmospheric_correction',
    version='0.6.1',
    description='Atmospheric Correction using 6S',
    author='Jonas Solvsteen',
    author_email='josl@dhi-gras.com',
    packages=find_packages(),
    include_package_data=True,
    entry_points="""
    [console_scripts]
    atmcorr=atmospheric_correction.scripts.atmcorr:cli
    """,
    install_requires=[
        'tqdm',
        'rasterio',
        'numpy>=1.11.1',
        'scipy>=0.17.1',
        'sensor_response_curves>=0.4'],
    dependency_links=[
        'https://github.com/DHI-GRAS/sensor_response_curves/archive/v0.4.tar.gz#egg=sensor_response_curves-0.4'])
