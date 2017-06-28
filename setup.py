from setuptools import setup, find_packages

setup(
    name='atmospheric_correction',
    version='0.4',
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
        'gdal>=2.1.0',
        'numpy>=1.11.1',
        'scipy>=0.17.1',
        'gdal_utils>=0.2'],
    dependency_links=[
        'https://github.com/DHI-GRAS/gdal_utils/archive/v0.2.tar.gz#egg=gdal_utils-0.2'])
