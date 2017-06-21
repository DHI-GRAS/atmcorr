from setuptools import setup

setup(
    name='atmospheric_correction',
    version='0.2',
    description='Atmospheric Correction using 6S',
    author='Jonas Solvsteen',
    author_email='josl@dhi-gras.com',
    find_packages=True,
    include_package_data=True,
    install_requires=[
        'tqdm',
        'gdal>=2.1.0',
        'numpy>=1.11.1',
        'scipy>=0.17.1',
        'gdal_utils==0.2'],
    dependency_links=[
        'https://github.com/DHI-GRAS/gdal_utils/archive/v0.2.tar.gz#egg=gdal_utils-0.2'])
