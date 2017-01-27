from setuptools import setup

setup(
    name='atmospheric_correction',
    version='0.1',
    description='Atmospheric Correction using 6S',
    author='Jonas Solvsteen',
    author_email='josl@dhi-gras.com',
    packages=['atmospheric_correction'],
    package_dir={'atmospheric_correction': ''},
    include_package_data=True,
    install_requires=['gdal_utils'],
    dependency_links=['https://github.com/DHI-GRAS/gdal_utils.git'])
