from setuptools import setup, find_packages

setup(
    name='atmcorr',
    version='0.14',
    description='Atmospheric Correction using 6S',
    author='Jonas Solvsteen',
    author_email='josl@dhi-gras.com',
    packages=find_packages(),
    entry_points="""
    [console_scripts]
    atmcorr=atmcorr.scripts.atmcorr:cli
    """,
    install_requires=[
        'dg_calibration',
        'sensor_response_curves',
        'satmeta'])
