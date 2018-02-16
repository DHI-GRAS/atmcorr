from setuptools import setup, find_packages
import versioneer

setup(
    name='atmcorr',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='Atmospheric Correction using 6S',
    author='Jonas Solvsteen',
    author_email='josl@dhigroup.com',
    packages=find_packages(),
    entry_points="""
    [console_scripts]
    atmcorr=atmcorr.scripts.atmcorr:cli
    """,
    install_requires=[
        'dg_calibration',
        'sensor_response_curves',
        'satmeta'])
