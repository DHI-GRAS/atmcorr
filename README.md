atmcorr
-------

[![Build Status](https://travis-ci.com/DHI-GRAS/atmcorr.svg?token=xFUvaoqNvYGuqLz5TEdJ&branch=master)](https://travis-ci.com/DHI-GRAS/atmcorr)

Wrapper for Py6S with automatic download of atmospheric
composition from MODIS.


## Testing

```
pip install -e .
pytest -v
```

Full testing requires that you install the [`atmcorr_testdata`](https://github.com/DHI-GRAS/atmcorr_testdata) package.


## Travis & Docker

```
docker build -t dhigras/atmcorr-env -f Dockerfile.environment .
docker push dhigras/atmcorr-env
```
