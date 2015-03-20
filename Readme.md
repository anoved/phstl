# phstl

Convert [GDAL](http://www.gdal.org/) rasters to STL mesh. Intended to produce landscape models from GeoTIFF heightmaps. Rewrite of [hmstl](https://github.com/anoved/hmstl).

## Usage

    usage: phstl.py [-h] [-x X] [-y Y] [-z Z] [-b BASE] [-m {surface,solid,box}]
                    [-c]
                    RASTER STL
    
    Convert a GDAL raster (like a GeoTIFF heightmap) to an STL terrain surface.
    
    positional arguments:
      RASTER                Input heightmap image
      STL                   Output terrain mesh
    
    optional arguments:
      -h, --help            show this help message and exit
      -x X                  Fit output x to extent (mm)
      -y Y                  Fit output y to extent (mm)
      -z Z                  Vertical scale factor
      -b BASE, --base BASE  Base height
      -m {surface,solid,box}, --mode {surface,solid,box}
                            Model mode. Base height must be >0 for solid or box.
      -c, --clip            Clip z to minimum elevation

Default `-z` scale factor is `1`. Default `--mode` is `surface`, yielding a partial mesh. Default `--base` height is `0` (must be nonzero for `solid` or `box` mode). Elevation `--clip` is disabled by default.

## Prerequisites

- [GDAL/OGR in Python](http://trac.osgeo.org/gdal/wiki/GdalOgrInPython) (`sudo easy_install gdal`)
- [Python-STL](https://github.com/apparentlymart/python-stl) (`sudo pip install stl`)


