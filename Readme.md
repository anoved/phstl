# phstl

Convert [GDAL](http://www.gdal.org/) rasters to partial STL mesh. Intended to produce landscape models from GeoTIFF heightmaps. Rewrite of [hmstl](https://github.com/anoved/hmstl).

## Usage

    usage: phstl.py [-h] [-x X] [-y Y] [-z Z] RASTER STL
    
    Convert a GDAL raster (like a GeoTIFF heightmap) to an STL terrain surface.
    
    positional arguments:
      RASTER      Input heightmap image
      STL         Output terrain mesh
    
    optional arguments:
      -h, --help  show this help message and exit
      -x X        Scale output to fit x extent
      -y Y        Scale output to fit y extent
      -z Z        Vertical scale factor

Only the upper surface of the landscape is generated, so the mesh is not manifold.


## Prerequisites

- [GDAL/OGR in Python](http://trac.osgeo.org/gdal/wiki/GdalOgrInPython) (`sudo easy_install gdal`)
- [Python-STL](https://github.com/apparentlymart/python-stl) (`sudo pip install stl`)


