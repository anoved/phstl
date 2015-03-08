# phstl

Convert [GDAL](http://www.gdal.org/) rasters to partial STL mesh. Intended to produce landscape models from GeoTIFF heightmaps. Rewrite of [hmstl](https://github.com/anoved/hmstl).

## Usage

Currently:

`./phstl.py /path/to/heightmap.tif /path/to/output/mesh.stl`

Only the upper surface of the landscape is generated, so the mesh is not manifold.


## Prerequisites

- [GDAL/OGR in Python](http://trac.osgeo.org/gdal/wiki/GdalOgrInPython) (`sudo easy_install gdal`)
- [Python-STL](https://github.com/apparentlymart/python-stl) (`sudo pip install stl`)


