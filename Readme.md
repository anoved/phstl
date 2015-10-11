# phstl

Convert [GDAL](http://www.gdal.org/) rasters to STL mesh. Intended to produce landscape models from GeoTIFF heightmaps.

## Usage

    usage: phstl.py [-h] [-x X] [-y Y] [-z Z] [-b BASE] [-c] [-v] RASTER [STL]
    
    Convert a GDAL raster (like a GeoTIFF heightmap) to an STL terrain surface.
    
    positional arguments:
      RASTER                Input heightmap image
      STL                   Output STL path (stdout)
    
    optional arguments:
      -h, --help            show this help message and exit
      -x X                  Fit output x to extent (mm)
      -y Y                  Fit output y to extent (mm)
      -z Z                  Vertical scale factor (1)
      -b BASE, --base BASE  Base height (0)
      -c, --clip            Clip z to minimum elevation
      -v, --verbose         Print log messages
      --band BAND           Raster data band (1)
      
## Tips

To make printable solids from the output surface, try using the Extrude and Bisect mesh edit tools in Blender.

## Prerequisites

- [GDAL/OGR in Python](http://trac.osgeo.org/gdal/wiki/GdalOgrInPython) (`sudo easy_install gdal`)

## License

This project is released under an open source MIT license.
