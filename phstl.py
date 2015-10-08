#!/usr/bin/python

from math import sqrt
import sys
import argparse
from collections import deque
from struct import unpack

import gdal
import stl

gdal.UseExceptions()

def fail(msg):
	print >> sys.stderr, msg
	exit(1)

def log(msg):
	if args.verbose:
		print msg

ap = argparse.ArgumentParser(description='Convert a GDAL raster (like a GeoTIFF heightmap) to an STL terrain surface.')
ap.add_argument('-x', action='store', default=0.0, type=float, help='Fit output x to extent (mm)')
ap.add_argument('-y', action='store', default=0.0, type=float, help='Fit output y to extent (mm)')
ap.add_argument('-z', action='store', default=1.0, type=float, help='Vertical scale factor')
ap.add_argument('-b', '--base', action='store', default=0.0, type=float, help='Base height')
ap.add_argument('-c', '--clip', action='store_true', default=False, help='Clip z to minimum elevation')
ap.add_argument('-v', '--verbose', action='store_true', default=False, help='Print log messages')
ap.add_argument('RASTER', help='Input heightmap image')
ap.add_argument('STL', help='Output terrain mesh')
args = ap.parse_args()

#
# XCoordinate
#
# Apply affine transformation to get output X coordinate from image coordinates.
#
# Parameters:
#  r, c: row, column image coordinates
#
# Accesses:
#  transform tuple
#
# Returns:
#  x value
#
def XCoordinate(r, c):
	return transform[0] + (c * transform[1]) + (r * transform[2])

#
# YCoordinate
#
# Apply affine transformation to get output Y coordinate from image coordinates.
#
# Parameters:
#  r, c: row, column image coordinates
#
# Accesses:
#  transform tuple
#
# Returns:
#  y value
#
def YCoordinate(r, c):
	return transform[3] + (c * transform[4]) + (r * transform[5])

#
# ZCoordinate
#
# Convert an elevation value to output Z units. Applies clipping to minimum
# elevation (zmin, if nonzero); z scaling (and exaggeration); and base offset.
#
# Parameters:
#  elevation
#
# Accesses:
#  zscale, zmin, and args.base
#
# Returns:
#  z value
#
def ZCoordinate(e):
	return (zscale * (float(e) - zmin)) + args.base

#
# NormalVector
#
# Calculate the normal vector of a triangle. (Unit vector perpendicular to
# triangle surface, pointing away from the "outer" face of the surface.)
#
# Parameters:
#  triangle vertices (nested x y z tuples)
#
# Returns:
#  normal vector (x y z tuple)
#
def NormalVector(t):
	
	(ax, ay, az) = t[0]
	(bx, by, bz) = t[1]
	(cx, cy, cz) = t[2]
	
	# first edge
	e1x = ax - bx
	e1y = ay - by
	e1z = az - bz
	
	# second edge
	e2x = bx - cx
	e2y = by - cy
	e2z = bz - cz
	
	# cross product
	cpx = e1y * e2z - e1z * e2y
	cpy = e1z * e2x - e1x * e2z
	cpz = e1x * e2y - e1y * e2x
	
	# return cross product vector normalized to unit length
	mag = sqrt((cpx * cpx) + (cpy * cpy) + (cpz * cpz))
	return (cpx/mag, cpy/mag, cpz/mag)

#
# AddTri
#
# a-c
# |/
# b
#
# Parameters:
#  m: the stl mesh to add the tri to
#  t: triangle vertices (a b c as nested x y z tuples)
#
# Results:
#  adds triangle facet to mesh
#
def AddTri(m, t):
	m.add_facet(NormalVector(t), t)

#
# AddQuad
#
# a-c
# |/|
# b-d
#
# Parameters:
#  m: the stl mesh to add the quad to
#  a, b, c, d: vertices (x y z tuples)
#
# Results:
#  adds two triangle facets to mesh
#
def AddQuad(m, a, b, c, d):
	AddTri(m, (a, b, c))
	AddTri(m, (d, c, b))

try:
	img = gdal.Open(args.RASTER)
except RuntimeError, e:
	fail(str(e).strip())

w = img.RasterXSize
h = img.RasterYSize

# get default transformation from image coordinates to world coordinates
transform = img.GetGeoTransform()

# save x pixel size if needed for scaling
xyres = abs(transform[1])

# initialize z scale to exaggeration factor, if any
zscale = args.z

# recalculate z scale and xy transform if different dimensions are requested
if args.x != 0.0 or args.y != 0.0:
	
	# recaculate xy scale based on requested x or y dimension
	# if both x and y dimension are set, select smaller scale
	if args.x != 0.0 and args.y != 0.0:
		pixel_scale = min(args.x / (w - 1), args.y / (h - 1))
	elif args.x != 0.0:
		pixel_scale = args.x / (w - 1)
	elif args.y != 0.0:
		pixel_scale = args.y / (h - 1)
	
	# adjust z scale to maintain proportions with new xy scale
	zscale *= pixel_scale / xyres
	
	# revise transformation matrix
	# image: 0,0 at top left corner of top left pixel (0.5,0.5 at pixel center)
	transform = (
			-pixel_scale * (w - 1) / 2.0, # 0 left edge of top left pixel
			 pixel_scale,                    # 1 pixel width
			 0,                              # 2
			 pixel_scale * (h - 1) / 2.0, # 3 top edge of top left pixel
			 0,                              # 4 
			-pixel_scale                     # 5 pixel height
	)

log(transform)

band = img.GetRasterBand(1)

# min, max, mean, sd; min used for z clipping
stats = band.GetStatistics(True, True)
log(stats)
if args.clip == True:
	zmin = stats[0]
else:
	zmin = 0

log('Initiating raster processing...')

mesh = stl.Solid(name="Surface")

# Space for two rows of image data is allocated. Extending the deque
# with a third (new) row of data automatically exposes the first (old).
pixels = deque(maxlen = (2 * w))

# not handling the data type flexibly here. should map raster datatype
# to an appropriate corresponding struct element specificier ('B', 'H', etc)
pixels.extend(unpack('H' * w, band.ReadRaster(0, 0, w, 1, w, 1, band.DataType)))

for y in range(h - 1):
	
	pixels.extend(unpack('H' * w, band.ReadRaster(0, y + 1, w, 1, w, 1, band.DataType)))
	
	for x in range(w - 1):
				
		ax = XCoordinate(y, x)
		ay = YCoordinate(y, x)
		az = ZCoordinate(pixels[x])
		
		bx = XCoordinate(y + 1, x)
		by = YCoordinate(y + 1, x)
		bz = ZCoordinate(pixels[w + x])
		
		cx = XCoordinate(y, x + 1)
		cy = YCoordinate(y, x + 1)
		cz = ZCoordinate(pixels[0 + x + 1])
		
		dx = XCoordinate(y + 1, x + 1)
		dy = YCoordinate(y + 1, x + 1)
		dz = ZCoordinate(pixels[w + x + 1])
		
		AddQuad(mesh, (ax, ay, az), (bx, by, bz), (cx, cy, cz), (dx, dy, dz))

log('Initiating mesh output...')
stl = open(args.STL, 'w')
mesh.write_binary(stl)
stl.close()
