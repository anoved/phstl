#!/usr/bin/python

from math import sqrt
import sys
import argparse
from collections import deque
from struct import pack, unpack

import gdal

gdal.UseExceptions()

class stlwriter():
	
	# facet_count: predicted number of facets
	# path: output binary stl file path (default: stdout)
	def __init__(self, facet_count, path=None):
		
		if path == None:
			self.f = sys.stdout
		else:
			self.f = open(path, 'w')
		
		# for future use: track number of facets actually written
		self.written = 0
		
		# write binary stl header with predicted facet count
		self.f.write('\0' * 80)
		self.f.write(pack('<I', facet_count))
	
	# n: normal vector to (a, b, c) facet
	# t: ((ax, ay, az), (bx, by, bz), (cx, cy, cz))
	def add_facet(self, n, t):
		self.f.write(pack('<3f', *n))
		for vertex in t:
			self.f.write(pack('<3f', *vertex))
		self.f.write('\0\0')
		self.written += 1
	
	def done(self):
		if self.f != sys.stdout:
			self.f.close()
	
	def __enter__(self):
		return self
	
	def __exit__(self, exc_type, exc_value, traceback):
		self.done()

def fail(msg):
	print >> sys.stderr, msg
	exit(1)

def log(msg):
	if args.verbose:
		print >> sys.stderr, msg

ap = argparse.ArgumentParser(description='Convert a GDAL raster (like a GeoTIFF heightmap) to an STL terrain surface.')
ap.add_argument('-x', action='store', default=0.0, type=float, help='Fit output x to extent (mm)')
ap.add_argument('-y', action='store', default=0.0, type=float, help='Fit output y to extent (mm)')
ap.add_argument('-z', action='store', default=1.0, type=float, help='Vertical scale factor')
ap.add_argument('-b', '--base', action='store', default=0.0, type=float, help='Base height')
ap.add_argument('-c', '--clip', action='store_true', default=False, help='Clip z to minimum elevation')
ap.add_argument('-v', '--verbose', action='store_true', default=False, help='Print log messages')
ap.add_argument('RASTER', help='Input heightmap image')
ap.add_argument('STL', nargs='?', default=None, help='Output STL path (stdout)')
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
# Could be bundled into stlwriter class.
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
	m.add_facet(NormalVector((a, b, c)), (a, b, c))
	m.add_facet(NormalVector((d, c, b)), (d, c, b))

try:
	img = gdal.Open(args.RASTER)
except RuntimeError, e:
	fail(str(e).strip())

# input raster dimensions
w = img.RasterXSize
h = img.RasterYSize

# output mesh dimensions are one row and column less than raster
mw = w - 1
mh = h - 1

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
		pixel_scale = min(args.x / mw, args.y / mh)
	elif args.x != 0.0:
		pixel_scale = args.x / mw
	elif args.y != 0.0:
		pixel_scale = args.y / mh
	
	# adjust z scale to maintain proportions with new xy scale
	zscale *= pixel_scale / xyres
	
	# revise transformation matrix
	# image: 0,0 at top left corner of top left pixel (0.5,0.5 at pixel center)
	transform = (
			-pixel_scale * mw / 2.0, # 0 left edge of top left pixel
			 pixel_scale,            # 1 pixel width
			 0,                      # 2
			 pixel_scale * mh / 2.0, # 3 top edge of top left pixel
			 0,                      # 4 
			-pixel_scale             # 5 pixel height
	)

log(transform)

band = img.GetRasterBand(1)

# map GDAL pixel data type to corresponding struct format character
typemap = {
	gdal.GDT_Byte:    'B',
	gdal.GDT_UInt16:  'H',
	gdal.GDT_Int16:   'h',
	gdal.GDT_UInt32:  'I',
	gdal.GDT_Int32:   'i',
	gdal.GDT_Float32: 'f',
	gdal.GDT_Float64: 'd'
}

if band.DataType not in typemap:
	fail('Unsupported data type: %s' % (gdal.GetDataTypeName(band.DataType)))
	
rowformat = typemap.get(band.DataType) * w

# min, max, mean, sd; min used for z clipping
stats = band.GetStatistics(True, True)
log(stats)
if args.clip == True:
	zmin = stats[0]
else:
	zmin = 0

log('Initiating raster processing...')

# Space for two rows of image data is allocated. Extending the deque
# with a third (new) row of data automatically exposes the first (old).
pixels = deque(maxlen = (2 * w))

# not handling the data type flexibly here. should map raster datatype
# to an appropriate corresponding struct element specificier ('B', 'H', etc)
pixels.extend(unpack(rowformat, band.ReadRaster(0, 0, w, 1, w, 1, band.DataType)))

with stlwriter(mw * mh * 2, args.STL) as mesh:

	for y in range(mh):
		
		pixels.extend(unpack(rowformat, band.ReadRaster(0, y + 1, w, 1, w, 1, band.DataType)))
		
		for x in range(mw):
					
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
