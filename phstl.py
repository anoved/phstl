#!/usr/bin/python

from math import sqrt
import sys
import argparse
from collections import deque
from struct import pack, unpack

from osgeo import gdal

gdal.UseExceptions()
gdal.TermProgress = gdal.TermProgress_nocb

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

class stlwriter():
	
	# path: output binary stl file path
	# facet_count: predicted number of facets
	def __init__(self, path, facet_count=0):
		
		self.f = open(path, 'wb')
		
		# track number of facets actually written
		self.written = 0
		
		# write binary stl header with predicted facet count
		self.f.write('\0' * 80)
		self.f.write(pack('<I', facet_count))
	
	# t: ((ax, ay, az), (bx, by, bz), (cx, cy, cz))
	def add_facet(self, t):
		self.f.write(pack('<3f', *NormalVector(t)))
		for vertex in t:
			self.f.write(pack('<3f', *vertex))
		self.f.write('\0\0')
		self.written += 1
	
	def done(self):
		# update final facet count in header before closing file
		self.f.seek(80)
		self.f.write(pack('<I', self.written))
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
ap.add_argument('-z', action='store', default=1.0, type=float, help='Vertical scale factor (1)')
ap.add_argument('-b', '--base', action='store', default=0.0, type=float, help='Base height (0)')
ap.add_argument('-c', '--clip', action='store_true', default=False, help='Clip z to minimum elevation')
ap.add_argument('-v', '--verbose', action='store_true', default=False, help='Print log messages')
ap.add_argument('--band', action='store', default=1, type=int, help='Raster data band (1)')
ap.add_argument('-m', '--minimum', action='store', default=None, type=float, help='Omit vertices below minimum elevation')
ap.add_argument('-M', '--maximum', action='store', default=None, type=float, help='Omit vertices above maximum elevation')
ap.add_argument('RASTER', help='Input heightmap image')
ap.add_argument('STL',  help='Output STL path')
args = ap.parse_args()

try:
	img = gdal.Open(args.RASTER)
except RuntimeError, e:
	fail(str(e).strip())

# input raster dimensions
w = img.RasterXSize
h = img.RasterYSize
log("raster dimensions = (%s, %s)" % (str(w), str(h)))

# output mesh dimensions are one row and column less than raster
mw = w - 1
mh = h - 1

# get default transformation from image coordinates to world coordinates
t = img.GetGeoTransform()

git = gdal.InvGeoTransform(t)[1]
log("geo->pixel transform: %s" % str(git))

# save x pixel size if needed for scaling
xyres = abs(t[1])

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
	t = (
			-pixel_scale * mw / 2.0, # 0 left edge of top left pixel
			 pixel_scale,            # 1 pixel width
			 0,                      # 2
			 pixel_scale * mh / 2.0, # 3 top edge of top left pixel
			 0,                      # 4 
			-pixel_scale             # 5 pixel height
	)

log("transform = %s" % str(t))

band = img.GetRasterBand(args.band)
nd = band.GetNoDataValue()

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

typeName = gdal.GetDataTypeName(band.DataType)
if band.DataType not in typemap:
	fail('Unsupported data type: %s' % typeName)

# rowformat is used to unpack a row of raw image data to numeric form
rowformat = typemap.get(band.DataType) * w
log("data type = %s" % typeName)
log("type format = %s" % typemap.get(band.DataType))

# min, max, mean, sd; min used for z clipping
stats = band.GetStatistics(True, True)
log("min, max, mean, sd = %s" % str(stats))

# zmin is subtracted from elevation values
if args.clip == True:
	zmin = stats[0]
else:
	zmin = 0
log("zmin = %s" % str(zmin))

# Rolling pixel buffer has space for two rows of image data.
# Old data is automatically discarded as new data is loaded.
pixels = deque(maxlen = (2 * w))
log("buffer size = %s" % str(pixels.maxlen))

# Initialize pixel buffer with first row of image data.
pixels.extend(unpack(rowformat, band.ReadRaster(0, 0, w, 1, w, 1, band.DataType)))

# precalculate output mesh size (STL is 50 bytes/facet + 84 byte header)
facetcount = mw * mh * 2
filesize = (facetcount * 50) + 84
log("predicted (max) facet count = %s" % str(facetcount))
log("predicted (max) STL file size = %s bytes" % str(filesize))

def skip(v):
	global nd
	global args
	if v == nd:
		return True
	if args.minimum != None and v < args.minimum:
		return True
	if args.maximum != None and v > args.maximum:
		return True
	return False

with stlwriter(args.STL, facetcount) as mesh:

	for y in range(mh):
		
		# Each row, extend pixel buffer with the next row of image data.
		pixels.extend(unpack(rowformat, band.ReadRaster(0, y + 1, w, 1, w, 1, band.DataType)))
		
		for x in range(mw):
			
			av = pixels[x]
			bv = pixels[w + x]
			cv = pixels[x + 1]
			dv = pixels[w + x + 1]
			
			# Apply transforms to obtain output mesh coordinates of the
			# four corners composed of raster points a (x, y), b, c,
			# and d (x + 1, y + 1):
			#
			# a-c   a-c     c
			# |/| = |/  +  /|
			# b-d   b     b-d
			
			# Points b and c are required for both facets, so if either
			# are unavailable, we can skip this pixel altogether.
			if skip(bv) or skip(cv):
				continue
			
			b = (
				t[0] + (x * t[1]) + ((y + 1) * t[2]),
				t[3] + (x * t[4]) + ((y + 1) * t[5]),
				(zscale * (float(bv) - zmin)) + args.base
			)
			
			c = (
				t[0] + ((x + 1) * t[1]) + (y * t[2]),
				t[3] + ((x + 1) * t[4]) + (y * t[5]),
				(zscale * (float(cv) - zmin)) + args.base
			)
			
			if not skip(av):
				a = (
					t[0] + (x * t[1]) + (y * t[2]),
					t[3] + (x * t[4]) + (y * t[5]),
					(zscale * (float(av) - zmin)) + args.base
				)
				mesh.add_facet((a, b, c))
			
			if not skip(dv):
				d = (
					t[0] + ((x + 1) * t[1]) + ((y + 1) * t[2]),
					t[3] + ((x + 1) * t[4]) + ((y + 1) * t[5]),
					(zscale * (float(dv) - zmin)) + args.base
				)
				mesh.add_facet((d, c, b))
		gdal.TermProgress(float(y + 1) / mh)

log("actual facet count: %s" % str(mesh.written))
