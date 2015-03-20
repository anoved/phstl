#!/usr/bin/python

# review cols width -1 stuff... I think it stretches it a wee bit much

from math import sqrt
import sys
import argparse

import gdal
import stl

gdal.UseExceptions()

def fail(msg):
	print >> sys.stderr, msg
	exit(1)

ap = argparse.ArgumentParser(description='Convert a GDAL raster (like a GeoTIFF heightmap) to an STL terrain surface.')
ap.add_argument('-x', action='store', default=0.0, type=float, help='Fit output x to extent (mm)')
ap.add_argument('-y', action='store', default=0.0, type=float, help='Fit output y to extent (mm)')
ap.add_argument('-z', action='store', default=1.0, type=float, help='Vertical scale factor')
ap.add_argument('-b', '--base', action='store', default=0.0, type=float, help='Base height')
ap.add_argument('-m', '--mode', action='store', default='surface', choices=['surface', 'solid', 'box'], help='Model mode. Base height must be >0 for solid or box.')
ap.add_argument('-c', '--clip', action='store_true', default=False, help='Clip z to minimum elevation')
ap.add_argument('RASTER', help='Input heightmap image')
ap.add_argument('STL', help='Output terrain mesh')
args = ap.parse_args()

if args.mode != 'surface' and args.base == 0:
	fail("Nonzero base height required for selected mode.")

#
#
#
def mapx(r, c):
	return transform[0] + (c * transform[1]) + (r * transform[2])

#
#
#
def mapy(r, c):
	return transform[3] + (c * transform[4]) + (r * transform[5])

#
# NormalVector
#
# Parameters:
#  tuple of triangle vertices (nested x y z tuples)
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
# ZElevation
#
# Convert an elevation value to output Z units. Applies clipping to minimum
# elevation (zmin, if nonzero); z scaling (and exaggeration); and base offset.
#
# Parameters:
#  elevation
#
# Returns;
#  z value
#
def ZElevation(e):
	return (zscale * (float(e) - zmin)) + args.base

#
# AddQuad
#
# a-b
# |/|
# c-d
#
# Parameters:
#  m: the stl mesh to add the quad to
#  a, b, c, d: vertices (x y z tuples)
#
# Results:
#  adds two triangle facets to mesh
#
def AddQuad(m, a, b, c, d):
	t1 = (a, b, c)
	t2 = (d, c, b)
	m.add_facet(NormalVector(t1), t1)
	m.add_facet(NormalVector(t2), t2)

try:
	img = gdal.Open(args.RASTER)
except RuntimeError, e:
	fail(str(e).strip())

cols = img.RasterXSize
rows = img.RasterYSize

transform = img.GetGeoTransform()
xyres = transform[1]
zscale = args.z

if args.x != 0.0 or args.y != 0.0:
	
	if args.x != 0.0 and args.y != 0.0:
		pixel_scale = min(args.x / (cols - 1), args.y / (rows - 1))
	elif args.x != 0.0:
		pixel_scale = args.x / (cols - 1)
	elif args.y != 0.0:
		pixel_scale = args.y / (rows - 1)
	
	zscale *= pixel_scale / xyres
	transform = (
			-pixel_scale * (cols - 1) / 2.0,
			pixel_scale,
			0,
			pixel_scale * (rows - 1) / 2.0,
			0,
			-pixel_scale
	)

print transform

band = img.GetRasterBand(1)

# min, max, mean, sd
# may use data min for z clipping
stats = band.GetStatistics(True, True)
if args.clip == True:
	zmin = stats[0]
else:
	zmin = 0

# use for masking (omit pixels with this value)
nodata = band.GetNoDataValue()

data = band.ReadAsArray()

mesh = stl.Solid(name="Surface")

for col in range(cols - 1):
	for row in range(rows - 1):

		ax = mapx(row, col)
		ay = mapy(row, col)
		ae = data[row, col]
		az = ZElevation(ae)

		bx = mapx(row + 1, col)
		by = mapy(row + 1, col)
		be = data[row + 1, col]
		if be == nodata:
			continue
		bz = ZElevation(be)

		cx = mapx(row, col + 1)
		cy = mapy(row, col + 1)
		ce = data[row, col + 1]
		if ce == nodata:
			continue
		cz = ZElevation(ce)

		dx = mapx(row + 1, col + 1)
		dy = mapy(row + 1, col + 1)
		de = data[row + 1, col + 1]
		dz = ZElevation(de)

		t1 = ((ax, ay, az), (bx, by, bz), (cx, cy, cz))
		t2 = ((dx, dy, dz), (cx, cy, cz), (bx, by, bz))

		if ae != nodata:
			mesh.add_facet(NormalVector(t1), t1)
		
		if de != nodata:
			mesh.add_facet(NormalVector(t2), t2)
		
		if args.mode != 'surface':
			
			faz = 0
			fbz = 0
			fcz = 0
			fdz = 0
			
			if args.mode == 'solid':
				faz = az - args.base
				fbz = bz - args.base
				fcz = cz - args.base
				fdz = dz - args.base
			
			# left wall
			if col == 0:
				AddQuad(mesh, (ax, ay, az), (ax, ay, faz), (bx, by, bz), (bx, by, fbz))
			
			# right wall
			if col == cols - 2:
				AddQuad(mesh, (dx, dy, dz), (dx, dy, fdz), (cx, cy, cz), (cx, cy, fcz))
			
			# top wall
			if row == 0:
				AddQuad(mesh, (cx, cy, cz), (cx, cy, fcz), (ax, ay, az), (ax, ay, faz))
			
			# bottom wall
			if row == rows - 2:
				AddQuad(mesh, (bx, by, bz), (bx, by, fbz), (dx, dy, dz), (dx, dy, fdz))
			
			# floor
			AddQuad(mesh, (ax, ay, faz), (cx, cy, fcz), (bx, by, fbz), (dx, dy, fdz))

stl = open(args.STL, 'w')
mesh.write_binary(stl)
stl.close()
