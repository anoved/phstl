#!/usr/bin/python

# todo:
# - validate elevation clip
# - elevation scaling
# - validate model scaling
# - parse cli options

from math import sqrt
import sys

import gdal
import stl

zscale = 0.1
modelw = 75.0
modelh = 90.0

def mapx(r, c):
	return transform[0] + (c * transform[1]) + (r * transform[2])

def mapy(r, c):
	return transform[3] + (c * transform[4]) + (r * transform[5])

# return normal vector of triangle surface
# t = (a, b, c) where a, b, and c are of form (x, y, z)
def norm(t):
	
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

def e2z(e):
	return zscale * (float(e) - zmin)
	# / 30.0

imgpath = sys.argv[1]

stlpath = sys.argv[2]

mesh = stl.Solid(name="Surface")

img = gdal.Open(imgpath)

cols = img.RasterXSize
rows = img.RasterYSize

# used to convert pixel coordinates to world coordinates 
#transform = img.GetGeoTransform()

# convert pixel coordinates to model coordinates (stretch to fit extent)
t = [0, 0, 0, 0, 0, 0]
t[1] = modelw / (cols)
t[5] = -modelh / (rows)
t[0] = -(modelw / 2) + (t[1] / 2.0)
t[3] = (modelh / 2) - (t[5] / 2.0)
transform = t

band = img.GetRasterBand(1)

# min, max, mean, sd
# may use data min for z clipping
stats = band.GetStatistics(True, True)
zmin = stats[0]
print zmin

# use for masking (omit pixels with this value)
nodata = band.GetNoDataValue()

data = band.ReadAsArray()

for col in range(cols - 1):
	for row in range(rows - 1):

		ax = mapx(row, col)
		ay = mapy(row, col)
		ae = data[row, col]
		az = e2z(ae)

		bx = mapx(row + 1, col)
		by = mapy(row + 1, col)
		be = data[row + 1, col]
		if be == nodata:
			continue
		bz = e2z(be)

		cx = mapx(row, col + 1)
		cy = mapy(row, col + 1)
		ce = data[row, col + 1]
		if ce == nodata:
			continue
		cz = e2z(ce)

		dx = mapx(row + 1, col + 1)
		dy = mapy(row + 1, col + 1)
		de = data[row + 1, col + 1]
		dz = e2z(de)

		t1 = ((ax, ay, az), (bx, by, bz), (cx, cy, cz))
		t2 = ((dx, dy, dz), (cx, cy, cz), (bx, by, bz))

		if ae != nodata:
			mesh.add_facet(norm(t1), t1)
		
		if de != nodata:
			mesh.add_facet(norm(t2), t2)

stl = open(stlpath, 'w')
mesh.write_binary(stl)
stl.close()

