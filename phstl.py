#!/usr/bin/python

from math import sqrt
import sys
import argparse

import gdal
import stl

ap = argparse.ArgumentParser(description='Placeholder phstl description')
ap.add_argument('--width', action='store', default=0.0, type=float, help='Width of model')
ap.add_argument('--height', action='store', default=0.0, type=float, help='Height of model')
ap.add_argument('infile', nargs='?')
ap.add_argument('outfile', nargs='?')

args = ap.parse_args()

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

# return normal vector of triangle surface
# t = (a, b, c) where a, b, and c are of form (x, y, z)
#
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

#
#
#
def e2z(e):
	return zscale * (float(e) - zmin)

img = gdal.Open(args.infile)
cols = img.RasterXSize
rows = img.RasterYSize

# Excise pixel output mode. It's only here because that's what hmstl did,
# but we're doing this whole rewrite because hmstl doesn't quite do as needed.
# Revise model mode to scale to fit rather than stretch. Then, we'll have a
# single scale value that can be used to adjust the z values proprotionally.
# Allow just one of width or height to be specified, in which case that
# dimension will be used as the scaling target. If both are specified, scale
# to fit (ie, select the smaller of the two implied scale factors).

#~ zscale = 1.0
#~ transform = [0, 0, 0, 0, 0, 0]
#~ 
 #~ args.mode == "model":
	#~ if args.width == 0 or args.height == 0:
		#~ print "must specify model --width and --height as well"
		#~ exit(1)
	#~ 
	#~ # with 75 x 90, getting 74.49 x 89.52 extent
	#~ 
	#~ transform[1] =  (args.width / cols)
	#~ transform[5] = -(args.height / rows)
	#~ transform[0] = -(args.width / 2.0) # + (transform[1] / 2)
	#~ transform[3] =  (args.height / 2.0) # - (transform[5] / 2)
	#~ 
	#~ zscale = (abs(transform[1]) + abs(transform[5])) / 2
	#~ 
#~ elif args.mode == "world":
	#~ # getting flat z
	#~ transform = img.GetGeoTransform()
	#~ # zscale 1 for world coords
	#~ zscale = 1.0

transform = img.GetGeoTransform()
zscale = 1.0

band = img.GetRasterBand(1)

# min, max, mean, sd
# may use data min for z clipping
stats = band.GetStatistics(True, True)
zmin = stats[0]

# use for masking (omit pixels with this value)
nodata = band.GetNoDataValue()

data = band.ReadAsArray()

mesh = stl.Solid(name="Surface")

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


stl = open(args.outfile, 'w')
mesh.write_binary(stl)
stl.close()

