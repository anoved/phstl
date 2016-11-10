# Float Operations

op-original.stl produced with pre-cd91616234c83a7fef30885c08b1703f66bd5591 phstl (using default 64-bit double NormalVector operations):

`../phstl.py -p 0 0 2 2 example.tif op-original.stl`

op-meshlab.stl is the result of the *Recompute Face Normals* Meshlab filter applied to the original output:

`meshlabserver -i op-original.stl -o op-meshlab.stl -s recompute-face-normals.mlx`

op-float32.stl produced with post-cd91616234c83a7fef30885c08b1703f66bd5591 phstl (using forced NumPy 32-bit float NormalVector operations):

`../phstl.py -p 0 0 2 2 example.tif op-float32.stl`

Apply hexdiff or inspect with hex editor to note op-float32.stl geometry data identical to op-meshlab.stl's (disregarding 80-byte header).

