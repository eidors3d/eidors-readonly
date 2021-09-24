all: box-4.1.msh box-4.0.msh box-2.2.msh \
	 cube-4.1.msh cube-4.0.msh cube-2.2.msh

box-4.1.msh: box.geo
	gmsh -2 -format msh41 -clmin 1 -o $@ $<
box-4.0.msh: box.geo
	gmsh -2 -format msh40 -clmin 1 -o $@ $<
box-2.2.msh: box.geo
	gmsh -2 -format msh22 -clmin 1 -o $@ $<

cube-4.1.msh: cube.geo
	gmsh -3 -format msh41 -clmin 1 -o $@ $<
cube-4.0.msh: cube.geo
	gmsh -3 -format msh40 -clmin 1 -o $@ $<
cube-2.2.msh: cube.geo
	gmsh -3 -format msh22 -clmin 1 -o $@ $<
