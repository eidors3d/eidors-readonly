		EIDORS 3D - Matlab base toolkit
		===============================
		  Version 2.01 - April 2003

This is a MATLAB toolkit for solving the forward and inverse problem
for Electrical Impedance Tomography in three dimensions. Of course
most EIT problems are three dimensional so it is long overdue that
people working in medical and industrial EIT collect 3D data and do 3D
forward modelling and reconstruction by default. We hope this
contribution will stimulate activity in 3D EIT research.


Contents
========
The EIDORS_3D package includes:

contents.m
fem_master_full.m
bld_master_full.m
bld_master.m
ref_master.m
set_3d_currents.m
set_multi_currents.m
forward_solver.m
get_3d_meas.m
get_multi_meas.m
m_3d_fields.m
figaro3d.m
slicer_plot.m
slicer_plot_n.m
solution_ext.m
jacobian_3d.m
adjoint_spin.m
integrofgrad.m
set_electrodes.m
set_inho.m
repaint_inho.m
laserbeam.m
paint_electrodes.m
potplot.m
db23d.m
find_boundary.m
iso_f_smooth.m
iso_s_smooth.m 
inverse_solver.m
demo_real.m
demo_complex.m
check_vols.m
delfix.m
triarea3d.m


datareal.mat
datacom.mat

mst-paper.pdf
readme.txt
gpl.html


Documentation
=============

At present the only documentation is a pre-print of a paper we
submitted to Measurement Science and Technology, "mst-paper.pdf"
accompaining this software and the MATLAB help files and descriptions
of the functions at the beginning of each function. There are also two
demos which illustrate how the package is used. They are imaginatively
entitled demo_real and demo_complex, run them as matlab scripts.

Licensing
=========

All the software within this package is covered by the Gnu General
Public License which is stated in the file "gpl.html". We would
appreciate if you publish results using this software that you refer
to our publications, in partuicular the MST paper included. Also if
you improve the software and fix bugs we would like to hear about it,
and we would like to include your fixes in future releases. We would
also welcome contributed demo-programs and test data.

Please remember that this is only a toolkit, and that the
demonstration programs are only simple worked examples. They are
not intended to be the state of the art in EIT reconstruction. By all
means publish your improved reconstruction algorithms and compare them
with the simple ones we demonstrate here. But please don't just run
the demo code on your own data, using the tiny demo mesh, a randomly
chosen regularization parameter and only one iteration of regularized
Newton's method then compare this with your own highly optimised code
and claim that your code is vastly superior than EIDORS' best shot(I
exagerate slightly but there have been such studies with EIDORS 2D).
This is a toolkit of routines for you to build a code which suits your
application.

Contact Adresses
================

This software  project is a part of a larger collaboration project called 
EIDORS (http://www.eidors.org). 

Please email any suggestions for improvements or bug reports at:

Bill.Lionheart@umist.ac.uk
Nicholas.Polydorides@umist.ac.uk

Have fun,


Nick Polydorides and Bill Lionheart
