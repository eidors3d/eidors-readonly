The files in this directory are provided by Marko Vauhkonen and
correspond to the EIDORS 2D toolkit in the paper:

    M Vauhkonen, W R B Lionheart, L M Heikkinen, P J
    Vauhkonen and J P Kaipio, (2001) "A MATLAB package for
    the EIDORS project to reconstruct two-dimensional EIT
    images"  Physiol. Meas. 22 107-111


Files were downloaded from Marko's web site on 14 Sept 2005,
and chequed into the EIDORS3D cvs.  http://venda.uku.fi/~mvauhkon/

The bug fix noted at http://venda.uku.fi/~mvauhkon/MatlabEIT2d.html
was implemented.

All files are licensed under the GPL as per the LICENSE.TXT file
provided with the code.

The remainder of this readme is the original file provided with EIDORS2D
--- Andy Adler
-------------------------------------------------------


		2D EIT Package for MATLAB
		=========================
		     Version 1.0
	             ===========

Installation
============

We suppose here that <PATH> stands for the path of the directory containing this README file.

- On Unix/Linux systems
  User should include the 2D EIT Package directories in his MATLAB path. This can be done in the
  startup file for automatic path setting or by using addpath in MATLAB in the following way

 >> addpath <PATH>/data <PATH>/demo <PATH>/forwardsol <PATH>/graphics <PATH>/inversesol <PATH>/mesh -begin                        


- On Windows systems
  User should include the 2D EIT Package directories in his MATLAB path. This can be done in the
  startup file for automatic path setting or by using addpath in MATLAB in the following way

 >> addpath <PATH>\data <PATH>\demo <PATH>\forwardsol <PATH>\graphics <PATH>\inversesol <PATH>\mesh -begin
 
In some m-files you will need the QMG mesh generator. This can be downloaded from the Vavasis' home pages
at http://www.cs.cornell.edu/Info/People/vavasis/qmg-home.html. However, it is not ncessary to have the QMG mesh
generator. Use , for example, the cirgrid_eit.m for simple circular meshes. See MakeMeshData.m for this.  

Contents
========

README.TXT	- This file.
LICENSE.TXT     - License text.
data		- Directory of measurement and simulation data.
demo		- Directory of demonstrations. 
forwardsol	- Directory of functions for solving the EIT forward problem in 2D.
graphics	- Directory of functions for plotting the results.
inversesol	- Directory of functions for solving the EIT inverse problem in 2D.
mesh		- Directory of functions for mesh handling.

Documentation
=============

The only documentation so far are the MATLAB help files and descriptions of the functions at the beginning of each function. 
The demonstration examples show how to use the Package.


License
=======

All the software within this package is covered by the license which is stated in the file "LICENSE.TXT".

Contact
=======

This Package has been developed in the University of Kuopio, Department of Applied Physics for EIT research work.
This software  project is a part of a larger collaboration project called EIDORS (http://www.ma.umist.ac.uk/bl/eidors).
If you have any suggestions for improvements or bug reports please let me know.

Best wishes,

Marko Vauhkonen
Department of Applied Physics
University of Kuopio, Finland
Marko.Vauhkonen@uku.fi 



