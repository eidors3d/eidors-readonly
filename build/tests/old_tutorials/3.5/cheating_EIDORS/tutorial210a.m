% Define Face Shapes: Small Face
% $Id$

clear p;
p.leye=   [78,97,98,117,118,141];
p.reye=   [66,82,83,102,103,123];
p.rsmile= [40:41, 53:55, 69:71, 86:88];
p.lsmile= [43:44, 57:59, 73:75, 91:93];
p.lsad =  [31,43,58,57,74,73,92,93,113,112,135];
p.rsad =  [28,40,53,54,69,87,70,107,88,108,129];
p.eyes= [p.leye,p.reye];
p.sad = [p.lsad,p.rsad];
p.smile= [p.rsmile, p.lsmile];
small_face = p;

% Simulate data for small face
imdl= mk_common_model('b2c');
e= size(imdl.fwd_model.elems,1);
simg= eidors_obj('image','small face');
simg.fwd_model= imdl.fwd_model;

% homogeneous data
simg.elem_data= ones(e,1);
vh= fwd_solve( simg);

% inhomogeneous data for sad face model
simg.elem_data(small_face.eyes)= 2;
simg.elem_data(small_face.sad)=  1.5;
vi= fwd_solve( simg);
