function paint_electrodes(sel,srf,vtx);
%function paint_electrodes(sel,srf,vtx);
%
%Auxilary function which plots the electrodes red at the boundaries.
%
%
%
% sel = The index of the electrode faces in the srf matrix
%       sel can be created by set_electrodes.m 
% srf = the boundary faces (triangles)
% vtx = The vertices matrix.


l = srf(sel,1); m = srf(sel,2); n = srf(sel,3);

Xs = [vtx(l,1);vtx(m,1);vtx(n,1)];
Ys = [vtx(l,2);vtx(m,2);vtx(n,2)];
Zs = [vtx(l,3);vtx(m,3);vtx(n,3)];

patch(Xs,Ys,Zs,'r');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%