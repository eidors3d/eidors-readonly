function [Reg] = iso_s_smooth(simp,vtx,deg,w);
%function [Reg] = iso_s_smooth(simp,vtx,deg,w);
%
%Calculates a second order discrete Gaussian smoothing operator of
%degree = 1. See help iso_f_smooth for more details.
%
%
%
%simp = The simplices matrix.
%vtx  = The vertices matrix.
%deg  = 1 for nodes, 2 for edges and 3 for faces
%w    = smoothing weight, w=1...k, default value = 1
%Reg  = The second order smoothing regulariser.


warning('EIDORS:deprecated','ISO_S_SMOOTH is deprecated as of 06-Jun-2012. ');

if nargin<2
   w=1;
end

if w<0
   error('Weight must be possitive');
end

[R_first] = iso_f_smooth(simp,vtx,deg,w);

Reg = R_first.'*R_first;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
