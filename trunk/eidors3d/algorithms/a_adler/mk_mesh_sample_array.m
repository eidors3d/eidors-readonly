function [eptr, pt_coord]= mk_mesh_sample_array( varargin )
% MK_MESH_SAMPLE_ARRAY - create rectangular coordinate map into mesh
% [eptr, pt_coord]= mk_mesh_sample_array( mdl1, mdl2, ... , npoints);
%
% pt_coord - n_points x dims matrix of sample points
% eptr -     element number containing the pt_coord
% mdl1, mdl2 - fwd_model objects
%
% for example: [ep, pt] = mk_coarse_fine_mapping( m1, m2 )
% ep= [x   y   z]  pt= [m1, m2 ]
%      1   1   1         1  2 % pt in el#1 in m1, el#2 in m2
%      1   1   2         2  2 % pt in el#2 in m1, el#2 in m2
%      1   2   1         2  0 % pt in el#2 in m1, pt not in m2
%
% npoints - total number of points for rectangular mapping

% (C) 2007 Andy Adler. License: GPL version 2 or version 3
% $Id: mk_mesh_sample_array.m,v 1.1 2007-09-21 14:02:12 aadler Exp $

[params, mdls] = proc_input (varargin);
