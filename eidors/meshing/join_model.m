function [fmdlo]= join_models(fmdl1, fmdl2, tol)
% JOIN_MODELS: Join two fmdl structures to create one
%
% [fmdlo]= crop_model(fmdl1, fmdl2, tol)
% fmdlo is the union of 

% (C) 2018 Andy Adler. License: GPL version 2 or version 3
% $Id$

if ischar(fmdl1) && strcmp(fmdl1,'UNIT_TEST'); do_unit_test; return; end



function do_unit_test
   imdl = mk_common_model('a2c0',8); fmdl= imdl.fwd_model;
   fmdl = crop_model(fmdl,inline('x<0','x','y','z'));
   unit_test_cmp('crop_model-a2c0-01',length(fmdl.electrode),5);
