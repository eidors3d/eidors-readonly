function voltage_sources(img)
% voltage_sources

% (C) 2013 Andy Adler.  License:  GPL v2 or v3
% $Id:$

if isstr(img) && strcmp(img,'UNIT_TEST'); do_unit_test; return; end

function do_unit_test
%%% EASY TEST CASE 

   fmdl = getfield( mk_common_model('a2c0',8), 'fwd_model');
   fmdl.stimulation(2:end) = [];
   fmdl.stimulation(1).stim_pattern = sparse([1,2],1,[1;-1],8,1);
   fmdl.stimulation(1).meas_pattern = sparse(1,[3,4],[1;-1],1,8);
   img = mk_image(fmdl,1);
   img.fwd_solve.get_all_meas = 1;
   vv = fwd_solve(img);
   imgv = mk_image(fmdl,vv.volt);
   show_fem(imgv);


   fmdl.stimulation(1)
