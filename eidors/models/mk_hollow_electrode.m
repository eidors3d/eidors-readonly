function fmdl = mk_hollow_electrode(fmdl, elec_idx)
% MK_HOLLOW_ELECTRODE: remove nodes with indicated electrdoes
% fmdl_out = mk_hollow_electrode(fmdl, elec_idx)
%   fmdl:     input model
%   elec_idx: electrodes for which we remove internal nodes
%             if not provided, remove for all electrodes
%   fmdl_out: output model
%
% This function is useful for internal electrodes created by
%   meshing software in which we don't want to calculate current
%   flow within the electrode

% (C) 2017 Andy Adler. License: GPL v2 or v3. $Id$

if ischar(fmdl) && strcmp(fmdl,'UNIT_TEST'); do_unit_test; return; end

function do_unit_test
    fmdl = unit_test_model1;
    subplot(221); show_fem(fmdl); view(0,60);

function fmdl = unit_test_model1
   R=0.65; Nel = 16;
   shape_str = [ ...
   'solid incyl  = cylinder (0,0,0; 0,0,1; 0.125) -maxh=0.015;\n', ...
   'solid farcyl = cylinder (0,0,0; 0,0,1;0.75) -maxh=0.55;\n' ...
   'solid pl1    =  plane(0,0,-0.09;0,0,-1);\n' ...
   'solid top    =  plane(0,0,0; 0,0,1) -maxh=0.13;\n' ...
   'solid mainobj= pl1 and top and farcyl and not incyl;\n'];
   elec_pos = zeros(Nel,6);
   Theta = 0:360/Nel:360-360/Nel;
   for i = 1:Nel
          elec_pos(i,1) = R*sind(Theta(i));
          elec_pos(i,2) = R*cosd(Theta(i));
          elec_pos(i,6) = 1;
          elec_obj(i) = {'top'};
   end
   elec_shape=[0.025];
   fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);
   fmdl = mdl2d_from3d(fmdl);
