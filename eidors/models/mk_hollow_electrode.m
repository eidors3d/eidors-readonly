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
%
% Limitations: currently works for 2D fems only
% Recommended replacement: mat_idx_to_electrode

% (C) 2017 Andy Adler and Bartek Grychtol. License: GPL v2 or v3. $Id$

if ischar(fmdl) && strcmp(fmdl,'UNIT_TEST'); do_unit_test; return; end

if nargin==1; elec_idx = 1:length(fmdl.electrode); end
elec_idx = elec_idx(:)'; % Row vector to index for loop

elidx=[];
for i = elec_idx
   elidx = [elidx, fmdl.electrode(i).nodes];
end

ELNODES = zeros(size(fmdl.nodes,1),1);
ELNODES(elidx) = 1;
% find all electrode elems
ELELEM = sum(ELNODES(fmdl.elems),2)==3;
% remove internal electrode elements
fmdl.elems = fmdl.elems(~ELELEM,:);
% used nodes
unodes = unique(fmdl.elems);
% remove internal electrode nodes
for i = elec_idx
   idx = ismember(fmdl.electrode(i).nodes,unodes);
   fmdl.electrode(i).nodes = fmdl.electrode(i).nodes(idx);
end

nodemap = zeros(size(fmdl.nodes,1),1);
nodemap(unodes) = 1:numel(unodes);
% remove unused nodes
fmdl.nodes = fmdl.nodes(unodes,:);
%remap elems
fmdl.elems = nodemap(fmdl.elems);
%remap elecs
for i = 1:length(fmdl.electrode)
   idx = ismember(fmdl.electrode(i).nodes,unodes);
   nidx =  nodemap( fmdl.electrode(i).nodes(idx) );
   if any(i==elec_idx);
      npts = fmdl.nodes(nidx,:);
      [npts_o, nnidx] = order_loop(npts);
      fmdl.electrode(i).nodes = nidx(nnidx);
   else % if not in specified index, don't reorder loop
      fmdl.electrode(i).nodes = nidx;
   end
end


fmdl.boundary = find_boundary(fmdl);

[~,idx] = min(sum(fmdl.nodes.^2,2));
fmdl.gnd_node = idx(1);


function do_unit_test
    fmdl = unit_test_model1;
    subplot(221); show_fem(fmdl); axis([-.1,.3,0.4,0.8])
    title 'original model - filled electrodes';
    fmdl1= mk_hollow_electrode(fmdl,[1,2,5]);
    
    subplot(222); show_fem(fmdl1); axis([-.1,.3,0.4,0.8])
    title 'original model - hollow electrode #1';

    fmdl2= mk_hollow_electrode(fmdl);
    subplot(223); show_fem(fmdl2); axis([-.1,.3,0.4,0.8])
    title 'original model - all hollow electrode';

    fmdl = unit_test_model2;
    fmdl4= mk_hollow_electrode(fmdl, length(fmdl.electrode));
    subplot(224); show_fem(fmdl4); axis([-1,1,-1,1])
    title 'original model - one electrode';
    

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

function fmdl = unit_test_model2
   extra={'ball', ...
         ['solid cyls= cylinder(0.2,0.2,0;0.2,0.2,1;0.2) -maxh=0.05;' ...
          'solid ball= cyls and orthobrick(-1,-1,0;1,1,0.5);']};
   fmdl= ng_mk_cyl_models(0,[6],[0.1,0,0.05],extra); 
   eln   = find(elem_select(fmdl, '(x-0.2).^2+(y-0.2).^2<0.2^2'));
   eln   = unique(fmdl.elems(eln,:));
   fmdl.electrode(end+1) = struct('nodes', eln(:)', 'z_contact', .01);
