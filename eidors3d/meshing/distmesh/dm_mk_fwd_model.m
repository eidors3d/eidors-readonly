function [fwd_mdl]= dm_mk_fwd_model( fd, h0, bbox, ...
                     elec_nodes, stim_pattern, z_contact, name)
% DM_MK_FWD_MODEL: create a fwd_model object using distmesh
% fwd_mdl= dm_mk_fwd_model( fd, fh, h0, bbox,...
%                          elec_nodes, stim_pattern, z_contact, name);
%
%  FD:        Distance function - 1 for space inside area, 0 outside
%  H0:        Initial edge length
%  BBOX:      Bounding box [xmin,ymin, {zmin}; xmax,ymax,{zmax}]
%
%  stim_pattern:      a stimulation pattern structure
%                     empty ([]) if stim_pattern is not available
%  elec_nodes:        cell of matrix N x [x,y,{z}] for each electrode
%  z_contact:         vector or scalar electrode contact impedance
%  name:              name for eidors object
%
%  fwd_mdl:           eidors format fwd_model

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id: dm_mk_fwd_model.m,v 1.4 2008-03-10 18:12:35 aadler Exp $

if nargin <8
   name = 'MDL from dm_mk_fwd_model';
end
fwd_mdl= create_refined_model(name, fd, h0, bbox, elec_nodes, ...
                        stim_pattern, z_contact); 

function fmdl= create_refined_model(name, fd, h0, bbox, elec_nodes, ...
                        stim_pattern, z_contact); 
   % fixed_node= [elec_nodes{:}]; - wish we could do it like this - matlab bug!!
   fixed_node= [];
   for i= 1:prod(size(elec_nodes)) 
      this_elec= [fixed_node; elec_nodes{i}];
      fixed_node= [fixed_node; elec_nodes{i}];
   end
   [vtx,simp] = call_distmesh(fd,h0,bbox,fixed_node);
   srf= find_boundary(simp);

   fmdl = construct_fwd_model(srf,vtx,simp, name, ...
                          stim_pattern, elec_nodes, z_contact);

% build fwd_model structure
function mdl= construct_fwd_model(srf,vtx,simp, name, ...
                       stim_pattern, elec_nodes, z_contact)
   mdl= eidors_obj('fwd_model', name);

   mdl.nodes    = vtx;
   mdl.elems    = simp;
   mdl.boundary = srf;
   mdl.gnd_node=           1;
   mdl.np_fwd_solve.perm_sym =     '{n}';
   mdl.name = name;

   % Model Stimulation
   if ~isempty(stim_pattern)
      mdl.stimulation= stim_pattern;
   end

   % Electrodes and z_contact

   % set the z_contact
   n_elec= prod(size(elec_nodes));
   z_contact= z_contact.*ones(n_elec,1);
   curr_e_node=0;
   for i= 1:n_elec
      this_elec= size(elec_nodes{i},1);

      electrodes(i).nodes    = curr_e_node + (1:this_elec);
      electrodes(i).z_contact= z_contact(i);

      curr_e_node= curr_e_node + this_elec;
   end


   mdl.electrode =     electrodes;
   mdl.solve=          'np_fwd_solve';
   mdl.jacobian=       'np_calc_jacobian';
   mdl.system_mat=     'np_calc_system_mat';

function [vtx,simp] = call_distmesh(fd,h0,bbox,fixed_node);
%  FH:        Scaled edge length function - decrease in refined areas
   fh= inline('ones(size(p,1),1)','p'); % initially uniform
   [vtx,simp] = distmeshnd(fd,fh,h0,bbox,fixed_node);
