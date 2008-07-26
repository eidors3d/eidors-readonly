function [fwd_mdl]= dm_mk_fwd_model( fd, fh, nnodes, bbox, elec_nodes, ...
                     refine_nodes, z_contact, name)
% DM_MK_FWD_MODEL: create a fwd_model object using distmesh
% fwd_mdl= dm_mk_fwd_model( fd, fh, h0, bbox,...
%                          elec_nodes, stim_pattern, z_contact, name);
%
%  fd:        Shape function - +ve for space inside area, -ve outside
%  fh:        Mesh refinement function - use [] for uniform
%  nnodes:    Estimate of number of nodes in model
%  bbox:      Bounding box [xmin,ymin, {zmin}; xmax,ymax,{zmax}]
%
%  elec_nodes:        cell of matrix N x [x,y,{z}] for each electrode
%  refine_nodes:      vector of fixed nodes to add to model to refine it
%  z_contact:         vector or scalar electrode contact impedance
%  name:              name for eidors object
%
%  fwd_mdl:           eidors format fwd_model

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id$

if nargin <7
   error('dm_mk_fwd_model requires 7 or 8 parameters');
else
   name = 'MDL from dm_mk_fwd_model';
end

if isempty(fh); fh= @huniform; end

h0= estimate_h0(bbox, nnodes);
fwd_mdl= create_refined_model(name, fd, fh, h0, bbox, elec_nodes, ...
                        refine_nodes, z_contact); 

% estimate initial edge length to get nnodes
function  h0= estimate_h0(bbox, nnodes);
   dims= size(bbox,2);
   area_est= prod(abs(diff(bbox,[],1)));
   h0 = (area_est/nnodes)^(1/dims);

function fmdl= create_refined_model(name, fd, fh, h0, bbox, elec_nodes, ...
                        refine_nodes, z_contact); 

   % fixed_node= [elec_nodes{:}]; - wish we could do it like this - matlab bug!!
   fixed_node= [];
   for i= 1:prod(size(elec_nodes)) 
      fixed_node= [fixed_node; elec_nodes{i}];
   end
   n_elec_nodes= size(fixed_node,1);
   fixed_node= [fixed_node; refine_nodes];
   [vtx,simp,srf] = call_distmesh(fd,fh, h0,bbox,fixed_node, n_elec_nodes);

   fmdl = construct_fwd_model(srf,vtx,simp, name, ...
                          elec_nodes, z_contact);


% build fwd_model structure
function mdl= construct_fwd_model(srf,vtx,simp, name, ...
                       elec_nodes, z_contact)
   mdl= eidors_obj('fwd_model', name);

   mdl.nodes    = vtx;
   mdl.elems    = simp;
   mdl.boundary = srf;
   mdl.gnd_node=           1;
   mdl.np_fwd_solve.perm_sym =     '{n}';
   mdl.name = name;

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

function [vtx,simp,srf] = call_distmesh(fd,fh,h0,bbox, ...
                                 fixed_node, n_elec_nodes);
   [vtx,simp] = distmeshnd(fd,fh,h0,bbox,fixed_node);
   srf= find_boundary(simp);
   % Test if distmesh puts extra unneeded nodes on boundary
   rmnode= [];
   % Get srf_nodes which aren't in electrodes
   srf_nodes= unique( srf(:));
   srf_nodes= srf_nodes( srf_nodes > n_elec_nodes);
   for nd = srf_nodes(:)'
       ff= find( any(srf == nd, 2) );
       this_bdy= srf(ff,:);
       this_bdy= unique(this_bdy(:));
       this_bdy( find(this_bdy==nd) ) = [];
       if all( this_bdy < n_elec_nodes);
           rmnode= [rmnode, nd];
%          disp([nd, this_bdy(:)']);
       end
   end

   if ~isempty(rmnode)
      vtx(rmnode,:) = [];
      simp = delaunayn( vtx);
      srf= find_boundary(simp);
   end

function h= huniform(p);
   h= ones(size(p,1),1);
