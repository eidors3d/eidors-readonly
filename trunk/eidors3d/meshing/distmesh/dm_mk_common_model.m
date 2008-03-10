function [fwd_mdl]= dm_mk_common_model( mdl_name, nnodes, elec_pattern, stim_pattern );
% DM_MK_COMMON_MODEL: create a prepared model with distmesh
% [fwd_mdl]= dm_mk_common_model( mdl_name, nnodes, elec_pattern, stim_pattern );
%
% mdl_name:     'c2'   - 2D circle
% nnodes:       estimate number of nodes in model
% elec_pattern:     [16, 1]  - 16 elecs per ring, one ring
% 
%  fwd_mdl:           eidors format fwd_model

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id: dm_mk_common_model.m,v 1.4 2008-03-10 18:12:35 aadler Exp $


if mdl_name == 'c2'
   fd=inline('sqrt(sum(p.^2,2))-1','p');
   bbox = [-1,-1;1,1];
else
   error(['Dont know what to do with mdl_name=',mdl_name]);
end

if elec_pattern == [16, 1]
   elecs_per_ring = elec_pattern(1);
   num_rings      = elec_pattern(2);
   stim_pattern = mk_stim_patterns(elecs_per_ring, num_rings, ...
                  '{ad}','{ad}',{},100);

   elec_width= .1;
   elec_nodes = mk_2d_circ_elec_nodes( elecs_per_ring, elec_width);
end

   h0= estimate_h0(bbox, nnodes);

   z_contact = 0;
   % Prevent singularity from zero contact impedance
   if any(z_contact)<.001
      z_contact=z_contact+0.001;
   end

   fwd_mdl= dm_mk_fwd_model( fd, h0, bbox, elec_nodes, stim_pattern, ...
                           z_contact, ['dm_mk_common_model= ',mdl_name]);


function  h0= estimate_h0(bbox, nnodes);
   dims= size(bbox,2);
   area_est= prod(abs(diff(dims,1)));
   h0 = (area_est/nnodes)^(1/dims);

function  elec_nodes = mk_2d_circ_elec_nodes( n_elecs, elec_width )
   radius=1;
   th_delta = elec_width/2/pi/radius;
% electrodes start top and go clockwise
   for i=1:n_elecs
      th= (i-1)*2*pi/n_elecs;
      th= th + th_delta*[-.5;0;.5];
      this_e_node= radius*[sin(th),cos(th)];
      elec_nodes{i} = this_e_node;
   end
