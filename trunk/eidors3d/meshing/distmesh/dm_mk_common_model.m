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
% $Id: dm_mk_common_model.m,v 1.2 2008-03-10 14:41:28 aadler Exp $


if mdl_name == 'c2'
   fd=inline('sqrt(sum(p.^2,2))-1','p');
   fh= inline('ones(size(p,1),1)','p'); % initially uniform
   bbox = [-1,-1;1,1];
else
   error(['Dont know what to do with mdl_name=',mdl_name]);
end

   h0= estimate_h0(bbox, nnodes);

   z_contact = 0;
   % Prevent singularity from zero contact impedance
   if any(z_contact)<.001
      z_contact=z_contact+0.001;
   end

   stim_pattern= [];
   f_mdl= dm_mk_fwd_model( fd, fh, h0, bbox, stim_pattern, ...
                           z_contact, ['dm_mk_common_model= ',mdl_name]);


function  h0= estimate_h0(bbox, nnodes);
   dims= size(bbox,2);
   area_est= prod(abs(diff(dims,1)));
   h0 = (area_est/nnodes)^(1/dims);
