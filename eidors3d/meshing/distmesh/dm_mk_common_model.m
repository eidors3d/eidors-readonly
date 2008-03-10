function [fwd_mdl]= dm_mk_common_model( mdl_name, elec_pattern );
% DM_MK_COMMON_MODEL: create a prepared model with distmesh
% [fwd_mdl]= dm_mk_common_model( mdl_name, nnodes, elec_pattern, stim_pattern );
%
% mdl_name:     'c2'   - 2D circle
% nnodes:       number of nodes in model
% elec_pattern:     [16, 1]  - 16 elecs per ring, one ring
% 
%  fwd_mdl:           eidors format fwd_model

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id: dm_mk_common_model.m,v 1.1 2008-03-10 14:16:19 aadler Exp $


      fd=inline('sqrt(sum(p.^2,2))-1','p');
      funi= inline('ones(size(p,1),1)','p');

