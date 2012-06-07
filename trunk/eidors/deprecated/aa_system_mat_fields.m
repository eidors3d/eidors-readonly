function FC= aa_system_mat_fields( fwd_model )
% AA_SYSTEM_MAT_FIELDS: fields (elem to nodes) fraction of system mat
% FC= aa_system_mat_fields( fwd_model )
% input: 
%   fwd_model = forward model
% output:
%   FC:        s_mat= C' * S * conduct * C = FC' * conduct * FC;

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id$

warning('EIDORS:deprecated','AA_SYSTEM_MAT_FIELDS is deprecated as of 07-Jun-2012. Use SYSTEM_MAT_FIELDS instead.');

FC = system_mat_fields( fwd_model );
