function [mapping, outside] = mk_coarse_fine_mapping(varargin)
%MK_COARSE_FINE_MAPPING: create a mapping matrix from coarse to fine FEM
% [c2f,out]= mk_coarse_fine_mapping( f_mdl, c_mdl, opt );
%
% C2F_ij is the fraction of f_mdl element i which is
%   contained in c_mdl element j. This is used to map
%   from data on the reconstruction model (c_mdl) to
%   the forward model f_mdl as 
%      elem_data_fine = Mapping*elem_data_coase
% OUT_i is the fraction of f_mdl element i which is not
%   contained in any c_mdl element.
%  
% Parameters:
%    c_mdl is coarse fwd_model
%    f_mdl is fine fwd_model
%    opt   is option structure (optional, see below)
%
% This function is a compatibility wrapper. The actual calculations are
% carried out by either MK_ANALYTIC_C2F or MK_APPROX_C2F (the older
% method). 
%
% To set a specific implementation use
%    eidors_default('set','mk_coarse_fine_mapping','mk_approx_c2f');
%
% To query current implementation use
%    eidors_default('get','mk_coarse_fine_mapping');
%
% NOTE that MK_ANALYTIC_C2F accepts an optional third argument whereas
% MK_APPROXIMATE_C2F does not.
%
% See also: EIDORS_DEFAULT, MK_ANALYTIC_C2F, MK_APPROX_C2F, CHECK_C2F_QUALITY

% (C) 2015 Bartlomiej Grychtol
% License: GPL version 2 or version 3
% $Id$


[mapping, outside] = eidors_default(varargin{:});
