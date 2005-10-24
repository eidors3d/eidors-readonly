function inv_mdl= mk_common_model( varargin )
% MK_COMMON_MODEL: make common EIT models
%
% Utility function to create common EIT FEM models,
% so that users do not need to re-write common code
%
% Usage: 
%   mk_common_model('ac',16) 
%
%   mk_common_model('dr',16)   - circular ring with 16 electrodes
%   mk_common_model('dr2',16)  - two circular rings with 16 electrodes
%
%   mk_common_model('dz',16)   - zigzag pattern electrodes
