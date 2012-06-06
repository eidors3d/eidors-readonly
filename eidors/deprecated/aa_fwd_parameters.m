function param = aa_fwd_parameters( fwd_model )
% AA_FWD_PARAMETERS: data= aa_fwd_solve( fwd_model, image)
% Extract parameters from a 'fwd_model' struct which are
% appropriate for Andy Adler's EIT code
%   param.n_elem     => number of elements
%   param.n_elec     => number of electrodes
%   param.n_node     => number of nodes (vertices)
%   param.n_stim     => number of current stimulation patterns
%   param.n_elec     => number of electrodes
%   param.n_dims     => dimentions (2= 2D, 3=3D)
%   param.n_meas     => number of measurements (total)
%   param.boundary   => FEM boundary
%   param.NODE       => vertex matrix
%   param.ELEM       => connection matrix
%   param.QQ         => Current into each NODE
%   param.VOLUME     => Volume (or area) of each element
%   param.normalize  => difference measurements normalized?
%   param.N2E        => Node to electrode converter
%
% If the stimulation patterns has a 'interior_sources' field,
%   the node current QQ, is set to this value for this stimulation.

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

warning('EIDORS:deprecated','AA_FWD_PARAMETERS is deprecated as of 06-Jun-2012. Use FWD_MODEL_PARAMETERS instead.');

param = fwd_model_parameters( fwd_model );
