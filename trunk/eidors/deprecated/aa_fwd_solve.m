function data =aa_fwd_solve(varargin)
% AA_FWD_SOLVE: data= aa_fwd_solve( fwd_model, img)
% Fwd solver for Andy Adler's EIT code
% Input:
%    fwd_model = forward model
%    img       = image struct
% Output:
%    data = measurements struct
% Options: (to return internal FEM information)
%    img.fwd_solve.get_all_meas = 1 (data.volt = all FEM nodes, but not CEM)
%    img.fwd_solve.get_all_nodes= 1 (data.volt = all nodes, including CEM)

% (C) 1995-2002 Andy Adler. License: GPL version 2 or version 3
% Ref: Adler & Guardo (1996) IEEE T. Med Imaging
% $Id$

% correct input paralemeters if function was called with only img
warning('EIDORS:deprecated','AA_FWD_SOLVE is deprecated as of 07-Jun-2012. Use FWD_SOLVE_1ST_ORDER instead.');

data = fwd_solve_1st_order( varargin{:} );
