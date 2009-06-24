function data= mv_fwd_solve( fwd_model, img)
% MV_FWD_SOLVE: data= mv_fwd_solve( fwd_model, img)
% Fwd solver for Marco Vauhkonen's EIDORS2D code
% data = measurements struct
% fwd_model = forward model
% img = image struct

% (C) 2009 Andy Adler. License: GPL version 2 or version 3
% $Id$

p= mv_fwd_parameters( fwd_model );

A= calc_system_mat( fwd_model, img );

[U,p,r]=ForwardSolution(p.n_node,p.n_node,A,p.C,p.T,[],'real');

% create a data structure to return
data.meas= U.Electrode(:);
data.time= -1; % unknown
data.name= 'solved by mv_fwd_solve';
