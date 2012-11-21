function data = fwd_solve_log_conductivity(jnk, img)
% FWD_SOLVE_LOG_CONDUCTIVITY solves voltage based on log_conductivity image
%
% Requires img.elem_data.log_conductivity OR img.node_data.log_conductivity
% Specific forward solver implementation for conductivity values can be set
% through img.fwd_model.fwd_solve_log_conductivity.conductivity_solver

% (C) 2012 Bartlomiej Grychtol
% License: GPL version 2 or 3
% $Id$

if isfield(img.log_conductivity,'elem_data')
    img.elem_data = exp(img.log_conductivity.elem_data);
else 
    img.node_data = exp(img.log_conductivity.node_data);
end

if isfield(img.fwd_model,'fwd_solve_log_conductivity')
    func = img.fwd_model.fwd_solve_log_conductivity.conductivity_solver;
else
    func = 'eidors_default';
end
img.fwd_model.solve = func;
data = fwd_solve(img.fwd_model, img);

