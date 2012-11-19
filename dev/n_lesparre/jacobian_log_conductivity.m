function J = jacobian_log_conductivity(jnk, img)
% JACOBIAN_LOG_CONDUCTIVITY - calculate Jacobian for log conductivity
% 
% By default, the calculation is based on JACOBIAN_ADJOINT, but can be
% changed through img.fwd_model.jacobian_log_conductivity.conductivity_jacobian.
%
% Requires img.elem_data.log_conductivity.
% Ignores img.fwd_model.normalize.

% (C) 2012 Nolwenn Lesparre and Bartlomiej Grychtol
% License: GPL version 2 or 3
% $Id$

if ~isfield(img.elem_data,'log_conductivity')
    error('img.elem_data.log_conductivity required');
end

logc = img.elem_data.log_conductivity;
img.elem_data = exp(logc);
try
    img.fwd_model.jacobian = ...
        img.fwd_model.jacobian_log_conductivity.conductivity_jacobian;
catch
    img.fwd_model.jacobian = 'eidors_default';
end

J = calc_jacobian(img);
if isfield(img.fwd_model,'coarse2fine') && ...
    size(img.fwd_model.coarse2fine,2) ~= length(img.elem_data);
    nc = size(img.fwd_model.coarse2fine,2);
    img.elem_data = mean(img.elem_data)*ones(nc,1);
end
dCond_dlogCond = img.elem_data; % d e^x /dx = e^x
J = J.*repmat((dCond_dlogCond),1,size(J,1))';



