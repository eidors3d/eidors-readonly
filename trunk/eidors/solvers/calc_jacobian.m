function J = calc_jacobian( fwd_model, img)
% CALC_JACOBIAN: calculate jacobian from an inv_model
% 
%  J = calc_jacobian( img )
%      calc Jacobian on img.fwd_model at conductivity given
%      in image (fwd_model is for forward and reconstruction)
%
% For reconstructions on dual meshes, the interpolation matrix
%    is defined as img.fwd_model.coarse2fine. This takes
%    coarse2fine * x_coarse = x_fine
%
% If the underlying jacobian calculator doesn't understand dual
%    meshes, then calc_jacobian will automatically postmultiply
%    by fwd_model.coarse2fine.
%
% img       is an image structure, with 'elem_data' or
%           'node_data' parameters

% (C) 2005-08 Andy Adler. License: GPL version 2 or version 3
% $Id$


if nargin == 1
   img= fwd_model;
else
   warning('EIDORS:DeprecatedInterface', ...
      ['Calling CALC_JACOBIAN with two arguments is deprecated and will cause' ...
       ' an error in a future version. First argument ignored.']);
end
ws = warning('query','EIDORS:DeprecatedInterface');
warning off EIDORS:DeprecatedInterface

try 
   fwd_model= img.fwd_model;
catch
   error('CALC_JACOBIAN requires an eidors image structure');
end

fwd_model_check(fwd_model);

cache_obj= jacobian_cache_params( fwd_model, img );

J= eidors_obj('get-cache', cache_obj, 'jacobian');
if ~isempty(J)
   eidors_msg('calc_jacobian: using cached value', 3);
   return
end

J= feval(fwd_model.jacobian, fwd_model, img);

if isfield(fwd_model,'coarse2fine')
   c2f= fwd_model.coarse2fine;
   if size(J,2)==size(c2f,1)
%     calc_jacobian did not take into account the coarse2fine
      J=J*c2f;
   end
end

warning(ws.state, 'EIDORS:DeprecatedInterface');

eidors_obj('set-cache', cache_obj, 'jacobian', J);
eidors_msg('calc_jacobian: setting cached value', 3);
        



% Make the Jacobian only depend on 
function cache_obj= jacobian_cache_params( fwd_model, img );
   img = physics_data_mapper(img);
   if isfield(img, 'elem_data')
      cache_obj = {fwd_model, img.elem_data, img.current_physics};
   elseif isfield(img, 'node_data')
      cache_obj = {fwd_model, img.node_data, img.current_physics};
   else
      error('calc_jacobian: execting elem_data or node_data in image');
   end

function fwd_model_check(fmdl)
pp = fwd_model_parameters(fmdl); % they cache, so no problem
if pp.n_elec == 0
    error('Cannot calculate Jacobian. No electrodes found.');
end
if pp.n_stim == 0
    error('Cannot calculate Jacobian. No stimulation found.');
end
if pp.n_meas == 0
    error('Cannot calculate Jacobian. No measurements found.');
end
