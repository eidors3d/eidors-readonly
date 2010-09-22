function data = fwd_solve( fwd_model, img)
% FWD_SOLVE: calculate data from a fwd_model object and an image
% 
% fwd_solve can be called as
%    data= fwd_solve( fwd_model, img)
% or
%    data= fwd_solve( img)
%
% in each case it will call the fwd_model.solve
%                        or img.fwd_model.solve method
%
% For reconstructions on dual meshes, the interpolation matrix
%    is defined as fwd_model.coarse2fine. If required, this takes
%    coarse2fine * x_coarse = x_fine
%
% data      is a measurement data structure
% fwd_model is a fwd_model structure
% img       is an img structure
%
% Options: (not available on all solvers)
%    img.fwd_solve.get_all_meas = 1 (data.volt = all FEM nodes, but not CEM)
%    img.fwd_solve.get_all_nodes= 1 (data.volt = all nodes, including CEM)

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

if nargin==1
   img= fwd_model;
   fwd_model= img.fwd_model;
end
fwd_model= eidors_model_params( fwd_model );

if isfield(fwd_model,'coarse2fine')
   c2f= fwd_model.coarse2fine;
   if size(img.elem_data,1)==size(c2f,2)
%     fwd_model data is provided on coarse mesh
      img.elem_data = c2f * img.elem_data; 
   end
end

if isfield(fwd_model,'background')
    img.elem_data = img.elem_data + fwd_model.background; 
end

data = feval( fwd_model.solve, fwd_model, img);
data= eidors_obj('data',data);  % create data object

eidors_obj('set-cache', img, 'fwd_solve_data', data);
eidors_msg('fwd_solve: setting cached value',3);
