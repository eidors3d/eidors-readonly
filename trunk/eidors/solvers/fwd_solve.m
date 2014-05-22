function data = fwd_solve(fwd_model, img)
% FWD_SOLVE: calculate data from a fwd_model object and an image
% 
% fwd_solve can be called as
%    data= fwd_solve( img)
% or (deprecated)
%    data= fwd_solve( fwd_model, img)
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
%
% Parameters:
%    fwd_model.background = constant conductivity offset added to elem_data
%    fwd_model.coarse2fine = linear mapping between img.elem_data and model parameters
%    img.params_mapping = function mapping img.elem_data to model parameters
%

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

if nargin == 1
   img= fwd_model;
else
   warning('EIDORS:DeprecatedInterface', ...
      ['Calling FWD_SOLVE with two arguments is deprecated and will cause' ...
       ' an error in a future version. First argument ignored.']);
end
ws = warning('query','EIDORS:DeprecatedInterface');
warning off EIDORS:DeprecatedInterface

fwd_model= img.fwd_model;


fwd_model= prepare_model( fwd_model );

% TODO: This should be handled by the data_mapper
if isfield(img,'params_mapping')
%     fwd_model data is provided using a mapping function
    mapping_function= img.params_mapping.function;
    img= feval(mapping_function,img);
end
if isfield(fwd_model,'coarse2fine') && isfield(img,'elem_data')
   c2f= fwd_model.coarse2fine;
   if size(img.elem_data,1)==size(c2f,2)
%     fwd_model data is provided on coarse mesh
      img.elem_data = c2f * img.elem_data; 

      if isfield(fwd_model,'background')
          img.elem_data = img.elem_data + fwd_model.background; 
      end
   end
end

if ~isfield(fwd_model, 'electrode')
   error('EIDORS: attempting to solve on model without electrodes');
end
if ~isfield(fwd_model, 'stimulation')
   error('EIDORS: attempting to solve on model without stimulation patterns');
end

data = feval( fwd_model.solve, fwd_model, img);
data= eidors_obj('data',data);  % create data object


if isa(fwd_model.solve,'function_handle')
    solver = func2str(fwd_model.solve);
else
    solver = fwd_model.solve;
end
if strcmp(solver,'eidors_default');
    solver = eidors_default('get','fwd_solve');
end
if isfield(fwd_model,'measured_quantity') && ~isfield(data,'measured_quantity')
   warning('EIDORS:MeasurementQuantityObliviousSolver',...
      ['The solver %s did not handle the requested measurement quantity properly.\n'...
       'The results may be incorrect. Please check the code to verify.'], ...
       solver);
elseif isfield(fwd_model,'measured_quantity') ... 
        && isfield(data,'measured_quantity') ...
        && ~strcmp(fwd_model.measured_quantity, data.measured_quantity)
   error('EIDORS:MeasurementQuantityDisagreement',...
       'The solver %s return measurements as %s, while %s was expected.',...
       solver, data.measured_quantity, fwd_model.measured_quantity);
end
    

warning on EIDORS:DeprecatedInterface

eidors_obj('set-cache', img, 'fwd_solve_data', data);
eidors_msg('fwd_solve: setting cached value',3);

function mdl = prepare_model( mdl )
mdl = mdl_normalize(mdl,mdl_normalize(mdl));
if ~isfield(mdl,'elems');
    return;
end

mdl.elems  = double(mdl.elems);
mdl.n_elem = size(mdl.elems,1);
mdl.n_node = size(mdl.nodes,1);
if isfield(mdl,'electrode');
    mdl.n_elec = length(mdl.electrode);
else
    mdl.n_elec = 0;
end
