function img_bkgnd = calc_jacobian_bkgnd( inv_model )
% CALC_JACOBIAN_BKGND: calculate background image around
%    which initial estimate of jacobian is calculated
% 
% img_bkgnd = calc_jacobian_bkgnd( inv_model )
% inv_model   is an EIDORS fwd_model 
% img_bkgnd   is an EIDORS struct
%
% Usage: 
%      The background for calc_jacobian may be specified
%      as an estimated value
%  inv_model.jacobian_bkgnd.value;  % scalar OR
%                                   % vector of Nx1
%
% Usage:
%      The background may be calculated by a function
%  inv_model.jacobian_bkgnd.func;

% (C) 2005-2012 Andy Adler and Bartlomiej Grychtol
% License: GPL version 2 or version 3
% $Id$

if ischar(inv_model) && strcmp(inv_model,'UNIT_TEST') do_unit_test; return;end

img_bkgnd = eidors_cache(@jacobian_bkgnd,{inv_model},'calc_jacobian_bkgnd');

function img_bkgnd = jacobian_bkgnd(inv_model)

%%% old interface %%%
if isfield(inv_model.jacobian_bkgnd,'func')
   img_bkgnd= feval( inv_model.jacobian_bkgnd.func, inv_model );
elseif isfield(inv_model.jacobian_bkgnd,'value')
   if has_params(inv_model.jacobian_bkgnd)
      warning('Ignoring parametrization-specific fields in Jacobian background')
   end
   % allow bkgnd to be scalar or vector
   fwd_model= inv_model.fwd_model;
   bkgnd = ones(size(fwd_model.elems,1),1);
   bkgnd(:)= inv_model.jacobian_bkgnd.value;

   img_bkgnd= eidors_obj('image', 'background image', ...
                         'elem_data', bkgnd, ...
                         'fwd_model', fwd_model );
else
%%% new interface with parametrizations (params) %%%
    fwd_model = inv_model.fwd_model;
    img_bkgnd = eidors_obj('image','background image','fwd_model',fwd_model);
    flds = fieldnames(inv_model.jacobian_bkgnd);
    % copy params from img
    for i = 1:length(flds)
       img_bkgnd.(flds{i}) = inv_model.jacobian_bkgnd.(flds{i});
    end
%     img_bkgnd = data_mapper(img_bkgnd);
end


function b = has_params(s)
b = false;
if isstruct(s)
   b = any(ismember(fieldnames(s),supported_params));
end

function do_unit_test
imdl = mk_common_model('d2c2');
test = calc_jacobian_bkgnd(imdl)
unit_test_cmp('t1:type',  test.type, 'image');
unit_test_cmp('t1:elem_data',  test.elem_data, ones(1024,1));

imdl.jacobian_bkgnd = rmfield(imdl.jacobian_bkgnd,'value');
imdl.jacobian_bkgnd.node_data  = 5;
test = calc_jacobian_bkgnd(imdl);
unit_test_cmp('t2:type',  test.type, 'image');
unit_test_cmp('t2:node_data',  test.node_data, 5);

imdl.jacobian_bkgnd = rmfield(imdl.jacobian_bkgnd,'node_data');
imdl.jacobian_bkgnd.node_data.val1  = 5;
imdl.jacobian_bkgnd.node_data.val2  = ones(length(imdl.fwd_model.nodes),1);
img = calc_jacobian_bkgnd(imdl);
unit_test_cmp('t3:node_data',  img.node_data, imdl.jacobian_bkgnd.node_data);

imdl.jacobian_bkgnd = rmfield(imdl.jacobian_bkgnd,'node_data');
imdl.jacobian_bkgnd.elem_data.val1  = 5;
imdl.jacobian_bkgnd.elem_data.val2  = ones(length(imdl.fwd_model.elems),1);
img = calc_jacobian_bkgnd(imdl);
unit_test_cmp('t4:elem_data',  img.elem_data, imdl.jacobian_bkgnd.elem_data);

imdl = rmfield(imdl,'jacobian_bkgnd');
imdl.jacobian_bkgnd.resistivity.elem_data = 3;
img = calc_jacobian_bkgnd(imdl);

unit_test_cmp('t5:resistivity',  img, struct( ...
    'type','image', 'name', 'background image', ...
    'fwd_model', imdl.fwd_model, 'resistivity', ...
     imdl.jacobian_bkgnd.resistivity ));
