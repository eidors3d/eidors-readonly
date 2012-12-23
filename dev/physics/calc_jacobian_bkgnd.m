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

if isstr(inv_model) && strcmp(inv_model,'UNIT_TEST') do_unit_test; return;end

img_bkgnd= eidors_obj('get-cache', inv_model, 'jacobian_bkgnd');
if ~isempty(img_bkgnd)
   eidors_msg('calc_jacobian_bkgnd: using cached value', 3);
   return
end
%%% old interface %%%
if isfield(inv_model.jacobian_bkgnd,'func')
   img_bkgnd= feval( inv_model.jacobian_bkgnd.func, inv_model );
elseif isfield(inv_model.jacobian_bkgnd,'value')
   % allow bkgnd to be scalar or vector
   fwd_model= inv_model.fwd_model;
   bkgnd = ones(size(fwd_model.elems,1),1);
   bkgnd(:)= inv_model.jacobian_bkgnd.value;

   img_bkgnd= eidors_obj('image', 'background image', ...
                         'elem_data', bkgnd, ...
                         'fwd_model', fwd_model );
else
%%% new interface with physics %%%
    fwd_model = inv_model.fwd_model;
    img_bkgnd = eidors_obj('image','background image','fwd_model',fwd_model);
    flds = fieldnames(inv_model.jacobian_bkgnd);
    % copy physics from to im
    for i = 1:length(flds)
       img_bkgnd.(flds{i}) = inv_model.jacobian_bkgnd.(flds{i});
    end
%     img_bkgnd = physics_data_mapper(img_bkgnd);
end


eidors_obj('set-cache', inv_model, 'jacobian_bkgnd', img_bkgnd);
eidors_msg('jacobian_bkgnd: setting cached value', 3);


function do_unit_test
imdl = mk_common_model('d2c2');
calc_jacobian_bkgnd(imdl)
imdl.jacobian_bkgnd = rmfield(imdl.jacobian_bkgnd,'value');
imdl.jacobian_bkgnd.node_data  = 5;
calc_jacobian_bkgnd(imdl)
imdl.jacobian_bkgnd.node_data.val1  = 5;
imdl.jacobian_bkgnd.node_data.val2  = ones(length(imdl.fwd_model.nodes),1);
img = calc_jacobian_bkgnd(imdl);
display(img.node_data);
imdl.jacobian_bkgnd = rmfield(imdl.jacobian_bkgnd,'node_data');
imdl.jacobian_bkgnd.elem_data.val1  = 5;
imdl.jacobian_bkgnd.elem_data.val2  = ones(length(imdl.fwd_model.elems),1);
img = calc_jacobian_bkgnd(imdl);
display(img.elem_data);
imdl = rmfield(imdl,'jacobian_bkgnd');
imdl.jacobian_bkgnd.resistivity.elem_data = 3;
img = calc_jacobian_bkgnd(imdl);


