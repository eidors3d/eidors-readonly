function img= mk_image(mdl, elem_data, params, name)
% MK_IMAGE: create eidors image object
%   img= mk_image(mdl, elem_data, name)
%
% Utility function to create an eidors_image object:
% Usage 1:
%   img = mk_image( inv_model )  -> uses jacobian backgnd for conductivity
%
% Usage 2: mdl can be a fwd_model or inv_model
%   img = mk_image( mdl, 1 ) -> uniform image with conductivity 1
%   img = mk_image( mdl, 2*ones(n_elems,1) ) -> uniform with c=2
%   img = mk_image( mdl, 2*ones(n_nodes,1) ) -> image with node data
%   img = mk_image( mdl, 1, 'This name') -> Specify a 'name' attribute
%   img = mk_image( mdl, ones(n_elems,3) );
%
% Usage 3: create image from previous image, override conductity
%  img = mk_image( other_image, 2 ) -> image with c=2
%
% Usage 4: create image with specific 'params' properties
%  img = mk_image( mdl, 3*ones(64,1),'resistivity' ); % resistivity image

% (C) 2008-12 Andy Adler and Bartlomiej Grychtol.
% Licenced under GPL version 2 or 3
% $Id$

if isstr(mdl) && strcmp(mdl,'UNIT_TEST'); do_unit_test; return; end

default_params = no_params; % later: conductivity ?
default_name    = 'Created by mk_image';
if nargin<3
    name = default_name;
    params = default_params;
elseif nargin < 4 % preserve old interface
    if ismember(params,supported_params)
        name = default_name;
    else
        name = params;
        params = default_params;
    end
end

% if nargin<2;
%    try
%      elem_data = mdl.jacobian_bkgnd.value;
%    catch
%      error('mk_image: for one parameter, needs a mdl.jacobian_bkgnd field');
%    end
% end

switch mdl.type
   case 'inv_model'
      if nargin == 1
         img = calc_jacobian_bkgnd(mdl);
         return
      else
         mdl = mdl.fwd_model;
         % warning('Ignoring inv_model.jacobian_bkgnd in favour of the supplied elem_data');
      end
   case 'fwd_model'
      mdl = mdl; % keep model
   case 'image'
      if nargin == 1
         img = data_mapper(mdl);
         return
      end
      mdl = mdl.fwd_model;
   otherwise; error('mk_image: no inv_model, fwd_model or image object');
end

img = eidors_obj('image',name);
img.fwd_model = mdl;
if isfield(mdl,'show_slices');
    img.show_slices = mdl.show_slices;
end
img = fill_in_data(img,elem_data,params);
% img.current_params = params;

function str = no_params
str = 'unspecified';

function img = fill_in_data(img,elem_data,params)
  % If image is presented as row vector, then transpose
  %  Shouldn't happen, but for matlab bugs
  if size(elem_data,1) == 1 && ...
     size(elem_data,2) >  1 
     elem_data = elem_data.';
  end
  switch size(elem_data,1)
     case 1
        sz = size(img.fwd_model.elems,1);
        elem_data_mat =    NaN*ones(sz,1);
        elem_data_mat(:) = elem_data;
        if strcmp(params,no_params);
            img.elem_data = elem_data_mat;
        else
            img.(params).elem_data = elem_data_mat;
        end

     case num_elems( img )
        if strcmp(params,no_params);
            img.elem_data = elem_data;
        else
            img.(params).elem_data = elem_data;
        end

     case size(img.fwd_model.coarse2fine,2)
        if strcmp(params,no_params);
            img.elem_data = elem_data;
        else
            img.(params).elem_data = elem_data;
        end

     case num_nodes( img )
        if strcmp(params,no_params);
            img.node_data = elem_data;
        else
            img.(params).node_data = elem_data;
        end

     otherwise
        error('Don''t understand number of elements.');
  end

% TESTS:
function do_unit_test
   imdl = mk_common_model('a2c2',8);
   im0 = mk_image( imdl );
   unit_test_cmp('im1',im0.elem_data, ones(64,1) );

   im0 = mk_image( imdl, 2 ); % warning
   unit_test_cmp('im2',im0.elem_data, 2*ones(64,1) );

   im0 = mk_image( imdl.fwd_model, 3*ones(64,1) );
   unit_test_cmp('im3',im0.elem_data, 3*ones(64,1) );

   im0 = mk_image(im0); % no change
   unit_test_cmp('im4',im0.elem_data, 3*ones(64,1) );

   im0.conductivity.elem_data = im0.elem_data;
   im0 = rmfield(im0, 'elem_data');
   im1 = mk_image(im0); % use data_mapper
   unit_test_cmp('im5',im1.elem_data, 3*ones(64,1) );

   im0.resistivity.node_data = 5;
   im0.data_mapper = 'resistivity';
   im1 = mk_image(im0); % use data_mapper
   unit_test_cmp('im6',im1.node_data(5), 5 );

   imdl.jacobian_bkgnd = rmfield(imdl.jacobian_bkgnd, 'value');
   imdl.jacobian_bkgnd.conductivity = im0.conductivity;
   im2 = mk_image(imdl);
   unit_test_cmp('im7',im2.conductivity.elem_data, 3*ones(64,1) );

   im0 = mk_image( imdl.fwd_model, 3*ones(64,1),'my_name' );
   unit_test_cmp('im8a',im0.elem_data, 3*ones(64,1) );

   im0 = mk_image( imdl.fwd_model, 3*ones(64,1),'resistivity' );
   unit_test_cmp('im8b',im0.resistivity.elem_data, 3*ones(64,1) );

   im0 = mk_image( imdl.fwd_model, 2*ones(64,3),'resistivity' );
   unit_test_cmp('im8c',im0.resistivity.elem_data, 2*ones(64,3) );

   im0 = mk_image( imdl.fwd_model, 3*ones(64,1),'resistivity','my name' );
   unit_test_cmp('im9a',im0.resistivity.elem_data, 3*ones(64,1) );
   unit_test_cmp('im9b',im0.name, 'my name' );
