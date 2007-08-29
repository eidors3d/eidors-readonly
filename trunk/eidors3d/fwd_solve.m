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
% data      is a measurement data structure
% fwd_model is a fwd_model structure
% img       is an img structure

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: fwd_solve.m,v 1.23 2007-08-29 09:20:57 aadler Exp $

if nargin==1
   img= fwd_model;
   fwd_model= img.fwd_model;
end
fwd_model= eidors_model_params( fwd_model );

% By default EIDORS doesn't cache images, so this doesn't normally help
data = eidors_obj('get-cache', img, 'fwd_solve_data');

if ~isempty(data)
   eidors_msg('fwd_solve: using cached value',3);
   return
end

data = feval( fwd_model.solve, fwd_model, img);
data= eidors_obj('data',data);  % create data object

eidors_obj('set-cache', img, 'fwd_solve_data', data);
eidors_msg('fwd_solve: setting cached value',3);
