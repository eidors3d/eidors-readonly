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

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: fwd_solve.m,v 1.10 2005-10-27 13:28:08 aadler Exp $

if nargin==1
   img= fwd_model;
   fwd_model= img.fwd_model;
end

% By default EIDORS doesn't cache images, so this doesn't normally help
data = eidors_obj('get-cache', img, 'fwd_solve_data');

if ~isempty(data)
   eidors_msg('fwd_solve: using cached value',2);
   return
end

data = feval( fwd_model.solve, fwd_model, img);
data= eidors_obj('data',data);  % create data object

eidors_obj('set-cache', img, 'fwd_solve_data', data);
eidors_msg('fwd_solve: setting cached value',2);
