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
%
% $Id: fwd_solve.m,v 1.3 2004-07-21 20:13:25 aadler Exp $

if nargin==1
   img= fwd_model;
   fwd_model= img.fwd_model;
end

data = eidors_obj('cache', img, 'fwd_solve_data');

if isempty(data)
   data = feval( fwd_model.solve, fwd_model, img);
   data= eidors_obj('data',data); 
%  disp('setting cached value');
   eidors_obj('cache', img, 'fwd_solve_data', data);
end

