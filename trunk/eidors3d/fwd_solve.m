function data = fwd_solve( fwd_model, image)
% FWD_SOLVE: calculate data from a fwd_model object and an image
% 
% fwd_solve can be called as
%    data= fwd_solve( fwd_model, image)
% or
%    data= fwd_solve( image)
%
% in each case it will call the fwd_model.solve
%                      or image.fwd_model.solve method
%
% data      is a measurement data structure
% fwd_model is a fwd_model structure
% image     is an image structure
%
% $Id: fwd_solve.m,v 1.2 2004-07-21 19:37:06 aadler Exp $

if nargin==1
   image= fwd_model;
   fwd_model= image.fwd_model;
end
data = feval( fwd_model.solve, fwd_model, image);

data= eidors_obj('data',data); 

