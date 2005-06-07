function img = inv_solve( inv_model, data1, data2)
% INV_SOLVE: calculate imag from an inv_model and data
% 
% inv_solve can be called as
%     img= inv_solve( inv_model, data1, data2)
%   if inv_model.type = 'differential'
% or
%     img= inv_solve( inv_model, data )
%   if inv_model.type = 'static'
%
% in each case it will call the inv_model.solve
%
% data      is a measurement data structure
% inv_model is a inv_model structure
% img       is an image structure
%
% $Id: inv_solve.m,v 1.6 2005-06-07 01:14:21 aadler Exp $

% TODO: does it make sense to cache solutions here?
%       if so, to where do they belong, data1 or data2?

eidors_msg('inv_solve',1);

% caching of images is disabled -
%  in order to implement image caching, it needs to be possible
%  to tag them with both data1 and data2

if     strcmp(inv_model.reconst_type,'static')
   if nargin~=2;
      error('only one data set is allowed for a static reconstruction');
   end
   img= feval( inv_model.solve, inv_model, data1);
elseif strcmp(inv_model.reconst_type,'differential')
   if nargin~=3;
      error('two data sets are required for a differential reconstruction');
   end
   img= feval( inv_model.solve, inv_model, data1, data2);
else
   error('inv_model.reconst_type not understood'); 
end

img= eidors_obj('image', img);
