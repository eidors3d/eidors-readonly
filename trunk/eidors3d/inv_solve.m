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
% For difference EIT:
% data1      => difference data at earlier time
% data2      => difference data at later time
%
% both data1 and data2 may be vectors of data
%  structures. In which case a vector of
%  reconstructed images will be returned
% if either data1 or data2 is a scalar structure, then it
%  is expanded to be the same size matrix
%
% $Id: inv_solve.m,v 1.7 2005-06-29 16:39:28 aadler Exp $

% COMMENT: There seems to be no general way to cache
%       inv_model parameters. Thus, each algorithm needs
%       to do it separately. It may be possible to cache
%       one-step inverse matrices, but that is not done. 

eidors_msg('inv_solve',1);

% caching of images is disabled -
%  in order to implement image caching, it needs to be possible
%  to tag them with both data1 and data2

if     strcmp(inv_model.reconst_type,'static')
   if nargin~=2;
      error('only one data set is allowed for a static reconstruction');
   end
   img= feval( inv_model.solve, inv_model, data1);
elseif strcmp(inv_model.reconst_type,'difference')
   if nargin~=3;
      error('two data sets are required for a difference reconstruction');
   end

   l_data1= length(data1);
   l_data2= length(data2);
   if l_data1 ~=1 & l_data2 ~=1 & ...
         l_data1 ~= l_data2
      error('inconsistent number of specified measurements');
   end 

   img= feval( inv_model.solve, inv_model, data1, data2);
else
   error('inv_model.reconst_type not understood'); 
end

img= eidors_obj('image', img);
