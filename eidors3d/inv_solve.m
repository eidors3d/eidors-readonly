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
% $Id: inv_solve.m,v 1.2 2004-07-21 21:09:14 aadler Exp $

if     strcmp(inv_model.reconst_type,'static')
   if nargin~=2;
      error('only one data set is allowed for a static reconstruction');
   end
   img= feval( inv_model.solve, inv_model, data1);
elseif strcmp(inv_model.reconst_type,'differential')
   if nargin~=3;
      error('two data set are required for a differential reconstruction');
   end
   img= feval( inv_model.solve, inv_model, data1, data2);
else
   error('inv_model.reconst_type not understood'); 
end
