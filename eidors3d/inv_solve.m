function image = inv_solve( inv_model, data1, data2)
% INV_SOLVE: calculate imag from an inv_model and data
% 
% inv_solve can be called as
%     image= inv_solve( inv_model, data1, data2)
%   if inv_model.type = 'differential'
% or
%     image= inv_solve( inv_model, data )
%   if inv_model.type = 'static'
%
% in each case it will call the inv_model.solve
%
% data      is a measurement data structure
% inv_model is a inv_model structure
% image     is an image structure
%
% $Id: inv_solve.m,v 1.1 2004-07-10 02:40:22 aadler Exp $

if     strcmp(inv_model.type,'static')
   if nargin~=2;
      error('only one data set is allowed for a static reconstruction');
   end
   image= feval( inv_model.solve, inv_model, data1);
elseif strcmp(inv_model.type,'differential')
   if nargin~=3;
      error('two data set are required for a differential reconstruction');
   end
   image= feval( inv_model.solve, inv_model, data1, data2);
end
