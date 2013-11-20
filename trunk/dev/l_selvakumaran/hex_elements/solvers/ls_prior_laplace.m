function Reg= ls_prior_laplace( inv_model );
% LS_PRIOR_LAPLACE calculate image prior for a mesh with quad/hex elements
% Reg= ls_prior_laplace( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct
%  or
% Reg= prior_laplace( fwd_model )
%
% This image prior is intended to be used as
%  R'*R, but may be used as R for as well.
%
% The Laplace prior is a 2nd order high pass filter.
% On a rectangular mesh, it is a convolution with
%   [-1,-1,-1;      [ 0;-1; 0
%    -1, 8,-1    or  -1; 4;-1
%    -1,-1,-1]        0;-1; 0]
%
% On a finite element mesh, we define the it as 
% -1 for each adjacent element, and 3 (in 2D) or 4 (in 3D)
% for the element itself

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: prior_laplace.m 3204 2012-06-29 17:45:03Z aadler $

if isstr(inv_model) && strcmp(inv_model,'UNIT_TEST'); do_unit_test; return; end

switch inv_model.type
  case 'inv_model'; fwd_model = inv_model.fwd_model;
  case 'fwd_model'; fwd_model = inv_model;
  otherwise; error('PRIOR_LAPLACE requires input type of inv_model or fwd_model');
end

pp= ls_fwd_model_parameters( fwd_model );

Reg = eidors_obj('get-cache', fwd_model , 'prior_laplace');
if ~isempty(Reg)
   eidors_msg('prior_laplace: using cached value', 3);
   return
end


Reg = speye( pp.n_elem );

   Iidx= [];
   Jidx= [];
   Vidx= [];
   for ii=1:pp.n_elem;
     el_adj = find_adjoin( ii, pp.ELEM );
     for jj=el_adj(:)';
         Iidx= [Iidx, ii, ii, jj, jj];
         Jidx= [Jidx, ii, jj, ii, jj];
         Vidx= [Vidx,  1, -1, -1,  1];
     end
   end

   Reg = sparse(Iidx,Jidx, Vidx, pp.n_elem, pp.n_elem );
   
eidors_obj('set-cache', fwd_model , 'prior_laplace', Reg);
eidors_msg('prior_laplace: setting cached value', 3);
   
   
% find elems which are connected to elems ee
function elems= find_adjoin(ee, ELEM)
   nn= ELEM(:,ee);
   [d,e]= size(ELEM);
   ss = false(size(ELEM));
   for i=1:d;
     ss= ss | ELEM==nn(i);
   end
   if d==4
       d1=2;
   elseif d==8
       d1=3;
   end
   elems= find(sum(ss,1)==d1);

