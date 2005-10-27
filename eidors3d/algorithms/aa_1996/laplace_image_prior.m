function Reg= laplace_image_prior( inv_model );
% LAPLACE_IMAGE_PRIOR calculate image prior
% Reg= laplace_image_prior( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct
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

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: laplace_image_prior.m,v 1.2 2005-10-27 13:28:08 aadler Exp $

pp= aa_fwd_parameters( inv_model.fwd_model );

Reg = speye( pp.n_elem );

   Iidx= [];
   Jidx= [];
   Vidx= [];
   for ii=1:pp.n_elem
     el_adj = find_adjoin( ii, pp.ELEM );
     for jj=el_adj(:)'
         Iidx= [Iidx, ii, ii, jj, jj];
         Jidx= [Jidx, ii, jj, ii, jj];
         Vidx= [Vidx,  1, -1, -1,  1];
     end
   end

   Reg = sparse(Iidx,Jidx, Vidx, pp.n_elem, pp.n_elem );

% find elems which are connected to elems ee
function elems= find_adjoin(ee, ELEM)
   nn= ELEM(:,ee);
   [d,e]= size(ELEM);
   ss= zeros(1,e);
   for i=1:d
     ss= ss+ any(ELEM==nn(i));
   end
   elems= find(ss==d-1);
