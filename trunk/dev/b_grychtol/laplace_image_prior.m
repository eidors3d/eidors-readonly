function Reg= prior_laplace( inv_model );
% PRIOR_LAPLACE calculate image prior
% Reg= prior_laplace( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct
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
% $Id$

pp= aa_fwd_parameters( inv_model.fwd_model );
pp = grid_or_elems( pp, inv_model);
Reg = speye( pp.n_elem );

   Iidx= [];
   Jidx= [];
   Vidx= [];
   for ii=1:pp.n_elem
     el_adj = find_adjoin( ii, pp.ELEM, pp.n_dims );
     for jj=el_adj(:)'
         Iidx= [Iidx, ii, ii, jj, jj];
         Jidx= [Jidx, ii, jj, ii, jj];
         Vidx= [Vidx,  1, -1, -1,  1];
     end
   end

   Reg = sparse(Iidx,Jidx, Vidx, pp.n_elem, pp.n_elem );

% find elems which are connected to elems ee
function elems= find_adjoin(ee, ELEM, n_dims)
   nn= ELEM(:,ee);
   [d,e]= size(ELEM);
   ss= zeros(1,e);
   for i=1:d
     ss= ss+ any(ELEM==nn(i));
   end
   elems= find(ss==n_dims);

   function pp = grid_or_elems( pp, inv_model)

    fwd_model = inv_model.fwd_model;
    if isfield(fwd_model,'grid')
        use_grid = 1;
        try if inv_model.prior_use_elems_not_grid== 1
                use_grid = 0;
        end; end
    
        if use_grid
            pp.ELEM = fwd_model.grid';
            pp.n_elem = size(pp.ELEM,2);
        end
    end

        