function fmdl= mk_grid_model(xvec, yvec, zvec);
% MK_GRID_MODEL: Create reconstruction model on pixelated grid 
%  fmdl= mk_grid_model(xvec, yvec, zvec);
%
%  fmdl - eidors forward model
%  xvec - x edges
%  yvec - y edges
%  zvec - z edges (optional - to create 3D model)

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id: mk_grid_model.m,v 1.3 2008-05-01 16:27:23 aadler Exp $

if nargin ==2
   fmdl = mk_2d_grid(xvec,yvec);
elseif nargin ==3
   fmdl = mk_3d_grid(xvec,yvec,zvec);
else
   error('check nargin');
end

function fmdl= mk_2d_grid(xvec, yvec);
   xlen = length(xvec);
   ylen = length(yvec);
   fmdl= eidors_obj('fwd_model', ...
            sprintf('Grid model %d x %d', xlen, ylen) );

   [x,y]= meshgrid( xvec, yvec);
   x=x';y=y';
   fmdl.nodes= [x(:),y(:)];
   k= 1:xlen-1;
   elem_frac = [ k;k+1;k+xlen; ...
                 k+1;k+xlen;k+xlen+1];
   elem_frac= reshape(elem_frac, 3,[])';
   fmdl.elems=  [];
   for j=0:ylen-2
      fmdl.elems=  [fmdl.elems; elem_frac + xlen*j];
   end

   fmdl.boundary = find_boundary( fmdl.elems);

% assign one single parameter to each square element
   e= size(fmdl.elems,1);
   params= ceil(( 1:e )/2);
   fmdl.coarse2fine = sparse(1:e,params,1,e,max(params));

