function [cmdl, c2f]= mk_grid_model(fmdl, xvec, yvec, zvec);
% MK_GRID_MODEL: Create reconstruction model on pixelated grid 
%  [cmdl,coarse2fine]= mk_grid_model(xvec, yvec, zvec, fmdl);
%
% Outputs:
%  cmdl - eidors reconstruction model (coarse model)
%  coarse2fine - c2f mapping to put onto fmdl (specify [] to not use)
%
% Inputs:
%  fmdl - fine model (forward model) to create coarse2fine mapping
%  xvec - x edges
%  yvec - y edges
%  zvec - z edges (optional - to create 3D model)

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id: mk_grid_model.m,v 1.4 2008-05-01 17:51:25 aadler Exp $

if nargin == 3
   cmdl = mk_2d_grid(xvec,yvec);
elseif nargin ==4
   cmdl = mk_3d_grid(xvec,yvec,zvec);
else
   error('check nargin');
end

if ~isempty( fmdl)
   c2f= mk_coarse_fine_mapping( fmdl, cmdl);
end

function cmdl= mk_2d_grid(xvec, yvec);
   xlen = length(xvec);
   ylen = length(yvec);
   cmdl= eidors_obj('fwd_model', ...
            sprintf('Grid model %d x %d', xlen, ylen) );

   [x,y]= meshgrid( xvec, yvec);
   x=x';y=y';
   cmdl.nodes= [x(:),y(:)];
   k= 1:xlen-1;
   elem_frac = [ k;k+1;k+xlen; ...
                 k+1;k+xlen;k+xlen+1];
   elem_frac= reshape(elem_frac, 3,[])';
   cmdl.elems=  [];
   for j=0:ylen-2
      cmdl.elems=  [cmdl.elems; elem_frac + xlen*j];
   end

   cmdl.boundary = find_boundary( cmdl.elems);

% assign one single parameter to each square element
   e= size(cmdl.elems,1);
   params= ceil(( 1:e )/2);
   cmdl.coarse2fine = sparse(1:e,params,1,e,max(params));

