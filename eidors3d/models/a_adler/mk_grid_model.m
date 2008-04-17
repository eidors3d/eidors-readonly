function fmdl= mk_grid_model(xvec, yvec);
% MK_GRID_MODEL: Create reconstruction model on pixelated grid 

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id: mk_grid_model.m,v 1.1 2008-04-17 19:44:47 aadler Exp $

xlen = length(xvec);
ylen = length(yvec);
fmdl= eidors_obj('fwd_model', ...
         sprintf('Grid model %d x %d', xlen, ylen) );

[x,y]= meshgrid( xvec, yvec);
    x=x';y=y';
    fmdl.nodes= [x(:),y(:)];
    k= 1:ylen-1;
    elem_frac = [ k;k+1;k+ylen; ...
                  k+1;k+ylen;k+ylen+1];
    elem_frac= reshape(elem_frac, 3,[])';
    fmdl.elems=  [];
    for j=0:ylen-2
       fmdl.elems=  [fmdl.elems; elem_frac + ylen*j];
    end

    fmdl.boundary = find_boundary( fmdl.elems);
