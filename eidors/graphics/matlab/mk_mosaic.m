function r_img = mk_mosaic(rimg, sep, vh, n_col)
%MK_MOSAIC Arrange multidimensional image matrix for display.
%  M = MK_MOSAIC(rimg, sep, vh, n_col)
%  Input:
%   rimg - an HxWxN or HxWxNxM array of pictures of size HxW
%   sep  - (optional) spacing between adjecent pictures (in pixels)
%          default: 0
%   vh   - an Mx2 array of positions for individual HxWxN blocks 
%          (N can be 1) default: []
%   n_col- force number of columns, otherwise adjusted to create a
%          roughly square output
%
% Output: A 2D array suitable for display with e.g. IMAGESC
% 
% This function is primarily used by SHOW_SLICES, but can also be called
% directly.
%
% See also: SHOW_SLICES

% (C) 2005-2012 Andy Adler and Bartlomiej Grychtol
% License: GPL v2 or v3.
% $Id$

if nargin==1 && ischar(rimg) && strcmp(rimg, 'UNIT_TEST'), do_unit_test, return, end

% jnk so that matab doesn't put larger dims in npy
[npy,npx,jnk] = size(rimg);
n_frames = size(rimg,3);
n_levels = size(rimg,4);

if nargin < 2
    sep = 0;
end
if nargin < 3
    vh = [];
end
if nargin < 4 
    n_col = 0;
end
vert_rows = 0;
if nargin > 2 && ~isempty(vh)
    img_cols = n_frames * max( vh(:,1) );
    img_rows = max( vh(:,2) );
else
    % vertical slices must be kept together
    % To nicely fill the image: img_cols ~ img_rows
    % Thus,  n_frames/vert_rows ~ vert_rows*n_levels;
    % or     vert_rows^2 ~ n_frames / n_levels
    vert_rows = ceil( sqrt(n_frames / n_levels) );
    if n_col > 0
        img_cols = n_col;
    else 
        img_cols = ceil( n_frames/vert_rows );
    end
    img_rows = ceil(n_frames*n_levels/img_cols);
    img_rows = ceil(img_rows/n_levels)*n_levels; % Ensure divisible by n_levels
end
% here include the separation
r_img = NaN*ones(img_rows*npy + (img_rows-1)*sep, ...
                 img_cols*npx + (img_cols-1)*sep );

idxx= (-npx:-1)+1;
idxy= (-npy:-1)+1;
imno= 1;
for img_no = 1:n_frames
   for lev_no = 1:n_levels
      if ~isempty(vh) %won't work for multiple image inputs
         i_col= n_frames*(vh( lev_no, 1)-1) + img_no;
         i_row= vh( lev_no, 2);
      else
         i_col= rem( img_no-1, img_cols) + 1;
         i_row= (ceil( img_no / img_cols) -1) * n_levels + lev_no ;
      end
%       disp([imno, vert_rows, img_cols, img_rows, img_no, lev_no, i_col, i_row]);
      r_img(i_row*npy + idxy + sep*(i_row-1), ...
            i_col*npx + idxx + sep*(i_col-1)) = rimg(:,:,img_no,lev_no);
      imno= imno + 1; 
   end
end

function do_unit_test
   img=calc_jacobian_bkgnd(mk_common_model('n3r2',[16,2]));
   img.calc_colours.npoints= 16;
   img.elem_data=toeplitz(1:size(img.fwd_model.elems,1),1);
   img.elem_data = img.elem_data(:,[1,1]);
   levels=[inf,inf,1,1,1;
           0,inf,inf,2,1;
           0,1,inf,  3,1];
   rimg= calc_slices( img, levels(:,1:3) );
   msk = mk_mosaic(rimg, 0, levels(:,4:5));
   subplot(3,1,1)
   imagesc(msk), axis equal tight
   
   msk = mk_mosaic(rimg);
   subplot(3,1,2)
   imagesc(msk), axis equal tight
   
