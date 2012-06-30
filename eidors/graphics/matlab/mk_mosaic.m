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

% jnk so that matab doesn't put larger dims in npy
[npy,npx,jnk] = size(rimg);
n_frames = size(rimg,3);
n_levels = size(rimg,4);

if nargin < 2
    sep = 0;
end
if nargin < 4 
    n_col = 0;
end
if nargin > 2 && ~isempty(vh)
    img_cols = max( vh(:,1) );
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
         i_col= vh( lev_no, 1) + img_no -1;
         i_row= vh( lev_no, 2);
      else
         i_col= rem( img_no-1, img_cols) + 1;
         i_row= (ceil( img_no / img_cols) -1) * n_levels + lev_no ;
      end
% disp([imno, vert_rows, img_cols, img_rows, img_no, lev_no, i_col, i_row]);
      r_img(i_row*npy + idxy + sep*(i_row-1), ...
            i_col*npx + idxx + sep*(i_col-1)) = rimg(:,:,img_no,lev_no);
      imno= imno + 1; 
   end
end
