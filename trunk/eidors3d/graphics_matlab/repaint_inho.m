function repaint_inho(mat,mat_ref,vtx,simp, thresh, clim);
%function repaint_inho(mat,mat_ref,vtx,simp, thresh);
%
%Repaints the simulated inhomogeneity according to the reference
%distribution. (Increase -> Red, Decrease -> Blue) 
%
%mat     = The simulated (targeted) distribution.
%mat_ref = The known initial (homogeneous) distribution.
%        = A default value of [] or 'auto' should scale reasonably
%vtx     = The vertices matrix.
%simp    = The simplices matrix.
%thresh  = Threshold to show imaged region (or [] for default)
%clim    = Colour limit

% (C) 2005 Andy Adler + Nick Polydorides. License: GPL version 2 or version 3
% $Id: repaint_inho.m,v 1.27 2007-08-30 03:37:32 aadler Exp $

if nargin<5
    thresh = [];
end
if nargin<6
    clim = [];
end

abs_inhomg= abs( scale_for_display( mat, mat_ref, clim) );
if isempty(thresh)
    thresh = 1/4;
end


ii= find( abs_inhomg > thresh);
   
this_x = simp(ii,:);
   
% looks best if eidors_colours.greylev < 0
colours= calc_colours( mat, clim, 0, mat_ref );
colours= colours(:,ii,:);
ELEM= vtx';

Xs=   zeros(3,length(ii));
Ys=   zeros(3,length(ii));
Zs=   zeros(3,length(ii));
for idx=[[1;2;3], ...
         [1;2;4], ...
         [1;3;4], ...
         [2;3;4]];
   Xs(:)=vtx(this_x(:,idx)',1);
   Ys(:)=vtx(this_x(:,idx)',2);
   Zs(:)=vtx(this_x(:,idx)',3);

   if size(colours,1)==1 & size(colours,2)==3
      % need to work around ^%$#%$# matlab bug which
      % forces an incorrect interpretation is colours of this size
      patch(Xs(:,[1:3,1]), ...
            Ys(:,[1:3,1]), ...
            Zs(:,[1:3,1]), ...
            colours(:,[1:3,1]), ...
            'EdgeColor','none', 'CdataMapping','Direct');
   else
      patch(Xs,Ys,Zs,colours, ...
            'EdgeColor','none', 'CdataMapping','Direct');
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
