function repaint_inho(mat,mat_ref,vtx,simp, thresh);
%function repaint_inho(mat,mat_ref,vtx,simp, thresh);
%
%Repaints the simulated inhomogeneity according to the reference
%distribution. (Increase -> Red, Decrease -> Blue) 
%
%mat     = The simulated (targeted) distribution.
%mat_ref = The known initial (homogeneous) distribution.
%vtx     = The vertices matrix.
%simp    = The simplices matrix.
%thresh  = Threshold to show imaged region
%
% $Id: repaint_inho.m,v 1.4 2005-10-14 15:55:35 aadler Exp $

inhomg= mat - mat_ref;

if nargin<5
    thresh = max(abs(inhomg(:)))/4;
end

ii= find( abs(inhomg(:)') > thresh);
   
this_x = simp(ii,:);
   
colour= calc_colours( inhomg(ii) );
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

   patch(Xs,Ys,Zs,colour,'EdgeColor','none');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
