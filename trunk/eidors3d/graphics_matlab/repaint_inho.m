function repaint_inho(mat,mat_ref,vtx,simp);
%function repaint_inho(mat,mat_ref,vtx,simp);
%
%Repaints the simulated inhomogeneity according to the reference
%distribution. (Increase -> Red, Decrease -> Blue) 
%
%mat     = The simulated (targeted) distribution.
%mat_ref = The known initial (homogeneous) distribution.
%vtx     = The vertices matrix.
%simp    = The simplices matrix.
%
% $Id: repaint_inho.m,v 1.2 2004-07-16 17:06:42 aadler Exp $

inhomg= mat - mat_ref;
hold('on');

for ii= find( inhomg(:)' )
   
   this_x = simp(ii,:);
   
   if inhomg(ii) > 0
      colour= 'r';
   elseif inhomg(ii) < 0;
      colour= 'b';
   else
      error('shouldn''t get here');
   end
      
      
      Xas = vtx(this_x(1:3),1);
      Yas = vtx(this_x(1:3),2);
      Zas = vtx(this_x(1:3),3);
      
      patch(Xas,Yas,Zas,colour,'EdgeColor','none');
      
      Xbs = vtx(this_x([1,2,4]),1);
      Ybs = vtx(this_x([1,2,4]),2);
      Zbs = vtx(this_x([1,2,4]),3);
      patch(Xbs,Ybs,Zbs,colour,'EdgeColor','none');
      
      
      Xcs = vtx(this_x([1,3,4]),1);
      Ycs = vtx(this_x([1,3,4]),2);
      Zcs = vtx(this_x([1,3,4]),3);
      patch(Xcs,Ycs,Zcs,colour,'EdgeColor','none');
      
      
      Xds = vtx(this_x(2:4),1);
      Yds = vtx(this_x(2:4),2);
      Zds = vtx(this_x(2:4),3);
      patch(Xds,Yds,Zds,colour,'EdgeColor','none');
      
end
hold('off')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
