function repaint_inho(mat,mat_ref,vtx,simp);
%function repaint_inho(mat,mat_ref,vtx,simp);
%
%Repaints the simulated inhomogeneity according to the reference
%distribution. (Increase -> Red, Decrease -> Blue) 
%
%
%
%mat     = The simulated (targeted) distribution.
%mat_ref = The known initial (homogeneous) distribution.
%vtx     = The vertices matrix.
%simp    = The simplices matrix.

for ii=1:length(mat)
       
   hold on;
   
   this_x = simp(ii,:);
   
   if abs(mat(ii)) > abs(mat_ref(ii))
      
      
      Xas = vtx(this_x(1:3),1);
      Yas = vtx(this_x(1:3),2);
      Zas = vtx(this_x(1:3),3);
      
      patch(Xas,Yas,Zas,'r','EdgeColor','none');
      
      Xbs = vtx(this_x([1,2,4]),1);
      Ybs = vtx(this_x([1,2,4]),2);
      Zbs = vtx(this_x([1,2,4]),3);
      patch(Xbs,Ybs,Zbs,'r','EdgeColor','none');
      
      
      Xcs = vtx(this_x([1,3,4]),1);
      Ycs = vtx(this_x([1,3,4]),2);
      Zcs = vtx(this_x([1,3,4]),3);
      patch(Xcs,Ycs,Zcs,'r','EdgeColor','none');
      
      
      Xds = vtx(this_x(2:4),1);
      Yds = vtx(this_x(2:4),2);
      Zds = vtx(this_x(2:4),3);
      patch(Xds,Yds,Zds,'r','EdgeColor','none');
   
   end
   
      if abs(mat(ii)) < abs(mat_ref(ii))
         
      Xas = vtx(this_x(1:3),1);
      Yas = vtx(this_x(1:3),2);
      Zas = vtx(this_x(1:3),3);
      patch(Xas,Yas,Zas,'b','EdgeColor','none');
      
      Xbs = vtx(this_x([1,2,4]),1);
      Ybs = vtx(this_x([1,2,4]),2);
      Zbs = vtx(this_x([1,2,4]),3);
      patch(Xbs,Ybs,Zbs,'b','EdgeColor','none');
      
      Xcs = vtx(this_x([1,3,4]),1);
      Ycs = vtx(this_x([1,3,4]),2);
      Zcs = vtx(this_x([1,3,4]),3);
      patch(Xcs,Ycs,Zcs,'b','EdgeColor','none');
      
      Xds = vtx(this_x(2:4),1);
      Yds = vtx(this_x(2:4),2);
      Zds = vtx(this_x(2:4),3);
      patch(Xds,Yds,Zds,'b','EdgeColor','none');

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