function [mat,grp] = set_inho(srf,simp,vtx,mat_ref,val);
%function [mat,grp] = set_inho(srf,simp,vtx,mat_ref,val);
%
%Auxiliary functions used to set small local inhomogeneities 
%near the boundary of the volume graphically.
%
%
%
%srf     = The boundary (surface) faces
%simp    = The simplices matrix
%vtx     = The vertices matrix
%mat_ref = The reference conductivity vector
%val     = The value to be asigned to the selected inhomogeneity
%mat     = The conductivity vector (updated) that contains the inhomogeneity
%grp     = Subset of the simplices matrix, the tetrahedral that have been selected. 



%figure;
trimesh(srf,vtx(:,1),vtx(:,2),vtx(:,3));
axis image;
set(gcf,'Colormap',[0 0 0]);
hidden off;

mat = mat_ref;

shifter = 1;

beg_wm = simp(1,:);


for dd=1:size(simp,1)
   
    
   if simp(dd,:) ~= beg_wm
      
      shifter = shifter + 1;
   end
   
   if simp(dd,:) == beg_wm
      
      shifter = dd-1;
      
      break;
   end
end


cnts = []; 

for i=1:size(srf,1)
   
   a = srf(i,1);
   b = srf(i,2);
   c = srf(i,3);
   
   ccnx = (vtx(a,1) + vtx(b,1) + vtx(c,1))/3;
   ccny = (vtx(a,2) + vtx(b,2) + vtx(c,2))/3;
   ccnz = (vtx(a,3) + vtx(b,3) + vtx(c,3))/3;
   
   ccn = [ccnx,ccny,ccnz];
   
   cnts = [cnts; ccn];
end


%Number of inho
   
disp('Click on the figure to locate the inhomogeneity');
pause;

   
   [sel] = laserbeam(vtx,srf,cnts); 
   
   %sel is the index of the triangular surface in the srf matrix
   %From that we extract the tet element.
   
   for rr=1:size(simp,1)
      
      if ismember(srf(sel,:),simp(rr,:)) == [1 1 1] 
         
         the_x = simp(rr+shifter,:);
         
         break;
         
      end
   
end

  
   the_nd = setdiff(the_x, srf(sel,:)); %The forth node
   
   grp = [];
   Ssimp = [];
  
   
   for qq=1:size(simp,1)
      
      if sum(ismember(simp(qq,:),the_nd)) == 1
         
         grp = [grp; simp(qq+shifter,:)];
         
         Ssimp = [Ssimp; (qq+shifter)];
         
         mat(qq+shifter) = val;
      end
   end
   
   
    
   for cc=1:size(grp,1)
      
   this_x = grp(cc,:);
   
   
   hold on;
   
   if val > mat_ref(Ssimp(cc))
      
      %Patching f_a
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
   
      if val < mat_ref(Ssimp(cc))
         
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
   
   
   view(3);
   camlight right; 
   lighting flat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%