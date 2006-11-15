function [Reg] = iso_f_smooth(simp,vtx,deg,w);
%function [Reg] = iso_f_smooth(simp,vtx,deg,w);
%
%Calculates a first order discrete Gaussian smoothing operator.
%
%
%
%simp = The simplices matrix.
%vtx  = The vertices matrix.
%deg  = 1 for nodes, 2 for edges and 3 for faces
%w    = smoothing weight, w=1...k, default value = 1
%Reg  = The first order smoothing regulariser. 

if nargin<2 
   w=1;
end
if w<0
   error('Weight must be possitive');
end

ndsrch = [simp,(1:size(simp,1))'];

Reg = spalloc(size(simp,1),size(simp,1),20*size(simp,1));

for i=1:size(ndsrch,1)
   
   t_id = ndsrch(i,1:size(simp,2));
   
   Xsimp = []; %The vector of simp indices that share deg nodes with simp(i)
   XX = [];
   
 for j = 1:length(t_id)
      
      t_nd = t_id(j);
      
      X1 = find(simp(:,1)==t_nd);
      X2 = find(simp(:,2)==t_nd);
      X3 = find(simp(:,3)==t_nd);
      X4 = find(simp(:,4)==t_nd);
      
      %The set of all indices containing the node t_nd
      XX = [XX;X1;X2;X3;X4];
      
  end %for j
      
      if deg == 1 %Neigbouring nodes
      Xu = unique(XX);
      Xsimp = [Xsimp;Xu];
      end
      
      if deg == 2 %Neihgouring edges
      Xs = sort(XX);
      Xd = diff(Xs);
      Xf = Xs(setdiff(1:length(Xd),find(Xd)));
      Xu = unique(Xf);
      Xsimp = [Xsimp;Xu];
      end
      
      if deg == 3 %Neihbouring faces
      Xs = sort(XX);
      Xd = diff(Xs);
      Xf = Xs(setdiff(1:length(Xd),find(Xd)));
      Xq = diff(Xf);
      Xz = Xf(setdiff(1:length(Xq),find(Xq)));
      Xu = unique(Xz);
      Xsimp = [Xsimp;Xu];
      end
      
      if deg > 3 | deg < 1 
          error('deg parameter can only be 1, 2 or 3')
      end
            
   %Remove the i'th simplex from the list of its neighbours
   Xdif = setdiff(Xsimp,i);
   
   if deg == 1
   Reg(i,Xdif) = -1/w;
   end

   if deg == 2
   for p=1:length(Xdif)
        Intr = intersect(t_id,simp(Xdif(p),:)); % 2 or 3 long vector
        dd=0;
       for h=1:length(Intr)-1
          [da] = db23d(vtx(Intr(h),1),vtx(Intr(h),2),vtx(Intr(h),3),...
                       vtx(Intr(h+1),1),vtx(Intr(h+1),2),vtx(Intr(h+1),3));
          dd = dd+da;
      end
   Reg(i,Xdif(p)) = -w/dd; 
   end
   end

   if deg == 3
    for p=1:length(Xdif)
       Intr = intersect(t_id,simp(Xdif(p),:)); % 3 long vector
       [ta] = triarea3d(vtx(Intr,:));
    end
   Reg(i,Xdif) = -w/ta;
   end
        
   Reg(i,i) = abs(sum(Reg(i,:)));
   
end %for i'th simplex



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

