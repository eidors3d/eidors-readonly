function [CC] = solution_ext(BB,vtx,simp);
%function [CC] = solution_ext(BB,vtx,simp);
%
%Auxiliary function that extracts a secondary NODE-wise solution CC 
%based on the calculated ELEMENT-wise solution BB. This function is 
%called from slicer_plot function.
%
%
%
%BB   = The calculated (element-wise) inverse solution 
%vtx  = The vertices matrix
%simp = The simplices matrix
 
%Initialise the destination matrix CC
CC = zeros(size(vtx,1),1);

for uu=1:length(BB)
   
   sim_val = BB(uu); % Solution value for simplex number uu
   
   % We aim to distribute this solution value to the nodes involved in the current simplex.
   
   sim_nodes = simp(uu,:);
   
   for oo=1:length(sim_nodes)
      
      s_nd = sim_nodes(oo); % Each of the simplex's nodes 
      
      CC(s_nd) = CC(s_nd) + sim_val;
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