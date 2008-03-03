function [V] = forward_solver(E,I,tol,pp,V);
%[V] = forward_solver(E,I,tol,pp,V);
%
%This function solves the forward problem using the Cholesky or LU method or 
%conjugate gradients. 
%
%
%
%E   = The full rank system matrix 
%I   = The currents matrix (RHS) 
%pp  = the column permutation vector 
%V   = The approximated nodal potential distribution
%tol = The tolerance in the forward solution, e.g. 1e-5


% d: number of current patterns
[n_nodes,d] = size(I);

% PCG solver is better if memory is tight. This limit, 100k nodes
% seems to be good for machines with ~ 1GB of memory
n_nodes_pcg = 1e5;

if nargin < 6
   V = zeros(size(E,1),d);
end

if isreal(E)==1

   if n_nodes < n_nodes_pcg
       V= E\I;
   else

       K = cholinc(E,tol*100);
       %flags needed for quiet output
       for i=1:d
           [V(:,i),flag] = pcg(E,I(:,i),tol*norm(I(:,i)),n_nodes,K',K,V(:,i));
       end

   end

   % Cholesky solver. Gives poor results matching others
   % so we no longer use it
   if 0 
       %Permute the rows and columns to make the factors sparser
       E = E(pp,pp);
       In = I(pp,:);
       rr(pp)=1:max(size(pp));  % this should be done only Once!
                                % actually much better just to do the
                                % renumbering when the mesh is generated!
       U = chol(E);
       q_c =  U' \ In;  
       Vn = U \ q_c;
       %De-permute the result for Cholesky
       V = Vn(rr,:);
   end

else % is complex
    
   if n_nodes < n_nodes_pcg
      V= E\I;
   else
      [L,U] = luinc(E,tol/10);
      
      for y=1:d
      
         V(:,y) = bicgstab(E,I(:,y),tol*norm(I(:,y)),1000,L,U);
      
      end 
   end 
   
end %is complex




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
