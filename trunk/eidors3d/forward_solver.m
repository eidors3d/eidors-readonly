function [V] = forward_solver(vtx,E,I,tol,pp,V);
%function [V] = forward_solver(vtx,E,I,tol,pp,V);
%
%This function solves the forward problem using the Cholesky or LU method or 
%conjugate gradients. 
%
%
%
%vtx = The vertices
%E   = The full rank system matrix 
%I   = The currents matrix (RHS) 
%pp  = the column permutation vector 
%V   = The approximated nodal potential distribution
%tol = The tolerance in the forward solution, e.g. 1e-5




% d: number of current patterns
[bla,d] = size(I);

if nargin < 6
V = zeros(size(E,1),d);
end

if isreal(E)==1

    if  pp ~= 1:size(I,1) %There is a colume permutation, hence Cholesky opted
        %Permute the rows and columns to make the factors sparser
        E = E(pp,pp);
        In = I(pp,:);
        rr(pp)=1:max(size(pp));
        U = cholinc(E,tol/10);
        q_c =  U' \ In;
        Vn = U \ q_c;
        %De-permute the result for Cholesky
        V = Vn(rr,:);
    else

        %Alternatively use pcg ********
        K = cholinc(E,tol*100);
        for i=1:d
            [V(:,i),flag,relres,iter,resvec] = pcg(E,I(:,i),tol*norm(I(:,i)),d^3,K',K,V(:,i));
        end

    end

end %is real

if isreal(E)==0
    
  [L,U] = luinc(E,tol/10);

  for y=1:d

  [V(:,y),flag,relres,iter,resvec] = bicgstab(E,I(:,y),tol*norm(I(:,y)),1000,L,U);

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