function [V] = forward_solver(E,I,tol,pp,V);
%[V] = forward_solver(E,I,tol,pp,V);
%
%This function solves the forward problem using matlab's \ sovler, or
%conjugate gradients (for large problems). 
%
%E   = The full rank system matrix 
%I   = The currents matrix (RHS) 
%tol = The tolerance in the forward solution, e.g. 1e-5
%pp  = UNUSED
%V   = The approximated nodal potential distribution (USED FOR PCG SOLN)

% (c) N. Polydorides 2003 % Copying permitted under terms of GNU GPL
% $Id$

[n_nodes,n_stims] = size(I);

t=cputime;
try
   V= E\I;
catch 
   [lasterr_str,lasterr_id]= lasterr;
   if lasterr_id ~= 'MATLAB:nomem'
      error(lasterr_str); % rethrow error
   end

   eidors_msg('Memory exhausted for inverse. Trying PCG',2);

   if nargin < 5
      sz= [size(E,1),n_stims];
      V = eidors_obj('get-cache', sz, 'forward_solver_V');
      if isempty(V); V= zeros(sz); end
   end

   if isreal(E)
      U = cholinc(E,tol*100); L = U'; 
      cgsolver = @pcg;
   else %Complex
      [L,U] = luinc(E,tol/10);
      cgsolver = @bicgstab;
   end

   for i=1:n_stims
      [V(:,i),flag] = feval( cgsolver, E,I(:,i), ...
               tol*norm(I(:,i)),n_nodes,L,U,V(:,i));
   end 
      eidors_obj('set-cache', sz, 'forward_solver_V', V);
end
disp(cputime-t);



%%% OLD CODE
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
