function [f,m,dcp] = mantis_solvecg(A, g, m_max, a, xi)
%MANTIS_SOLVECG inner iteration of the REGINN algorithm. This cg iteration
% approximates the solution to Af=g, where A is the Jacobian
% (cf. Rieder 2003, Lechleiter and Rieder 2006, Winkler and Rieder 2015)
%
% (C) 2015 Robert Winkler. License: GPL version 2 or version 3

%%  norm definitions
if size(a,1)==size(a,2) % we got a full (possibly non-diagonal) weight matrix
  fullMat = true;
  as = sqrtm(full(a));
  nX = @(x) norm(as*x);
  nY = @(y) norm(y);
  AH = a\transpose(A);
else % weight matrix is given by its diagonal entries (faster)
  fullMat = false;
  as = sqrt(a(:));
  nX = @(x) norm(as.*x);                       % weighted norm in X
  nY = @(y) norm(y);                           % standard norm in Y
  %AH = repmat(1./a,1,size(A,1)).*transpose(A); % adjoint matrix in weighted space
  ai = 1./a;
end

%%  cg initialization, cf. Rieder 2003, p. 126
f = zeros(size(A,2),1); % initial guess
r = g;                  % r = g-Af; initial guess f = 0!!
if fullMat; d = AH*r; else d = transpose((transpose(r)*A)).*ai; end
p = d;
m = 0;

nd = nX(d);     % norm of the defect
nr0 = nY(r);    % norm of the resitual at the beginning
nr = nr0;       % norm of the iteration residual
scond = nr0*xi; % stopping residual

%%  start cg iteration (cf. Rieder 2003, p. 126).
%This is the repeat-loop of REGINN
while (nr >= scond && m < m_max)
  m = m+1;
  q = A*p;
  
  alpha = (nd/nY(q))^2;
  
  f = f+alpha*p;
  r = r-alpha*q;
  
  if fullMat; d = AH*r; else d = transpose((transpose(r)*A)).*ai; end
    
  nd_old = nd;
  
  nd = nX(d);
  nr = nY(r);
  
  beta = (nd/nd_old)^2;
  
  p = d+beta*p;
  
end % end of while loop

dcp = nr/nr0; % relative discrepancy decrease that was actually achieved in this step


end
