function S = spdiag(V,K)
%SPDIAG Sparse diagonal matrices and diagonals of a matrix.
%
% This is a sparse version of matlab's diag command, see the
% help of that command for its use.
%
% See also: DIAG, SPDIAGS
%
% (C) 2013 Andy Adler. License GPL v2 or v3

% Matlab's spdiags is really annoying. This function provides
%  an interface like G-d intended

if nargin == 1; K=0; end

N = length(V) + K;
V = [zeros(K,1);V(:)];
S = spdiags(V(:), K, N, N);








