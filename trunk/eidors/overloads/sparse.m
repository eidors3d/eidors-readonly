function S = sparse(X, varargin)
%SPARSE Create sparse matrix (EIDORS overload).
%   S = SPARSE(X) converts a sparse or full matrix to sparse form by
%   squeezing out any zero elements.
%
%  See also SPARFUN/SPARSE

% (C) 2011 Bartlomiej Grychtol. License: GPL v2 or v3.
% $Id$

S = builtin('sparse',double(X),varargin{:});