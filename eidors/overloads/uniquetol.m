function out = uniquetol(in, tol, varargin)
% C = uniquetol(A,TOL):  unique values in A using tolerance TOL.
% For recent versions (>=2015a) call built in function. Otherwise
% find a alternate solve, using either different matlab builtin
% functions, or the uniquetol provided by code (C) Siyi Deng 2010.
%
% Only provides the functions of 
%   C = uniquetol(A,tol,'ByRows',true,'DataScale',1);

% (C) Andy Adler 2015 Licenced under GPL v2 or v3
% $Id$ 

if isstr(in) && strcmp(in,'UNIT_TEST'); do_unit_test; return; end
DEBUG = 0; % IF == ## then only this test

if (~DEBUG && exist('uniquetol','builtin')) || (DEBUG == 1)
   if DEBUG;   disp('using builtin uniquetol'); end
   out = builtin('uniquetol',in, tol, varargin{:});
   return;
end

% Now we provide a backup, but only if byrows and Datascale
if length(varargin)<4
   error('ByRows and DataScale must be provided')
end
for i=1:2:length(varargin);
   switch varargin{i}
     case 'ByRows'
       if ~varargin{i+1};   error('Only support the ByRows option'); end
     case 'DataScale'
       if varargin{i+1}~=1; error('DataScale must by 1'); end
     otherwise 
       error('Option %s not supported',varargin{i});
   end
end

if (~DEBUG && exist('_mergesimpts','builtin')) || (DEBUG == 2)
   if DEBUG;   disp('using _mergesimpts'); end
   out = builtin('_mergesimpts',in,tol,'first');
   return;
end

if DEBUG;   disp('using uniquetol_repl'); end
out = uniquetol_repl(in,tol,'rows','first');
%   UNIQUETOL(...,'ROWS')

function [z,ii,jj] = uniquetol_repl(x,tol,varargin)
%UNIQUETOL Unique element within a tolerance.
%   [Y,I,J] = UNIQUETOL(X,TOL) is very similar to UNIQUE, but allows an
%   additional tolerance input, TOL. TOL can be taken as the total absolute
%   difference between similar elements. TOL must be a none negative
%   scalar. If not provided, TOL is assumed to be 0, which makes UNIQUETOL
%   identical to UNIQUE.
%
%   UNIQUETOL(...,'ROWS')
%   UNIQUETOL(...,'FIRST')
%   UNIQUETOL(...,'LAST')
%   These expressions are identical to the UNIQUE counterparts.
%
%   See also UNIQUE.

% Siyi Deng; 03-19-2010; 05-15-2010; 10-29-2010;

% Licence:
% Copyright (c) 2010, Siyi Deng
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

if size(x,1) == 1, x = x(:); end
if nargin < 2 || isempty(tol) || tol == 0
    [z,ii,jj] = unique(x,varargin{:});
    return;
end
[y,ii,jj] = unique(x,varargin{:});
if size(x,2) > 1
    [~,ord] = sort(sum(x.^2,1),2,'descend');
    [y,io] = sortrows(y,ord);
    [~,jo] = sort(io);
    ii = ii(io);
    jj = jo(jj);
end
d = sum(abs(diff(y,1,1)),2);

isTol = [true;d > tol];
z = y(isTol,:);
bin = cumsum(isTol); % [n,bin] = histc(y,z);
jj = bin(jj);
ii = ii(isTol);

% UNIQUETOL;
   
function do_unit_test
   x = [1 2; 1.06 2; 1.1 2; 1.1 2.03];

   for tol = logspace(-4,1,5);
      uu = uniquetol(x,tol,'ByRows',true,'DataScale',1);
      um = builtin('_mergesimpts',x,tol,'first');
      ur = uniquetol_repl(x,tol,'rows','first');
      unit_test_cmp(sprintf('um=uu %f',tol),uu,um);
      unit_test_cmp(sprintf('um=ur %f',tol),uu,um);
   end
