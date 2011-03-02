function unit_test_cmp(txt,a,b,tol)
% UNIT_TEST_CMP: compare matrices in eidors output
% unit_test_cmp(txt,a,b,tol)
% if a==b print ok, otherwise print fail

% License GPL v2 or v3: $Id$

   if nargin < 4; tol = 0; end
   fprintf('TEST: %20s = ',txt);
   ok='fail';
   try; if isnan(a) == isnan(b); a(isnan(a))=0; b(isnan(b))=0; end; end
   try; if all(abs(a - b) <= tol);  ok='ok'; end; catch; end
   disp(ok)
   if ~strcmp(ok,'ok'); keyboard; end
