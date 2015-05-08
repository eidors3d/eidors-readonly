function unit_test_cmp(txt,a,b,tol)
% UNIT_TEST_CMP: compare matrices in eidors output
% unit_test_cmp(txt,a,b,tol)
% if tol = -negative it is an expected fail
% if tol = -inf then expented fail with no tolerance
% if a==b print ok, otherwise print fail

% License GPL v2 or v3: $Id$

persistent ntotal;
persistent npass;
if strcmp(txt,'RESET_COUNTER'); ntotal=0; npass=0; return; end
if strcmp(txt,'SHOW_COUNTER');
  eidors_msg('%s: pass %d/%d',a, npass, ntotal,0); return;
end

if strcmp(txt,'UNIT_TEST'); do_unit_test; return; end


if nargin < 4; tol = 0; end
if tol<0;
   expect_fail = 1;
   if tol==-inf; tol= 0; end
else
   expect_fail= 0;
end
tolstr='';

   fprintf('TEST: %20s = ',txt);
   ok='fail';
   if (isnumeric(a) || islogical(a)) && ...
      (isnumeric(b) || islogical(b))
      if isnan(a) == isnan(b);
          a(isnan(a))=0; b(isnan(b))=0;
       end;
      if all(abs(double(a) - double(b)) <= tol);
         ok='OK';
      end;

      if abs(tol)>0
         tolstr= sprintf('(%1.2f x tol)', full(max(abs(a(:)-b(:))/tol))); 
      end
   else
      if strcmp(eidors_var_id(a), eidors_var_id(b)); ok='OK';end
   end

   if expect_fail
      ok = 'OK (fail as expected)';
   end

   fprintf('%4s %s\n', ok, tolstr);
   if strcmp(ok(1:2),'OK'); npass= npass+1; end
   ntotal= ntotal+1;

function do_unit_test
   unit_test_cmp('Expect OK'  ,1,1);
   unit_test_cmp('Expect Fail',1,1.01, -inf);
   unit_test_cmp('Expect OK'  ,1,.99,.02);
   unit_test_cmp('Expect Fail',1,.99,-.002);

   a= rand(10); b = a;
   unit_test_cmp('Expect OK'  ,a,b);
   unit_test_cmp('Expect Fail',a,b+.001, -inf);
   unit_test_cmp('Expect OK  ',a,b+.001, .002);

   a(1,1) = NaN; b=a;
   unit_test_cmp('Expect OK'  ,a,b);
   unit_test_cmp('Expect Fail',a,b+.001, -inf);
   unit_test_cmp('Expect OK  ',a,b+.001, .002);

   unit_test_cmp('Expect Fail',a,'boo', -inf);
   unit_test_cmp('Expect OK','boo','boo');

   t.a= a; s.a=b;
   unit_test_cmp('Expect OK'  ,t,s);
   s2.b= b; 
   unit_test_cmp('Expect Fail'  ,t,s2, -inf);
