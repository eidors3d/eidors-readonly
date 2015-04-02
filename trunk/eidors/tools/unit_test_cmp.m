function unit_test_cmp(txt,a,b,tol)
% UNIT_TEST_CMP: compare matrices in eidors output
% unit_test_cmp(txt,a,b,tol)
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
tolstr='';

   fprintf('TEST: %20s = ',txt);
   ok='fail';
   if (isnumeric(a) || islogical(a)) && ...
      (isnumeric(b) || islogical(b))
      if isnan(a) == isnan(b);
          a(isnan(a))=0; b(isnan(b))=0;
       end;
      if all(abs(double(a) - double(b)) <= tol);
         ok='ok';
      end;

      if tol>0
         tolstr= sprintf('(%1.2f x tol)', full(max(abs(a(:)-b(:))/tol))); 
      end
   else
      if strcmp(eidors_var_id(a), eidors_var_id(b)); ok='ok';end
   end

   fprintf('%4s %s\n', ok, tolstr);
   if strcmp(ok,'ok'); npass= npass+1; end
   ntotal= ntotal+1;

function do_unit_test
   unit_test_cmp('Expect OK'  ,1,1);
   unit_test_cmp('Expect Fail',1,1.01);
   unit_test_cmp('Expect OK'  ,1,.99,.02);
   unit_test_cmp('Expect Fail',1,.99,.002);

   a= rand(10); b = a;
   unit_test_cmp('Expect OK'  ,a,b);
   unit_test_cmp('Expect Fail',a,b+.001);
   unit_test_cmp('Expect OK  ',a,b+.001, .002);

   a(1,1) = NaN; b=a;
   unit_test_cmp('Expect OK'  ,a,b);
   unit_test_cmp('Expect Fail',a,b+.001);
   unit_test_cmp('Expect OK  ',a,b+.001, .002);

   unit_test_cmp('Expect Fail',a,'boo');
   unit_test_cmp('Expect OK','boo','boo');

   t.a= a; s.a=b;
   unit_test_cmp('Expect OK'  ,t,s);
   s2.b= b; 
   unit_test_cmp('Expect Fail'  ,t,s2);
