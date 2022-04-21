function res = unit_test_cmp(txt,a,b,tol)
% UNIT_TEST_CMP: compare matrices in eidors output
% unit_test_cmp(txt,a,b,tol)
% if tol = -negative it is an expected fail
% if tol = -inf then expented fail with no tolerance
% if a==b print ok, otherwise print fail
%
% To run unit tests and keep statistics
%   unit_test_cmp UNIT_TEST_FCN fun_name

% License GPL v2 or v3: $Id$

persistent ntotal;
persistent test_start_time;
persistent npass;
if strcmp(txt,'RESET_COUNTER');
    ntotal=0; npass=0; test_start_time=cputime(); return;
end
if strcmp(txt,'SHOW_COUNTER');
  if ntotal == 0;
     eidors_msg('%s: pass %d/%d (t=%6.2fs)',a, ...
         npass, ntotal, cputime() - test_start_time, 0 );
  else
     eidors_msg('%s: pass %d/%d = %5.1f%% (t=%6.2fs)', a, ...
         npass, ntotal, 100*npass/ntotal, ...
         cputime() - test_start_time, 0 );
  end
  return;
end
if strcmp(txt,'UNIT_TEST_FCN');
  unit_test_cmp('RESET_COUNTER');
  feval(a,'UNIT_TEST');
  unit_test_cmp( 'SHOW_COUNTER',a);
  return
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
   ok='Fail';
   if (isnumeric(a) || islogical(a)) && ...
      (isnumeric(b) || islogical(b))
      sza = size(a); szb= size(b);
      eqsz= isequal( size(a), size(b));
      sza1 = all(sza==1); szb1 = all(szb==1);
      if ~eqsz && ~sza1 && ~szb1
         ok='Fail (size change)';
      elseif any(isnan(a(:)) ~= isnan(b(:)))
             ok='Fail (NaNs change)';
      else
         rel_err = full(max(abs(double(a(:))-double(b(:)))));
         if rel_err <= tol;
            ok='OK';
         end

         if abs(rel_err)>0
            tolstr= sprintf('(%1.3f x tol)', rel_err/tol); 
         end
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
   
   if nargout>0
       res = strcmp(ok(1:2),'OK');
   end

function do_unit_test
   unit_test_cmp('Expect OK'  ,3,3);

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

   unit_test_cmp('Expect Fail', ones(3,3), ones(3,3,3), -inf);
   unit_test_cmp('Expect Fail', ones(3,1), ones(1,3), -inf);
   unit_test_cmp('Expect OK'  ,3,[3,3,3,3]);
   unit_test_cmp('Expect OK'  ,3,[3,3,3,3]);
