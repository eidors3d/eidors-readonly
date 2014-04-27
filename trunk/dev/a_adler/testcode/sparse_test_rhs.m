function sparse_test()
% this tests the speed of sparse matrix multiplies with a high sparsity level

for sz = [16 32 64 128 256 512 1024 2048 4096]
rhs = rand(10e4,1e3);
lhs = sparse(sz,10e4);
lhs(:,end-(sz-1):end) = eye(sz);
%lhs(1,end-(sz/2):end) = -1.5;
fprintf('size(lhs)= %dx%d\n',size(lhs));
lhs=lhs';
rhs=rhs';

tic; for i=1:100
   ans = (rhs * lhs)';              % MATLAB MULTIPLY
end;
tbase = toc;
fprintf('MATLAB MULTIPLY  %fs x1.0\n', tbase);

%tic; for i=1:100
%   ans = rhs(end-(sz-1):end, :)';     % DIRECT SELECT
%end;
%t = toc;
%fprintf('DIRECT SELECT    %fs x%0.1f\n', t, tbase/t);

tic; for i=1:100
   idx= any(lhs);
   ans = (rhs(:,idx)*lhs(idx,:))';  % SELECT WITH IDX
end;
t = toc;
fprintf('SELECT WITH IDX  %fs x%0.1f\n', t, tbase/t);

tic; for i=1:100
   idx= find(any(lhs));
   ans = (rhs(:,idx)*lhs(idx,:))';  % SELECT WITH FIND
end;
t = toc;
fprintf('SELECT WITH FIND %fs x%0.1f\n', t, tbase/t);
end

%TIMES ARE:
%Elapsed time is 0.866917 seconds.
%Elapsed time is 0.005294 seconds.
%Elapsed time is 0.014516 seconds.
%Elapsed time is 0.014434 seconds.
