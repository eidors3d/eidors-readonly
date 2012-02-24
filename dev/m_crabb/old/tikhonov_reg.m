function [x]= tikhonov_reg(A,b,alpha)
%Tikhnonov regularised solution of Ax=b with alpha the regular parameter
%x=A'*(A*A'+a^2*I)^-1*b

B = A*A' + alpha^2*eye(size(A,1)); %(A*A'+a^2*I)
C=B\b; %(A*A'+a^2*I)^-1*b

x=A'*C; %x=A'*(A*A'+a^2*I)^-1*b
