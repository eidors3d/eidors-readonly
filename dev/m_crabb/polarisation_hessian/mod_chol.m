function [ L, D ] = mod_chol( A, delta, beta )
%MOD_CHOL
%[ L, D ] = mod_chol( A )
% Modified cholesky factorisation
% Produces a lower triangular matrix L and diagonal matrix D such that and
% LDL^T = A + E for symmetric matrices A, where E = 0 for positive definite
% A but is 
%
% Methodby Gill, Murray and Wright
%
% F Watson 2016

if size(A,1)~=size(A,2)
    error('mod_chol for square matrices only');
end

if norm(A - A.') > 1e-12
    error('mod_chol accepts symmetric matrices only')
end

n = size(A,1);

if n==1
    L = sqrt(A);
    D = 1;
    return
end

% Max off-diagonal
Adash = A - diag(diag(A));
xi = max(Adash(:));

% Max-on diagonal
eta = max(abs(diag(A)));

% Tolerances
switch nargin
    case 1
        beta = max([eta, xi/n, eps]);
        delta = max(1e-15, 1e-15 * norm(A,'inf'));
    case 2
        beta = max([eta, xi/n, eps]);
    case 3
        
    otherwise
        error('invalid number of  inputs for mod_chol')
end

% Initialise
d = zeros(1,size(A,1));
L = eye(size(A));
C = zeros(size(A));

% Modified cholesky factorisation
for jj=1:size(A,1)
    
    C(jj,jj) = A(jj,jj) - sum(d(1:jj-1).*(L(jj,1:jj-1).^2)); % Is the square correct?
    
    theta_i = max(abs(C(jj:end,jj)));
    
    d(jj) = max( [abs(C(jj,jj)), (theta_i/beta)^2, delta] );
    
    for ii=jj+1:size(A,1)
        
        C(ii,jj) = A(ii,jj) - sum(d(1:jj-1).*L(ii,1:jj-1).*L(jj,1:jj-1));
        
        L(ii,jj) = C(ii,jj)/d(jj);
        
        
    end
    

end

% % Modified cholesky factorisation -- tried other ordering for timing but
% % no differences/slightly slower. Not verified modification is correct in
% % this version (e.g. theta_i along rows?) but correctly produces
% % non-modified upper triangular factorisation
% for ii=1:size(A,1)
%     
%     C(ii,ii) = A(ii,ii) - sum(d(1:ii-1).'.*(L(1:ii-1,ii).^2));
%     
%     theta_i = max(abs(C(ii,ii:end)));
%     
%     d(ii) = max( [abs(C(ii,ii)), (theta_i/beta)^2, delta] );
%     
%     for jj=ii+1:size(A,1)
%         
%         C(ii,jj) = A(ii,jj) - sum(d(1:ii-1).'.*L(1:ii-1,jj).*L(1:ii-1,ii));
%         
%         L(ii,jj) = C(ii,jj)/d(ii);
%         
%         
%     end
%     
% 
% end


D = diag(d);

end
