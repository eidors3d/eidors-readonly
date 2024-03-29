function [V] = left_divide(E,I,tol,pp,V)
%[V] = left_divide(E,I,tol,pp,V);
%
% Implements left division and overcomes many small inefficiencies
% of matlab's left division. Also uses conjugate gradients (for large problems).
%
%E   = The full rank system matrix
%I   = The currents matrix (RHS)
%tol = The tolerance in the forward solution, e.g. 1e-5
%
% pp,V are old options from previous solver

% (c) N. Polydorides 2003 % Copying permitted under terms of GNU GPL
% $Id$

if ~exist('tol'); tol = 1e-8; end

[n_nodes,n_stims] = size(I);

try
    % V= E\I;
    % This takes MUCH longer when you have  more vectors in I,
    %  even if they are repeated. There must be some way to simplify
    %  this to speed it up. Matlab's sparse operators really should
    %  do this for you.
    
    % TODO, Q&R should be cached somewhere
    if ~isequal(E,E.')
        
        F=abs(E-E.');
        
        % is the matrix symmetric up to tollerance
        issym=all(F(F>0)<tol*(abs(E(F>0))));
        
        if issym
            % symmeterise
            F=0.5*(E'+E);
            
            [Q,R] = qr(I,0);
            rnotzeros = any(R~=0,2);
            Q= Q(:,rnotzeros);
            R= R(rnotzeros,:);
            
            V= (F \ Q)*R;
        else
            clearvars F issym
            
            [Q,R] = qr(I,0);
            rnotzeros = any(R~=0,2);
            Q= Q(:,rnotzeros);
            R= R(rnotzeros,:);
            V= (E \ Q)*R;
        end
    else
        
        [Q,R] = qr(I,0);
        rnotzeros = any(R~=0,2);
        Q= Q(:,rnotzeros);
        R= R(rnotzeros,:);
        V= (E \ Q)*R;
        
    end
    % TODO: Iteratively refine
    %  From GH Scott: "once we have
    %   computed the approximate solution x, we perform one step
    %   of iterative refinement by computing the residual: r = Ax - b
    %   and then recalling the solve routine to solve
    %   Adx = r for the correction dx.
    % However, we don't want to repeat the '\', so we implement
    %   the underlying algorithm:
    %   If A is sparse, then MATLAB software uses CHOLMOD (after 7.2) to compute X.
    %    The computations result in  P'*A*P = R'*R
    %   where P is a permutation matrix generated by amd, and R is
    %   an upper triangular matrix. In this case, X = P*(R\(R'\(P'*B)))
    %
    % See also:
    % http://www.cs.berkeley.edu/~wkahan/MxMulEps.pdf
    % especially page 15 where it discusses the value of iterative refinement
    %  without extra precision bits.  ALso, we need to enable
    
    
catch
    [lasterr_str,lasterr_id]= lasterr;
    if ~strcmp(lasterr_id , 'MATLAB:nomem')
        error(lasterr_str); % rethrow error
    end
    
    eidors_msg('Memory exhausted for inverse. Trying PCG',2);
    
    if nargin < 5
        sz= [size(E,1),n_stims];
        V = eidors_obj('get-cache', sz, 'left_divide_V');
        if isempty(V); V= zeros(sz); end
    end
    
    ver = eidors_obj('interpreter_version'); % Matlab2013 renamed cholinc -> ichol
    if isreal(E)
        opts.droptol = tol*100;
        opt.type = 'ict';
        if ver.isoctave || ver.ver < 7.012
            U = cholinc(E, opt.droptol);
        else
            U = ichol(E, opt);
        end
        L = U';
        cgsolver = @pcg;
    else %Complex
        opts.droptol = tol/10;
        if ver.isoctave || ver.ver < 7.012 % Matlab2007 introduced ilu, luinc has now been dropped
            [L,U] = luinc(E, opt.droptol);
        else
            [L,U] = ilu(E, opt);
        end
        cgsolver = @bicgstab;
    end
    
    for i=1:n_stims
        [V(:,i),flag] = feval( cgsolver, E,I(:,i), ...
            tol*norm(I(:,i)),n_nodes,L,U,V(:,i));
    end
    eidors_obj('set-cache', sz, 'left_divide_V', V);
end


%%% OLD CODE
% Cholesky solver. Gives poor results matching others
% so we no longer use it
if 0
    %Permute the rows and columns to make the factors sparser
    E = E(pp,pp);
    In = I(pp,:);
    rr(pp)=1:max(size(pp));  % this should be done only Once!
    % actually much better just to do the
    % renumbering when the mesh is generated!
    U = chol(E);
    q_c =  U' \ In;
    Vn = U \ q_c;
    %De-permute the result for Cholesky
    V = Vn(rr,:);
end
