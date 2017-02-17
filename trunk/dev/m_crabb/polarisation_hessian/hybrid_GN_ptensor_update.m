function [ dx ] = hybrid_GN_ptensor_update( fmdl, img, DU0, J, W, hp2RtR, dv, de, opt )
%HYBRID_GN_PTENSOR_UPDATE
% NB fix fmdl, img and DU0 in function handle to follow Eidors format

img.elem_data = de+1;

if isfield(opt,'update_U0')

end
if isfield(opt,'neumann')
    neumann = opts.neumann;
else
    neumann = 'freespace';
end

[~,~,CC] = calc_phessian_obj(fmdl,img,DU0,dv,neumann);

if ~any(strcmp(opt.update_method, {'lu','pcg'}))
    error(['unsupported update_method: ',opt.update_method]);
end
if strcmp(opt.update_method, 'lu')
    try
        % the actual update
        warning('needs a modified Cholesky +ve def decomp')
        h_min = 1e-9;
        GN = J'*W*J;
        max_C = diag(GN) - h_min;
        max_C(max_C <= 0) = 0;
        max_C =max_C./CC;
        max_C(isinf(max_C) | max_C < 0 | isnan(max_C)) = inf;
        max_C = min([max_C;1]);
        
        dx = -(GN + hp2RtR + max_C*diag(CC))\(J'*W*dv + hp2RtR*de); % LU decomp
    catch ME % boom
        fprintf('      LU decomp failed: ');
        disp(ME.message)
        fprintf('      try opt.update_method = ''pcg'' on future runs\n');
        opt.update_method = 'pcg';
    end
end
if strcmp(opt.update_method, 'pcg')
    warning('This Hessian may not be +ve definite')
    tol = 1e-6; % default 1e-6
    maxit = []; % default [] --> min(n,20)
    M = []; % default [] --> no preconditioner
    x0 = []; % default [] --> zeros(n,1)
    
    % try Preconditioned Conjugate Gradient: A x = b, solve for x
    % avoids J'*J for n x m matrix with large number of m cols --> J'*J becomes an m x m dense matrix
    LHS = @(x) (J'*(W*(J*x)) + hp2RtR*x + CC(:).*x(:));
    RHS = -(J'*W*dv + hp2RtR*de);
    
    tol=100*eps*size(J,2)^2; % rough estimate based on multiply-accumulates
    %      maxit = 10;
    
    [dx, flag, relres, iter, resvec] = pcg(LHS, RHS, tol, maxit, M', M', x0);
    % TODO if verbose...
    switch flag
        case 0
            if opt.verbose > 1
                fprintf('      PCG: relres=%g < tol=%g @ iter#%d\n',relres,tol,iter);
            end
        case 1
            if opt.verbose > 1
                fprintf('      PCG: relres=%g > tol=%g @ iter#%d (max#%d) [max iter]\n',relres,tol,iter,maxit);
            end
        case 2
            error('error: PCG ill-conditioned preconditioner M');
        case 3
            if opt.verbose > 1
                fprintf('      PCG: relres=%g > tol=%g @ iter#%d (max#%d) [stagnated]\n',relres,tol,iter,maxit);
            end
        case 4
            error('error: PCG a scalar quantity became too large or small to continue');
        otherwise
            error(sprintf('error: PCG unrecognized flag=%d',flag));
    end

end

