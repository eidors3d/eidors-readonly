function rs=primaldual_tvrecon_lsearch(inv_mdl, vmeas, ...
              maxiter,alpha1,alpha2, epsilon, beta)
% rs=tvrecon(msh,c,vmeas,maxiter,alpha1,alpha2)

% 11/02/01 By Andrea Borsic
% 01/02/02 Modified: new steplength search on dual variable
% function rs=tvrecon(msh,c,vmeas,maxiter,alpha1,alpha2)
%
% ex. rs=pd_tvrecon(msh,c,vmeas,10,5e-3,1e-9); for simulated data
% ex. rs=pd_tvrecon(msh,c,vmeas,10,5e-3,5e-9); for tank data
%
% PARAMETERS
%   alpha1 - hyperparameter for the initial Tikhonov step
%   alpha2 - hyperparameter for the TV update
%       A good alpha11 is 2e-5, a good alpha2 is 1e-9
%   epsilon - termination tolerance (recommend 1e-6);
%   beta    - initial value of approx: abs ~ sqrt(sqr+beta)
%
% PRIMAL DUAL IMPLEMENTATION
%
% Total variation reconstruction with norm-2 fitting of the measurements
%
% As reference see.. "A non linear primal-dual method for total variation-based image restoration" T.F. Chan, G.H. Golub, P. Mulet
% Note about notation: A is A' for us, and y is s, C is B, beta=mu^2
%

% (C) 2002-2006 Andrea Borsic. License: GPL version 2 or version 3
% $Id: primaldual_tvrecon_lsearch.m,v 1.15 2008-03-13 19:28:04 aadler Exp $

% Initialisation
fwd_model= inv_mdl.fwd_model;

msh.TC = fwd_model.elems';
msh.PC = fwd_model.nodes';



decay_beta=0.7;

% this is used for the line search procedure,
% the last element, and the biggest must be one
len=([0,1e-4,1e-3,1e-2,0.1,0.2,0.5,0.8,1]);

% Inizialisation
A=calc_R_prior( inv_mdl);

n=size(A,1); % num_rows_L
m=size(A,2); % num_elem

% Create homogeneous model
IM= eidors_obj('image','');
IM.fwd_model= fwd_model;

if 0 % static EIT - this code doesn't work yet
   scaling=vmeas\vsim;
   s=s*scaling;
   %u=potentials(msh,s,c);
   %vsim=measures(msh,u);
   vsim= sim_measures( IM, s);
   de_v=vmeas-vsim;
   %J=jacobian(msh,u);
   IM.s= s;
   J= calc_jacobian( IM );
   de_s=[J;alpha1*A]\[de_v;zeros(n,1)];
   s=s+de_s;
else
   s= zeros(m,1);
   x=zeros(n,1);
   J= calc_jacobian( calc_jacobian_bkgnd(inv_mdl) );
   IM.elem_data= s;
   IM.reconst_type = 'difference';
   IM.J = J;
   vsim= sim_measures( IM, s);
   de_v=vmeas-vsim;
   de_s=[J;alpha1*A]\[de_v;zeros(n,1)];
   s=s+de_s;
   scaling= 1;
end

%dispmsh(msh,s); colorbar; drawnow;

% Iterative procedure

terminate=0;
iter=1;
ind=1;
    rs(:,iter)=s;

while (~terminate)&(iter<maxiter)
    
%   u=potentials(msh,s,c);
%   vsim=measures(msh,u); 
%   J=jacobian(msh,u);  - Jacobian is same in difference EIT
    vsim= sim_measures( IM, s);
%   J= calc_jacobian( IM );
%   plot([vsim, vmeas]);% pause

    z=A*s;	% This is the dual variable
    
    eta=sqrt(z.^2+beta);
    
    grad=J'*(vsim-vmeas);
    
    for i=1:n
        
        grad=grad+alpha2*A(i,:)'*A(i,:)*s/eta(i);
        
    end % for
    
    primal=sum(abs(z)); % we don't care here about 0.5*norm(de_v)
    dual=sum(x.*z);
    
    eidors_msg('PDIPM: %2d & %1.3e & %1.3e & %1.3e & %1.3e & %1.3e & %1.3e & %1.3e & ',iter,primal,dual,primal-dual,norm(vsim-vmeas),beta,len(ind),norm(grad),3);
        
    E=spdiags(eta,0,n,n);
    F=spdiags(ones(n,1)-(1./eta).*x.*z,0,n,n);
    
    B=(J'*J+alpha2*A'*inv(E)*F*A);
    
    %B=0.5*(B+B');
    
    de_s=-B\(alpha2*A'*inv(E)*z+J'*(vsim-vmeas));
    
    ang=acos((dot(de_s,-grad)/(norm(de_s)*norm(-grad))))*(360/(2*pi));
    
    fprintf('+'); eidors_msg('PDIPM angle=%+3.1f deg',ang,3);
    
    % line search
    
    for k=1:length(len)
%       meas_k = measures(msh,s+len(k)*de_s,c);
        meas_k= sim_measures( IM, s+len(k)*de_s);
        func_val(k)=0.5*norm( meas_k - vmeas )^2 + ...
                    alpha2*sum(abs(A*(s+len(k)*de_s)));
    end % for
    
    [temp,ind]=min(func_val);% disp(len(ind));
    
    % conductivity update
        
    s=s+len(ind)*de_s;
               
    de_x=-x+inv(E)*z+inv(E)*F*A*de_s;
    
    % dual step length rule
    
    lims=sign(de_x); % this are the limits (+1 or -1) torward x is pushed is de_x is added to it
    
    clearance=lims-x; % this is the signed distance to the limits
    
 % this protects against division, other values wil dominate, it doesn't affect the algorithm
    de_x(de_x==0)=1e-6;
        
 % stemps that will make on compunent of x exceed the limits of step*de_x is applied
    steps=clearance./de_x;
    
    idx= steps==0;
    steps(idx)=de_x(idx);
    
    % we need to pick up the smallest, and have some safety room
       
    x=x+min(1,0.99*min(steps))*de_x;
        
    % Upper and lower limits enforcement
    s( s<0.01*scaling )=0.01*scaling;
    % and upper bounds, dynamic range=1e4
    s( s>100*scaling )=100*scaling;
    
    beta=beta*decay_beta;       % beta is reduced
    decay_beta=decay_beta*0.8;  % the rate at wich beta is reduced is also adjusted
    if beta<1e-12
        beta=1e-12;
    end
    
%    if norm(A*x)>gap(x,z) error('Rounding errors are spoiling the calculation, stopping.'); end % if
    
    % The primal-dual gap has been reduced and measures match
    if (sum(abs(z)-x.*z)<epsilon)&(norm(vsim-vmeas)<epsilon)
        terminate=1;
    end % if
    
    iter=iter+1;

    rs(:,iter)=s;
    
end % while

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function vsim= sim_measures( IM, s);
   if IM.reconst_type == 'difference'
      vsim = IM.J*s;
   else
      vh= fwd_solve( IM );
      IM.elem_data= s;
      vi= fwd_solve( IM );

      vsim = calc_difference_data( vh, vi, IM.fwd_model);
   end
