function [imgC,m_v,xi_v,sd,m] = mantis_solve(imgC,VM,el,reco)
%MANTIS_SOLVEMANTIS Inexact Gauss-Newton method, see Winkler & Rieder,
% IOP Inverse Problems 31 (2015)
%
% imgC    - EIDORS image model
% VM      - measured potentials/voltages
% el      - struct containing electrode information:
%             el.Npp (number of electrodes per plane),
%             el.P   (number of electrode planes),
%             el.Sizes (electrode surface areas).
% reco    - struct containing information/constants for "MANTIS".
%           If it is empty or fields are missing, default values are used.
% verbose - level of verbosity: 0 = 0ff,
%                               1 = timing,
%                               2 = important info + timing,
%                               3 = detailed info.
%
% (C) 2015 Robert Winkler. License: GPL version 2 or version 3

if nargin < 3; error('Not enough arguments provided.'); end
if nargin < 4 || ~isstruct(reco); clear reco; reco.dummy = true; end

if ~isfield(el,'P'); el.P = 1; end
N = el.Npp*el.P;                                                          % number of electrodes

% some MANTIS constants (default initializations) %%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(reco,'ct');      reco.ct       = 'log';            end        % conductivity transformation as in [W & R 2015]
if ~isfield(reco,'weights'); reco.weights  = 'Ssigma';         end        % specify type of weights as in [W & R 2015]
if ~isfield(reco,'R');       reco.R        = 1.1;              end        % stopping criterion for Morozov's discrepancy principle (no need to modify in general. increase if convergence issues ocurr)
if ~isfield(reco,'sigma');   reco.sigma    = 0.99;             end        % reduction parameter (no need to modify)
if ~isfield(reco,'ximax');   reco.ximax    = sqrt(reco.sigma); end        % REGINN parameter (no need to modify)
if ~isfield(reco,'out_max'); reco.out_max  = 9999;             end        % maximal number of outer iterations (while-loop, no need to modify)
if ~isfield(reco,'inn_max'); reco.inn_max  = (N*(N-1))/2;      end        % maximal number of inner iterations (repeat-loop, no need to modify)
if ~isfield(reco,'cmin');    reco.cmin     = 1e-4;             end        % minimal allowed conductivity
if ~isfield(reco,'cmax');    reco.cmax     = 1/reco.cmin;      end        % maximal allowed conductivity
if ~isfield(reco,'zmin');    reco.zmin     = 1e-5;             end        % minimal allowed contact impedance
if ~isfield(reco,'verbose'); reco.verbose  = 0;                end        % verbosity level (0=quiet ... 3=detailed)
verbose = reco.verbose;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% step 1: estimate noise level from data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORTANT: This ONLY works for a FULL set of measurements per plane,
%            including potentials on driving electrodes. Otherwise, some
%            other redundancy (e.g. reciprocity of single measurements)
%            must be used to estimate the noise level!
noisetime = tic;
reco.noiseEst = mantis_initnoise(imgC.fwd_model.stimulation,VM,el.Npp,...
                                                             el.P,verbose);
noisetime = toc(noisetime);

if verbose >= 2
  fprintf('Estimated RELATIVE noise: %3.2e.\n',reco.noiseEst/norm(VM,'fro'));% rel. noise estimated from data
  fprintf('Estimated ABSOLUTE noise: %3.2e.\n',reco.noiseEst);            % abs. noise estimated from data
end

%% step 2: estimate contact impedances and background conductivity %%%%%%%%
dataC = fwd_solve(imgC);                                                  % simulate reference data

% Extract reference contact impedances as vector out of this freaking struct array.
% Note: These are the impedances used for reference computations, not the true contact impedances.
zref = struct2cell(imgC.fwd_model.electrode);
zref = cell2mat(permute(zref(2,1,:),[3 1 2]));

if isfield(imgC.fwd_model,'electrode') && isfield(imgC.fwd_model.electrode,'nodes')
  nd = Inf;
  for k=1:length(imgC.fwd_model.electrode)
    nd = min(nd,numel(imgC.fwd_model.electrode(k).nodes));
  end
  if nd == Inf || nd <= 10
    warning('Poor electrode discretization (#nodes per el. <= 10). Skipping estimation of contact impedances (inaccurate!).');
    reco.estZ = false;
  end
end

% compute initial guesses for c0 and z
c0ztime = tic;
[c0,z] = mantis_initc0z(VM,dataC.meas,imgC.fwd_model.stimulation,el,zref,reco,verbose);
c0ztime = toc(c0ztime);

if verbose >= 2
  fprintf('Estimated background: %3.2e\n',c0);
  fprintf('Estimated  z_contact: %3.2e\n',z{1});
end

imgC.elem_data = c0*ones(size(imgC.elem_data));                           % set estimated background
[imgC.fwd_model.electrode.z_contact] = deal(z{:});                        % set estimated contact impedances

%% step 3: initial guess to start Gauss-Newton iteration %%%%%%%%%%%%%%%%%%
% some variables to track progress
dataC = fwd_solve(imgC);                                                  % simulate GN initial guess
residual = VM-dataC.meas;                                                 % residual

xi_v  = zeros(reco.out_max,1);                                            % track evolution of regularization parameters
m_v   = zeros(reco.out_max,1);                                            % track evolution of inner iterations
sd    = zeros(reco.out_max+1,1);                                          % track evolution of residual
sd(1) = norm(residual);                                                   % track residual norm (we stop by discrepancy principle!)

inittime = now;

n = 0;                                                                    % outer iterations counter
cgtime = 0;                                                               % track time consumed for cg iterations
jctime = 0;                                                               % track time for computing Jacobians

if verbose >= 2
  fprintf('== REGINN: Inexact GN iteration ==\n');
  fprintf('[GNit]   |residual|   cg-par   #cg\n');
  fprintf('----------------------------------\n');
end

% See also: Rieder 2005, DOI:10.1137/040604029 %

%% step 4: start REGularized INexact Newton (REGINN) method of MANTIS %%%%%
while 1                                                                   % outer REGINN loop
  n = n+1;
  jcstart = tic;
  Jac = calc_jacobian(imgC);                                              % jacobian
  jctime = jctime + toc(jcstart);
  
  %% Apply conductivity transformation before Newton step %%%%%%%%%%%%%%%%%
  ctUpd = mantis_ctts(imgC.elem_data,reco.ct);                            % conductivity transformation
  ww    = mantis_ctdst(ctUpd,reco.ct);                                    % transformation weights (chain rule)
  Jac   = repmat(ww',size(Jac,1),1).*Jac;                                 % overwrite Jacobian with transformed Jacobian
  wUpd  = mantis_ctweights(reco.weights,Jac,imgC.elem_data,ww);           % weights for update

  %% automatically compute cg parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if n == 1
    inmx = 1; xi = 1-eps;                                                 % 1. iteration: force 1 inner iteration (steepest descent, "save route")
  elseif n == 2
    inmx = 2; xi = dcp; reco.ximin = 0.5*dcp;                             % 2. iteration: force 1--2 inner iterations (# inner iterations is regularization. allow only moderate increase each GN step)
  else
    inmx = min(m_v(n-2)+m_v(n-1),reco.inn_max);                           % compute max. number of inner loop iterations to avoid single "spikes"
    if m_v(n-2) < m_v(n-1); xi = 1-(m_v(n-2)/m_v(n-1))*(1-xi_v(n-1));     % increase tolerance
    else xi = reco.sigma*xi_v(n-1); end                                   % decrease tolerance
    xi = max(xi*reco.ximax,reco.ximin);                                   % limit xi (just for convergence theory in Banach spaces...)
  end
  xi_v(n) = xi;                                                           % track cg regularization parameters

  %% call (inner) weighted cg loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  cgstart = tic;
  [s,m_v(n),dcp] = mantis_solvecg(Jac, residual, inmx, wUpd, xi_v(n));
  cgtime = cgtime + toc(cgstart);
  
  %% update conductivity and recompute forward problem %%%%%%%%%%%%%%%%%%%%
  imgC.elem_data = ...                                                    % transform back to conductivity space + update conductivity
    min(max(mantis_ctst(s+ctUpd,reco.ct),reco.cmin),reco.cmax);
  
  eidors_cache('clear_new',inittime);

  dataC = fwd_solve(imgC);                                                % simulate GN initial guess
  residual = VM-dataC.meas;                                               % residual
  
  %% track residual %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  sd(n+1) = norm(residual);
  
  if verbose >= 2
    fprintf('[%4i]    %3.2e     %4.3f   %3i\n',n,sd(n+1),xi_v(n),m_v(n));
  end
  
  %% Check for convergence %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if sd(n+1) <= reco.R*reco.noiseEst
    if verbose >= 2
      fprintf('==================================\n');
      fprintf('Tolerance %3.2e reached.\n',reco.R*reco.noiseEst);
    end
    break;
  elseif n >= reco.out_max
    if verbose >= 2
      fprintf('==================================\n');
    end
    warning('Max. number of outer iterations reached.\n');
    break;
  end  
end
% END REGINN loop of MANTIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% some analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_v  = m_v(1:n);                                                          % crop inner iterations vector
xi_v = xi_v(1:n);                                                         % crop cg parameter vector
sd = sd(1:n+1);                                                           % crop residual norm vector
m = sum(m_v);                                                             % total #cg iterations

if verbose >= 1
  fprintf('MANTIS converged after %i Newton steps, %i overall cg iterations.\n',n,m);
  fprintf('Speed: %6.2fs for computing %3i Jacobians of size %ix%i.\n',jctime,n,size(Jac,1),size(Jac,2));
  fprintf('       %6.2fs for computing %3i cg iterations.\n',cgtime,m);
  fprintf('       %6.2fs for estimating noise level.\n',noisetime);
  fprintf('       %6.2fs for estimating background and contacts.\n',c0ztime);
end
  
end

