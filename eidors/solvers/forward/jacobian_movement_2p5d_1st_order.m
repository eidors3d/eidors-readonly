function J = jacobian_movement_2p5d_1st_order( fwd_model, img)
% JACOBIAN_MOVEMENT_2P5D: J = jacobian_movement_2p5d_1st_order( img )
% Calculate Jacobian Matrix for current stimulation EIT
% J = Jacobian matrix
% img.fwd_model = forward model
% img.elem_data = background for jacobian calculations
%
% fwd_model.normalize_measurements if param exists, calculate
%                                  a Jacobian for normalized
%                                  difference measurements
%
%    img.fwd_solve_2p5d_1st_order.k = [ a .. b ]
%        solve, integrating over the range k = a .. b      (default: [0 Inf])
%        - provide a single k to get a point solution
%        - solve over a reduced range a = 0, b = 3 for a faster solution
%    img.fwd_solve_2p5d_1st_order.method = 'name'
%        perform numerical integration using the selected method (default: 'quadv')
%        'trapz' - trapezoidal integration across the listed points k
%        'quadv' - adaptive quadrature (vectorized), from a to b
%        'integral' - adaptive quadrature (matlab2012+), from a to b


% (C) 2016 A Boyle
% License: GPL version 2 or version 3

% correct input paralemeters if function was called with only img
if ischar(fwd_model) && strcmp(fwd_model,'UNIT_TEST'); do_unit_test; return; end

if nargin == 1
   img = fwd_model;
elseif  strcmp(getfield(warning('query','EIDORS:DeprecatedInterface'),'state'),'on')
   warning('EIDORS:DeprecatedInterface', ...
      ['Calling JACOBIAN_ADJOINT with two arguments is deprecated and will cause' ...
       ' an error in a future version. First argument ignored.']);
end

img = data_mapper(img);
if ~ismember(img.current_params, supported_params)
    error('EIDORS:PhysicsNotSupported', '%s does not support %s', ...
    'JACOBIAN_ADJOINT',img.current_params);
end

org_params = img.current_params;
% all calcs use conductivity
img = convert_img_units(img, 'conductivity');

img.elem_data = check_elem_data_and_apply_c2f(img);
% TODO we don't really want to expand out the c2f

fwd_model = img.fwd_model;

pp = fwd_model_parameters( fwd_model );
gnd = img.fwd_model.gnd_node;

img.fwd_model.system_mat = @system_mat_2p5d_1st_order;
gnd = fwd_model.gnd_node;
img.fwd_model.system_mat_2p5d_1st_order.k = 0;
img.fwd_model.system_mat_2p5d_1st_order.factory = 1;

k = [0 Inf]; try k = img.fwd_model.jacobian_movement_2p5d_1st_order.k; end
method = 'quadv'; try method = img.fwd_model.jacobian_movement_2p5d_1st_order.method; end

pp.Sf = system_mat_2p5d_1st_order( img ); % returns a function Sf(k) that builds a system matrix for 'k'
pp.CC = connectivity_matrix( pp );
% build sparse DD * CC: conductivity * connectivity
lCC = size(pp.CC,1);
DD = kron(img.elem_data,ones(pp.n_dims,1));
DD(end+1:lCC) = 1; % CEM
pp.DD = spdiags(DD, 0, lCC, lCC); % sparse diagonal matrix
pp.DC = pp.DD * pp.CC;
clear DD;
% shape derivatives
[pp.dSx, pp.dTx] = shape_derivatives(pp, img);

% numerical integration
if length(k) == 1 % singleton k
   J = 2*jacobian_k(k, pp, gnd, img.fwd_model.stimulation);
else
   switch method
      case 'trapz'
         % less accurate: trapz
         trace = 0;
         if trace; fprintf('%8s %12s %16s %16s %16s\n','fcnt','a','b-a','||Q||','d||Q||'); end
         n = 0; kil = k(1); nJl=0;
         tol = 1e-8;
         k(isinf(k)) = tol^(-1/6); % 1/k^2 ^3
         Jf = zeros(pp.n_meas, pp.n_elec * pp.n_dims, length(k)); % voltages under electrodes elec x stim, frequency domain
         for ki = k
            n = n + 1;
            Jf(:,:,n) = jacobian_k(ki, pp, gnd, img.fwd_model.stimulation);
            if trace
               nJ = norm(Jf(:,:,n));
               fprintf('%8d     %12e %16e %16e %16e\n',n,ki,ki-kil,nJ,nJ-nJl);
               kil = ki; nJl = nJ;
            end
         end
         J = 2/pi*trapz(k,Jf,3);
         if 0 % draw J(k) to check we are integrating over a large enough range
            Jff = squeeze(reshape(Jf,pp.n_meas*pp.n_elec*pp.n_dims,1,length(k)));
            slope = (Jff(:,3)-Jff(:,2))./max(abs(Jff),[],2);
            clf; plot(slope); drawnow; pause(1);
            Jff_normalized = bsxfun(@rdivide, Jff, Jff(:,1)); % Jff ./ J_0
            clf;semilogx(repmat(k,pp.n_meas*pp.n_elec*pp.n_dims,1)',Jff_normalized');
            xlabel('k'); ylabel('J_k'); title('J_k/J_0')
            drawnow; pause(1);
         end
      case 'quadv'
         % more accurate: adaptive gaussian quadrature
         trace = 0;
         reltol = 1e-4;
         tol = norm(jacobian_k(0, pp, gnd, img.fwd_model.stimulation))*reltol;
         kend = min(tol^(-1/6), max(k)); % don't go too far... k = Inf is a singular matrix, stop adjacent to numeric singularity
         % quadv is scheduled to be removed from matlab eventually... but it is
         % WAY faster than integral with any tolerance configuration I could identify
         % disp('       9     0.0000000000     2.71580000e+15 3286857082845.8784179688');
         if trace; fprintf('%8s %12s %16s %16s\n','fcnt','a','b-a','Q(1)'); end
         J = 2/pi*quadv(@(kk) jacobian_k(kk, pp, gnd, img.fwd_model.stimulation), k(1), kend, tol, trace);
      case 'integral'
         reltol = 1e-4;
         tol = norm(jacobian_k(0, pp, gnd, img.fwd_model.stimulation))*reltol;
         kend = min(tol^(-1/6), max(k)); % don't go too far... k = Inf is a singular matrix, stop adjacent to numeric singularity
         opts = {'ArrayValued', true,
                 'AbsTol', tol, % default: 1e-10
                 'RelTol', reltol}; % default:  1e-6
         opts = opts';
         % the integral solution is about 10x slower (5.47 seconds vs. 0.60 seconds for UNIT_TEST)
         % ... I played with AbsTol and RelTol but wasn't able to affect the outcome
         J = 2/pi*integral(@(kk) jacobian_k(kk, pp, gnd, img.fwd_model.stimulation), k(1), kend, opts{:});
      otherwise
         error(['unrecognized method: ' method]);
   end
end

if ~strcmp(org_params,'conductivity')
    J = apply_chain_rule(J, img, org_params);
    if isfield(fwd_model, 'coarse2fine') && ...
          size(img.elem_data,1) == size(fwd_model.coarse2fine,1)
            J = J*fwd_model.coarse2fine;
            nparam = size(fwd_model.coarse2fine,2);
    end
end

%restore img to original condition
if ~strcmp(org_params,'conductivity') || isfield(img, org_params)
    img = rmfield(img,'elem_data');
    img.current_params = [];
end

% calculate normalized Jacobian
if pp.normalize
   data = fwd_solve( img );
   J = J ./ (data.meas(:)*ones(1,nparam));
end

function J = jacobian_k(k, pp, gnd, stim)
   SS = pp.Sf(k);
   idx = 1:size(SS,1);
   idx( gnd ) = [];
   
   sv = potential_k(SS, pp, gnd);
   zi2E = -(pp.N2E(:,idx)/ SS(idx,idx)) * pp.CC(:,idx)';
   SE = jacobian_calc(pp, zi2E, pp.dSx, k, pp.dTx, sv);
   
   idx = 0;
   J = zeros( pp.n_meas, pp.n_elec * pp.n_dims );
   for j = 1:pp.n_stim
      meas_pat = stim(j).meas_pattern;
      m = size(meas_pat,1); % new measurements added
      SEj = reshape( SE(:,j,:), pp.n_elec, pp.n_elec * pp.n_dims);
      J( idx+(1:m),: ) = meas_pat * SEj;
      idx = idx + m;
   end

function v = potential_k(S, pp, gnd)
   idx = 1:size(S,1);
   idx( gnd ) = [];
   idx2 = find(any(pp.N2E));
   v = zeros(pp.n_node, pp.n_stim); % voltages at all nodes x number of stim, frequency domain
   v(idx,:) = left_divide( S(idx,idx), 1/2*pp.QQ(idx,:) );

% SE_{i,j,k} is dV_i,j / dS_k
%  where V_i is change in voltage on electrode i for
%        stimulation pattern j
%        S_k is k=n*ne+e for change in position on electrode e of ne, for movement axis n=(0,1)
function SE = jacobian_calc(pp, zi2E, dSx, k, dTx, sv)
   N = pp.n_elec*pp.n_dims;
   SE = zeros(pp.n_elec, pp.n_stim, N);
   for n = 1:N;
      dSxk2dTx = dSx{n} + k^2 * dTx{n};
      SE(:,:,n) = (zi2E * (dSxk2dTx * pp.DC)) * sv;
   end

function J = apply_chain_rule(J, img, org_params)
   switch(org_params)
      case 'resistivity'
         dCond_dPhys = -img.elem_data.^2;
      case 'log_resistivity'
         dCond_dPhys = -img.elem_data;
      case 'log_conductivity'
         dCond_dPhys = img.elem_data;
      otherwise
         error('not implemented yet')
   end
   J = J.*repmat(dCond_dPhys ,1,size(J,1))';

function elem_data = check_elem_data_and_apply_c2f(img)
   elem_data = img.elem_data;
   sz_elem_data = size(elem_data);
   if sz_elem_data(2) ~= 1;
      error('jacobian_adjoin: can only solve one image (sz_elem_data = %)', ...
            sz_elem_data);
   end

   if isfield(img.fwd_model, 'coarse2fine');
      c2f = img.fwd_model.coarse2fine;
      sz_c2f = size(c2f);
      switch sz_elem_data(1)
         case sz_c2f(1); % Ok
         case sz_c2f(2); elem_data = c2f * elem_data;
         otherwise; error(['jacobian_adjoint: provided elem_data ' ...
              ' (sz = %d) does not match c2f (sz = %d %d)'], sz_elem_data(1), sz_c2f);
      end
   else
      if sz_elem_data(1) ~= num_elems(img.fwd_model)
         error(['jacobian_adjoint: provided elem_data (sz = %d) does ' ...
            ' not match fwd_model (sz = %d)'], sz_elem_data(1), num_elems(sz_c2f));
      end
   end

function str = supported_params
    str = {'conductivity'
           'resistivity'
           'log_conductivity'
           'log_resistivity'};

% Define the element connectivity matrix Ce
function Ce = connectivity_matrix( pp, img );
   m = (pp.n_dims+1)*pp.n_elem;
   n = size(pp.N2E,2);
   ii = 1:m; % row indices
   jj = pp.ELEM(:); % col indices
   ss = ones(m,1); % values = 1
   Ce = sparse(ii, jj, ss, m, n); % m x n sparse matrix

function [dS, dT] = shape_derivatives( pp, img )
   d1 = pp.n_dims+1; % 2d --> always 3
   sz = d1 * pp.n_elem; % dS and dT matrices should be 'sz' rows x cols
   Ec = cell(pp.n_elem,1);
   for elec = 1:length(img.fwd_model.electrode)
      ii = []; jj = []; ss = []; tt = [];
      for node = img.fwd_model.electrode(elec).nodes(:)';
         % find elements touching that node
         [ idx, elem ] = find(pp.ELEM == node);
         for j = 1:length(elem(:))
            e = elem(j);
            % build shape matrices
            if ~any(size(Ec{e})) % already calculated?
               Ec{e} = inv([ ones(d1,1), pp.NODE( :, pp.ELEM(:,e) )' ]); %%% (12) [Boyle2016]
            end
            E = Ec{e};
            E1 = E(2:end,:); % E_/1
            absdetE = 0.5/pp.VOLUME(e); % | det E | = 1 / (n_dim! * element volume) for n_dim=2
            for d = 1:pp.n_dims
               u = zeros(d1,1); v = u; % init
               u(idx(j)) = 1; % element node# to perturb
               v(d+1) = 1; % direction for perturbation: x,z?
               B = (v'*E*u) / absdetE; %%% (13) [Boyle2016]
               dE1 = -E*u*v'*E; dE1 = dE1(2:end,:); %%% (14) [Boyle2016]
               dSe = (v'*E*u*E1'*E1 + dE1'*E1 + E1'*dE1)*pp.VOLUME(e); %%% eqn (11) [Boyle2016]
%               dTe = -absdetE * (v'*E*u) * (ones(d1)+eye(d1))/12; %%% eqn (20) [Boyle2016]
               dTe = (v'*E*u)/(2*absdetE) * (ones(d1)+eye(d1))/12; %%% eqn (20) [Boyle2016]
               dSs(:,d) = dSe(:); % square matrix to column
               dTs(:,d) = dTe(:);
               if 0 % test
                  dfact = factorial(pp.n_dims);
                  Ae = pp.NODE(:,pp.ELEM(:, e))'; Ae = [ones(pp.n_dims+1,1), Ae]; Ae = inv(Ae);
                  Be = Ae(2:pp.n_dims+1,:);
                  absdetAe = abs(det(Ae));
                  u = zeros(pp.n_dims+1,1); v = u; u(idx(j)) = 1; v(d+1) = 1;
                  delBe = -Ae*u*v'*Ae; delBe = delBe(2:pp.n_dims+1,:);
                  dSe_old = 1/dfact*(1/absdetAe*v'*Ae*u*Be'*Be + delBe'*Be/absdetAe + Be'*delBe/absdetAe);
                  %%%
                  Ae = pp.NODE(:,pp.ELEM(:, e))'; Ae = [ones(pp.n_dims+1,1), Ae]; Ae = inv(Ae);
                  Be = Ae(2:pp.n_dims+1,:);
                  Se = Be'*Be/dfact/abs(det(Ae));
                  Te = 1/(2*abs(det(Ae)))*(ones(d1)+eye(d1))/12;
                  delta = 1e-8; d_NODE = pp.NODE; d_NODE(d,node) = d_NODE(d,node) + delta;
                  Ae = d_NODE(:,pp.ELEM(:, e))'; Ae = [ones(pp.n_dims+1,1), Ae]; Ae = inv(Ae);
                  Be = Ae(2:pp.n_dims+1,:);
                  Se_delta = Be'*Be/dfact/abs(det(Ae));
                  Te_delta = 1/(2*abs(det(Ae)))*(ones(d1)+eye(d1))/12;
                  dSe_pert = (Se_delta - Se) / delta;
                  dTe_pert = (Te_delta - Te) / delta;

                  if norm(dSe_pert - dSe,1)/norm(dSe_pert) > 5e-6
                     eidors_msg(sprintf('@@@ dSe: elec#%d, elem#%d (%d of %d), %s-axis: calc wrong', ...
                                        elec, e, j, length(elem(:)), 'w'+d), ...
                                1);
                     error('stop');
                  end
                  if norm(dTe_pert - dTe)/norm(dTe_pert) > 5e-6
                     eidors_msg(sprintf('@@@ dTe: elec#%d, elem#%d (%d of %d), %s-axis: calc wrong', ...
                                        elec, e, j, length(elem(:)), 'w'+d), ...
                                1);
                     error('stop');
                  end
               end
            end
            se_idx = (1:d1)+(e-1)*d1;
            [xx,yy] = meshgrid(se_idx, se_idx);
            ii(end+1:end+d1^2) = xx(:);
            jj(end+1:end+d1^2) = yy(:);
            ss(end+1:end+d1^2,:) = dSs;
            tt(end+1:end+d1^2,:) = dTs;
            % note that for any element, we usually have more than one node
            % that is perturbed; all ss (and, separately, tt) for that element
            % are summed by the sparse() operation so that we have a sparse
            % block-wise matrix with 3x3 blocks on the diagonal in 2D
         end
      end
      % return a dS and dT cell array with a cell per dimension
      for d = 1:pp.n_dims
         idx = (d-1)*pp.n_elec + elec;
         dS{idx} = sparse(ii,jj,ss(:,d),sz,sz);
         dT{idx} = sparse(ii,jj,tt(:,d),sz,sz);
         if 0 % test dS against perturbation
            elec_nodes = img.fwd_model.electrode(elec).nodes(:);
            dS_old = calc_delVm( elec_nodes, pp, img.fwd_model, img, d, 1, 1);
            clf; spy(pp.CC'*dS_old*pp.CC,'go'); hold on; spy(pp.CC'*dS{idx}*pp.CC,'rx'); hold off;
            if norm(dS{idx} - dS_old, 1) > 2e-5
               eidors_msg(sprintf('@@@ dS{%d}: elec#%d, %s-axis movement: calc wrong', ...
                                  idx, ...
                                  mod(idx-1,pp.n_elec)+1, ...
                                  'x'+floor((idx-1)/pp.n_elec)), ...
                          1);
               error('stop');
            end
         end
      end
   end

%%% TODO rm start (2d movement Jacobian)
function delVm=  calc_delVm( elec_nodes_array, pp, fwd_model, img, colidx, Re_Ce, cond_Ce_Vc)
   I = []; J=[]; S= [];
   pp.dfact = factorial(pp.n_dims);
   for elec_nodes= elec_nodes_array(:)';
      [rowidx, elemidx] = find(pp.ELEM == elec_nodes);
      % Define the system sensitivity matrix to movement delSm
      sz= (pp.n_dims+1)*pp.n_elem;
      %delSm = sparse(sz,sz);
      % For each touching element, calculate the perturbation
      jcount = 1;
      for j = elemidx'
          % Extract the coordinates of the element's four nodes
          Ae = pp.NODE(:,pp.ELEM(:, j))';
          % Define the invertible matrix P: augment Ae by adding a
          % column of ones to invert
          P = [ones(pp.n_dims+1,1), Ae];
          Ae = inv(P);
          absdetAe = abs(det(Ae));
          % Define Be as the matrix Ae with row 1 deleted
          Be = Ae(2:pp.n_dims+1,:);
          % For this coordinate, perturb P by [rowidx,colidx], which are
          % our paper's perturbation vectors [a,b]
          a = zeros(pp.n_dims+1,1);
          b = a;
          a(rowidx(jcount)) = 1;
          jcount = jcount + 1;
          b(colidx+1) = 1;
          % Calculate the system submatrix subSm for the element j by
          % asymmetric perturbation of the electrode node k
          deldetAe =   1/absdetAe*b'*Ae*a;
          delBe = -Ae*a*b'*Ae;
          delBe = delBe(2:pp.n_dims+1,:);
          subSm = 1/pp.dfact*(...
              deldetAe*Be'*Be + ...
              delBe'*Be/absdetAe + ...
              Be'*delBe/absdetAe);
      
          % Embed subSm into delSm such that subSm(1,1) is the
          % (4j+1,4j+1) element of delSm
          se_idx= (pp.n_dims+1)*j+(-pp.n_dims : 0);
          switch pp.n_dims
             case 2
                Iidx = vertcat(se_idx,se_idx,se_idx);
                I = [I Iidx(:)];
                J = [J, se_idx,se_idx,se_idx];
                
             case 3
                Iidx = vertcat(se_idx,se_idx,se_idx,se_idx);
                I = [I Iidx(:)];
                J = [J, se_idx,se_idx,se_idx,se_idx];
          end
          S = [S subSm(:)];
      end
   end
   delSm = sparse(I,J,S,sz,sz);
   delVm = Re_Ce * delSm * cond_Ce_Vc;

function do_unit_test()
   nn = 16; % number of electrodes
   imdl2 = mk_geophysics_model('h22c',nn);
   img2 = mk_image(imdl2,1);
   nm = size(stim_meas_list(img2.fwd_model.stimulation),1); % number of measurements
   ne = size(imdl2.rec_model.elems,1); % number of elements in coarse model
   img2.fwd_model.conductivity_jacobian = zeros(nm,ne);

   % for the 3d model, we throw out the rec_model and inject the
   % imdl2.rec_model, then recalculate the c2f so that we can compare apples to
   % apples when we take ||Jxxx - J3||
   imdl3 = mk_geophysics_model('h32c',nn);
   for s = {'nodes', 'elems', 'boundary', 'name','electrode'}
      imdl3.rec_model.(s{:}) = imdl2.rec_model.(s{:});
   end
   disp('recalculating 3D c2f, so that coarse meshes agree (2D vs 3D)');
   [c2f,bkgnd] = mk_approx_c2f(imdl3.fwd_model, imdl3.rec_model);
   imdl3.fwd_model.coarse2fine = c2f;
   imdl3.fwd_model.background = bkgnd;
   img3 = mk_image(imdl3,1);
   img3.fwd_model.conductivity_jacobian = zeros(nm,ne);

   % confirm c2f is sane
   ctr = interp_mesh(imdl2.rec_model);
   sel = ((ctr(:,1)-100).^2 < 80) & ((ctr(:,2)+50).^2<80);
   img2.elem_data = imdl2.fwd_model.coarse2fine*(1 + sel*9) + -10*imdl2.fwd_model.background;
   img3.elem_data = imdl3.fwd_model.coarse2fine*(1 + sel*9) + -10*imdl3.fwd_model.background;
   clf;
   subplot(121);show_fem(img2,1);
   subplot(122);show_fem(img3,1);
   drawnow;
   % set data to homogeneous
   img2.elem_data = imdl2.fwd_model.background*0+1;
   img3.elem_data = imdl3.fwd_model.background*0+1;

   t = tic;
   img2.fwd_model.nodes(1,:) = img2.fwd_model.nodes(1,:) + rand(1,2)*1e-8; % defeat cache
   J2p = jacobian_movement_perturb(img2);
   J2p(:,1:(end-2*nn)) = [];
   fprintf(' 2D perturb = %.2f sec\n', toc(t));
   t = tic;
   img2.fwd_model.nodes(1,:) = img2.fwd_model.nodes(1,:) + rand(1,2)*1e-8; % defeat cache
   J2a = jacobian_movement(img2);
   J2a(:,1:(end-2*nn)) = [];
   fprintf(' 2D adjoint = %.2f sec\n', toc(t));

   t = tic;
   img2.fwd_model.nodes(1,:) = img2.fwd_model.nodes(1,:) + rand(1,2)*1e-8; % defeat cache
   img2.fwd_model.jacobian = @jacobian_movement_2p5d_1st_order;
   img2.fwd_model.jacobian_movement_2p5d_1st_order.k = 0;
   J2p50 = jacobian_movement_2p5d_1st_order(img2);
   J2p50(:,1:(end-2*nn)) = [];
   fprintf(' 2.5D (k = 0) = %.2f sec\n', toc(t));

   ke = [ 0 0 0 ];
   t = tic;
   img2.fwd_model.nodes(1,:) = img2.fwd_model.nodes(1,:) + rand(1,2)*1e-8; % defeat cache
   img2.fwd_model.jacobian = @jacobian_movement_2p5d_1st_order;
   img2.fwd_model.jacobian_movement_2p5d_1st_order.k = [0 logspace(-3,+1,20)];
   img2.fwd_model.jacobian_movement_2p5d_1st_order.method = 'trapz';
   J2p5kt = jacobian_movement_2p5d_1st_order(img2);
   J2p5kt(:,1:(end-2*nn)) = [];
   ke(1) = img2.fwd_model.jacobian_movement_2p5d_1st_order.k(end);
   fprintf(' 2.5D (k = 0..%.1f, trapz) = %.2f sec\n', ke(1), toc(t));
   t = tic;
   img2.fwd_model.nodes(1,:) = img2.fwd_model.nodes(1,:) + rand(1,2)*1e-8; % defeat cache
   img2.fwd_model.jacobian = @jacobian_movement_2p5d_1st_order;
   img2.fwd_model.jacobian_movement_2p5d_1st_order.k = [0 Inf];
   img2.fwd_model.jacobian_movement_2p5d_1st_order.method = 'quadv';
   J2p5kq = jacobian_movement_2p5d_1st_order(img2);
   J2p5kq(:,1:(end-2*nn)) = [];
   ke(2) = img2.fwd_model.jacobian_movement_2p5d_1st_order.k(end);
   fprintf(' 2.5D (k = 0..%.1f, quadv) = %.2f sec\n', ke(2), toc(t));
   t = tic;
   img2.fwd_model.nodes(1,:) = img2.fwd_model.nodes(1,:) + rand(1,2)*1e-8; % defeat cache
   img2.fwd_model.jacobian = @jacobian_movement_2p5d_1st_order;
   img2.fwd_model.jacobian_movement_2p5d_1st_order.k = [0 Inf];
   img2.fwd_model.jacobian_movement_2p5d_1st_order.method = 'integral';
   J2p5ki = jacobian_movement_2p5d_1st_order(img2);
   J2p5ki(:,1:(end-2*nn)) = [];
   ke(3) = img2.fwd_model.jacobian_movement_2p5d_1st_order.k(end);
   fprintf(' 2.5D (k = 0..%.1f, integral) = %.2f sec\n', ke(3), toc(t));

   t = tic;
   img3.fwd_model.nodes(1,:) = img3.fwd_model.nodes(1,:) + rand(1,3)*1e-8; % defeat cache
   J3a = jacobian_movement(img3);
   J3a(:,1:(end-3*nn)) = []; J3a(:,(1*nn+1):(2*nn)) = []; % delete y-axis
   fprintf(' 3D adjoint = %.2f sec\n', toc(t));
   if 0 % skip really slow calculations
      t = tic;
      img3.fwd_model.nodes(1,:) = img3.fwd_model.nodes(1,:) + rand(1,3)*1e-8; % defeat cache
      J3p = jacobian_movement_perturb(img3);
      J3p(:,1:(end-3*nn)) = []; J3p(:,(1*nn+1):(2*nn)) = []; % delete y-axis
      fprintf(' 3D perturb = %.2f sec\n', toc(t));
   else
      J3p = J3a;
      fprintf(' 3D perturb = <SKIP>\n');
   end


   tol = 1e-8;
   reltol = 1e2;
   fprintf('tol = %0.1e, reltol = %0.1e, ||J_0|| = %0.1e, ||J|| = %0.1e\n', tol, reltol,norm(J2a),norm(J3a));
   unit_test_cmp('2D perturb                 vs 3D', max(max(abs(J2p-J3p)))>tol,1);
   unit_test_cmp('2D adjoint                 vs 3D', max(max(abs(J2a-J3a)))>tol,1);
   unit_test_cmp('2D adjoint vs 2D perturb        ', J2a, J2p,  tol*10);
   unit_test_cmp('3D adjoint vs 3D perturb        ', J3a, J3p,  tol);
   unit_test_cmp('2.5D (k = 0)                  vs 2D',  J2p50, J2a, tol);
   unit_test_cmp('2.5D (k = 0) matches          == 2D',  max(max(abs(J2p50-J2a)))<tol,1);
   unit_test_cmp('2.5D (k = 0) does not match   != 3D',  max(max(abs(J2p50-J3a)))>tol,1);
   unit_test_cmp(sprintf('2.5D (k = 0..%-4.1f) (trapz)    != 2D',ke(1)), max(max(abs(J2a./J2p5kt-1)))>reltol, 1);
   unit_test_cmp(sprintf('2.5D (k = 0..%-4.1f) (trapz)    == 3D',ke(1)), max(max(abs(J3a./J2p5kt-1)))<reltol, 1);
   unit_test_cmp(sprintf('2.5D (k = 0..%-4.1f) (trapz)    vs 3D',ke(1)), J3a./J2p5kt, 1, reltol);
   unit_test_cmp(sprintf('2.5D (k = 0..%-4.1f) (quadv)    vs 3D',ke(2)), J3a./J2p5kq, 1, reltol);
   unit_test_cmp(sprintf('2.5D (k = 0..%-4.1f) (integral) vs 3D',ke(3)), J3a./J2p5ki, 1, reltol);
