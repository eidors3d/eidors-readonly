function J= jacobian_adjoint_2p5d_1st_order( fwd_model, img)
% JACOBIAN_ADJOINT_2P5D: J= jacobian_adjoint_2p5d_1st_order( img )
% Calculate Jacobian Matrix for current stimulation EIT
% J         = Jacobian matrix
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
%        - solve over a reduced range a=0, b=3 for a faster solution
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
   img= fwd_model;
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

img.elem_data = check_elem_data(img);

fwd_model= img.fwd_model;

pp= fwd_model_parameters( fwd_model, 'skip_VOLUME' );
gnd = img.fwd_model.gnd_node;

img.fwd_model.system_mat = @system_mat_2p5d_1st_order;
gnd = fwd_model.gnd_node;
img.fwd_model.system_mat_2p5d_1st_order.k = 0;
img.fwd_model.system_mat_2p5d_1st_order.factory = 1;

pp.Sf = system_mat_2p5d_1st_order( img ); % returns a function Sf(k) that builds a system matrix for 'k'
pp.FC= system_mat_fields( fwd_model );
pp.FT= system_mat_2p5d_fields( fwd_model );

c2f = speye(pp.n_elem); try c2f = img.fwd_model.coarse2fine; end % make code uniform by assigning 1:1 c2f if no c2f provided
k = [0 Inf]; try k = img.fwd_model.jacobian_adjoint_2p5d_1st_order.k; end
method = 'quadv'; try method = img.fwd_model.jacobian_adjoint_2p5d_1st_order.method; end

if length(k) == 1 % singleton k
   J = 2*jacobian_k(k, pp, gnd, img.fwd_model.stimulation, c2f);
else
   switch method
      case 'trapz'
         % less accurate: trapz
         trace = 0;
         if trace; fprintf('%8s %12s %16s %16s\n','fcnt','a','b-a','||Q||'); end
         n = 0; kil = k(1);
         tol = 1e-8;
         k(isinf(k)) = tol^(-1/6);
         Jf = zeros(pp.n_meas, size(c2f,2), length(k)); % voltages under electrodes elec x stim, frequency domain
         for ki = k
            n = n + 1;
            Jf(:,:,n) = jacobian_k(ki, pp, gnd, img.fwd_model.stimulation, c2f);
            if trace; fprintf('%8d     %12e %16e %16e\n',n,ki,ki-kil,norm(Jf(1,1,n))); kil = ki; end
         end
         J=2/pi*trapz(k,Jf,3);
         % check
         Jff = squeeze(reshape(Jf,pp.n_meas*size(c2f,2),1,length(k)));
         assert(max(abs(Jff(:,end)) < tol), sprintf('trapz k=%e truncated too early as k->Inf',k(end)));
         slope = (Jff(:,3)-Jff(:,2))./Jff(:,2); slope(Jff(:,2) == 0) = nan;
         %clf; plot(slope); drawnow;
         %[idx]=find(abs(slope) > median(abs(slope))+sqrt(tol)); clf;semilogx(k,Jff(idx,:)); drawnow;
         assert(median(abs(slope)) < sqrt(tol)*10, sprintf('trapz k->0 (tol=%e), slope=%e != 0',tol,median(abs(slope))));
         if 0 % draw J(k) to check we are integrating over a large enough range
            clf;semilogx(repmat(k,pp.n_meas*size(c2f,2),1)',Jff');
            xlabel('k'); ylabel('J_k'); title('J_k')
            drawnow; pause(1);
         end
         clear Jff;
         clear slope;
      case 'quadv'
         % more accurate: adaptive gaussian quadrature
         trace = 0;
         tol = norm(jacobian_k(0, pp, gnd, img.fwd_model.stimulation, 1))*1e-2;
         kend = min(tol^(-1/6), max(k)); % don't go too far... k=Inf is a singular matrix, stop adjacent to numeric singularity
         % quadv is scheduled to be removed from matlab eventually... but it is
         % WAY faster than integral with any tolerance configuration I could identify
         if trace; fprintf('%8s %12s %16s %16s\n','fcnt','a','b-a','Q(1)'); end
         J=2/pi*quadv(@(kk) jacobian_k(kk, pp, gnd, img.fwd_model.stimulation, c2f), k(1), kend, tol, trace);
      case 'integral'
         reltol = 1e-8;
         tol = norm(jacobian_k(0, pp, gnd, img.fwd_model.stimulation, 1))*reltol;
         kend = min(tol^(-1/6), max(k)); % don't go too far... k=Inf is a singular matrix, stop adjacent to numeric singularity
         opts = {'ArrayValued', true,
                 'AbsTol', tol, % default: 1e-10
                 'RelTol',reltol}; % default:  1e-6
         opts = opts';
         % the integral solution is about 10x slower (5.47 seconds vs. 0.60 seconds for UNIT_TEST)
         % ... I played with AbsTol and RelTol but wasn't able to affect the outcome
         J=2/pi*integral(@(kk) jacobian_k(kk, pp, gnd, img.fwd_model.stimulation, c2f), k(1), kend, 'ArrayVAlued', true);
      otherwise
         error(['unrecognized method: ' method]);
   end
end

if ~strcmp(org_params,'conductivity')
    J = apply_chain_rule(J, img, org_params);
    if isfield(fwd_model, 'coarse2fine') && ...
          size(img.elem_data,1)==size(fwd_model.coarse2fine,1)
            J=J*fwd_model.coarse2fine;
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
   data= fwd_solve( img );
   J= J ./ (data.meas(:)*ones(1,nparam));

end

function J = jacobian_k(k, pp, gnd, stim, c2f)
   SS = pp.Sf(k);
   idx = 1:size(SS,1);
   idx( gnd ) = [];
   
   sv = potential_k(SS, pp, gnd);
   
   zi2E= zeros(pp.n_elec, pp.n_node);
   zi2E(:, idx)= -pp.N2E(:,idx)/ SS(idx,idx);
   
   % NOTE: this is k*FT here, not k^2*FT: FT is squared in jacobian_calc (FT*DE*FT)
   DE = jacobian_calc(pp, zi2E, pp.FC, k*pp.FT, sv, c2f);
   
   J = zeros( pp.n_meas, size(c2f,2) );
   idx=0;
   for j= 1:pp.n_stim
      meas_pat= stim(j).meas_pattern;
      n_meas  = size(meas_pat,1);
      DEj = reshape( DE(:,j,:), pp.n_elec, size(c2f,2) );
      J( idx+(1:n_meas),: ) = meas_pat*DEj;
      idx= idx+ n_meas;
   end

function v = potential_k(S, pp, gnd)
   idx = 1:size(S,1);
   idx( gnd ) = [];
   idx2 = find(any(pp.N2E));
   v = zeros(pp.n_node, pp.n_stim); % voltages at all nodes x number of stim, frequency domain
   v(idx,:) = left_divide( S(idx,idx), 1/2*pp.QQ(idx,:) );

% DE_{i,j,k} is dV_i,j / dS_k
%  where V_i is change in voltage on electrode i for
%        stimulation pattern j
%        S_k is change in conductivity on element k
function DE = jacobian_calc(pp, zi2E, FC, FT, sv, c2f)
   d0 = pp.n_dims + 0;
   d1 = pp.n_dims + 1;
   
   zi2E_FCt = zi2E * FC';
   zi2E_FTt = zi2E * FT';
   
   FC_sv = FC * sv;
   FT_sv = FT * sv;

   DE= zeros(pp.n_elec, pp.n_stim, size(c2f,2) );
   de0= pp.n_elem * d0;
   de1= pp.n_elem * d1;
   if all(all(speye(size(c2f)) == c2f))
      for k= 1:size(c2f,2);
         idx = d0*(k-1)+1 : d0*k;
         dq1= zi2E_FCt(:,idx) * FC_sv(idx,:);
         dq2= zi2E_FTt(:,idx) * FT_sv(idx,:);
         DE(:,:,k)= dq1 + dq2;
      end
   elseif 0 % Code is slower
      for k= 1:size(c2f,2);
         chg_col = kron( c2f(:,k), ones(d0,1));
         dDD_dEj = spdiags(chg_col, 0, de0, de0);
         dq1 = zi2E_FCt * dDD_dEj * FC_sv;
         chg_col = kron( c2f(:,k), ones(d1,1));
         dDD_dEj = spdiags(chg_col, 0, de1, de1);
         dq2 = zi2E_FTt * dDD_dEj * FT_sv;
         DE(:,:,k)= dq1 + dq2;
      end
   else
      for k= 1:size(c2f,2);
         ff = find( c2f(:,k) > 1e-8 );
         lff= length(ff)*d0;
         ff1= ones(d0,1) * ff(:)';
         ffd= d0*ff1 + (-(d0-1):0)'*ones(1,length(ff));
         dDD_dEj = spdiags(c2f(ff1,k), 0, lff, lff);
         dq1= zi2E_FCt(:,ffd) * dDD_dEj * FC_sv(ffd,:);
         lff= length(ff)*d1;
         ff1= ones(d1,1) * ff(:)';
         ffd= d1*ff1 + (-(d1-1):0)'*ones(1,length(ff));
         dDD_dEj = spdiags(c2f(ff1,k), 0, lff, lff);
         dq2= zi2E_FTt(:,ffd) * dDD_dEj * FT_sv(ffd,:);
         DE(:,:,k)= dq1 + dq2;
      end
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

function elem_data = check_elem_data(img)
   elem_data = img.elem_data;
   sz_elem_data = size(elem_data);
   if sz_elem_data(2) ~= 1;
      error('jacobian_adjoin: can only solve one image (sz_elem_data=%)', ...
            sz_elem_data);
   end

   if isfield(img.fwd_model, 'coarse2fine');
      c2f = img.fwd_model.coarse2fine;
      sz_c2f = size(c2f);
      switch sz_elem_data(1)
         case sz_c2f(1); % Ok
         case sz_c2f(2); elem_data = c2f * elem_data;
         otherwise; error(['jacobian_adjoint: provided elem_data ' ...
              ' (sz=%d) does not match c2f (sz=%d %d)'], sz_elem_data(1), sz_c2f);
      end
   else
      if sz_elem_data(1) ~= num_elems(img.fwd_model)
         error(['jacobian_adjoint: provided elem_data (sz=%d) does ' ...
            ' not match fwd_model (sz=%d)'], sz_elem_data(1), num_elems(sz_c2f));
      end
   end

function str = supported_params
    str = {'conductivity'
           'resistivity'
           'log_conductivity'
           'log_resistivity'};

function do_unit_test()
   imdl2 = mk_geophysics_model('h22c',16);
   assert(length(imdl2.rec_model.elems) > 20, 'expect sufficient rec_model density');
   img2 = mk_image(imdl2,1);
   if 0 % check model c2f
      clf; subplot(121); show_fem(imdl2.fwd_model); subplot(122); show_fem(imdl2.rec_model);
      c2f=imdl2.fwd_model.coarse2fine; img2r = mk_image(imdl2); img2r.fwd_model = imdl2.rec_model; img2r.elem_data = zeros(size(c2f,2),1);
      for i=1:size(c2f,2)
         subplot(121);
         img2.elem_data = c2f(:,i); show_fem(img2); title(sprintf('%d',i));
         subplot(122);
         img2r.elem_data(:) = 0; img2r.elem_data(i)=1; show_fem(img2r);  title(sprintf('%d',i));
         drawnow; pause(0.25);
      end
   end

   % for the 3d model, we throw out the rec_model and inject the
   % imdl2.rec_model, then recalculate the c2f so that we can compare apples to
   % apples when we take ||Jxxx - J3||
   imdl3 = mk_geophysics_model('h32c',16);
   for s = {'nodes', 'elems', 'boundary', 'name','electrode'}
      imdl3.rec_model.(s{:}) = imdl2.rec_model.(s{:});
   end
   disp('recalculating 3D c2f, so that coarse meshes agree (2D vs 3D)');
   [c2f,bkgnd] = mk_approx_c2f(imdl3.fwd_model, imdl3.rec_model);
   imdl3.fwd_model.coarse2fine = c2f;
   imdl3.fwd_model.background = bkgnd;
   img3 = mk_image(imdl3,1);
   if 0 % check model c2f
      clf; subplot(121); show_fem(imdl3.fwd_model); subplot(122); show_fem(imdl3.rec_model);
      c2f=imdl3.fwd_model.coarse2fine; img3r = mk_image(imdl3); img3r.fwd_model = imdl3.rec_model; img3r.elem_data = zeros(size(c2f,2),1);
      for i=1:size(c2f,2)
         subplot(121);
         img3.elem_data = c2f(:,i); show_fem(img3); title(sprintf('%d',i));
         subplot(122);
         img3r.elem_data(:) = 0; img3r.elem_data(i)=1; show_fem(img3r);  title(sprintf('%d',i));
         drawnow; pause(0.25);
      end
   end

   t = tic;
   img2.fwd_model.nodes(1,:) = img2.fwd_model.nodes(1,:) + rand(1,2)*1e-8; % defeat cache
   J2 = calc_jacobian(img2);
   fprintf(' 2D                         = %.2f sec\n', toc(t));

   t = tic;
   img2.fwd_model.nodes(1,:) = img2.fwd_model.nodes(1,:) + rand(1,2)*1e-8; % defeat cache
   img2.fwd_model.jacobian = @jacobian_adjoint_2p5d_1st_order;
   img2.fwd_model.jacobian_adjoint_2p5d_1st_order.k = 0;
   J2p50 = calc_jacobian(img2);
   fprintf(' 2.5D (k=0)                 = %.2f sec\n', toc(t));

   ke = [0 0 0];
   t = tic;
   img2.fwd_model.nodes(1,:) = img2.fwd_model.nodes(1,:) + rand(1,2)*1e-8; % defeat cache
   img2.fwd_model.jacobian = @jacobian_adjoint_2p5d_1st_order;
   img2.fwd_model.jacobian_adjoint_2p5d_1st_order.k = [0 logspace(-4,1,100)]; % capture all the effects to (k^2 T)!
   img2.fwd_model.jacobian_adjoint_2p5d_1st_order.method = 'trapz';
   J2p5kt = calc_jacobian(img2);
   ke(1) = img2.fwd_model.jacobian_adjoint_2p5d_1st_order.k(end);
   fprintf(' 2.5D (k=0..%.1f, trapz)    = %.2f sec\n', ke(1), toc(t));
   t = tic;
   img2.fwd_model.nodes(1,:) = img2.fwd_model.nodes(1,:) + rand(1,2)*1e-8; % defeat cache
   img2.fwd_model.jacobian = @jacobian_adjoint_2p5d_1st_order;
   img2.fwd_model.jacobian_adjoint_2p5d_1st_order.k = [0 Inf];
   img2.fwd_model.jacobian_adjoint_2p5d_1st_order.method = 'quadv';
   J2p5kq = calc_jacobian(img2);
   ke(2) = img2.fwd_model.jacobian_adjoint_2p5d_1st_order.k(end);
   fprintf(' 2.5D (k=0..%.1f, quadv)    = %.2f sec\n', ke(2), toc(t));
   t = tic;
   img2.fwd_model.nodes(1,:) = img2.fwd_model.nodes(1,:) + rand(1,2)*1e-8; % defeat cache
   img2.fwd_model.jacobian = @jacobian_adjoint_2p5d_1st_order;
   img2.fwd_model.jacobian_adjoint_2p5d_1st_order.k = [0 Inf];
   img2.fwd_model.jacobian_adjoint_2p5d_1st_order.method = 'integral';
   J2p5ki = calc_jacobian(img2);
   ke(3) = img2.fwd_model.jacobian_adjoint_2p5d_1st_order.k(end);
   fprintf(' 2.5D (k=0..%.1f, integral) = %.2f sec\n', ke(3), toc(t));

   t = tic;
   img3.fwd_model.nodes(1,:) = img3.fwd_model.nodes(1,:) + rand(1,3)*1e-8; % defeat cache
   J3 = calc_jacobian(img3);
   fprintf(' 3D                         = %.2f sec\n', toc(t));


   tol = 1e-8;
   reltol = norm(J3)*2e-2;
   unit_test_cmp('2D                          vs 3D', J2, J3, -tol);
   unit_test_cmp('2.5D (k=0)                  vs 2D', J2p50, J2, tol);
   unit_test_cmp(sprintf('2.5D (k=0..%-4.1f) (trapz)    vs 3D',ke(1)), J2p5kt, J3, 2*reltol);
   unit_test_cmp(sprintf('2.5D (k=0..%-4.1f) (quadv)    vs 3D',ke(2)), J2p5kq, J3, reltol);
   unit_test_cmp(sprintf('2.5D (k=0..%-4.1f) (integral) vs 3D',ke(3)), J2p5ki, J3, reltol);

   imgr = img2;
   imgr.fwd_model = imdl2.rec_model;
   clf; subplot(221); imgr.elem_data = sens(J3); show_fem(imgr,1); title('3D [log_{10}]');
        subplot(222); imgr.elem_data = sens(J2); show_fem(imgr,1); title('2D [log_{10}]');
        subplot(223); imgr.elem_data = sens(J2p5kq); show_fem(imgr,1); title('2.5D (k=0..3.0) [log_{10}]');
        subplot(224); imgr.elem_data = sens(J2p50); show_fem(imgr,1); title('2.5D (k=0) [log_{10}]');
   if 1
        subplot(222); imgr.elem_data = sens(J2p5kt); show_fem(imgr,1); title('2.5D (k=0..3.0) trapz [log_{10}]');
        subplot(223); imgr.elem_data = sens(J2p5kq); show_fem(imgr,1); title('2.5D (k=0..Inf) quadv [log_{10}]');
        subplot(224); imgr.elem_data = sens(J2p5ki); show_fem(imgr,1); title('2.5D (k=0..Inf) integral [log_{10}]');
   end

function S = sens(J)
   S = log10(sum(J.^2,1));
