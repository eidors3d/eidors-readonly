function data =fwd_solve_2p5d_1st_order(fwd_model, img)
% FWD_SOLVE_2P5D_1ST_ORDER: data= fwd_solve_2p5d_1st_order( img)
% Fwd solver for Andy Adler's EIT code
% Input:
%    img       = image struct
% Output:
%    data = measurements struct
% Options: (to return internal FEM information)
%    img.fwd_solve.get_all_meas = 1 (data.volt = all FEM nodes, but not CEM)
%    img.fwd_solve.get_all_nodes= 1 (data.volt = all nodes, including CEM)
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
end
fwd_model= img.fwd_model;

img = data_mapper(img);
if ~ismember(img.current_params, supported_params)
    error('EIDORS:PhysicsNotSupported', '%s does not support %s', ...
    'FWD_SOLVE_2P5D_1ST_ORDER',img.current_params);
end
% all calcs use conductivity
img = convert_img_units(img, 'conductivity');

pp= fwd_model_parameters( fwd_model, 'skip_VOLUME' );

k = [0 Inf]; try k = img.fwd_model.fwd_solve_2p5d_1st_order.k; end % default to 8 full waves: k pi h / 2 z = 8 * 2 pi
method = 'quadv'; try method = img.fwd_model.fwd_solve_2p5d_1st_order.method; end % default to 8 full waves: k pi h / 2 z = 8 * 2 pi

img.fwd_model.system_mat = @system_mat_2p5d_1st_order;
gnd = fwd_model.gnd_node;
img.fwd_model.system_mat_2p5d_1st_order.k = 0;
img.fwd_model.system_mat_2p5d_1st_order.factory = 1;
Sf = system_mat_2p5d_1st_order( img ); % returns a function Sf(k) that builds a system matrix for 'k'
if length(k) == 1 % singleton k
   v=2*potential_k(Sf(k), pp, gnd); % integral of a delta function
else
   switch method
      case 'trapz'
         % less accurate: trapz
         n = 0;
         tol = 1e-8;
         k(isinf(k)) = 1/sqrt(tol);
         vf = zeros(pp.n_elec*pp.n_stim,length(k)); % voltages under electrodes elec x stim, frequency domain
         for ki = k
            n = n + 1;
            vf(:,n) = potential_k(Sf(ki), pp, gnd);
         end
         v=2/pi*trapz(k,vf,2);
      case 'quadv'
         % more accurate: adaptive gaussian quadrature
         tol = 1e-8;
         kend = min(1/sqrt(tol), max(k)); % don't go too far... k=Inf is a singular matrix, stop adjacent to numeric singularity
         % quadv is scheduled to be removed from matlab eventually... but it is
         % WAY faster than integral with any tolerance configuration I could identify
         v=2/pi*quadv(@(kk) potential_k(Sf(kk), pp, gnd), k(1), kend);
      case 'integral'
         tol = 1e-8;
         kend = min(1/sqrt(tol), max(k)); % don't go too far... k=Inf is a singular matrix, stop adjacent to numeric singularity
         % the integral solution is about 10x slower (5.47 seconds vs. 0.60 seconds for UNIT_TEST)
         % ... I played with AbsTol and RelTol but wasn't able to affect the outcome
         v=2/pi*integral(@(kk) potential_k(Sf(kk), pp, gnd), k(1), kend, 'ArrayVAlued', true);
      otherwise
         error(['unrecognized method: ' method]);
   end
end
% up to now, we've been dealing with voltages at the electrodes,
% now we transform to dipole measurements v=(m+ - m-)
v2meas = get_v2meas(pp.n_elec, pp.n_stim, img.fwd_model.stimulation);
data.meas = v2meas' * v;
data.time= NaN; % unknown
data.name= 'solved by fwd_solve_2p5d_1st_order';
try; if img.fwd_solve.get_all_meas == 1
   data.volt = v(1:pp.n_node,:); % but not on CEM nodes
end; end
try; if img.fwd_solve.get_all_nodes== 1
   data.volt = v;                % all, including CEM nodes
end; end

% v0: a zero initialized array n x m for n nodes and m stimulus
function v = potential_k(S, pp, gnd)
   vf_nodes = zeros(pp.n_node, pp.n_stim); % voltages at all nodes x number of stim, frequency domain
   idx = 1:size(S,1); % pp.n_node + #CEM elec
   idx( gnd ) = [];
   idx2 = find(any(pp.N2E));
   vf_nodes(idx,:) = left_divide( S(idx,idx), 1/2*pp.QQ(idx,:));
   vf_elecs = pp.N2E(:,idx2) * vf_nodes(idx2,:); % electrode voltages, frequency domain
   v = vf_elecs(:); 

function v2meas = get_v2meas(n_elec,n_stim,stim)
    v2meas = sparse(n_elec*n_stim,0);
    for i=1:n_stim
        meas_pat= stim(i).meas_pattern;
        n_meas  = size(meas_pat,1);
        v2meas((i-1)*n_elec + 1: i*n_elec,end+(1:n_meas)) = meas_pat';
    end

function do_unit_test
   ne = 16;
   isOctave= exist('OCTAVE_VERSION');
   % 2D
   imdl2 = mk_geophysics_model('h2a', ne);
   imdl2.fwd_model.stimulation = stim_pattern_geophys(ne, 'Wenner');
   img2 = mk_image(imdl2.fwd_model, 1);
   % we get a 2D solution for klim=0; z = 1/2; h = 1;
   img2.fwd_model.fwd_solve_2p5d_1st_order.k = [0 logspace(-2,0,12) 3]; % default
   tic
   img2.fwd_model.nodes(1,:) = img2.fwd_model.nodes(1,:) + rand(1,2)*1e-8; % defeat cache
   v2 = fwd_solve(img2);
   fprintf('fwd_solve 2D t= %f s\n', toc);
   % 2.5D
   tic
   img2.fwd_model.fwd_solve_2p5d_1st_order.method = 'trapz';
   img2.fwd_model.nodes(1,:) = img2.fwd_model.nodes(1,:) + rand(1,2)*1e-8; % defeat cache
   v2p5t = fwd_solve_2p5d_1st_order(img2);
   fprintf('fwd_solve_2p5d_1st_order 2D t= %f s (%d k''s, trapz)\n', toc, length(img2.fwd_model.fwd_solve_2p5d_1st_order.k));
   tic
   img2.fwd_model.fwd_solve_2p5d_1st_order.method = 'quadv';
   img2.fwd_model.nodes(1,:) = img2.fwd_model.nodes(1,:) + rand(1,2)*1e-8; % defeat cache
   v2p5q = fwd_solve_2p5d_1st_order(img2);
   fprintf('fwd_solve_2p5d_1st_order 2D t= %f s (%d k''s, quadv)\n', toc, length(img2.fwd_model.fwd_solve_2p5d_1st_order.k));
   if isOctave
      fprintf('fwd_solve_2p5d_1st_order 2D t= --- s (-- k''s, integral) SKIPPED [not supported in Octave]\n');
   else
      tic
      img2.fwd_model.fwd_solve_2p5d_1st_order.method = 'integral';
      img2.fwd_model.nodes(1,:) = img2.fwd_model.nodes(1,:) + rand(1,2)*1e-8; % defeat cache
      v2p5i = fwd_solve_2p5d_1st_order(img2);
      fprintf('fwd_solve_2p5d_1st_order 2D t= %f s (%d k''s, integral)\n', toc, length(img2.fwd_model.fwd_solve_2p5d_1st_order.k));
   end
   tic
   img2.fwd_model.fwd_solve_2p5d_1st_order.k = 0;
   img2.fwd_model.nodes(1,:) = img2.fwd_model.nodes(1,:) + rand(1,2)*1e-8; % defeat cache
   v2p5k0 = fwd_solve_2p5d_1st_order(img2);
   fprintf('fwd_solve_2p5d_1st_order 2D t= %f s (%d k''s)\n', toc, length(img2.fwd_model.fwd_solve_2p5d_1st_order.k));
   % 3D
   imdl3 = mk_geophysics_model('h3a', ne);
   imdl3.fwd_model.stimulation = stim_pattern_geophys(ne, 'Wenner');
   img3 = mk_image(imdl3.fwd_model, 1);
   tic
   img3.fwd_model.nodes(1,:) = img3.fwd_model.nodes(1,:) + rand(1,3)*1e-8; % defeat cache
   v3 = fwd_solve(img3);
   fprintf('fwd_solve 3D t= %f s\n', toc);
   % analytic half-space (dipole)
   tic
   img2.fwd_model.nodes(1,:) = img2.fwd_model.nodes(1,:) + rand(1,2)*1e-8; % defeat cache
   vh = fwd_solve_halfspace(img2);
   fprintf('fwd_solve_halfspace 3D t= %f s\n', toc);
   scale=1e3; uvh=unique(round(vh.meas*scale)/scale);
   unit_test_cmp('h2a 2p5d values TEST', uvh, ...
   [ 0.0060; 0.0080; 0.0110; 0.0160; 0.0320 ], 1/scale);
   tol = norm(vh.meas)*0.025; % 2.5% error tolerance
   unit_test_cmp('2D   (h2a)                    vs analytic TEST', norm(vh.meas - v2.meas),     0, -tol);
   unit_test_cmp('2.5D (h2a + 2.5D, k=0)        vs 2D       TEST', norm(v2.meas - v2p5k0.meas), 0, tol);
   unit_test_cmp('2.5D (h2a + 2.5D trapz,2*tol) vs analytic TEST', norm(vh.meas - v2p5t.meas),   0, 2*tol);
   unit_test_cmp('2.5D (h2a + 2.5D quadv**)     vs analytic TEST', norm(vh.meas - v2p5q.meas),   0, tol);
   if ~isOctave
      unit_test_cmp('2.5D (h2a + 2.5D integral)    vs analytic TEST', norm(vh.meas - v2p5i.meas),   0, tol);
   end
   unit_test_cmp('3D   (h3b)                    vs analytic TEST', norm(vh.meas - v3.meas),     0, 2*tol);
pause(2);
clf; h=plot([vh.meas, v2p5q.meas, v3.meas, v2.meas],':');
     legend('analytic', 'FEM 2.5D', 'FEM 3D', 'FEM 2D', 'Location','Best');
     set(gca,'box','off');
     set(h,'LineWidth',2);
     set(h(1),'Marker','o','MarkerSize',7);
     set(h(2),'Marker','square','MarkerSize',7);
     set(h(3),'Marker','diamond','MarkerSize',3);
     set(h(4),'Marker','+');
pause(2);
clf; h=plot([vh.meas, v2p5q.meas, v3.meas],':');
     legend('analytic', 'FEM 2.5D', 'FEM 3D','Location','Best');
     set(gca,'box','off');
     set(h,'LineWidth',2);
     set(h(1),'Marker','o','MarkerSize',7);
     set(h(2),'Marker','square','MarkerSize',7);
     set(h(3),'Marker','diamond','MarkerSize',3);
