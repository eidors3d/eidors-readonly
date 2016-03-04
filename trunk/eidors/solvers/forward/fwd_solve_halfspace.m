function data = fwd_solve_halfspace(fwd_model, img)
% FWD_SOLVE_HALFSPACE: data = fwd_solve_halfspace(img)
% Fwd solver for an analytic halfspace
% Input:
%    img       = image struct
% Output:
%    data = measurements struct
% Options: (to return internal node information)
%    img.fwd_solve.get_all_meas = 1 (data.volt = all FEM nodes, but not CEM)
%    img.fwd_solve.get_all_nodes= 1 (data.volt = all nodes, including CEM)
%
% *** Only supports point electrode models (PEM), a CEM will be approximated as
%     a PEM at the center of the electrode geometry.
%
% *** Surface topology is NOT modelled. It is assumed that the electrodes are
%     approximately coplanar and on the surface.
%
% *** The domain is assumed to be homogeneous. The mean conductivity will be used.
%
% *** Contact impedances z_c are ignored. (TODO)
%
% (C) 2016 A. Boyle. License: GPL version 2 or version 3
% dipole model: [Loke , "Tutorial: 2-D and 3-D electrical imaging surveys", 2004]
% Ref: TODO

% correct input paralemeters if function was called with only img
if ischar(fwd_model) && strcmp(fwd_model,'UNIT_TEST'); do_unit_test; return; end

if nargin == 1
   img= fwd_model;
elseif  strcmp(getfield(warning('query','EIDORS:DeprecatedInterface'),'state'),'on')
   error('EIDORS:DeprecatedInterface', ...
      'Calling FWD_SOLVE_HALFSPACE with two arguments is deprecated');
end
fwd_model= img.fwd_model;

img = data_mapper(img);
if ~ismember(img.current_params, supported_params)
    error('EIDORS:PhysicsNotSupported', '%s does not support %s', ...
    'FWD_SOLVE_HALFSPACE',img.current_params);
end
pp= fwd_model_parameters( fwd_model, 'skip_VOLUME' );

% all calcs use resistivity
img = convert_img_units(img, 'conductivity');
rho = 1/mean(img.elem_data);

% create a data structure to return
data.meas= analytic_meas(fwd_model, rho);
data.time= NaN; % unknown
data.name= 'solved by fwd_solve_halfspace';
if isfield(img, 'fwd_solve') && isfield(img.fwd_solve, 'get_all_meas') && img.fwd_solve.get_all_meas
   v = analytic_nodal(fwd_model, rho);
   data.volt = v(1:pp.n_node,:); % but not on CEM nodes
end
try; if img.fwd_solve.get_all_nodes
   error('does not support CEM: can''t get_all_nodes; try get_all_meas');
end; end


function do_unit_test
   ne = 16;
   imdl = mk_geophysics_model('h2p5a', ne);
   imdl.fwd_model.stimulation = stim_pattern_geophys(ne, 'Wenner');
   img = mk_image(imdl.fwd_model, 1);
   img.fwd_solve.get_all_meas = 1;
   vh = fwd_solve_halfspace(img);
   vd = fwd_solve(img);
   scale=1e3; uvh=unique(round(vh.meas*scale)/scale);
   unit_test_cmp('h2a halfspace values TEST', uvh, ...
  [ 0.0060; 0.0080; 0.0110; 0.0160; 0.0320 ], 1/scale);
   unit_test_cmp('h2a halfspace vs default TEST', norm(vh.meas - vd.meas), 0, 4e-3);
clf; h=plot([vh.meas vd.meas],'o--'); legend('analytic','FEM'); set(gca,'box','off'); set(h,'LineWidth',2);
%   img.fwd_solve.get_all_meas = 1;
%   vh = fwd_solve_halfspace(img);
%   gnd_node = img.fwd_model.nodes(img.fwd_model.gnd_node,:);
%clf; show_fem(img); hold on; plot3(gnd_node(1),gnd_node(2),gnd_node(3),'o');
%clf; plot(vh.volt);

% returns the potential difference
% half-space with PEM electrodes
% [Loke , "Tutorial: 2-D and 3-D electrical imaging surveys", 2004]
function dv = analytic_meas(mdl, rho)
   ne = length(mdl.electrode); % number of electrodes
   nd = size(mdl.nodes, 2); % dimensions in model
   ep = zeros(ne,nd); % init electrode position
   for i = 1:ne
      en = mdl.electrode(i).nodes;
      ep(i,:) = mean(mdl.nodes(en), 2);
   end
   list = stim_meas_list(mdl.stimulation);
   a = list(:,1); % s+
   b = list(:,2); % s-
   m = list(:,3); % m+
   n = list(:,4); % m-
   AM = dist(ep(a), ep(m));
   BM = dist(ep(b), ep(m));
   AN = dist(ep(a), ep(n));
   BN = dist(ep(b), ep(n));
% TODO  zc = [mdl.electrode(:).z_contact]';
% TODO  MNzc = zc(m) + zc(n);

   k = geometric_factor(AM, BM, AN, BN);
   I = stim_current(mdl);
   dv = (rho*I)./k; % TODO AB missing zc

% returns the potential everywhere relative to the ground node
% half-space with PEM electrodes
% [Loke , "Tutorial: 2-D and 3-D electrical imaging surveys", 2004]
function v = analytic_nodal(mdl, rho)
   ne = length(mdl.electrode); % number of electrodes
   ns = length(mdl.stimulation); % number of stim patterns
   nn = size(mdl.nodes, 1); % nodes in the model
   nd = size(mdl.nodes, 2); % dimensions in model
   ep = zeros(ne,nd); % init electrode position
   for i = 1:ne
      en = mdl.electrode(i).nodes;
      ep(i,:) = mean(mdl.nodes(en), 2);
   end
   a = ones(ns,1)*nan; b = a; I = a;
   for i=1:ns % same order as pp.QQ
      s = mdl.stimulation(i).stim_pattern;
      assert_sane_stim(mdl.stimulation(i));
      a(i) = find(s > 0);
      b(i) = find(s < 0);
      I(i) = sum(s(s > 0));
   end
   epm = mdl.nodes(mdl.gnd_node); % ground node
   epn = mdl.nodes; % all nodes
   AM = dist(ep(a), epm);
   BM = dist(ep(b), epm);
   abN = dist(kron(ep,ones(nn,1)), repmat(epn,ne,1)); % every electrode to every node
   abN = reshape(abN, nn, ne); % cols = electrode#, rows = node#
% TODO  zc = [mdl.electrode(:).z_contact]';
% TODO  MNzc = zc(m) + zc(n);

   v = ones(nn,length(I));
   for i=1:ne
      ai = a(i); bi = b(i);
      k = geometric_factor(AM(ai), BM(bi), abN(:,ai), abN(:,bi));
      v(:,i) = rho*I(i)./k; % TODO AB missing zc
   end

% transfer resistance:" stimulus current over measured voltage
function [I, a, b] = stim_current(mdl)
   stim = mdl.stimulation;
   I = [];
   N = 0;
   for i = 1:length(stim)
      s = stim(i).stim_pattern;
      assert_sane_stim(stim(i));
      n = size(stim(i).meas_pattern, 1);
      I(N+[1:n]) = sum(s(s > 0));
      N = N + n;
   end
   I=I';

function assert_sane_stim(stim)
   assert(length(stim) == 1); % only one stim pattern per "assert_sane" check
   s = stim.stim_pattern;
   assert(strcmp(stim.stimulation, 'Amp')); % scaling?
   assert(length(find(abs(s) > 0)) == 2); % pairwise stimulus only
   assert(sum(s) == 0); % dipole
   assert(size(s,2) == 1); % single dipole per stimulation(i)

% distance from A to B
function AB = dist(A, B)
   AB = sqrt(sum((B - A).^2,2));

% analytic model, for half-space with point electrodes (PEM)
function k = geometric_factor(AM, BM, AN, BN)
   k = 2*pi./(1./AM-1./BM-1./AN+1./BN);
