function J=jacobian_tall(fwd_model, img)
% JACOBIAN_TALL(img)
% 
% A wrapper function around the selected Jacobian calculator
% (calc_jacobian_XXX), which helps when there are many
% measurements (typically > 1000) where memory consumption and
% calculation time become excessive.
%
%  fwd_model.jacobian_tall.jacobian    [default:'eidors_default']
%     Jacobian calculator to use in inner loops
%
%  fwd_model.jacobian_tall.threshold   [default: 250]
%     Number of measurements per block to throw at the inner
%     Jacobian calculator. Results are concatentated, as if the
%     inner Jacobian calculator were given the entire stimulus in
%     one shot.
%
% (C) 2016, A. Boyle
% License: GPL version 2 or version 3

% UNIT_TEST?
J = 0;
if nargin == 1
   img = fwd_model;
else
   warning('EIDORS:DeprecatedInterface', ...
      ['Calling CALC_JACOBIAN with two arguments is deprecated and will cause' ...
       ' an error in a future version. First argument ignored.']);
end
if isstr(img) && strcmp(img,'UNIT_TEST'); do_unit_test; return; end

opt.jacobian = 'eidors_default'; try; opt.jacobian = img.fwd_model.jacobian_tall.jacobian; end
opt.thres = 250; try; opt.thres = img.fwd_model.jacobian_tall.threshold; end
assert(opt.thres > 0, 'threshold must be > 0');

img.fwd_model.jacobian = opt.jacobian;

% note that stim_meas_list conversions are too slow to be efficient,
% here, we split approximately based on the stimulation(#) and we get some arbitrary number of measurements
stim = img.fwd_model.stimulation;
n = count_meas(stim);
m=size(img.fwd_model.elems,1); try; m=size(img.fwd_model.coarse2fine,2); end
J=zeros(n,m);
if 0 % bypass
   J = calc_jacobian(img);
elseif 0 % slow
   % we assume the same stimulus and measurement gains
   stim_amp = max(max(stim(1).stim_pattern));
   meas_amp = max(max(stim(1).meas_pattern));
   ne = size(stim(1).stim_pattern,1);
   stim_list = stim_meas_list(stim);
   for i=0:ceil(n/opt.thres)-1
      m0 = i*opt.thres +1; m1 = min(n, (i+1)*opt.thres);
      img.fwd_model.stimulation = stim_meas_list(stim_list(m0:m1,:), ne, stim_amp, meas_amp);
      J(m0:m1,:) = calc_jacobian(img);
   end
else % fast
   n_stim = length(img.fwd_model.stimulation);
   n1 = 0; j = 0; % init
   while n1 < n
      n0 = n1 + 1;
      i = j+1;
      while (j < n_stim) && (n1-n0+1 < opt.thres)
         j = j + 1;
         n1 = n1 + count_meas(stim(j));
      end
      img.fwd_model.stimulation = stim(i:j);
      J(n0:n1,:) = calc_jacobian(img);
   end
end

function n = count_meas(stim)
   n = 0;
   for i=1:length(stim)
      n = n + size(stim(i).meas_pattern,1);
   end

function do_unit_test
ne=64;
for mdl={'h2a','h2p5a'}
   stim=mk_stim_patterns(ne,1,[0 2],[0 1]);
   imdl=mk_geophysics_model(mdl{:}, ne);
   fprintf('model = %s, # stim = %d, fwd elems = %d, inv elems = %d\n', ...
      mdl{:}, count_meas(stim), size(imdl.fwd_model.elems,1), size(imdl.rec_model.elems,1));
   nd = size(imdl.fwd_model.nodes,2);
   % we twiddle a model node to defeat the top level of eidors_cache-ing
   % while still allowing caching at the lower levels
   imdl.fwd_model.nodes(1,:) = imdl.fwd_model.nodes(1,:) + rand(1,nd)*eps*1e3;
   imdl.fwd_model.stimulation = stim; 
   img=mk_image(imdl.fwd_model,1);    
   tic; Jr=calc_jacobian(img); tr=toc; fprintf('run time calc_jacobian:                    %f s\n',tr);
   imdl.fwd_model.jacobian = @jacobian_tall;
   imdl.fwd_model.jacobian_tall.jacobian = 'eidors_default';
   tol = eps;
   for thres=[250 500 1000 5000];
      imdl.fwd_model.nodes(1,:) = imdl.fwd_model.nodes(1,:) + rand(1,nd)*eps*1e3;
      imdl.fwd_model.jacobian_tall.threshold = thres;
      img=mk_image(imdl.fwd_model,1);    
      tic; Jt=calc_jacobian(img); tt=toc; fprintf('run time jacobian_tall (thres=%6.0d) %f s (%0.3fx)\n',thres,tt,tt/tr);
      unit_test_cmp(sprintf('tall (thres=%d)',thres),Jr,Jt,tol);
   %   disp([norm(Jr) norm(Jt) norm(Jr - Jt)/norm(Jr)])
   end
end
