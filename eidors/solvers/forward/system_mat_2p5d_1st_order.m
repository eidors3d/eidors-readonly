function s_mat= system_mat_2p5d_1st_order( fwd_model, img)
% SYSTEM_MAT_2P5D_1ST_ORDER: 2.5D system matrix
% SS= SYSTEM_MAT_2P5D_1ST_ORDER( FWD_MODEL, IMG)
% fwd_model = forward model
% img       = image background for system matrix calc
% s_mat = CC' * (SS + k^2 * TT) * conductivites * CC;
% where:
%   SS  = Unconnected system Matrix
%   TT  = Unconnected Fourier correction
%   CC  = Connectivity Matrix
%
% Parameters
%   fwd_model.system_mat_2_5D_1st_order.k
%        the k-th harmonic 0..inf [default 0]
%   fwd_model.system_mat_2_5D_1st_order.factory
%        return a factory in s_mat.assemble(s_mat, k) to build
%        the system matrix based on argument k [default 0]
%
% NOTE: This code is adapted from dev/a_adler/system_mat_2_5D_1st_order.m

% (C) 2005, 2013, 2016 A. Adler, A. Boyle. Licensed under the GPL Version 2

if ischar(fwd_model) && strcmp(fwd_model,'UNIT_TEST'); do_unit_test; return; end

if nargin == 1
   img= fwd_model;
elseif  strcmp(getfield(warning('query','EIDORS:DeprecatedInterface'),'state'),'on')
   warning('EIDORS:DeprecatedInterface', ...
      ['Calling SYSTEM_MAT_1ST_ORDER with two arguments is deprecated and will cause' ...
       ' an error in a future version. First argument ignored.']);
end
fwd_model= img.fwd_model;

% check parametrizations
if isfield(img,'current_params') && ...
     (~ischar(img.current_params) || ~strcmp(img.current_params,'conductivity'))
    if ischar(img.current_params)
       error('system_mat_2p5d_1st_order does not work for %s',img.current_params);
    else
       img.current_params
       error('system_mat_2p5d_1st_order does not work for the above');
    end
end

factory = 0; try factory = fwd_model.system_mat_2p5d_1st_order.factory; end
k = 0; try k = fwd_model.system_mat_2p5d_1st_order.k; end

FS = system_mat_fields(img.fwd_model);      % 2D solution
FT = system_mat_2p5d_fields(img.fwd_model); % Fourier domain correction
[ES, ET] = calc_sigma(img, FS, FT); % conductivity

if factory
   s_mat = @(k) assemble(k,FS,FT,ES,ET);
else
   s_mat.E = assemble(k,FS,FT,ES,ET);
end


function SS = assemble(k,FS,FT,ES,ET)
   assert(all(size(k) == [1 1]), 'expected scalar k');
   SS = FS'*ES*FS + k^2*FT'*ET*FT;

function [ES,ET] = calc_sigma(img, FS, FT)
   lFS= size(FS,1);
   lFT= size(FT,1);
   ED = elem_dim(img.fwd_model);
   lNE= ED*num_elems(img.fwd_model);
   elem_data = check_elem_data(img.fwd_model, img);
   if size(elem_data,3) == 1
   % Scalar conductivity == isotropic
      elem_sigma = kron( elem_data, ones(ED,1) );
      elem_sigma(end+1:lFS) = 1; % add ones for CEM
      ES= spdiags(elem_sigma,0,lFS,lFS);

      elem_sigma = kron( elem_data, ones(ED+1,1) );
      elem_sigma(end+1:lFT) = 1; % add ones for CEM
      ET= spdiags(elem_sigma,0,lFT,lFT);
   else
%      switch elem_dim(img.fwd_model)
%        case 2;
%          idx = 1:2:lNE;
%          ES= sparse([idx,idx+1,idx,idx+1]', ...
%                     [idx,idx,idx+1,idx+1]', elem_data(:), lFS,lFS);
%
%          ES(lNE+1:lFS,lNE+1:lFS) = speye(lFS-lNE);
%        case 3;
%          idx = 1:3:lNE;
%          ES= sparse([idx,idx+1,idx+2,idx,idx+1,idx+2,idx,idx+1,idx+2]', ...
%                     [idx,idx,idx,idx+1,idx+1,idx+1,idx+2,idx+2,idx+2]', ...
%                      elem_data(:), lFS,lFS);
%
%          ES(lNE+1:lFS,lNE+1:lFS) = speye(lFS-lNE);
%        otherwise;
          error('%d D anisotropic elements not implemented', elem_dim(fwd_model));
%      end
   end

function check_2d(img)
   assert(size(img.fwd_model.nodes,2) == 2, 'expecting 2D model for 2.5D reconstruction');
   % check that all electrodes are actually part of the mesh
   nn = unique(img.fwd_model.elems);
   for i=1:length(img.fwd_model.electrode)
      ne = img.fwd_model.electrode(i).nodes;
      for j=1:length(ne)
         assert(length(find(ne(j) == nn))>0, sprintf('elec#%d, node#%d: not connected to mesh',i,j));
      end
   end

function elem_data = check_elem_data(fwd_model, img);
   elem_data = img.elem_data;
   if any(size(elem_data) == [1 1])
      elem_data = elem_data(:);
   end
   sz_elem_data = size(elem_data);
   if sz_elem_data(2) ~= 1;
      error('system_mat_2p5d_1st_order: can only solve one image (sz_elem_data=%)', ...
            sz_elem_data);
   end

   if isfield(fwd_model, 'coarse2fine');
     c2f = fwd_model.coarse2fine;
     sz_c2f = size(c2f);
     switch sz_elem_data(1)
       case sz_c2f(1); % Ok
       case sz_c2f(2); elem_data = c2f * elem_data;
          if isfield(fwd_model, 'background')
              elem_data = elem_data + fwd_model.background;
          end

       otherwise; error(['system_mat_2p5d_1st_order: provided elem_data ' ...
            ' (sz=%d) does not match c2f (sz=%d %d)'], sz_elem_data(1), sz_c2f);
     end
   else
     if sz_elem_data(1) ~= num_elems(fwd_model)
       error(['system_mat_2p5d_1st_order: provided elem_data (sz=%d) does ' ...
          ' not match fwd_model (sz=%d)'], sz_elem_data(1), num_elems(fwd_model));
     end
   end

function do_unit_test()
   ne = 16;
   k=0; % for k=0, the result is a 2D reconstruction
   % 2D
   imdl2 = mk_geophysics_model('h2a', ne);
   imdl2.fwd_model.stimulation = stim_pattern_geophys(ne, 'Wenner');
   img2 = mk_image(imdl2.fwd_model, 1);
   check_2d(img2);
   v2 = fwd_solve(img2);
   % 2.5D
   img2.fwd_model.system_mat = @system_mat_2p5d_1st_order;
   img2.fwd_model.system_mat_2p5d_1st_order.k = k;
   v2p5 = fwd_solve(img2);
   % 3D
   imdl3 = mk_geophysics_model('h3a', ne);
   imdl3.fwd_model.stimulation = stim_pattern_geophys(ne, 'Wenner');
   img3 = mk_image(imdl3.fwd_model, 1);
   v3 = fwd_solve(img3);
   % analytic half-space (dipole)
   vh = fwd_solve_halfspace(img2);
   scale=1e3; uvh=unique(round(vh.meas*scale)/scale);
   unit_test_cmp('h2a 2p5d values TEST', uvh, ...
   [ 0.0060; 0.0080; 0.0110; 0.0160; 0.0320 ], 1/scale);
   tol = norm(vh.meas)*0.050; % 5% error tolerance
   unit_test_cmp('2D   (h2a)                vs analytic TEST', norm(vh.meas - v2.meas),   0, -tol);
   unit_test_cmp('2.5D (h2a + 2.5D sys_mat) vs 2D       TEST', norm(v2.meas - v2p5.meas), 0, tol);
   unit_test_cmp('3D   (h3a)                vs analytic TEST', norm(vh.meas - v3.meas),   0, tol);
clf; h=plot([vh.meas, v2.meas, v2p5.meas, v3.meas],':');
     legend('analytic','FEM 2D', sprintf('FEM 2.5D k=%0.1g',k), 'FEM 3D','Location','Best');
     set(gca,'box','off');
     set(h,'LineWidth',2);
     set(h(1),'Marker','o','MarkerSize',7);
     set(h(2),'Marker','+');
     set(h(3),'Marker','square','MarkerSize',7);
     set(h(4),'Marker','diamond','MarkerSize',3);
