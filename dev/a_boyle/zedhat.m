function out = zedhat(in, filename, fopt)
% function img = zedhat('img.zh')          % load
% function zedhat(img, 'mdl.zhm', [fopt])  % store
% function zedhat(imdl, 'mdl.zhm', [fopt]) % store
% function zedhat(fmdl, 'mdl.zhm', [fopt]) % store
% function zedhat({ data }, 'mdl.zhm', [fopt]) % store
%
% img/inv_model/fwd_model: the mesh, stim/meas, conductivities
%
% fopt = 'a' for "append to file", 'w' for create or truncate to zero [default]
%
% NOTE: load only returns the *first* data in the file
%
% ex:
%   % create model in EIDORS
%   imdl = mk_common_model('d2d4c',16);
%   zedhat(imdl, 'imdl.zh');
%   % generate difference data
%   fmdl = mk_common_model('f2d4c',16);
%   imgh = mk_image(fmdl,1);
%   vh = fwd_solve(imgh);
%   imgi = create_inclusion(imgh, [+0.3 +0.3], 0.2, 0.9);
%   imgi = create_inclusion(imgi, [-0.4 -0.1], 0.1, 1.1);
%   vi = fwd_solve(imgi);
%   zedhat({vh, vi}, 'imdl.zh', 'a');
%   % generate expected solution for comparison
%   imge = inv_solve(imdl, vh, vi);
%   % solve in zedhat ...
%   % load solution
%   imgz = zedhat('img.zh');
%   clf; subplot(221); show_fem(imgi,1); xlabel('fwd'); subplot(223); show_fem(imge,1); xlabel('eidors'); subplot(224); show_fem(imgz,1); xlabel('zedhat');
%   assert(norm(imgz.elem_data(:) - imge.elem_data(:)) < eps*size(imge.elem_data,1), 'ERROR: zedhat-eidors img comparison failed');
%
% (C) 2019, A. Boyle
if isstr(in) && strcmp(in,'UNIT_TEST'); out=do_unit_test(); return; end

if nargin == 1
   filename = in;
   fopt = 'r';
elseif nargin < 3
   fopt = 'w';
end

FILE = fopen(filename, fopt);

if nargin == 1
   out = deserialize_img(FILE);
else
   out = in;
   if iscell(in)
      in = cell2mat(in);
      in = [ in(:).meas ];
   end
   if isnumeric(in)
      serialize_data(FILE, in);
   else
      assert(isstruct(in), 'ERROR: zedhat export: expected struct or raw data');
      assert(isfield(in, 'type'), 'ERROR: zedhat export: expected ''type'' field');
      switch in.type
      case 'data'
         serialize_data(FILE, in.meas);
      case 'image'
         serialize_img(FILE, in);
      case 'inv_model'
         serialize_imdl(FILE, in);
      case 'fwd_model'
         serialize_fmdl(FILE, in);
      otherwise
         error('unrecognized in.type');
      end
   end
end

fclose(FILE);
end

function serialize_data(FILE, data)
   serialize_array(FILE, 'data', '%24.16e', data, '# measurements [V]');
end

function serialize_imdl(FILE, imdl)
   serialize_header(FILE);
   serialize_int(FILE, 'format', 1);
   serialize_mdl_type(FILE, 'forward');
   serialize_mdl(FILE, imdl.fwd_model);
   if isfield(imdl, 'rec_model')
      serialize_mdl_type(FILE, 'reconstruction');
      serialize_mdl(FILE, imdl.rec_model);
   end
   if isfield(imdl, 'hyperparameter') & isfield(imdl.hyperparameter, 'value')
      serialize_block(FILE, 'hyperparameter', '%24.16e', [imdl.hyperparameter.value], '# λ_conductivity');
   end
   if isfield(imdl, 'elem_data')
      parameters_check(imdl);
      serialize_array(FILE, 'parameters', '%24.16e', [imdl.elem_data], '# conductivity [S/m])');
   end
   if isfield(imdl, 'meas')
      serialize_array(FILE, 'data', '%24.16e', imdl.meas, '# measurements [V]');
   end
end

function parameters_check(imdl)
   fwd_elems = size(imdl.fwd_model.elems,1);
   if isfield(imdl.fwd_model, 'coarse2fine')
      fwd_elems = size(imdl.fwd_model.coarse2fine, 2);
   end
   assert(size(imdl.elem_data,1) == fwd_elems);
   if isfield(imdl, 'rec_model')
      rec_elems = size(imdl.rec_model.elems,1);
      if isfield(imdl.rec_model, 'coarse2fine')
         rec_elems = size(imdl.rec_model.coarse2fine, 2);
      end
      assert(size(imdl.elem_data,1) == rec_elems);
   end
end

function serialize_img(FILE, img)
   serialize_imdl(FILE, img);
end

function serialize_fmdl(FILE, fmdl)
   serialize_header(FILE);
   serialize_int(FILE, 'format', 1, '# file format version');
   serialize_mdl_type(FILE, 'forward');
   serialize_mdl(FILE, fmdl);
end

function img = deserialize_img(FILE)
   if deserialize_test(FILE, 'hyperparameter')
      img.hyperparameter.value = deserialize_block(FILE, 'hyperparameter', '%e', 1);
   end
   if deserialize_test(FILE, 'parameters')
      img.elem_data = deserialize_array(FILE, 'parameters', '%e');
   end
   if deserialize_test(FILE, 'data')
      img.meas = deserialize_array(FILE, 'data', '%e');
   end
   img.fwd_model = deserialize_mdl(FILE);
   % TODO c2f in fwd_model and rec_model

   % make the imported data happy in eidors
   img.name = 'zedhat';
   if isfield(img, 'elem_data')
      img.type = 'image';
      img.name = 'zedhat';
   elseif isfield(img, 'hyperparameter')
      img.type = 'inv_model';
      img.solve = 'eidors_default';
      img.RtR_prior = 'eidors_default';
      img.jacobian_bkgnd.value = 1;
      img.reconst_type = 'difference';
   else
      img = img.fwd_model;
   end
end

function serialize_mdl_type(FILE, type)
   block = 'modeltype';
   fprintf(FILE,[ block '\n']);
   fprintf(FILE, '%s\n\n', type);
end

function serialize_mdl(FILE, mdl)
   dim = size(mdl.nodes,2);
   mdl.boundary = find_boundary(mdl.elems);

   surfnr = 1;
   bcnr = 0;
   domin = 1;
   domout = 0;
   xx = ones(size(mdl.boundary,1),1); % expand scalar to number of rows in boundary list
   TODO = xx;
   % find electrodes, using surface elements
   bcnr = bcnr*xx;
   thres = size(mdl.boundary,2);
   pem = [];
   for e = 1:length(mdl.electrode)
      nn = mdl.electrode(e).nodes;
      if length(nn) == 1 % PEM
         pem = [ pem; surfnr e domin domout dim nn zeros(1,dim-1) ]; % TODO or better as [ ... 1 nn ];
         continue;
      end
      srfcnt = 0*xx;
      for n = nn(:)'
         srfcnt = srfcnt + sum(n == mdl.boundary,2);
      end
      bcnr(srfcnt >= thres) = e;
   end
   surfaceelements = [surfnr*TODO bcnr domin*TODO domout*TODO dim*xx mdl.boundary; pem];

   matnr = 1; % TODO
   xx = ones(size(mdl.elems,1),1); % expand scalar to number of rows in boundary list
   TODO = xx;
   volumeelements = [matnr*TODO (dim+1)*xx mdl.elems];

   stimmeas = stim_meas_list(mdl.stimulation);
   assert(size(stimmeas,2) == 4);
   amp = [];
   for i = 1:length(mdl.stimulation)
       stim = mdl.stimulation(i);
       stim_amp = full(max(stim.stim_pattern(:)));
       meas_gain =  full(max(stim.meas_pattern,[],2));
       amp = [ amp; stim_amp * meas_gain ];
   end
   stimmeas = [stimmeas amp];

   % netgen
   points_desc = '#          X             Y';
   volumeelements_desc = '#  matnr      np      p1      p2      p3';
   surfaceelems_desc = '# surfnr    bcnr   domin  domout      np      p1      p2';
   parametermap_desc = '# elementnr parameternr overlapfraction [row sum normally = 1.0]';
   zc_desc = '# bcnr zc [Ω/m]';
   if(dim == 3)
      points_desc = [points_desc '             Z'];
      volumeelements_desc = [volumeelements_desc '      p4'];
      surfaceelems_desc = [surfaceelems_desc '      p3'];
      zc_desc = '# bcnr zc [Ω/m²]';
   end
   serialize_int(FILE, 'dimension', dim);
   serialize_int(FILE, 'geomtype', 0);
   serialize_block(FILE, 'surfaceelements', '%d', surfaceelements, surfaceelems_desc);

   serialize_block(FILE, 'points', '%24.16e', mdl.nodes, points_desc);
   serialize_block(FILE, 'volumeelements', '%d', volumeelements, volumeelements_desc);
   if isfield(mdl, 'electrode')
      ne = length(mdl.electrode);
      serialize_block(FILE, 'contactimpedances', {'%d','%24.16e'}, [ [1:ne]' [mdl.electrode(:).z_contact]'], zc_desc);
   end
   % zedhat
   stimmeas_desc = '#  A  B  M  N  amp*gain';
   serialize_block(FILE, 'stimmeas', {'%d','%d','%d','%d','%24.16e'}, stimmeas, stimmeas_desc);
   if isfield(mdl,'coarse2fine')
      serialize_sparse(FILE, 'parametermap', mdl.coarse2fine, parametermap_desc);
   end
end

function mdl = deserialize_mdl(FILE)
   % netgen
   dim = deserialize_int(FILE, 'dimension');
   assert(dim > 1 & dim < 4, 'ERROR: zedhat import: bad dimension in model dim=2..3');
   geomtype = deserialize_int(FILE, 'geomtype');
   assert(geomtype == 0, 'ERROR: zedhat import: bad geomtype in model != 0');
   surfaceelements = deserialize_block(FILE, 'surfaceelements', '%d', 5+dim);
   points = deserialize_block(FILE, 'points', '%e', dim);
   volumeelements = deserialize_block(FILE, 'volumeelements', '%d', 2+(dim+1));
   % zedhat
   stimmeas = deserialize_block(FILE, 'stimmeas', {'%d', '%d', '%d', '%d', '%e'}, 5);
   contactimpedances = deserialize_block(FILE, 'contactimpedances', {'%d', '%e'}, 2);
   parametermap = [];
   if deserialize_test(FILE, 'parametermap')
      parametermap = deserialize_block(FILE, 'parametermap',{'%d','%d','%e'}, 3);
   end

   % EIDORS fwd_model
   mdl.nodes = points;
   mdl.elems = uint32(volumeelements(:, (end-dim):end));
   mdl.type = 'fwd_model';
   mdl.name = 'zedhat';
   n = 0;
   bc_list = sort(unique(surfaceelements(:,2)))';
   for bc = bc_list(2:end)
      n = n + 1;
      idx = find(surfaceelements(:,2) == bc);
      nodes = surfaceelements(idx,(end-(dim-1)):end);
      nodes = nodes(nodes > 0); % remove any PEM dummy nodes
      mdl.electrode(n).nodes = unique(sort(nodes(:)));
      mdl.electrode(n).z_contact = 0.01;
   end
   mdl.gnd_node = 1;
   mdl.normalize_measurements = 0;
   %remove any PEM dummy nodes from surfaceelements
   [i,~] = find(surfaceelements(:,6:end) == 0);
   surfaceelements(i,:) = [];
   mdl.boundary = uint32(surfaceelements(:, (end-(dim-1)):end));
   if size(parametermap,2) > 0
      ii = parametermap(:,1);
      jj = parametermap(:,2);
      vv = parametermap(:,3);
      mdl.coarse2fine = sparse(ii, jj, vv, size(mdl.elems,1), max(jj));
   end

   num_elec_mesh = length(bc_list) - 1;
   stimmeas_abmn = stimmeas(:,1:4);
   num_elec_stimmeas = length(unique(sort(stimmeas_abmn(:))));
   num_elec = max(num_elec_mesh, num_elec_stimmeas);
   mdl.stimulation = stim_meas_list(stimmeas_abmn, num_elec, 0.010, 1);
   n_meas = 0;
   idx = 1;
   while n_meas < size(stimmeas,1)
      stim = mdl.stimulation(idx);
      incr = size(stim.meas_pattern,1);
      gain = stimmeas(n_meas+1:n_meas+incr,5);
      mdl.stimulation(idx).meas_pattern = diag(gain) * stim.meas_pattern;
      n_meas = n_meas + incr;
      idx = idx + 1;
   end
   for i = 1:size(contactimpedances,1)
      e = contactimpedances(i,1);
      zc = contactimpedances(i,2);
      mdl.electrode(e).z_contact = zc;
   end

   mdl.solve = 'eidors_default';
   mdl.system_mat = 'eidors_default';
   mdl.jacobian = 'eidors_default';
end


function ret = deserialize_test(FILE, block)
   fseek(FILE, 0, 'bof');
   line = '';
   while ~strcmp(sprintf('%s\n', block), line) && ~feof(FILE)
      line = fgets(FILE);
   end
   ret = ~feof(FILE);
end

function out = deserialize_int(FILE, block)
   out = deserialize_type(FILE, block, '%d');
end

function blk = deserialize_type(FILE, block, type)
   ret = deserialize_test(FILE, block);
   if ~iscell(type)
      type = {type};
   end
   assert(ret, ['ERROR: zedhat import: missing ' block]);
   n = length(type);
   if n > 1
      [blk, cnt] = fscanf(FILE, [type{1} ' ' type{2} '\n'], n);
      assert(cnt == n, ['ERROR: zedhat import: missing ' block ' count']);
   else
      [blk, cnt] = fscanf(FILE, [type{1} '\n'], n);
      assert(cnt == n, ['ERROR: zedhat import: missing ' block ' count']);
   end
end

function blk = deserialize_block(FILE, block, type, cols)
   blocktype = '%d';
   if cols < 0
      blocktype = {'%d','%d'};
   end
   rows = deserialize_type(FILE, block, blocktype);
   if length(rows) > 1
      cols = rows(2);
   end
   if ~iscell(type)
      type = {type};
   end
   for i=length(type):cols
      type{end+1} = type{end};
   end
   blk = zeros(rows(1),cols);
   for i = 1:rows
      for j = 1:cols-1
         [blk(i, j), cnt] = fscanf(FILE, [type{j} ' '], 1);
         assert(cnt == 1);
      end
      [blk(i, cols), cnt] = fscanf(FILE, [type{end} '\n'], 1);
      assert(cnt == 1);
   end
end

function blk = deserialize_array(FILE, block, type)
   blk = deserialize_block(FILE, block, type, -1);
end

function blk = serialize_header(FILE)
   fprintf(FILE, 'zedhat\n\n');
end

function blk = serialize_int(FILE, block, in, comment)
   if nargin > 3
      fprintf(FILE, [comment '\n']);
   end
   fprintf(FILE,[ block '\n']);
   fprintf(FILE, '%d\n\n', in);
end

function serialize_block_or_array(FILE, block, rows, type, in, comment)
   if nargin > 5
      fprintf(FILE, [comment '\n']);
   end
   if ~iscell(type)
      type = {type};
   end
   fprintf(FILE,[ block '\n']);

   if length(rows) == 1
      fprintf(FILE, '%d\n', rows);
      cols = size(in,2);
   else
      fprintf(FILE, '%d %d\n', rows(1), rows(2));
      cols = rows(2);
      rows = rows(1);
   end
   for i=length(type):cols
      type{end+1} = type{end};
   end
   if ~iscell(in)
      in = num2cell(in);
   end
   for i = 1:rows
      for j = 1:cols-1
          fprintf(FILE, [type{j} ' '], in{i,j});
      end
      fprintf(FILE, [type{end} '\n'], in{i,end});
   end
   fprintf(FILE, '\n');
end

function serialize_block(FILE, block, type, in, comment)
   rows = size(in,1);
   serialize_block_or_array(FILE, block, rows, type, in, comment);
end

function serialize_array(FILE, block, type, in, comment)
   [rows, cols] = size(in);
   serialize_block_or_array(FILE, block, [rows cols], type, in, comment);
end

function serialize_sparse(FILE, block, in, comment)
   [i, j, v] = find(in);
   in = [i j v];
   serialize_block(FILE, block, {'%d','%d','%24.16e'}, in, comment);
end


function x = unit_test_svd_inv_solve(imdl, b, hp)
   eta = 1;
   J = calc_jacobian(imdl);
   [U,S,V] = svd(J, 'econ');
   Vt = V';
   sv = diag(S);
   Utb = U'*b;
   sv_f = (sv.^2)./((sv.^2) + hp);
   x = V * (sv_f .* (Utb ./ (sv*eta + eps) ));
end

function out = do_unit_test()
   ENABLE_C2F = 0;
   eidors_cache clear;
   for sel = 1:4
      switch sel
      case 4 % complex case
         NE = 16;
         imdl_sel = 'c2c';
         fmdl_sel = 'd2c';
         disp(' --- 16x CEM + large mesh ---');
      case 3 % complex case
         NE = 16;
         imdl_sel = 'd2d4c';
         fmdl_sel = 'f2d4c';
         disp(' --- 16x CEM + large mesh ---');
      case 2 % simpler model, for debug (CEM)
         NE = 4;
         imdl_sel = 'a2C';
         fmdl_sel = 'b2C';
         disp(' --- 4x CEM ---');
      case 1 % simpler model, for debug (PEM)
         NE = 4;
         imdl_sel = 'a2c';
         fmdl_sel = 'b2c';
         disp(' --- 4x PEM ---');
      otherwise
         error('oops');
      end
      imdl = mk_common_model(imdl_sel, NE);
      imdl.fwd_model.boundary = find_boundary(imdl.fwd_model);
      if ENABLE_C2F
         fmdl = mk_common_model(fmdl_sel, NE);
         fmdl.fwd_model.boundary = find_boundary(fmdl.fwd_model);
         imdl.rec_model = imdl.fwd_model;
         imdl.fwd_model = fmdl.fwd_model;
         c2f = mk_coarse_fine_mapping( imdl.fwd_model, imdl.rec_model);
         f2c = mk_coarse_fine_mapping( imdl.rec_model, imdl.fwd_model);
         c2f = c2f ./ sum(c2f,2);
         imdl.fwd_model.coarse2fine = c2f;
      end
      clf; show_fem(imdl.fwd_model)
      imgh = mk_image(imdl,1);
      imgi = create_inclusion(imgh, [+0.3 +0.3], 0.2, 0.9);
      imgi = create_inclusion(imgi, [-0.4 -0.1], 0.1, 1.1); 
      clf; show_fem(imgi,1);
      vh = fwd_solve(imgh);
      vi = fwd_solve(imgi);
      imdl.RtR_prior = 'prior_tikhonov';
      imdl.hyperparameter.value = 1e-5;
      t=tic;
      J=calc_jacobian(imgh);
      fprintf('calc_jacobian in %0.2f s\n', toc(t));
      img = inv_solve(imdl, vh, vi);
      clf;subplot(121); show_fem(imgi,1); subplot(122); show_fem(img,1);
      vz = fwd_solve(img);
      imge = imdl;
      imge.elem_data = [imgh.elem_data imgi.elem_data img.elem_data];
      zedhat(imge, sprintf('m%d.zh',sel));
      zedhat({vh, vi, vz}, sprintf('m%d.zh',sel), 'a');

%      % create model in EIDORS
%      disp('test: test.zh write');
%      imdl = mk_common_model(imdl_sel, NE);
%      imdl.fwd_model.boundary = find_boundary(imdl.fwd_model);
%      % generate difference data
%      %%fmdl = mk_common_model(fmdl_sel, NE);
%      %%fmdl.fwd_model.boundary = find_boundary(fmdl.fwd_model);
%      %%imdl.rec_model = imdl.fwd_model;
%      %%imdl.fwd_model = fmdl.fwd_model;
%      %%c2f = mk_coarse_fine_mapping( imdl.fwd_model, imdl.rec_model);
%      %%f2c = mk_coarse_fine_mapping( imdl.rec_model, imdl.fwd_model);
%      %%c2f = c2f ./ sum(c2f,2);
%      %%imdl.fwd_model.coarse2fine = c2f;
%      imgh = mk_image(imdl,1);
%      imgi = create_inclusion(imgh, [+0.3 +0.3], 0.2, 0.9);
%      imgi = create_inclusion(imgi, [-0.4 -0.1], 0.1, 1.1);
%      %%imgh.elem_data = f2c*imgh.elem_data;
%      %%imgi.elem_data = f2c*imgi.elem_data;
%      t=tic;
%      J=calc_jacobian(imgh);
%      fprintf('calc_jacobian in %0.2f s\n', toc(t));
%      vh = fwd_solve(imgh);
%      vi = fwd_solve(imgi);
%      % generate expected solution for comparison
%      imdl.RtR_prior = 'prior_tikhonov';
%      imdl.hyperparameter.value = 1e-6;
%      imdl.elem_data = imgh.elem_data;
%      imge = inv_solve(imdl, vh, vi);
%      img_eidors_inv_solve = imge;
%      x = unit_test_svd_inv_solve(imdl, vi.meas - vh.meas, imdl.hyperparameter.value);
%      fprintf('||SVD - inv_solve()||=%g\n', norm(x - imge.elem_data));
%      imdl.elem_data = [imgh.elem_data imgi.elem_data x];
%      imgt = imge; imgt.elem_data = x; ve = fwd_solve(imgt);
%      zedhat(imdl, sprintf('m%d.zh',sel));
%      disp('test: test.zh write + append');
%      zedhat({vh, vi, ve}, sprintf('m%d.zh',sel), 'a');
%      % we'd solve in zedhat, but instead we just check it all loads back safely
%      % load solution
%      disp('test: test.zh read');

      % test read-back of file
      imgz = zedhat(sprintf('m%d.zh',sel));
      clf; subplot(221); show_fem(imgi,1); xlabel('fwd');
           if isfield(imge, 'rec_model')
              imgr = imgi; imgr.fwd_model = imge.rec_model;
               subplot(222); show_fem(imgr,1); xlabel('rec');
           else
              imgr = imgh; imgr.elem_data = imge.elem_data(:,3);
              subplot(222); show_fem(imgr,1); xlabel('soln');
           end
           imge2 = imgi; imge2.elem_data = imge.elem_data(:,2);
           imgz2 = imgi; imgz2.elem_data = imgz.elem_data(:,2);
           subplot(223); show_fem(imge2,1); xlabel('eidors'); subplot(224); show_fem(imgz2,1); xlabel('zedhat');
      assert(norm(imgz.elem_data(:) - imge.elem_data(:)) < eps*size(imge.elem_data,1), 'ERROR: zedhat-eidors img.elem_data comparison failed');
      assert(norm(imgz.fwd_model.nodes(:) - imge.fwd_model.nodes(:)) < eps*size(imge.fwd_model.nodes(:),1), 'ERROR: zedhat-eidors img.fmdl.nodes comparison failed');
      assert(norm(double(imgz.fwd_model.elems(:)) - imge.fwd_model.elems(:)) < eps*size(imge.fwd_model.elems(:),1), 'ERROR: zedhat-eidors img.fmdl.elems comparison failed');
      assert(norm(double(imgz.fwd_model.boundary(:)) - double(imge.fwd_model.boundary(:))) < eps*size(imge.fwd_model.boundary(:),1), 'ERROR: zedhat-eidors img.fmdl.boundary comparison failed');
      assert(norm(imgz.hyperparameter.value - imge.hyperparameter.value) < eps*size(imge.fwd_model.boundary(:),1), 'ERROR: zedhat-eidors img.hyperparameter.value comparison failed');
      assert(norm(stim_meas_list(imgz.fwd_model.stimulation) - stim_meas_list(imge.fwd_model.stimulation)) < eps*length(imge.fwd_model.stimulation)*100, 'ERROR: zedhat-eidors img.fmdl.stimulation comparison failed');
      for e = 1:length(imge.fwd_model.electrode)
         assert(norm(imgz.fwd_model.electrode(e).nodes - imgz.fwd_model.electrode(e).nodes) < eps*10, 'ERROR: zedhat-eidors img.fmdl.electrode comparison failed');
      end
      out(sel).imgz = imgz;
      out(sel).imge = imge;
      %out(sel).imge = img_eidors_inv_solve;

      vd = fwd_solve(imgz);
      %vd.meas
      %subplot(222); plot(vd.meas); xlabel('vd');
   end
end
