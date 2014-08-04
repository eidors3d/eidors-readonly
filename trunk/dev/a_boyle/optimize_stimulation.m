function [stim, vsel] = optimize_stimulation(stim, verbose);
% function [stim, vsel] = optimize_stimulation(stim);
% This function optimizes the stimulation and measurement
% patterns of a stim = fwd_model.stimulation(:) to be more
% computationally efficient in a MATLAB environment where
% for loops are poorly optimized and large matrix operations
% can be performed efficiently in parallel.
%
% The optimized stimulation is returned in 'stim'. The
% optimization most likely results in reordering the
% measurement vector returned from a fwd_solve() using the
% new patterns. A matrix 'vsel' is provided to reverse this
% reordering.
%   old_meas = fwd_solve(fmdl); % SLOW (~45min?)
%   [fmdl.stimulation, vsel] = optimize_stimulations(fmdl.stimulation(:));
%   new_meas = fwd_solve(fmdl); % FAST (~1min?)
%   new_meas.meas = vsel * new_meas.meas;
%   % where new_meas == old_meas
% Times were observed on a 32-electrode linear array using
% simple 4-electrode patterns across 32-electrodes with
% electrode movements.
%
% The optimizations are
%  (1) reduce length(stimulation) by combining two
%      stim_patterns that are the same, the respective
%      meas_patterns are merged
%  (2) reduce the total measurements calculated by combining
%      duplicate meas_patterns within a single stimulation
%  (3) reduce length(stimulation) by merging two
%      stim_patterns that have the same meas_pattern, the
%      respective meas_patterns are combined, optionally
%      some computational waste can be introduced by merging
%      nearly matching meas_patterns and dropping the
%      unnecessary voltages from vsel. (Currently DISABLED.)
% stim is the optimized stimulation structure, vsel is a matrix to multiply vi.meas
% by to get what would have been the expected voltage sequence prior to this
% optimization
% (C) 2013 Alistair Boyle. License: GPL version 2 or version 3

  % comment these out to disable each of the optimizations
  MERGE_STIM_PATTERNS=1;
  MERGE_SELF_MEAS_PATTERNS=1;
%  MERGE_MEAS_PATTERNS=1;

  if nargin < 2
    verbose = 1;
  end

  % branch off to do unit test if requested
  if ischar(stim) && strcmp(stim, 'UNIT_TEST');
    vsel = 0;
    stim = do_unit_test(); % pass = f()
    return;
  elseif ischar(stim) && strcmp(stim, 'SINGLE_TEST');
    vsel = 0;
    stim = do_single_test(verbose);
    return;
  end

  % stim/meas_pattern near match threshold
  % (optimization #3, MERGE_MEAS_PATTERNS)
  % 0 = 0% - only merge stim_patterns that have exactly the same meas_pattern
  % 1 = 100% - merge everything into a single stimulation
  waste_thres = 0;

  % count measurements in we will get from fwd_solve
  nst = length(stim);
  nvt = 0;
  for i = 1:length(stim);
    nvt = nvt + size(stim(i).stim_pattern, 2) * size(stim(i).meas_pattern, 1);
  end
  % default - no changes
  vsel = speye(nvt);
  % V = vsel*VV;
  % where V is the original measurements and VV are the optimizated measurements
  % with V and VV as column vectors

  if verbose
    fprintf('  stim/meas pattern optimization\n');
    fprintf('    %d stim_patterns, %d measurements\n', nst, nvt);
  end

  % stimulation patterns are by row
  % (1,1) = 1
  % (2,1) = -1
  % measurement patterns are by column
  % (1,3) = -1
  % (1,4) = 1
  %
  % we combine stimulation pattterns first since they are the
  % most expensive (non-matrix) loop in the forward solution
  %
  % NOTE -- this is kind of a sh*ty bubble sort but it does the trick for now
  if exist('MERGE_STIM_PATTERNS')
    if verbose
      fprintf('    1. stim merge\n');
      fprintf('      stim');
    end
    len_stim_start = length(stim);
    len_meas_start = size(vsel,2);
    na = 1;
    i = 1;
    while i <= length(stim)
      if verbose
        fprintf(' %d', i);
      end
      j = i+1;
      % count current measurements, to track insertion point
      ni = size(stim(i).stim_pattern, 2) * size(stim(i).meas_pattern, 1);
      na = na + ni; % last entry that is fixed
      % and reset deletion point
      nd = na; % first entry of the deletion/move-to-the-front list
      while j <= length(stim)
        % number of measurements for this stim pattern
        nj = size(stim(j).stim_pattern, 2) * size(stim(j).meas_pattern, 1);
        % stim patterns match --> remove redundant stim and merge measurement patterns
        stim_match = all(size(stim(i).stim_pattern) == size(stim(j).stim_pattern)) && ...
                     all(all(stim(i).stim_pattern == stim(j).stim_pattern));
        stim_match_inv = all(size(stim(i).stim_pattern) == size(stim(j).stim_pattern)) && ...
                         all(all(stim(i).stim_pattern == -stim(j).stim_pattern));
        % TODO possibly fix fwd_model_parameter: it and the
        % functions downstream from there can't handle
        % stimulations with more than one column even though
        % it concatenates them together. For now we avoid the
        % case where we would optimize this by combining
        % stim_pattern columns
        if (~stim_match) && (~stim_match_inv)
          meas_match = 0; % HACK, see TODO text above
        else
          meas_match = all(size(stim(i).meas_pattern) == size(stim(j).meas_pattern)) && ...
                       all(all(stim(i).meas_pattern == stim(j).meas_pattern));
          meas_match_inv = all(size(stim(i).meas_pattern) == size(stim(j).meas_pattern)) && ...
                           all(all(stim(i).meas_pattern == -stim(j).meas_pattern));
        end

        if (stim_match && meas_match) || ...
           (stim_match_inv && meas_match_inv)
          if (verbose>= 2)
            fprintf('x%d',j);
          end
          % merge vsel and drop the stim(j)
          [vsel, na, nd] = vsel_merge_col(vsel, na, nd, nj);
          % delete stim pattern, don't increment j
          stim(j) = [];
        elseif stim_match || stim_match_inv
          if (verbose>= 2)
            fprintf('s%d',j);
          end
          % append measurement pattern, columns
          new_meas = [1:size(stim(j).meas_pattern, 1)] + size(stim(i).meas_pattern,1);
          stim(i).meas_pattern(new_meas,:) = stim(j).meas_pattern;
          % update voltage selector
          inv = +1;
          if stim_match_inv
            inv = -1; % invert the measured voltages if the stimulus is inverted
          end
          [vsel, na, nd] = vsel_move_col(vsel, na, nd, nj, inv);
          % delete stim pattern, don't increment j
          stim(j) = [];
        elseif meas_match
          if (verbose>= 2)
            fprintf('m%d',j);
          end
          error('meas_match code is broken FIXME');
          % append stimulation pattern, rows
          new_stim = [1:size(stim(j).stim_pattern, 2)] + size(stim(i).stim_pattern,2);
          stim(i).stim_pattern(:,new_stim) = stim(j).stim_pattern;
          [vsel, na, nd] = vsel_move_col(vsel, na, nd, nj);
          % delete stim pattern, don't increment j
          stim(j) = [];
        else % nothing changed
          j = j +1; % nothing found: next
          nd = nd + nj; % incr. deletion point
        end
      end
      i = i +1; % nothing found: next
      if verbose >= 2
        fprintf('\n');
      end
    end
    len_stim_end = length(stim);
    len_meas_end = size(vsel, 2);
    if verbose
      fprintf('\n');
      fprintf('      stim pattern reduction %d -> %d\n', len_stim_start, len_stim_end);
      fprintf('      meas pattern reduction %d -> %d\n', len_meas_start, len_meas_end);
    end
  end

  % next we look for duplicate measurement patterns within any stim
  if exist('MERGE_SELF_MEAS_PATTERNS')
    if verbose
      fprintf('    2. meas merge - within a stim pattern\n');
      fprintf('      stim');
    end
    len_meas_start = size(vsel, 2);
    na = 0;
    i=1;
    while i <= length(stim)
      if verbose
        fprintf(' %d', i);
      end
      % remove duplicates within each measurement pattern
      j=1;
      while j <= size(stim(i).meas_pattern, 1)
        k=j+1;
        while k <= size(stim(i).meas_pattern, 1)
          % check for a pattern match
          if all(stim(i).meas_pattern(j, :) == stim(i).meas_pattern(k, :))
            stim(i).meas_pattern(k, :) = []; % delete pattern
            % update vsel, merge duplicated measurement
            vsel(:,na+j) = vsel(:,na+j) + vsel(:,na+k);
            vsel(:,na+k) = [];
          elseif all(stim(i).meas_pattern(j, :) == -stim(i).meas_pattern(k, :))
            stim(i).meas_pattern(k, :) = []; % delete pattern
            % update vsel, merge duplicated measurement
            vsel(:,na+j) = vsel(:,na+j) - vsel(:,na+k);
            vsel(:,na+k) = [];
          else
            k=k+1; % try the next pattern
          end
        end
        j=j+1;
      end
      ni = size(stim(i).stim_pattern, 2) * size(stim(i).meas_pattern, 1);
      na = na + ni; % count finalized patterns
      i = i +1; % nothing found: next
    end
    len_meas_end = size(vsel, 2);
    if verbose
      fprintf('\n');
      fprintf('      meas pattern reduction %d -> %d\n', len_meas_start, len_meas_end);
    end
  end

% TODO see comment about fwd_model_parameters above -- it
% can't handle multiple column stim_pattern, which leaves
% this optimization unavailable until it is fixed

  % check for matching measurement patterns across all stimulation patterns
  if exist('MERGE_MEAS_PATTERNS')
    if verbose
      fprintf('    3. stim/meas merge - merge across stim patterns: %d%% waste threshold\n', waste_thres*100);
      fprintf('      stim');
    end
    len_meas_start = size(vsel, 2);
    len_stim_start = length(stim);
    na = 1;
    i=1;
    while i <= length(stim)
      if verbose
        fprintf(' %d', i);
      end
      stim_i = size(stim(i).stim_pattern, 2);
      meas_i = size(stim(i).meas_pattern, 1);
      ni = stim_i * meas_i;
      na = na + ni; % count finalized patterns
      nd = na; % start of next measurements
      j=i+1;
      while j <= length(stim)
        stim_j = size(stim(j).stim_pattern, 2);
        meas_j = size(stim(j).meas_pattern, 1);
        nj = stim_j * meas_j;
        % and reset deletion point
        % count matches between patterns
        meas1 = stim(i).meas_pattern;
        meas2 = stim(j).meas_pattern;
        [meas, meas_j_sel, waste_i, waste_j] = merge_meas(meas1, meas2);
        if all(size(meas1) == size(meas2)) && ...
           all(all(meas1 == meas2)) % full match
          % TODO a full match, including ordering...
        end
        % merge stim patterns if its 'useful'
        % if we merge, approximately how many useless measurements will we calculate?
        % from this merger, it would be:
        waste = waste_i*stim_i + waste_j*stim_j;
        calcs = (meas_i + meas_j)*(stim_i + stim_j);
        if waste/calcs <= waste_thres % is that more than our threshold?
          % commit our merged pattern
          % [meas, meas2_sel, waste1, waste2] = merge_meas(meas1, meas2);
          stim(i).meas_pattern = meas;
          % and update vsel to match
          si = stim_i;
          mi = meas_i;
          mi_p = size(meas, 2);
          sj = stim_j;
          mj = meas_j;
          vsel = merge_vsel(vsel, na, nd, si, mi, mi_p, sj, mj, meas_j_sel);
          % then merge the stimulation patterns
          stim(i).stim_pattern(:,end+1:end+stim_j) = stim(j).stim_pattern;
          % and finally, delete the old stim
          stim(j) = [];
          na = na + nj; % count finalized patterns
          % and don't increment the j counter since we just vaporized the old j entry
        else
          j=j+1;
        end
        nd = nd + nj; % count finalized patterns
      end
      i = i +1; % nothing found: next
    end
    len_meas_end = size(vsel, 2);
    len_stim_end = length(stim);
    len_meas_waste = full(sum(sum(abs(vsel)) == 0)); % count all zero columns
    if(verbose)
      fprintf('\n');
      fprintf('      meas pattern reduction %d -> %d\n', len_meas_start, len_meas_end);
      fprintf('      meas pattern waste 0 -> %d\n', len_meas_waste);
      fprintf('      stim pattern reduction %d -> %d\n', len_stim_start, len_stim_end);
    end
  end

  if verbose
    fprintf('    SUMMARY stimulation loops       %d -> %d (%d%%)\n', nst, length(stim), floor((1-length(stim)/nst)*100));
    fprintf('    SUMMARY voltage meas calculated %d -> %d (%d%%)\n', size(vsel,1), size(vsel,2), floor((1-size(vsel,2)/size(vsel,1))*100));
  end

  % TODO look for inversions of stim or meas pattern (vsel(x,y)=-1;)
  % TODO look for inversions of both stim and meas pattern (vsel(x,y)=+1;)

% insert nj entries starting at nd prior to na
% mult allows for inverted patterns
function [vsel, na , nd] = vsel_move_col(vsel, na, nd, nj, mult)
  if nargin < 5
    mult = 1;
  end
  if abs(mult) ~= 1
    error('expected mult to be +/-1');
  end
  if (nj == 0)
    return;
  end
  if nd < na
    error('vsel_move_col() with nd prior to na');
  end
  if nj < 0
    error('vsel_move_col() with nj < 0');
  end
  nvt = size(vsel,2);
  if (nd + nj -1) > nvt % expand vsel with empty columns if required
    e = (nd + nj -1) - nvt;
    r = size(vsel, 1);
    vsel(:,end+1:end+e) = zeros(r,e);
    nvt = size(vsel,2);
  end
  if mult == -1
    % flip measurements if its a reciprocal measurement
    vsel(:, [nd:nd+nj-1]) = -vsel(:, [nd:nd+nj-1]);
  end
  % The KEY moment:
  % Reorder to {MOVED, OLD_START, OLD_END}
  % where MOVED are elements nd:nd+nj-1.
  new = [nd:nd+nj-1 na:nd-1 nd+nj:nvt]; % new order
  old = [na:nvt]; % old order
  vsel(:, old) = vsel(:, new);

  % adjust insert pointer for new patterns
  na = na + nj;
  % update deletion point
  nd = nd + nj;

% merge the nj columns starting at nd into the columns prior to na
function [vsel, na, nd] = vsel_merge_col(vsel, na, nd, nj)
  if (nj == 0)
    return;
  end
  if nd < na
    error('vsel_merge_col() with nd prior to na');
  end
  if nj < 0
    error('vsel_merge_col() with nj < 0');
  end
  d = nd : nd+nj-1; % delete these
  m = na - nj : na -1; % after merging with these
  if length(m) ~= length(d)
    error('vsel_merge_col() with m != d');
  end
  if max(d) > size(vsel, 2)
    fprintf('vsel_merge_col() vsel[%d,%d] - max delete[%d], max merge[%d]',...
            size(vsel,1), size(vsel, 2), ...
            max(d), max(m));
    error('vsel_merge_col() index beyond max range');
  end
  vsel(:,m) = vsel(:,m) + vsel(:,d);
  vsel(:,d) = [];

% append non-matching entries from meas2 to meas1
% also returns a selector to get meas2 entries from the returned result
% note that meas1_sel would be meas1_sel=1:size(meas1, 1)
function [meas, meas2_sel, waste1, waste2] = merge_meas(meas1, meas2);
  % meas2 row select from meas
  meas2_sel = zeros(1,size(meas2, 1)); % preallocate
  meas = meas1;
  for i = 1:size(meas2, 1)
    match = 0;
    for j = 1:size(meas1, 1)
      % check for a pattern match
      if all(meas2(i, :) == meas1(j, :))
        match = j;
      end
    end
    if match == 0 % no match
      % append to meas1
      meas(end+1, :) = meas2(i, :);
      meas2_sel(i) = size(meas, 1); % last entry
    else % match
      meas2_sel(i) = match; % remember which one it is that matches
    end
  end
  waste1 = size(meas,1) - size(meas1,1);
  waste2 = size(meas,1) - size(meas2,1);

% na is the first entry of the original measurements in vsel
% mi is the original number of measurements
% nd is the first entry of the merge measurements
% mj is the number of meas2 measurements prior to the merger
% meas2_sel is the selection matrix for nd starting at na, after the new patterns are appended
% si is the original stimulation pattern count
% sj is the new stimulation pattern count
function vsel = UNUSED_merge_vsel(vsel, na, nd, si, mi, mi_p, sj, mj, meas2_sel)
  ni = si*mi; % former total number of measurements in meas1
  nj = sj*mj; % former total number of measurements in meas2

  % how many extra measurements are calculated per stimulation pattern?
  waste_i = mi_p - mi;
  waste_j = mi_p - mj;

  % first pad out vsel where any new waste will occur
  % this will happen for each stim pattern where there are new voltage
  % patterns that weren't required before
  for i=0:si-1
    vsel = vsel_move_col(vsel, na+mi_p*i+mi, size(vsel,2), waste_i);
  end
  na = na + (mi_p*si); % new start of meas2 patterns
  nd = nd + (waste_i*si);
  % move the stim2*meas2 columns to the correct place
  vsel = vsel_move_col(vsel, na, nd, sj*mj);
  % pad out the meas2 columns per stim2 for waste
  % first find which columns of meas2 will become waste
  tmp = ones(1,mi_p);
  tmp(meas2_sel) = 0;
  meas2_pad = find(tmp == 1);
  % then for each stim pattern
  for i=0:sj-1
    % add zero columns to drop the waste
    vsel = vsel_move_col(vsel, na+mi_p*i+mj, size(vsel,2), waste_j);
    % and reorder the columns so that they match the measurements
    t = (na+mi_p*i) -1;
    vsel(:, t + [meas2_sel meas2_pad]) = vsel(:, t + [1:mi_p]);
  end

function pass = do_unit_test();
  pass = 1;

  % move column 3&4 prior to column 2
  vsel_expected = speye(4);
  vsel_expected = vsel_expected(:,[1 3 4 2]);
  [vsel, na, nd] = vsel_move_col(speye(4),2,3,2);
  if any(any(vsel ~= vsel_expected))
    fprintf('error: vsel_move_col() not working as expected, columns in wrong order');
    pass = 0;
  else
    fprintf('vsel_move_col() gave the correct column order\n');
  end
  if na ~= 4
    fprintf('error: vsel_move_col() not working as expected, updated na=%d\n',na);
    pass = 0;
  else
    fprintf('vsel_move_col() gave the correct updated na\n');
  end
  if nd ~= 5
    fprintf('error: vsel_move_col() not working as expected, updated nd=%d\n',nd);
    pass = 0;
  else
    fprintf('vsel_move_col() gave the correct updated nd\n');
  end
  % merge column 3&4 prior to column 3
  vsel_expected = speye(4);
  vsel_expected = vsel_expected(:,[1 2]) + vsel_expected(:,[3 4]);
  [vsel, na, nd] = vsel_merge_col(speye(4),3,3,2);
  if any(size(vsel_expected) ~= size(vsel)) || any(any(vsel ~= vsel_expected))
    fprintf('error: vsel_merge_col() not working as expected, columns in wrong order');
    vsel=vsel
    vsel_expected=vsel_expected
    pass = 0;
  else
    fprintf('vsel_merge_col() gave the correct column order\n');
  end
  if na ~= 3
    fprintf('error: vsel_merge_col() not working as expected, updated na=%d\n',na);
    pass = 0;
  else
    fprintf('vsel_merge_col() gave the correct updated na\n');
  end
  if nd ~= 3
    fprintf('error: vsel_merge_col() not working as expected, updated nd=%d\n',nd);
    pass = 0;
  else
    fprintf('vsel_merge_col() gave the correct updated nd\n');
  end

  % some pathological cases for debugging
  imdl = mk_common_model('n3r2',[16,2]);
  img = mk_image(imdl);
  select_fcn = inline('(x-0.2).^2+(y-0.5).^2<0.1^2','x','y','z');
  img.elem_data = 1 + elem_select(img.fwd_model, select_fcn);
  for i=1:4
    s(i).stimulation = 'Amp';
    s(i).stim_pattern = zeros(32,1);
    s(i).stim_pattern(1) =  1;
    s(i).stim_pattern(8) = -1;
    s(i).meas_pattern = zeros(1,32);
    s(i).meas_pattern(2)   =  1;
    s(i).meas_pattern(i+9) = -1;
  end
  img.fwd_model.stimulation = s;
  pass = unit_test_run(pass, img, 'n3r2 3D, 16x2 electrodes, all the same stim_pattern');
  for i=1:4
    s(i).meas_pattern = zeros(1,32);
    s(i).meas_pattern(2)  =  1;
    s(i).meas_pattern(12) = -1;
  end
  img.fwd_model.stimulation = s;
  pass = unit_test_run(pass, img, 'n3r2 3D, 16x2 electrodes, single repeated pattern');
  so = s;
  for i=1:4
    s = so;
    s(i).stim_pattern(3) =  1;
    s(i).stim_pattern(7) = -1;
    img.fwd_model.stimulation = s;
    pass = unit_test_run(pass, img, sprintf('n3r2 3D, 16x2 electrodes, almost the same stim_pattern except %d',i));
  end
  s = so; % reciprocal stimulus
  for i=3:4
    s(i).stim_pattern(3) =  1;
    s(i).stim_pattern(7) = -1;
    img.fwd_model.stimulation = s;
    pass = unit_test_run(pass, img, sprintf('n3r2 3D, 16x2 electrodes, almost the same stim_pattern except %d',i));
  end
  for i=1:4
    s = so; % reciprocal measurements
    s(i).stim_pattern = -s(i).stim_pattern;
    img.fwd_model.stimulation = s;
    pass = unit_test_run(pass, img, sprintf('n3r2 3D, 16x2 electrodes, almost the same stim_pattern except reciprocal stim %d',i));
  end
  for i=1:4
    s = so; % reciprocal measurements
    s(i).meas_pattern = -s(i).meas_pattern;
    img.fwd_model.stimulation = s;
    pass = unit_test_run(pass, img, sprintf('n3r2 3D, 16x2 electrodes, almost the same stim_pattern except reciprocal meas %d',i));
  end
  for i=1:4
    s = so;
    s(i).meas_pattern(5)  =  1;
    s(i).meas_pattern(12) = -1;
    img.fwd_model.stimulation = s;
    pass = unit_test_run(pass, img, sprintf('n3r2 3D, 16x2 electrodes, almost the same meas_pattern except %d',i));
  end
  for i=1:4
    s(i).stim_pattern = zeros(32,1);
    s(i).stim_pattern(1) =  1;
    s(i).stim_pattern(i+6) = -1;
    s(i).meas_pattern = zeros(1,32);
    s(i).meas_pattern(2)   =  1;
    s(i).meas_pattern(i+9) = -1;
  end
  img.fwd_model.stimulation = s;
  pass = unit_test_run(pass, img, 'n3r2 3D, 16x2 electrodes, all different');

  % try a generic model and stim/meas patttern
  % run it through the optimization and compare the fwd_solve results
  imdl= mk_common_model('d2d1c',19);
  img = mk_image(imdl);
  select_fcn = inline('(x-0.2).^2+(y-0.5).^2<0.1^2','x','y','z');
  img.elem_data = 1 + elem_select(img.fwd_model, select_fcn);
  img.fwd_model.stimulation = mk_stim_patterns(19,1,[0,1],[0,1],{},1);
  pass = unit_test_run(pass, img, 'd2d1c 2D, 19 electrodes');

  % another simple 2D circular model
  imdl = mk_common_model('c2C0',16);
  img = mk_image(imdl);
  select_fcn = inline('(x-0.2).^2+(y-0.5).^2<0.1^2','x','y','z');
  img.elem_data = 1 + elem_select(img.fwd_model, select_fcn);
  pass = unit_test_run(pass, img, 'c2C0 2D, 16 electrodes');

  % another simple 2D circular model
  imdl = mk_common_model('c2C0',16);
  img = mk_image(imdl);
  select_fcn = inline('(x-0.2).^2+(y-0.5).^2<0.1^2','x','y','z');
  img.elem_data = 1 + elem_select(img.fwd_model, select_fcn);
  s = mk_stim_patterns(16, 1, '{ad}','{op}', {}, 1);
  img.fwd_model.stimulation = s;
  pass = unit_test_run(pass, img, 'c2C0 2D, 16 electrodes w/ adjacent stim/meas');


  % a 3D model
  imdl = mk_common_model('n3r2',[16,2]);
  img = mk_image(imdl);
  select_fcn = inline('(x-0.2).^2+(y-0.5).^2<0.1^2','x','y','z');
  img.elem_data = 1 + elem_select(img.fwd_model, select_fcn);
  pass = unit_test_run(pass, img, 'n3r2 3D, 16x2 electrodes');

  % geophysics patterns
  % this stim pattern doesn't make much sense on a cylindrical 3D model
  % (its supposed to be a surface linear array) but it will do for this test
  % Wenner
  s = stim_pattern_geophys(32, 'Wenner', {'spacings', 1:32});
  imdl = mk_common_model('n3r2',[16,2]);
  img = mk_image(imdl);
  select_fcn = inline('(x-0.2).^2+(y-0.5).^2<0.1^2','x','y','z');
  img.elem_data = 1 + elem_select(img.fwd_model, select_fcn);
  img.fwd_model.stimulation = s;
  pass = unit_test_run(pass, img, 'n3r2 3D, 16x2 electrodes, Wenner stimulation');
  % Schlumberger
  spacing= [1 1 1 2 3 4 6 8 8 11 12 14 17]; % a
  multiples= [1 2 3 2 5/3 6/4 7/6 1 10/8 1 13/12 15/14 1];% n
  s = stim_pattern_geophys(32, 'Schlumberger', {'spacings', spacing,'multiples',multiples});
  imdl = mk_common_model('n3r2',[16,2]);
  img = mk_image(imdl);
  select_fcn = inline('(x-0.2).^2+(y-0.5).^2<0.1^2','x','y','z');
  img.elem_data = 1 + elem_select(img.fwd_model, select_fcn);
  img.fwd_model.stimulation = s;
  pass = unit_test_run(pass, img, 'n3r2 3D, 16x2 electrodes, Schumberger stimulation');
  % Dipole-dipole
  spacing= [1 1 1 2 3 3 4 4 5 6 6 7 8 8 9 10 10 11 12 12 13 14 14 15 16 17]; % a
  multiples= [1 2 3 2 1 5/3 1 2  1 1 7/6 1 1 10/8 1 1 12/10 1 1 13/12 1 1 15/14 1 1 1]; % n
  s = stim_pattern_geophys(32, 'DipoleDipole', {'spacings', spacing,'multiples',multiples});
  imdl = mk_common_model('n3r2',[16,2]);
  img = mk_image(imdl);
  select_fcn = inline('(x-0.2).^2+(y-0.5).^2<0.1^2','x','y','z');
  img.elem_data = 1 + elem_select(img.fwd_model, select_fcn);
  img.fwd_model.stimulation = s;
  pass = unit_test_run(pass, img, 'n3r2 3D, 16x2 electrodes, Dipole-dipole stimulation');

  % a nastier one that does electrode movement
  x = ones(1,8);
  S = [1*x 2*x 3*x 4*x];
  M = repmat(1:8,1,4);
  s=stim_pattern_geophys(32, 'DipoleDipole', {'spacings', S, 'multiples', M, 'reciprocal_meas',1});
  imdl = mk_common_model('n3r2',[16,2]);
  img = mk_image(imdl);
  select_fcn = inline('(x-0.2).^2+(y-0.5).^2<0.1^2','x','y','z');
  img.elem_data = 1 + elem_select(img.fwd_model, select_fcn);
  img.fwd_model.stimulation = s;
  pass = unit_test_run(pass, img, 'n3r2 3D, 32 electrodes w/ Dipole-dipole');
  % then add some electrode movements
  x = ones(1,8);
  S = [1*x 2*x 3*x 4*x];
  M = repmat(1:8,1,4);
  s=stim_pattern_geophys(16, 'DipoleDipole', {'spacings', S, 'multiples', M, 'reciprocal_meas',1});
  imdl = mk_common_model('n3r2',[16,2]);
  img = mk_image(imdl);
  select_fcn = inline('(x-0.2).^2+(y-0.5).^2<0.1^2','x','y','z');
  img.elem_data = 1 + elem_select(img.fwd_model, select_fcn);
  s = unit_test_expand_for_electrode_movement(s, 1);
  img.fwd_model.stimulation = s;
  pass = unit_test_run(pass, img, 'n3r2 2D, 16 electrodes w/ Dipole-dipole & movement in 1D');

  if pass == 1
    fprintf('\nPASSED optimize_stimulation() unit_test\n');
  else
    fprintf('\nFAILED optimize_stimulation() unit_test\n');
  end

function pass = unit_test_run(pass, img, desc);
  fprintf('%s\n',desc);
  % original solution
  t = tic;
  vo = fwd_solve(img);
  t_orig = toc(t);

  % room for error: numeric precision
  thres = eps(max(vo.meas))*length(vo.meas)*1e2;

  % and now rerun the fwd_solve after an optimize_stimulation()
  t = tic;
  [s, vsel] = optimize_stimulation(img.fwd_model.stimulation(:));
  img.fwd_model.stimulation = s;
  ts = tic;
  vn = fwd_solve(img);
  vn.meas = vsel*vn.meas; % correct for optimizations
  t_optim_solve = toc(ts);
  t_optim = toc(t);

  % calculate runtimes
  t_orig_minutes = floor(t_orig/60);
  t_orig_seconds = floor(t_orig - t_orig_minutes*60);
  t_orig_ms = floor((t_orig - t_orig_seconds - t_orig_minutes*60)*100);
  t_optim_minutes = floor(t_optim/60);
  t_optim_seconds = floor(t_optim - t_optim_minutes*60);
  t_optim_ms = floor((t_optim - t_optim_seconds - t_optim_minutes*60)*100);
  t_optim_solve_minutes = floor(t_optim_solve/60);
  t_optim_solve_seconds = floor(t_optim_solve - t_optim_solve_minutes*60);
  t_optim_solve_ms = floor((t_optim_solve - t_optim_solve_seconds - t_optim_solve_minutes*60)*100);
  spd = (1 - t_optim/t_orig)*100;
  if spd >= 0
    spd_str = 'faster';
  else
    spd_str = 'slower';
  end
  spd_solve = (1 - t_optim_solve/t_orig)*100;
  if spd_solve >= 0
    spd_solve_str = 'faster';
  else
    spd_solve_str = 'slower';
  end
  fprintf('  runtime: old=%dm:%ds:%dms, new=%dm:%ds:%dms (%d%% %s), solve only %dm:%ds:%dms (%d%% %s)\n', ...
          t_orig_minutes, t_orig_seconds, t_orig_ms, ...
          t_optim_minutes, t_optim_seconds, t_optim_ms, ...
          abs(round(spd)), spd_str, ...
          t_optim_solve_minutes, t_optim_solve_seconds, t_optim_solve_ms, ...
          abs(round(spd_solve)), spd_solve_str);

  dd = vo.meas - vn.meas;
  d = norm(dd); % difference in results
  if d > thres
    fprintf('FAIL: norm(diff)=%g > %g - %s\n', d, thres, desc);
    i = find(abs(dd) > thres, 10);
    if length(dd) > 5
      fprintf('  first %g mismatched measurements:\n',length(i));
      orig=vo.meas(i)
      err=dd(i)
    else
      fprintf('  all measurements, %g mismatches:\n',length(i));
      orig_opt_err=[vo.meas vn.meas dd]
      vsel=vsel
    end
    pass = 0;
  end

function s = unit_test_expand_for_electrode_movement(s, N);
  fprintf('  adjust patterns to estimate movement\n');
  % expand initial stim/meas matrices
  % te = total number of electrodes, after delta electrodes are added
  n_elec = size(s(1).stim_pattern,1);
  te = n_elec*(N+1);
  for k = 1:length(s)
    s(k).stim_pattern(te,:) = zeros(1,size(s(k).stim_pattern,2));
    s(k).meas_pattern(:,te) = zeros(size(s(k).meas_pattern,1),1);
  end
  so = s; % original stimulation(:)
  sn = length(s);
  for i = 1:N;
    fprintf('   ');
    % TODO we could do this in a smarter way by reducing the additional calculations
    % TODO to only those that have moved BUT first we'll do it the dumb way
    for j = 1:n_elec
      fprintf(' %s%d', 'U'+i-1, j);
      eio = j; % electrode index (original)
      ein = i*n_elec + j; % electrode index (new)
      % stimulation patterns are by row
      % (1,1) = 1
      % (2,1) = -1
      % measurement patterns are by column
      % (1,3) = -1
      % (1,4) = 1
      si1 = ((i-1)*n_elec+j)*sn +1; % stimulation index 1
      sie = si1 + sn -1; % stimulation index end
      s(si1:sie) = so; % duplicate all stimulations
      for k = si1:sie % for each new stimulation
        % move stimulation electrodes
        s(k).stim_pattern(ein,:) = s(k).stim_pattern(eio, :);
        s(k).stim_pattern(eio,:) = zeros(1,size(s(k).stim_pattern,2));
        % move measurement electrodes
        s(k).meas_pattern(:,ein) = s(k).meas_pattern(:,eio);
        s(k).meas_pattern(:,eio) = zeros(size(s(k).meas_pattern,1),1);
      end
    end
    fprintf('\n');
  end

function pass = do_single_test(stim);
  imdl = mk_common_model('n3r2',[16,2]);
  img = mk_image(imdl);
  img.fwd_model.stimulation = stim;
  img.fwd_model.electrode(:) = [];
  for i = 1:size(stim(1).stim_pattern,1)
    img.fwd_model.electrode(i).z_contact = 100;
    img.fwd_model.electrode(i).nodes = 10 + 3*i;
  end
%    select_fcn = inline('(x-0.2).^2+(y-0.5).^2<0.1^2','x','y','z');
%    img.elem_data = 1 + elem_select(img.fwd_model, select_fcn);

  pass = unit_test_run(1, img, 'SINGLE_TEST'); % pass = f(...)

  N = 1;
  s = img.fwd_model.stimulation;
  s = unit_test_expand_for_electrode_movement(s, N);
  img.fwd_model.stimulation = s;
  for i=length(img.fwd_model.electrode)+1:...
        length(img.fwd_model.electrode)*(N+1)
    img.fwd_model.electrode(i).z_contact = 100;
    img.fwd_model.electrode(i).nodes = 10 + 3*i;
  end
  pass = unit_test_run(1, img, 'SINGLE_TEST + electrode movements'); % pass = f(...)


  if pass == 1
    fprintf('\nPASSED optimize_stimulation() single_test\n');
  else
    fprintf('\nFAILED optimize_stimulation() single_test\n');
  end
