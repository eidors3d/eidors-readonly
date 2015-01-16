function list = stim_meas_to_list(stim)
% function list = stim_meas_to_list(stim)
% Converts an EIDORS stimulation structure (stimulus and measurement selection
% matrices) into a stimulus and measurement list along with the measurement
% gain for each.
%
%  [ stim+  stim-  meas+  meas-  stim_current  meas_gain ]
%
% For each stimulus and measurement, a pair of electrodes is assumed. Patterns
% that do not perform pair drive and measurement will result in an error.
%
% (C) A. Boyle 2014

if ~isstruct(stim)
  error('expecting EIDORS stimulus strucutre');
end
if ~isfield(stim, 'stim_pattern')
  error('expecting EIDORS stimulus structure, missing stim_pattern');
end
if ~isfield(stim, 'meas_pattern')
  error('expecting EIDORS stimulus structure, missing meas_pattern');
end
if ~isfield(stim, 'stimulation')
  error('expecting EIDORS stimulus structure, missing stimulation');
end

ne = 0;
nm = 0;
for i = 1:length(stim)
  nel = size(stim(i).stim_pattern,1);
  ne = max([ne nel]);
  nm = nm + size(stim(i).meas_pattern,1) * size(stim(i).stim_pattern,2);
  if i == 1
    stim_units = stim(i).stimulation;
  else
    if ~strcmp(stim_units, stim(i).stimulation)
      error('stimulation units mismatch');
    end
  end
end

list = zeros(nm, 6);
nm = 1;
for i = 1:length(stim)
  nml = size(stim(i).meas_pattern,1);
  [mp,tmp] = find(stim(i).meas_pattern' > 0);
  [mm,tmp] = find(stim(i).meas_pattern' < 0);
  if length(mp) ~= nml
    error('meas_pattern not pairwise measurements');
  end
  if length(mm) ~= length(mp)
    error('meas_pattern not pairwise stimulations');
  end
  nml = size(stim(i).stim_pattern,2);
  [sp,tmp] = find(stim(i).stim_pattern > 0);
  [sm,tmp] = find(stim(i).stim_pattern < 0);
  if length(sp) ~= nml
    error('stim_pattern not pairwise stimulations');
  end
  if length(sm) ~= length(sp)
    error('stim_pattern not pairwise stimulations');
  end
  for j = 1:size(stim(i).stim_pattern,2);
    for k = 1:size(stim(i).meas_pattern,1);
      list(nm,:) = [sp(j) sm(j) mp(k) mm(k) stim(i).stim_pattern(sp(j),j) stim(i).meas_pattern(k,mp(k)) ];
      nm = nm +1;
    end
  end
end
