function idxr = reciprocity_idx( fmdl);
% RECIPROCITY_IDX: find indices of stim, meas pairs that are recirocal
%     ie. stimulation/measurement is same as measurement/stimulation on other
% usage:  idx = reciprocity_idx( fwd_model);
%  fmdl: an eidors fwd_model structure
%  idx:  index of corresponding reciprocal pairs
%   ie measurement #3 is reciprocal to idx(3)
%
%  if a measurement has no reciprocal pair. idx(m) = NaN;
%
% example 
%     imdl= mk_common_model('a2c2',8);
%     idx = reciprocity_idx(imdl.fwd_model);

% (C) 2010 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isstr(fmdl) && strcmp(fmdl,'UNIT_TEST'); do_unit_test; return; end

sm_0= []; sm_r = [];
for stim = fmdl.stimulation(:)';
  sp = stim.stim_pattern';
  nsp = sum(abs(sp))/2; %norm so max unidirectional current is 1
  mp = stim.meas_pattern;
  sp = ones(size(mp,1),1)*sp;
  sm_0= [sm_0; [sp/nsp, mp*nsp]];
  sm_r= [sm_r; [mp, sp]];
end
idxr = ones(size(sm_0,1),1);
for i=1:size(idxr)
  sm0i= ones(size(sm_0,1),1)*sm_0(i,:);
  mm= all( abs(sm0i - sm_r) < 1e-10, 2);
  mp= all( abs(sm0i + sm_r) < 1e-10, 2);
  m = [find(mm), find(mp)];
  if length(m)>1;
      meas= sprintf('%d,',m); 
      error('More than one reciprocal measure %d=>(%s). Giving up', i,meas);
  elseif length(m)==0;
      idxr(i) = NaN;
  else
      idxr(i) = m;
  end
end

function do_unit_test

%     [01] [12] [23] [34] [45] [56 [67] [70]
% [01] X    X    15   19   23   27  31   X
% [12] X    X    X    20   24   28  32   36
% [23] 1    X    X    X    25   29  33   37
% [34] 2    6    X    X    X    30  34   38
% [45] 3    7    11   X    X    X   35   39
% [56] 4    8    12   16   X    X   X    40
% [67] 5    9    13   17   21   X   X    X
% [70] X    10   14   18   22   26  X    X
tst.stimulation = mk_stim_patterns(8,1,[0,1],[0,1],{'rotate_meas'},1);
idx = reciprocity_idx( tst );
reshape(idx,5,8) 
