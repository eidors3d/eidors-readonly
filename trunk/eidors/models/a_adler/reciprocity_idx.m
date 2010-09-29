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
idx = reciprocity_idx( tst ); idx = reshape(idx,5,8);
do_indiv_test('8-[01]-[01]-rotate',idx(:,[1,5]), ...
    [15,19,23,27,31;35,39,3,7,11]');

tst.stimulation = mk_stim_patterns(8,1,[0,1],[0,1],{'rotate_meas'},10);
idx = reciprocity_idx( tst ); idx = reshape(idx,5,8);
do_indiv_test('8(10)-[01]-[01]-rotate',idx(:,[1,5]), ...
    [15,19,23,27,31;35,39,3,7,11]');

%     [01]   [12]   [23]   [34]   [45]   [50] 
% [01] X      X      7  9   10 11  13 13  X   
% [12] X      X      X      11 12  14 14  16 16 
% [23] 1  1   X      X      X      15 15  17 17 
% [34] 2  2   4  4   X      X      X      18 18 
% [45] 3  3   5  5   8  7   X      X      X   
% [50] X      6  6   9  8   12 10  X      X   
tst.stimulation = mk_stim_patterns(6,1,[0,1],[0,1],{'rotate_meas'},1);
idx = reciprocity_idx( tst ); idx = reshape(idx,3,6); 
do_indiv_test('6-[01]-[01]-rotate',idx(:,[1,4]), [9,11,13;18 2 4]');

tst.stimulation = mk_stim_patterns(6,1,[0,1],[0,1],{'no_rotate_meas'},1);
idx = reciprocity_idx( tst ); idx = reshape(idx,3,6); 
do_indiv_test('6-[01]-[01]-no_rotate',idx(:,[1,4]), [7,10,13;2 4 18]');
    
%     [02]   [13]   [24]   [35]   [40]   [51] 
% [02] X      4   6  X      10 11  X      16 16
% [13] 1  1   X      7  9   X      13 14  X
% [24] X      5   4  X      11 12  X      17 17
% [35] 2  2   X      8  7   X      14 15  X
% [40] X      6   5  X      12 10  X      18 18
% [51] 3  3   X      9  8   X      15 13  X
tst.stimulation = mk_stim_patterns(6,1,[0,2],[0,2],{'rotate_meas'},1);
idx = reciprocity_idx( tst ); idx = reshape(idx,3,6);
do_indiv_test('6-[02]-[02]-rotate',idx(:,[1,4]), [6,11,16;15,2,7]');

tst.stimulation = mk_stim_patterns(6,1,[0,2],[0,2],{'no_rotate_meas'},1);
idx = reciprocity_idx( tst ); idx = reshape(idx,3,6); 
do_indiv_test('6-[02]-[02]-no_rotate',idx(:,[1,4]), [4,10,16;2,8,14]');

%     [02]   [13]   [24]   [35]   [40]   [51] 
% [02] 1  1   7  12  13 17  19 22  25 27  31 32
% [13] 2  2   8   7  14 18  20 23  26 28  32 33
% [24] 3  3   9   8  15 13  21 24  27 29  33 34
% [35] 4  4   10  9  16 14  22 19  28 30  34 35
% [40] 5  5   11 10  17 15  23 20  29 25  35 36
% [51] 6  6   12 11  18 16  24 21  30 26  36 31
tst.stimulation = mk_stim_patterns(6,1,[0,2],[0,2],{'meas_current','rotate_meas'},1);
idx = reciprocity_idx( tst ); idx = reshape(idx,6,6);
do_indiv_test('6-[02]-[02]-mc-rotate',idx(:,[1,4]),  ...
            [ 1,12,17,22,27,32;19,30,35,4,9,14]');

tst.stimulation = mk_stim_patterns(6,1,[0,2],[0,2],{'meas_current','no_rotate_meas'},1);
idx = reciprocity_idx( tst ); idx = reshape(idx,6,6); 
do_indiv_test('6-[02]-[02]-mc-no_rotate',idx(:,[1,4]),  ...
            [ 1,7,13,19,25,31; 4,10,16,22,28,34]');

tst.stimulation = mk_stim_patterns(6,1,[0,4],[0,4],{'meas_current','no_rotate_meas'},1);
idx = reciprocity_idx( tst ); idx = reshape(idx,6,6); 
do_indiv_test('6-[04]-[04]-mc-no_rotate',idx(:,[1,4]),  ...
            [ 1,7,13,19,25,31; 4,10,16,22,28,34]');


function do_indiv_test(txt,a,b,tol)
   if nargin < 4; tol = 0; end
   fprintf('%25s = ',txt);
   ok='fail';
   try; if isnan(a) == isnan(b); a(isnan(a))=0; b(isnan(b))=0; end; end
   try; if all(abs(a - b) <= tol);  ok='ok'; end; catch; end
   disp(ok)
   if ~strcmp(ok,'ok'); keyboard; end
