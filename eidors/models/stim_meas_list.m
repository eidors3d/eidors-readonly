function [stim, meas_sel]= stim_meas_list( sp_mp , Nelec, current, gain);
%STIM_MEAS_LIST: mk stimulation pattern from list of electrodes
% [stim, meas_sel]= stim_meas_list( sp_mp , Nelec, amplitude, gain);
% [sp_mp]         = stim_meas_list( stim);
%
% stim =    EIDORS stimulation structure
%     use: fwd_model.stimulation = stim;
% meas_sel =EIDORS meas_select (select all values specified)
%                to form part of a fwd_model object
%
% sp_mp = matrix [stim+, stim-, meas+, meas-] x N_patterns
%
% Nelec  = Num of electrodes         DEFAULT = max sp_mp
% current= drive current levels,     DEFAULT = .010 Amp
% gain   = gain on voltage channels  DEFAULT = 1
%
% example: stim_meas_list([1,2,3,4; 1,2,4,5])
%    create stim pattern in to elecs 1,2 with differential
%    measurements on electrodes 3,4 and 4,5;
%    then convert the stim struct back to a list of stim/meas pairs
%  s=stim_meas_list([1,2,3,4; 1,2,4,5]);
%  stim_meas_list(s)
%  ans=[1,2,3,4; 1,2,4,5]

% (C) 2010,2015 Andy Adler, Alistair Boyle. License: GPL version 2 or version 3
% $Id$

if ischar(sp_mp) && strcmp(sp_mp,'UNIT_TEST'); do_unit_test; return; end

if isstruct(sp_mp)
  stim = sp_mp;
  meas_sel = [];
  nst = length(stim);
  nvt = 0;
  for i = 1:nst;
    nvt = nvt + size(stim(i).stim_pattern, 2) * size(stim(i).meas_pattern, 1);
  end
  stim = flatten_stim(stim, nst, nvt);
  stim = stim(:,[2 1 4 3]);
  return
end

if nargin <2; Nelec = max(sp_mp(:)); end
if nargin <3; current = 1;           end
if nargin <4; gain    = 1;           end

if any(sp_mp(:) > Nelec);
    error('Electrode patterns require more electrodes than Nelec');
end
stim = struct([]);
Npat = size(sp_mp,1);
is = 0;
current = current.*ones(Npat,1); % Extend if required
for i=1:Npat
   cur = sp_mp(i,1:2);
   new_stim = sparse( cur, 1, current(i)*[-1,1], Nelec,1);
   % create a new stim if it isn't the same as the last one
   if (is < 1) || any(any(stim(is).stim_pattern ~= new_stim))
     is = is + 1;
   end
   stim(is).stimulation = 'Amp';
   stim(is).stim_pattern = new_stim;
   mes = sp_mp(i,3:4);
   if isfield(stim(is),'meas_pattern') % append pattern if required
     stim(is).meas_pattern = [ stim(is).meas_pattern; sparse( 1, mes, gain *  [-1,1], 1, Nelec)];
   else
     stim(is).meas_pattern = sparse( 1, mes, gain *  [-1,1], 1, Nelec);
   end
end

% take in a stim/meas struct
% return a matrix of stim/meas pairs, per row with drive current 'i' and measurement gain 'g' calculated
% [ +s -s +m -m i g ]
function stim_flat = flatten_stim(stim, nst, nvt)
  stim_flat = zeros(nvt, 6);
  idx = 1;
  % TODO calculate 'gain' when it matched (m+ == - m-)
  % TODO calculate 'gain' when it is unmatched (m+ ~= m-)
  % TODO calculate 'current' when it matched (s+ == - s-)
  % TODO calculate 'current' when it is unmatched (s+ ~= s-)
  for i = 1:nst
      nmp= size(stim(i).meas_pattern, 1); % number of measurement patterns for this stim pair
      [sp, jnk, spv]= find(stim(i).stim_pattern>0);
      [sn, jnk, snv]= find(stim(i).stim_pattern<0);      
      [order, mp, mpv]= find(stim(i).meas_pattern>0);  
      [~, idx2sort] = sort(order); mp = mp(idx2sort); mpv = mpv(idx2sort);
      [order, mn, mnv]= find(stim(i).meas_pattern<0);    
      [~, idx2sort] = sort(order); mn = mn(idx2sort); mnv = mnv(idx2sort);
      % expand s+/s- to match the size of m+/m-
      sp  = zeros(nmp,1)+sp;
      sn  = zeros(nmp,1)+sn;
      spv = zeros(nmp,1)+spv;
      snv = zeros(nmp,1)+snv;
      stim_flat(idx:idx+nmp-1,:) = ...
        [ sp sn ... % stim pairs
            mp mn spv mpv];  % meas pairs
      idx = idx + nmp;
  end

function  do_unit_test
   imdl = mk_common_model('a2c0',16);
   img = mk_image(imdl);
   list_in = [1,2,3,4;1,2,4,5];
   img.fwd_model.stimulation = stim_meas_list(list_in,16);
   list_out = stim_meas_list(img.fwd_model.stimulation);
   unit_test_cmp('pattern#1', list_in, list_out);


   vh = fwd_solve(img);
   list_in = [6,7,3,4;1,2,4,5];
   img.fwd_model.stimulation = stim_meas_list(list_in,16);
   list_out = stim_meas_list(img.fwd_model.stimulation);
   unit_test_cmp('pattern#2', list_in, list_out);

   vh = fwd_solve(img);
   
   
   nElecs = 16;
   stim = mk_stim_patterns( nElecs, 1, '{ad}', '{ad}');
   sp_mp = stim_meas_list(stim);
   inj_diff = mod(sp_mp(:,2) - sp_mp(:,1), nElecs);
   meas_diff = mod(sp_mp(:,3) - sp_mp(:,4), nElecs);
   unit_test_cmp('pattern#3', true, all(inj_diff == 1) && all(meas_diff == 1));
   
   stim = mk_stim_patterns( nElecs, 1, '{ad}', '{ad}', {'no_meas_current', 'no_rotate_meas'});
   sp_mp = stim_meas_list(stim);
   inj_diff = mod(sp_mp(:,2) - sp_mp(:,1), nElecs);
   meas_diff = mod(sp_mp(:,3) - sp_mp(:,4), nElecs);
   unit_test_cmp('pattern#4', true, all(inj_diff == 1) && all(meas_diff == 1));
   
   stim = mk_stim_patterns( nElecs, 1, '{ad}', '{ad}', {'no_meas_current'});
   sp_mp = stim_meas_list(stim);
   inj_diff = mod(sp_mp(:,2) - sp_mp(:,1), nElecs);
   meas_diff = mod(sp_mp(:,3) - sp_mp(:,4), nElecs);
   unit_test_cmp('pattern#5', true, all(inj_diff == 1) && all(meas_diff == 1));
   
   Skip = 3;
   stim = mk_stim_patterns( nElecs, 1, [0 1+Skip], [0 1+Skip], {'no_meas_current', 'no_rotate_meas'});
   sp_mp = stim_meas_list(stim);
   inj_diff = mod(sp_mp(:,2) - sp_mp(:,1), nElecs);
   meas_diff = mod(sp_mp(:,3) - sp_mp(:,4), nElecs);
   unit_test_cmp('pattern#6', true, all(inj_diff == Skip+1) && all(meas_diff == Skip+1));
   
   nElecs = 32;
   Skip = 7;
   stim = mk_stim_patterns( nElecs, 1, [0 1+Skip], [0 1+Skip], {'meas_current', 'rotate_meas'});
   sp_mp = stim_meas_list(stim);
   inj_diff = mod(sp_mp(:,2) - sp_mp(:,1), nElecs);
   meas_diff = mod(sp_mp(:,3) - sp_mp(:,4), nElecs);


   unit_test_cmp('pattern#7', true, all(inj_diff == Skip+1) && all(meas_diff == Skip+1));
   
   
   

