function [stim, meas_sel]= stim_meas_list( sp_mp , Nelec, current, gain);
%STIM_MEAS_LIST: mk stimulation pattern from list of electrodes
% [stim, meas_sel]= stim_meas_list( sp_mp , Nelec, amplitude, gain);
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
%    measurements on electrodes 3,4 and 4,5

% (C) 2010 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isstr(sp_mp) && strcmp(sp_mp,'UNIT_TEST'); do_unit_test; return; end

if nargin <2; Nelec = max(sp_mp(:)); end
if nargin <3; current = 1;           end
if nargin <4; gain    = 1;           end

if any(sp_mp(:) > Nelec);
    error('Electrode patterns require more electrodes than Nelec');
end
stim = struct([]);
Npat = size(sp_mp,1);
is = 0;
for i=1:Npat
   cur = sp_mp(i,1:2);
   new_stim = sparse( cur, 1, current*[-1,1], Nelec,1);
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

function  do_unit_test
   imdl = mk_common_model('a2c0',16);
   img = mk_image(imdl);
   img.fwd_model.stimulation = stim_meas_list([1,2,3,4;1,2,4,5],16);
   vh = fwd_solve(img);
   img.fwd_model.stimulation = stim_meas_list([6,7,3,4;1,2,4,5],16);
   vh = fwd_solve(img);
