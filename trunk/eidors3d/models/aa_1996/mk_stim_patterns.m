function stim= mk_stim_pattern( n_elec, n_rings, inj, meas, options)
%MK_STIM_PATTER: create an EIDORS3D stimulation pattern structure
%                to form part of a fwd_model object
% stim= mk_stim_pattern( n_elec, n_rings, inj, meas, options)
%
% where
% stim(#).stimulation = 'mA'
%     (#).stim_pattern= [vector n_elec*n_rings x 1 ]
%     (#).meas_pattern= [matrix n_elec*n_rings x n_meas_patterns]
%
% for example, for an adjacent pattern for 4 electrodes, with 0.5mA
%   if all electrodes are used for measurement
% stim(1).stim_pattern= [0.5;-0.5;0;0]
% stim(1).meas_pattern= [1, 0, 0,-1 
%                       -1, 1, 0, 0 
%                        0,-1, 1, 0 
%                        0, 0,-1, 1]
%
% PARAMETERS:
%   n_elec:   number of electrodes per ring
%   n_rings:  number of electrode rings (1 for 2D)
%
%   inj: injection pattern
%      '{ad}' or 'adjacent' -> adjacent drive: equivalent to [0 1]
%      '{op}' or 'opposite' -> opposite drive: equivalent to [0, n_elec/2]
%      [x y]: First pattern is [x,y] next is [x+1,y+1] 
%
%   meas: measurement pattern
%      '{ad}' or 'adjacent' -> adjacent measurement
%      No other measurement patterns supported
%
%   options: cell array of options, eg {'no_meas_current'}
%     if contradictory options are specified, only the last applies
%      'no_meas_current'  -> don't make measurements on current
%                            carrying electrodes
%      'meas_current'  -> do make measurements on current
%                            carrying electrodes

if isstr(inj)
   if      strcmp(inj,'{ad}') | strcmp(inj,'adjacent')
      inj= [0, 1];
   elseif  strcmp(inj,'{op}') | strcmp(inj,'opposite')
      inj= [0, floor(n_elec/2)];
   else
      error(['parameter inj=',inj,' not understood']);
   end

end

% iterate through the options cell array
for opt = options
   if     strcmp(opt, 'no_meas_current')
   elseif strcmp(opt, 'meas_current')
   else
      error(['option parameter opt=',opt,' not understood']);
   end
end

Ne = n_elec * n_rings;
for i=1:Ne
   stim(i).stimulation = 'mA';
end
