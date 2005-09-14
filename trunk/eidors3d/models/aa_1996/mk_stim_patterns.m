function [stim, meas_sel]= mk_stim_patterns( ...
            n_elec, n_rings, inj, meas, options, amplitude)
%MK_STIM_PATTER: create an EIDORS3D stimulation pattern structure
%                to form part of a fwd_model object
% [stim, meas_sel] = mk_stim_patterns( n_elec, n_rings, ...
%                                      inj, meas, options, amplitude)
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
% meas_sel: when not using data from current injection electrodes,
%           it is common to be given a full measurement set.
%           For example 16 electrodes gives 208 measures, but 256
%           measure sets are common. 'meas_sel' indicates which
%           electrodes are used
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
%      'rotate_meas'     -> rotate measurements with stimulation
%      'no_rotate_meas'  -> don't rotate measurements with stimulation
%
%   amplitude: drive current levels, DEFAULT = 1mA

if nargin<6; amplitude= 1; end
if nargin<5; options= {};  end
v = process_args(n_elec, n_rings, inj, meas, options, amplitude );

curr_pat = v.inj(:) * ones(1,n_elec);
meas_pat = v.meas(:) * ones(1,n_elec);
offset   = [1;1]*(0:n_elec-1);

i=1;
for ring = 0:v.n_rings-1
   for elec= 0:v.n_elec-1
       stim(i).stimulation = 'mA';
       stim(i).stim_pattern= mk_stim_pat(v, elec, ring );
       stim(i).meas_pattern= mk_meas_pat(v, elec, ring );

       i=i+1;
   end
end

meas_sel= meas_select( n_elec, v.inj);

% when not using data from current injection electrodes,
% it is common to be given a full measurement set.
% For example 16 electrodes gives 208 measures, but 256
% measure sets are common.
%
% This function calculates a selector matrix to remove the extra
function meas_sel= meas_select( n_elec, inj)
  n2_elec= n_elec^2;
  e_idx= 0:n2_elec-1;
  ELS=rem(   rem(e_idx,n_elec)- ...
           floor(e_idx/n_elec)+n_elec , ...
           n_elec)';
  injx= n_elec+[-1 0 [-1 0]+inj*[-1;1] ];
  ELS= rem( injx ,n_elec)' * ones(1,n2_elec) ==ones(4,1)*ELS';
  meas_sel= ~any( ELS )';



function stim_pat = mk_stim_pat(v, elec, ring, amplitude)
   stim_idx = rem( v.inj + elec, v.n_elec) + 1 + v.n_elec*ring;
   stim_pat = sparse(v.tn_elec, 1);
   stim_pat( stim_idx ) = v.amplitude*[-1;1];

% Measurement config can stay static, or can rotate with
% the stim pattern. This code keeps measurements static
function meas = mk_meas_pat(v, elec, ring, amplitude)
   meas= sparse(v.tn_elec, v.tn_elec);

   if v.rotate_meas
      ofs = elec;
   else
      ofs = 0;
   end 

   within_ring = rem(0:v.tn_elec-1, v.n_elec);
   ouside_ring = floor( (0:v.tn_elec-1)/ v.n_elec) * v.n_elec;
   meas_seq    = (0:v.tn_elec-1)*v.tn_elec + 1;

   meas_pat = rem( v.meas(1) + within_ring + ofs, v.n_elec ) + ...
                    ouside_ring + meas_seq;
   meas(meas_pat) =1;

   meas_pat = rem( v.meas(2) + within_ring + ofs, v.n_elec ) + ...
                    ouside_ring + meas_seq;
   meas(meas_pat) = -1;

   if v.use_meas_current == 0
       stim_idx = rem( v.inj + elec, v.n_elec) + 1 + v.n_elec*ring;
   % each column of meas is a measurement pattern
   % Test whether each col has contribution from stim
       elim= any(meas(stim_idx,:));
       meas(:,elim) = [];
   end


function v = process_args(n_elec, n_rings, inj, meas, options, amplitude )

% Stimulation (injection) pattern
% This currently does not handle complex injection patterns
if isstr(inj)
   if      strcmp(inj,'{ad}') | strcmp(inj,'adjacent')
      inj= [0, 1];
   elseif  strcmp(inj,'{op}') | strcmp(inj,'opposite')
      inj= [0, floor(n_elec/2)];
   else
      error(['parameter inj=',inj,' not understood']);
   end
elseif prod(size(inj))==2
% OK
else
      error(['parameter inj not understood']);
end

v.inj= inj;

% Measurement configuration. 
% All systems I know of use adjacent measurement,
% are there actually any others?
if isstr(meas)
   if      strcmp(meas,'{ad}') | strcmp(meas,'adjacent')
      meas= [0, 1];
   else
      error(['parameter meas=',meas,' not understood']);
   end
elseif prod(size(meas))==2
% OK
else
      error(['parameter meas not understood']);
end

v.meas= meas;

v.use_meas_current = 0;
v.rotate_meas = 0;

% iterate through the options cell array
for opt = options
   if     strcmp(opt, 'no_meas_current')
      v.use_meas_current = 0;
   elseif strcmp(opt, 'meas_current')
      v.use_meas_current = 1;
   elseif strcmp(opt, 'rotate_meas')
      v.rotate_meas = 1;
   elseif strcmp(opt, 'no_rotate_meas')
      v.rotate_meas = 0;
   else
      error(['option parameter opt=',opt,' not understood']);
   end
end

v.n_elec = n_elec;
v.n_rings= n_rings;
v.tn_elec= n_rings * n_elec;
v.amplitude = amplitude;
