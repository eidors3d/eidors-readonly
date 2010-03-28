function [stim, meas_sel]= mk_stim_patterns( ...
            n_elec, n_rings, inj, meas, options, amplitude)
%MK_STIM_PATTERNS: create an EIDORS stimulation pattern structure
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
%      '{ad}'        -> adjacent drive: equivalent to [0 1]
%      '{op}'        -> opposite drive: equivalent to [0, n_elec/2]
%      '{trig}'      -> trigonometric drive [cos,sin,cos,sin ...]
%                       '{trig}' implies the 'meas_current' option.
%      '{trigcscs}'  -> trigonometric drive [cos,sin,cos,sin ...]
%      '{trigccss}'  -> trigonometric drive [cos,cos, ... sin,sin, ...]
%      '{mono}'      -> Drive via each elec, current leaves by ground
%      Bi-polar injection patterns:
%        [x y]: First pattern is [x,y] next is [x+1,y+1] 
%      Mono-polar injection patterns:
%        [x]:   First pattern is [x]   next is [x+1] 
%
%   meas: measurement pattern
%      '{ad}'        -> adjacent measurement
%      '{op}'        -> opposite drive: equivalent to [0, n_elec/2]
%      '{trig}'      -> trigonometric drive [sin,cos,sin,cos ...]
%      '{mono}'      -> Meas at each elec
%      Bi-polar measurement patterns:
%        [x y]: First pattern is [x,y] next is [x+1,y+1] 
%      Mono-polar measurement patterns:
%        [x]:   First pattern is [x]   next is [x+1] 
%
%   options: cell array of options, eg {'no_meas_current'}
%     if contradictory options are specified, only the last applies
%      'no_meas_current' / 'meas_current'
%         -> do / don't make mesurements on current carrying electrodes
%      'rotate_meas' / 'no_rotate_meas'
%         -> do / don't rotate measurements with stimulation pattern
%      'do_redundant' / 'no_redundant'
%         -> do / don't make reciprocally redundant measures
%      'balance_inj' / 'no_balance_inj'
%         -> do / don't draw current from all electrodes so total
%            injection is zero (useful for mono patterns)
%      'balance_meas' / 'no_balance_meas'
%         -> do / don't subtrant measurement from all electrodes so total
%            average measurement is zero (useful for mono patterns)
%
%   amplitude: drive current levels, DEFAULT = 1mA

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isstr(n_elec) && strcmp(n_elec,'UNIT_TEST'); do_unit_test; return; end

if nargin<6; amplitude= 1; end
if nargin<5; options= {};  end
v = process_args(n_elec, n_rings, inj, meas, options, amplitude );

curr_pat = v.inj(:) * ones(1,n_elec);
meas_pat = v.meas(:) * ones(1,n_elec);
offset   = [1;1]*(0:n_elec-1);

i=1;
for ring = 0:v.n_rings-1
   seen_patterns= struct;
   for elec= 0:v.n_elec-1
       if v.trig_inj && elec == v.n_elec-1 ; continue; end % No indep patterns
       s_pat= mk_stim_pat(v, elec, ring );
       m_pat= mk_meas_pat(v, elec, ring );

       if v.do_redundant == 0 % elim redudant
          [m_pat, seen_patterns] = elim_redundant(m_pat, s_pat, seen_patterns);
       end

       if ~isempty(m_pat) 
           stim(i).stimulation = 'mA';
           stim(i).stim_pattern= sparse(s_pat);
           stim(i).meas_pattern= sparse(m_pat);
           i=i+1;
       end

   end
end


meas_sel= meas_select( n_elec, v.inj, v);

% when not using data from current injection electrodes,
% it is common to be given a full measurement set.
% For example 16 electrodes gives 208 measures, but 256
% measure sets are common.
% 
% meas_select only makes sense for bipolar drive and measurement
% scenarios. Otherwise it is set to []
%
% This function calculates a selector matrix to remove the extra
%
% reshape(meas_select( 6, [0,1], v),6,6) 0 0 1 1 1 0
%                                        0 0 0 1 1 1
%                                        1 0 0 0 1 1
%                                        1 1 0 0 0 1
%                                        1 1 1 0 0 0
%                                        0 1 1 1 0 0
%                                        
% reshape(meas_select( 6, [0,2], v),6,6) 0 0 1 1 0 0
%                                        0 0 0 1 1 0
%                                        0 0 0 0 1 1
%                                        1 0 0 0 0 1
%                                        1 1 0 0 0 0
%                                        0 1 1 0 0 0
%                                        
% reshape(meas_select( 6, [0,3], v),6,6) 0 0 1 0 0 1
%                                        1 0 0 1 0 0
%                                        0 1 0 0 1 0
%                                        0 0 1 0 0 1
%                                        1 0 0 1 0 0
%                                        0 1 0 0 1 0

function meas_sel= meas_select( n_elec, inj, v)
  if prod(size(inj))~=2 || prod(size(v.meas))~=2
     meas_sel = [];
     return;
  end
 
  n2_elec= n_elec^2;
  e_idx= 0:n2_elec-1;
% ELS rotates the measurement around for each stim pattern
%  0 5 4 3 2 1
%  1 0 5 4 3 2
%  2 1 0 5 4 3
%  3 2 1 0 5 4
%  4 3 2 1 0 5
%  5 4 3 2 1 0
  ELS=rem(   rem(e_idx,n_elec)- ...
           floor(e_idx/n_elec)+n_elec , ...
           n_elec)';

  injx= n_elec+[-1 0 [-1 0]+inj*[-1;1] ];
% MEAS or STIM electrodes are equal?
%  5 5 5 5 5 5 5 5 5 5 5   0 1 2 3 4 5 5 0 1 2 3 4
%  0 0 0 0 0 0 0 0 0 0 0 = 0 1 2 3 4 5 5 0 1 2 3 4
%  2 2 2 2 2 2 2 2 2 2 2   0 1 2 3 4 5 5 0 1 2 3 4
%  3 3 3 3 3 3 3 3 3 3 3   0 1 2 3 4 5 5 0 1 2 3 4
  ELS= rem( injx ,n_elec)' * ones(1,n2_elec) ==ones(4,1)*ELS';
  inj_meas_sel= ~any( ELS )';

% reshape(rem(e_idx,n_elec)- floor(e_idx/n_elec),n_elec,n_elec)
%  0 -1 -2 -3 -4 -5
%  1  0 -1 -2 -3 -4
%  2  1  0 -1 -2 -3
%  3  2  1  0 -1 -2  -> Test if >0
%  4  3  2  1  0 -1
%  5  4  3  2  1  0
  if v.do_redundant==0
     ELS = rem(e_idx,n_elec)- floor(e_idx/n_elec);
     inj_meas_sel = inj_meas_sel & (ELS'>0);
  end
  
  % Insert electrode indices for multiple ring measurements
  meas_sel = [];
  inj_meas_sel = reshape(inj_meas_sel,n_elec,n_elec);
  for i = 1:v.n_rings
      meas_sel = blkdiag(meas_sel, ~inj_meas_sel); 
  end
  meas_sel = ~logical(meas_sel(:));
  

function stim_pat = mk_stim_pat(v, elec, ring );
   stim_pat = sparse(v.tn_elec, 1);
   if v.balance_inj
      stim_pat= stim_pat - sum(v.i_factor)/ (v.tn_elec-1);
   elseif v.trig_inj
      stim_pat = trig_pat( elec, v.tn_elec, v.trig_inj);
      return;
   end

   stim_idx = rem( v.inj + elec, v.n_elec) + 1 + v.n_elec*ring;
   stim_pat( stim_idx ) = v.amplitude*v.i_factor;

% Measurement config can stay static, or can rotate with
% the stim pattern. This code keeps measurements static
function meas = mk_meas_pat(v, elec, ring );
   meas= sparse(v.tn_elec, v.tn_elec);
   if v.balance_meas
      meas= meas - sum(v.m_factor)/ (v.tn_elec-1);
   elseif v.trig_meas
      meas= trig_pat( 1:v.tn_elec, v.tn_elec)';
      return;
   end

   if v.rotate_meas
      ofs = elec;
   else
      ofs = 0;
   end 

   mseq= 0:v.tn_elec-1;
   within_ring = rem(v.n_elec +  mseq , v.n_elec);
   ouside_ring = floor( mseq / v.n_elec) * v.n_elec;
   meas_seq    = mseq *v.tn_elec + 1;

   for i=1:length(v.meas)
      meas_pat = rem( v.meas(i) + within_ring + ofs, v.n_elec ) + ...
                       ouside_ring + meas_seq;
      meas(meas_pat) = v.m_factor(i);
   end

   if v.use_meas_current == 0
       stim_idx = rem( v.inj + elec, v.n_elec) + 1 + v.n_elec*ring;
   % each column of meas is a measurement pattern
   % Test whether each col has contribution from stim
       elim= any(meas(stim_idx,:));
       meas(:,elim) = [];
   end

   meas= meas';


function v = process_args(n_elec, n_rings, inj, meas, options, amplitude )

% SET DEFAULTS
   v.trig_meas= 0;
   v.trig_inj = 0;

   v.use_meas_current = 0;
   v.rotate_meas = 0;
   v.do_redundant = 1;
   v.balance_inj = 0;
   v.balance_meas= 1;

% Stimulation (injection) pattern
% This currently does not handle complicated injection patterns
if isstr(inj)
   if      strcmp(inj,'{ad}')
      inj= [0, 1];
      rel_ampl= [-1;1];
   elseif  strcmp(inj,'{op}')
      inj= [0, floor(n_elec/2)];
      rel_ampl= [-1;1];
   elseif  strcmp(inj,'{trig}')
      v.trig_inj = 1;
      v.use_meas_current = 1; % We need to measure on the electrodes
      rel_ampl= [];
   elseif  strcmp(inj,'{trigcscs}')
      v.trig_inj = 1;
      v.use_meas_current = 1; % We need to measure on the electrodes
      rel_ampl= [];
   elseif  strcmp(inj,'{trigccss}')
      v.trig_inj = 2;
      v.use_meas_current = 1; % We need to measure on the electrodes
      rel_ampl= [];
   elseif  strcmp(inj,'{mono}')
      inj= [0];
      rel_ampl= [1];
   else
      error(['parameter inj=',inj,' not understood']);
   end
elseif prod(size(inj))==1
      rel_ampl= [1];
elseif prod(size(inj))==2
      rel_ampl= [-1;1];
else
      error(['parameter inj not understood']);
end

v.inj= inj;
v.i_factor=      rel_ampl;

% Measurement configuration. 
% All systems I know of use adjacent measurement,
% are there actually any others?
if isstr(meas)
   if      strcmp(meas,'{ad}')
      meas=     [0, 1];
      rel_ampl= [ 1;-1];
   elseif  strcmp(meas,'{op}')
      meas= [0, floor(n_elec/2)];
      rel_ampl= [ 1;-1];
   elseif  strcmp(meas,'{trig}')
      v.trig_meas= 1;
      rel_ampl= [ 1;-1];
   elseif  strcmp(meas,'{mono}')
      meas= [0];
      rel_ampl= [-1];
   else
      error(['parameter meas=',meas,' not understood']);
   end
elseif prod(size(meas))==1
      rel_ampl= [1];
elseif prod(size(meas))==2
      rel_ampl= [ 1;-1];
else
      error(['parameter meas not understood']);
end

v.meas=          meas;
v.m_factor=      rel_ampl;

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
   elseif strcmp(opt, 'do_redundant')
      v.do_redundant = 1;
   elseif strcmp(opt, 'no_redundant')
      v.do_redundant = 0;
   elseif strcmp(opt, 'balance_inj')
      v.balance_inj = 1;
   elseif strcmp(opt, 'no_balance_inj')
      v.balance_inj = 0;
   elseif strcmp(opt, 'balance_meas')
      v.balance_meas= 1;
   elseif strcmp(opt, 'no_balance_meas')
      v.balance_meas= 0;
   else
      error(['option parameter opt=',opt,' not understood']);
   end
end

v.n_elec = n_elec;
v.n_rings= n_rings;
v.tn_elec= n_rings * n_elec;
v.amplitude = amplitude;

function [m_pat, seen_patterns] = elim_redundant(m_pat, s_pat, seen_patterns);
   m_pat_new= sparse([]);
   s_pat_str= ['s',sprintf('%d_', find(s_pat) ),'m'];
   for j=1:size(m_pat,1);
      this_m_pat= m_pat(j,:);
      pat_str= [s_pat_str, sprintf('%d_', find(this_m_pat))];
      if ~isfield(seen_patterns,pat_str);
         m_pat_new= [m_pat_new;this_m_pat];
         % we've seen this pattern
         seen_patterns.(pat_str)= 1;
         % and it's dual by reciprocity
         pat_str= ['s',sprintf('%d_', find(this_m_pat) ), ...
                   'm',sprintf('%d_', find(s_pat))];
         seen_patterns.(pat_str)= 1;
      end
   end
   m_pat= m_pat_new;

% create trig patterns.
%
% n_elecs is total number of electrodes
% elec    is the electrodes selected (can be multiple)
%         (elec goes from 0 to n_elecs-1
% if sel = 1 -> cos|sin|cos|sin (default)
% if sel = 2 -> cos|cos|sin|sin
% 
function pat= trig_pat( elec_sel, n_elecs, sel);
    if nargin<3; sel=1; end
    idx= linspace(0,2*pi,n_elecs+1)'; idx(end)= [];
    omega= idx*[1:n_elecs/2];
    meas_pat= [cos(omega), sin(omega) ];
    if sel==1;
       % reorder so we get cos|sin|cos|sin
       order = reshape(1:n_elecs,[],2)';
       meas_pat= meas_pat(:,order(:));
    end
    meas_pat= meas_pat(:,1:end-1); % only n_elecs-1 independent patterns
    pat  = meas_pat(:, elec_sel+1);

function do_unit_test
   stim = mk_stim_patterns(4,1,[0,1],[0,1],{},1);
   do_indiv_test('t1',stim(1).stim_pattern, [-1;1;0;0]);
   do_indiv_test('t2',stim(4).stim_pattern, [1;0;0;-1]);
   do_indiv_test('t3',stim(1).meas_pattern, [0,0,1,-1]);
   do_indiv_test('t4',stim(4).meas_pattern, [0,1,-1,0]);

%      'no_meas_current' / 'meas_current'
%         -> do / don't make mesurements on current carrying electrodes
   stim = mk_stim_patterns(4,1,[0,1],[0,1],{'meas_current'},1);
   do_indiv_test('meas_current: t1',stim(1).meas_pattern,  ...
         [1,-1,0,0; 0,1,-1,0; 0,0,1,-1; -1,0,0,1]);
   do_indiv_test('meas_current: t2',stim(4).stim_pattern, [1;0;0;-1]);
   stim = mk_stim_patterns(4,1,[0,1],[0,1],{'no_meas_current'},1);
   do_indiv_test('meas_current: t3',stim(1).meas_pattern,  [0,0,1,-1]);
   do_indiv_test('meas_current: t2',stim(4).stim_pattern, [1;0;0;-1]);

%      'rotate_meas' / 'no_rotate_meas'
%         -> do / don't rotate measurements with stimulation pattern

   stim = mk_stim_patterns(6,1,[0,1],[0,1],{'no_rotate_meas'},1);
   do_indiv_test('no_rotate_meas: t1',stim(2).stim_pattern, [0;-1;1;0;0;0]);
   do_indiv_test('no_rotate_meas: t2',stim(2).meas_pattern, ...
         [0,0,0,1,-1,0; 0,0,0,0,1,-1; -1,0,0,0,0,1]);
   do_indiv_test('no_rotate_meas: t3',stim(3).stim_pattern, [0;0;-1;1;0;0]);
   do_indiv_test('no_rotate_meas: t4',stim(3).meas_pattern, ...
         [1,-1,0,0,0,0;0,0,0,0,1,-1; -1,0,0,0,0,1]);

   stim = mk_stim_patterns(6,1,[0,1],[0,1],{'rotate_meas'},1);
   do_indiv_test('rotate_meas: t1',stim(2).stim_pattern, [0;-1;1;0;0;0]);
   do_indiv_test('rotate_meas: t2',stim(2).meas_pattern, ...
         [0,0,0,1,-1,0; 0,0,0,0,1,-1; -1,0,0,0,0,1]);
   do_indiv_test('rotate_meas: t3',stim(3).stim_pattern, [0;0;-1;1;0;0]);
   do_indiv_test('rotate_meas: t4',stim(3).meas_pattern, ...
         [0,0,0,0,1,-1; -1,0,0,0,0,1; 1,-1,0,0,0,0]);

%      'do_redundant' / 'no_redundant'
%         -> do / don't make reciprocally redundant measures
   stim = mk_stim_patterns(6,1,[0,1],[0,1],{'do_redundant'},1);

   do_indiv_test('do_redundant: t0',length(stim), 6);
   do_indiv_test('do_redundant: t1',stim(2).stim_pattern, [0;-1;1;0;0;0]);
   do_indiv_test('do_redundant: t2',stim(2).meas_pattern, ...
         [0,0,0,1,-1,0; 0,0,0,0,1,-1; -1,0,0,0,0,1]);
   do_indiv_test('do_redundant: t3',stim(3).stim_pattern, [0;0;-1;1;0;0]);
   do_indiv_test('do_redundant: t4',stim(3).meas_pattern, ...
         [1,-1,0,0,0,0;0,0,0,0,1,-1; -1,0,0,0,0,1]);

   stim = mk_stim_patterns(6,1,[0,1],[0,1],{'no_redundant'},1);
   do_indiv_test('no_redundant: t0',length(stim), 4);
   do_indiv_test('no_redundant: t1',stim(2).stim_pattern, [0;-1;1;0;0;0]);
   do_indiv_test('no_redundant: t2',stim(2).meas_pattern, ...
         [0,0,0,1,-1,0; 0,0,0,0,1,-1; -1,0,0,0,0,1]);
   do_indiv_test('no_redundant: t3',stim(3).stim_pattern, [0;0;-1;1;0;0]);
   do_indiv_test('no_redundant: t4',stim(3).meas_pattern, ...
         [0,0,0,0,1,-1;-1,0,0,0,0,1]);
   do_indiv_test('no_redundant: t5',stim(4).meas_pattern, ...
         [-1,0,0,0,0,1]);

%      'balance_inj' / 'no_balance_inj'
%         -> do / don't draw current from all electrodes so total
%            injection is zero (useful for mono patterns)
   stim = mk_stim_patterns(4,1,'{mono}',[0,1],{'balance_inj'},1);
   do_indiv_test('balance_inj: t0',length(stim), 4);
   do_indiv_test('balance_inj: t1',stim(2).stim_pattern, -[1;-3;1;1]/3);
   do_indiv_test('balance_inj: t2',stim(2).meas_pattern, ...
         [0,0,0,1,-1,0; 0,0,0,0,1,-1; -1,0,0,0,0,1]);
% WRONG - FIx

   stim = mk_stim_patterns(4,1,'{mono}',[0,1],{'no_balance_inj'},1);
   do_indiv_test('no_balance_inj: t0',length(stim), 4);
   do_indiv_test('no_balance_inj: t1',stim(2).stim_pattern, [0;1;0;0]);

   stim = mk_stim_patterns(4,1,'{mono}',[0,1],{},1);
   do_indiv_test('no_balance_inj: t0',length(stim), 4);
   do_indiv_test('no_balance_inj: t1',stim(2).stim_pattern, [0;1;0;0]);

function do_indiv_test(txt,a,b,tol)
   if nargin < 4; tol = 0; end
   fprintf('%20s = ',txt);
   ok='fail';
   try; if isnan(a) == isnan(b); a(isnan(a))=0; b(isnan(b))=0; end; end
   try; if all(abs(a - b) <= tol);  ok='ok'; end; end
   disp(ok)
   if ~strcmp(ok,'ok'); keyboard; end
