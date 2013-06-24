function [stim, meas_sel]= mk_stim_patterns( ...
            n_elec, n_rings, inj, meas, options, amplitude)
%MK_STIM_PATTERNS: create an EIDORS stimulation pattern structure
%                to form part of a fwd_model object
% [stim, meas_sel] = mk_stim_patterns( n_elec, n_rings, ...
%                                      inj, meas, options, amplitude)
%
% where
% stim(#).stimulation = 'Amp'
%     (#).stim_pattern= [vector n_elec*n_rings x 1 ]
%     (#).meas_pattern= [matrix n_elec*n_rings x n_meas_patterns]
%
% for example, for an adjacent pattern for 4 electrodes, with 0.5 Amp
%   if all electrodes are used for measurement
% stim(1).stim_pattern= [0.5;-0.5;0;0]
% stim(1).meas_pattern= [1,-1, 0, 0 
%                        0, 1,-1, 0 
%                        0, 0, 1,-1 
%                       -1, 0, 0, 1]
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
%      '{mono}'      -> Meas at each elec (Reminder, you may want to set meas_current);
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
%   amplitude: drive current levels, DEFAULT = 0.010 Amp

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isstr(n_elec) && strcmp(n_elec,'UNIT_TEST'); do_unit_test; return; end

if nargin<6; amplitude= .01; end
if nargin<5; options= {};  end
v = process_args(n_elec, n_rings, inj, meas, options, amplitude );

[stim,mpat ] = calc_stim(v, n_elec, n_rings);

% Meas_sel selects meas which would exist if we measured all of them
v.do_redundant = 1; v.use_meas_current = 1; 
[jnk ,mpat0] = calc_stim(v, n_elec, n_rings);

%meas_sel= meas_select_old( n_elec, v.inj, v);
 meas_sel= meas_select( mpat, mpat0);

function [stim,mpat] = calc_stim(v, n_elec, n_rings)
   curr_pat = v.inj(:) * ones(1,n_elec);
   meas_pat = v.meas(:) * ones(1,n_elec);
   offset   = [1;1]*(0:n_elec-1);

   stim= struct([]);
   mpat= struct([]);
   i=1; j=1;
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
              stim(i).stimulation = 'Amp';
              stim(i).stim_pattern= sparse(s_pat);
              stim(i).meas_pattern= sparse(m_pat);
              i=i+1;
          end
              mpat(j).meas_pattern= sparse(m_pat);
              j=j+1;
      end
   end

% when not using data from current injection electrodes,
% it is common to be given a full measurement set.
% For example 16 electrodes gives 208 measures, but 256
% measure sets are common.
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
function meas_sel= meas_select( mpat, mpat0);
   meas_sel = [];
   for i=1:length(mpat);
      [mset,  err ] = mk_meas_set( mpat(i).meas_pattern );
      [mset0, err0] = mk_meas_set( mpat0(i).meas_pattern );
      if err || err0; msel_i = 1 + 0*mset; % Set to 1 for error
      else            msel_i = ismember(mset0,mset);
      end
      meas_sel = [meas_sel; msel_i];
   end

   meas_sel = logical(meas_sel(:));

function [mset,err] = mk_meas_set( meas_pat )
   mpats= size(~meas_pat,1);
   mset = zeros(mpats,1);
   err=0;
   for i=1:size(meas_pat,1);
      mpat = meas_pat(i,:);
      fp = find(mpat>0); if length(fp)==0; fp = 0; end
      fn = find(mpat<0); if length(fn)==0; fn = 0; end
      if length(fp)>1 || length(fn)>1 ; err=1; return; end
      mset(i) =  fp*1e7 + fn;
   end

function stim_pat = mk_stim_pat(v, elec, ring );
   stim_pat = sparse(v.tn_elec, 1);
   if v.balance_inj
      stim_pat= stim_pat - sum(v.i_factor)/ (v.tn_elec-1);
   elseif v.trig_inj
      stim_pat = trig_pat( elec, v.tn_elec, v.trig_inj) * v.amplitude;
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
       elim= any(meas(stim_idx,:),1);
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
   v.balance_meas= 0;

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
      rel_ampl= [1];
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

v.n_elec = n_elec;
v.n_rings= n_rings;
v.tn_elec= n_rings * n_elec;
v.amplitude = amplitude;

v= parse_options(v, options);

function v= parse_options(v, options);
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

function trig_tests;
   stim = mk_stim_patterns(8,1,'{trig}',[0,1],{},2);
   t= linspace(0,2*pi,8+1)'; t(end)= [];
   unit_test_cmp('trig: t1',[stim(1:4).stim_pattern], ...
         2*[cos(t),sin(t),cos(2*t),sin(2*t)], 1e-10);

   stim = mk_stim_patterns(8,1,'{trigcscs}',[0,1],{},2);
   unit_test_cmp('trig: t2',[stim(1:4).stim_pattern], ...
         2*[cos(t),sin(t),cos(2*t),sin(2*t)], 1e-10);

   stim = mk_stim_patterns(8,1,'{trigccss}',[0,1],{},2);
   unit_test_cmp('trig: t2',[stim([1,2,5,6]).stim_pattern], ...
         2*[cos(t),sin(t),cos(2*t),sin(2*t)], 1e-10);
%keyboard

function do_unit_test
   trig_tests;
   stim = mk_stim_patterns(4,1,[0,1],[0,1],{},1);
   unit_test_cmp('t1',stim(1).stim_pattern, [-1;1;0;0]);
   unit_test_cmp('t2',stim(4).stim_pattern, [1;0;0;-1]);
   unit_test_cmp('t3',stim(1).meas_pattern, [0,0,1,-1]);
   unit_test_cmp('t4',stim(4).meas_pattern, [0,1,-1,0]);

%      'no_meas_current' / 'meas_current'
%         -> do / don't make mesurements on current carrying electrodes
   stim = mk_stim_patterns(4,1,[0,1],[0,1],{'meas_current'},1);
   unit_test_cmp('meas_current: t1',stim(1).meas_pattern,  ...
         [1,-1,0,0; 0,1,-1,0; 0,0,1,-1; -1,0,0,1]);
   unit_test_cmp('meas_current: t2',stim(4).stim_pattern, [1;0;0;-1]);
   stim = mk_stim_patterns(4,1,[0,1],[0,1],{'no_meas_current'},1);
   unit_test_cmp('meas_current: t3',stim(1).meas_pattern,  [0,0,1,-1]);
   unit_test_cmp('meas_current: t2',stim(4).stim_pattern, [1;0;0;-1]);

%      'rotate_meas' / 'no_rotate_meas'
%         -> do / don't rotate measurements with stimulation pattern

   stim = mk_stim_patterns(6,1,[0,1],[0,1],{'no_rotate_meas'},1);
   unit_test_cmp('no_rotate_meas: t1',stim(2).stim_pattern, [0;-1;1;0;0;0]);
   unit_test_cmp('no_rotate_meas: t2',stim(2).meas_pattern, ...
         [0,0,0,1,-1,0; 0,0,0,0,1,-1; -1,0,0,0,0,1]);
   unit_test_cmp('no_rotate_meas: t3',stim(3).stim_pattern, [0;0;-1;1;0;0]);
   unit_test_cmp('no_rotate_meas: t4',stim(3).meas_pattern, ...
         [1,-1,0,0,0,0;0,0,0,0,1,-1; -1,0,0,0,0,1]);

   stim = mk_stim_patterns(6,1,[0,1],[0,1],{'rotate_meas'},1);
   unit_test_cmp('rotate_meas: t1',stim(2).stim_pattern, [0;-1;1;0;0;0]);
   unit_test_cmp('rotate_meas: t2',stim(2).meas_pattern, ...
         [0,0,0,1,-1,0; 0,0,0,0,1,-1; -1,0,0,0,0,1]);
   unit_test_cmp('rotate_meas: t3',stim(3).stim_pattern, [0;0;-1;1;0;0]);
   unit_test_cmp('rotate_meas: t4',stim(3).meas_pattern, ...
         [0,0,0,0,1,-1; -1,0,0,0,0,1; 1,-1,0,0,0,0]);

%      'do_redundant' / 'no_redundant'
%         -> do / don't make reciprocally redundant measures
   stim = mk_stim_patterns(6,1,[0,1],[0,1],{'do_redundant'},1);

   unit_test_cmp('do_redundant: t0',length(stim), 6);
   unit_test_cmp('do_redundant: t1',stim(2).stim_pattern, [0;-1;1;0;0;0]);
   unit_test_cmp('do_redundant: t2',stim(2).meas_pattern, ...
         [0,0,0,1,-1,0; 0,0,0,0,1,-1; -1,0,0,0,0,1]);
   unit_test_cmp('do_redundant: t3',stim(3).stim_pattern, [0;0;-1;1;0;0]);
   unit_test_cmp('do_redundant: t4',stim(3).meas_pattern, ...
         [1,-1,0,0,0,0;0,0,0,0,1,-1; -1,0,0,0,0,1]);

   stim = mk_stim_patterns(6,1,[0,1],[0,1],{'no_redundant'},1);
   unit_test_cmp('no_redundant: t0',length(stim), 4);
   unit_test_cmp('no_redundant: t1',stim(2).stim_pattern, [0;-1;1;0;0;0]);
   unit_test_cmp('no_redundant: t2',stim(2).meas_pattern, ...
         [0,0,0,1,-1,0; 0,0,0,0,1,-1; -1,0,0,0,0,1]);
   unit_test_cmp('no_redundant: t3',stim(3).stim_pattern, [0;0;-1;1;0;0]);
   unit_test_cmp('no_redundant: t4',stim(3).meas_pattern, ...
         [0,0,0,0,1,-1;-1,0,0,0,0,1]);
   unit_test_cmp('no_redundant: t5',stim(4).meas_pattern, ...
         [-1,0,0,0,0,1]);

%      'balance_inj' / 'no_balance_inj'
%         -> do / don't draw current from all electrodes so total
%            injection is zero (useful for mono patterns)
   stim = mk_stim_patterns(4,1,'{mono}',[0,1],{'balance_inj','meas_current'},1);
   unit_test_cmp('balance_inj: t0',length(stim), 4,1);
   unit_test_cmp('balance_inj: t1',stim(2).stim_pattern, -[1;-3;1;1]/3);
   unit_test_cmp('balance_inj: t2',stim(2).meas_pattern, ...
         [1,-1,0,0;0,1,-1,0;0,0,1,-1;-1,0,0,1]);

   stim = mk_stim_patterns(4,1,'{mono}',[0,1],{'no_balance_inj','no_meas_current'},1);
   unit_test_cmp('no_balance_inj: t0',length(stim), 4,1);
   unit_test_cmp('no_balance_inj: t1',stim(2).stim_pattern, [0;1;0;0]);
   unit_test_cmp('no_balance_inj: t2',stim(2).meas_pattern, ...
         [0,0,1,-1;-1,0,0,1]);

   stim = mk_stim_patterns(4,1,'{mono}',[0,1],{},1);
   unit_test_cmp('no_balance_inj: t0',length(stim), 4,1);
   unit_test_cmp('no_balance_inj: t1',stim(2).stim_pattern, [0;1;0;0]);

%      'balance_meas' / 'no_balance_meas'
%         -> do / don't subtrant measurement from all electrodes so total
%            average measurement is zero (useful for mono patterns)
   stim = mk_stim_patterns(4,1,[0,1],'{mono}',{'no_balance_meas','meas_current'},1);
   unit_test_cmp('no_balance_meas: t0',length(stim), 4,1);
   unit_test_cmp('no_balance_meas: t1',stim(2).stim_pattern, [0;-1;1;0]);
   unit_test_cmp('no_balance_meas: t1',stim(2).meas_pattern, eye(4));

   stim = mk_stim_patterns(4,1,[0,1],'{mono}',{'meas_current'},1);
   unit_test_cmp('no_balance_meas: t0',length(stim), 4,1);
   unit_test_cmp('no_balance_meas: t1',stim(2).stim_pattern, [0;-1;1;0]);
   unit_test_cmp('no_balance_meas: t1',stim(2).meas_pattern, eye(4));

   stim = mk_stim_patterns(4,1,[0,1],'{mono}',{'no_meas_current'},1);
   unit_test_cmp('no_balance_meas: t0',length(stim), 4,1);
   unit_test_cmp('no_balance_meas: t1',stim(2).stim_pattern, [0;-1;1;0]);
   unit_test_cmp('no_balance_meas: t1',stim(2).meas_pattern, [1,0,0,0;0,0,0,1]);

   stim = mk_stim_patterns(4,1,[0,1],'{mono}',{},1); % DO WE WANT THIS AS DEFAULT??
   unit_test_cmp('no_balance_meas: t0',length(stim), 4,1);
   unit_test_cmp('no_balance_meas: t1',stim(2).stim_pattern, [0;-1;1;0]);
   unit_test_cmp('no_balance_meas: t1',stim(2).meas_pattern, [1,0,0,0;0,0,0,1]);

   stim = mk_stim_patterns(4,1,[0,1],'{mono}',{'balance_meas','meas_current'},1);
   unit_test_cmp('balance_meas: t0',length(stim), 4);
   unit_test_cmp('balance_meas: t1',stim(2).stim_pattern, [0;-1;1;0]);
   unit_test_cmp('balance_meas: t1',stim(2).meas_pattern, (4*eye(4)-ones(4))/3);

   stim = mk_stim_patterns(4,1,[0,1],[0,1],{},2);
   unit_test_cmp('amplitude: t1',stim(2).stim_pattern, [0;-2;2;0]);
   stim = mk_stim_patterns(4,1,'{ad}',[0,1],{},2);
   unit_test_cmp('amplitude: t2',stim(2).stim_pattern, [0;-2;2;0]);
   stim = mk_stim_patterns(4,1,'{mono}',[0,1],{'no_balance_inj'},2);
   unit_test_cmp('amplitude: t3',stim(2).stim_pattern, [0;2;0;0]);
   stim = mk_stim_patterns(4,1,'{mono}',[0,1],{},2);
   unit_test_cmp('amplitude: t4',stim(2).stim_pattern, [0;2;0;0]);
  
   [stim,msel] = mk_stim_patterns(6,1,[0,2],[0,1],{'meas_current','no_rotate_meas'},2);
   msel = reshape(msel, 6, 6);
   unit_test_cmp('meas_sel: t1',msel(:,6), [1;1;1;1;1;1]);

   [stim,msel] = mk_stim_patterns(6,1,[0,2],[0,1],{'meas_current','rotate_meas'},2);
   msel = reshape(msel, 6, 6);
   unit_test_cmp('meas_sel: t1',msel(:,6), [1;1;1;1;1;1]);

   [stim,msel] = mk_stim_patterns(6,1,[0,1],[0,1],{'no_meas_current','no_rotate_meas'},2);
   msel = reshape(msel, 6, 6);
   unit_test_cmp('meas_sel: t1',msel(:,6), [0;1;1;1;0;0]);

   [stim,msel] = mk_stim_patterns(6,1,[0,2],[0,1],{'meas_current','no_rotate_meas'},2);
   msel = reshape(msel, 6, 6);
   unit_test_cmp('meas_sel: t1',msel(:,6), [1;1;1;1;1;1]);

   [stim,msel] = mk_stim_patterns(6,1,[0,1],[0,1],{'no_meas_current','no_rotate_meas'},2);
   msel = reshape(msel, 6, 6);
   unit_test_cmp('meas_sel: nnp01',msel(:,[4,5]), [1,1;1,1;0,1;0,0;0,0;1,0]);

   [stim,msel] = mk_stim_patterns(6,1,[0,2],[0,1],{'no_meas_current','no_rotate_meas'},2);
   msel = reshape(msel, 6, 6);
   unit_test_cmp('meas_sel: nnp02',msel(:,[4,5]), [1,0;1,1;0,1;0,0;0,0;0,0]);

   [stim,msel] = mk_stim_patterns(6,1,[0,3],[0,1],{'no_meas_current','no_rotate_meas'},2);
   msel = reshape(msel, 6, 6);
   unit_test_cmp('meas_sel: nnp03',msel(:,[4,5]), [0,0;1,0;0,1;0,0;1,0;0,1]);

   [stim,msel] = mk_stim_patterns(6,1,[1,2],[0,1],{'no_meas_current','no_rotate_meas'},2);
   msel = reshape(msel, 6, 6);
   unit_test_cmp('meas_sel: nnp12',msel(:,[4,5]), [1,0;1,1;1,1;0,1;0,0;0,0]);

   [stim,msel] = mk_stim_patterns(6,1,[2,4],[0,1],{'no_meas_current','no_rotate_meas'},2);
   msel = reshape(msel, 6, 6);
   unit_test_cmp('meas_sel: nnp24',msel(:,[4,5]), [0,0;0,0;1,0;1,1;0,1;0,0]);

   [stim,msel] = mk_stim_patterns(6,1,[2,4],[3,6],{'no_meas_current','no_rotate_meas'},2);
   msel = reshape(msel, 6, 6);
   unit_test_cmp('meas_sel: nnp2436',msel(:,[4,5]), [1,0;0,1;0,0;1,0;0,1;0,0]);

   [stim,msel] = mk_stim_patterns(6,1,[0,1],[0,1],{'no_meas_current','rotate_meas'},2);
   msel = reshape(msel, 6, 6);
   unit_test_cmp('meas_sel: nrp01',msel(:,6), [0;0;1;1;1;0]);

   [stim,msel] = mk_stim_patterns(6,1,[0,2],[0,1],{'no_meas_current','rotate_meas'},2);
   msel = reshape(msel, 6, 6);
   unit_test_cmp('meas_sel: nrp02',msel(:,6), [0;0;0;1;1;0]);

   [stim,msel] = mk_stim_patterns(6,1,[0,2],[3,4],{'no_meas_current','rotate_meas'},2);
   msel = reshape(msel, 6, 6);
   unit_test_cmp('meas_sel: nrp0234',msel(:,6), [1;1;0;0;0;0]);

   [stim,msel] = mk_stim_patterns(6,1,[0,1],[0,1],{'no_meas_current','no_rotate_meas','no_redundant'},2);
   msel = reshape(msel, 6, 6);
   unit_test_cmp('meas_sel: nnp01',msel(:,[2,5]), [0,0;0,0;0,0;1,0;1,0;1,0]);

   [stim,msel] = mk_stim_patterns(6,1,[0,2],'{mono}',{'no_meas_current','no_rotate_meas'},2);
   msel = reshape(msel, 6, 6);
   unit_test_cmp('meas_sel: nnp2436',msel(:,[4,5]), [1,0;1,1;1,1;0,1;1,0;0,1]);

% TESTS FROM OLD mk_stim_patterns_test CODE
   pat= mk_stim_patterns(16,1,'{ad}','{ad}',{}, 1);
   test_adj(pat);

   options= {'no_rotate_meas'};
   pat= mk_stim_patterns(16,1,'{ad}','{ad}', options,1);
   test_adj(pat);

   options= {'no_rotate_meas', 'no_meas_current'};
   pat= mk_stim_patterns(16,1,'{ad}','{ad}', options,1);
   test_adj(pat);

   options= {'no_rotate_meas', 'meas_current'};
   pat= mk_stim_patterns(16,1,'{ad}','{ad}', options,1);
   test_adj_full(pat);

   options= {'meas_current'};
   pat= mk_stim_patterns(16,1,'{ad}','{ad}', options,1);
   test_adj_full(pat);

   options= {'rotate_meas'};
   pat= mk_stim_patterns(16,1,'{ad}','{ad}', options,1);
   test_adj_rotate(pat);

   options= {'rotate_meas', 'no_meas_current'};
   pat= mk_stim_patterns(16,1,'{ad}','{ad}', options,1);
   test_adj_rotate(pat);

   options= {'rotate_meas','no_redundant', 'no_meas_current'};
   pat= mk_stim_patterns(16,1,'{ad}','{ad}', options,1);
   test_adj_no_redund(pat);

function ok= test_adj(pat)
   %%%  test adjacent current pattern

   unit_test_cmp('pt#01', length(pat), 16);
   unit_test_cmp('pt#02', pat(1).stimulation, 'Amp');
   unit_test_cmp('pt#03', pat(1).stim_pattern, [-1;1;zeros(14,1)]); % Stim pattern # 1

   meas= pat(1).meas_pattern;
   unit_test_cmp('pt#04', size(meas), [13 16] );
   unit_test_cmp('pt#05', meas(1,:), [0,0,1,-1,zeros(1,12)] );
   unit_test_cmp('pt#06', meas(13,:), [zeros(1,14),1,-1] );
   unit_test_cmp('pt#07', pat(10).stim_pattern , [zeros(9,1);-1;1;zeros(5,1)] );

   meas= pat(10).meas_pattern;
   unit_test_cmp('pt#08', size(meas), [13 16] );
   unit_test_cmp('pt#09', meas(1,:), [1,-1,zeros(1,14)] );
   unit_test_cmp('pt#10', meas(13,:), [-1,zeros(1,14),1] );

function ok= test_adj_full(pat)
   %%% test adjacent current pattern (full)

   unit_test_cmp('pt#11', length(pat), 16);
   unit_test_cmp('pt#12', pat(1).stimulation, 'Amp');
   unit_test_cmp('pt#13', pat(1).stim_pattern, [-1;1;zeros(14,1)]);

   meas= pat(1).meas_pattern;
   unit_test_cmp('pt#14', size(meas), [16 16] );
   unit_test_cmp('pt#15', meas(1,:), [1,-1,zeros(1,14)] );
   unit_test_cmp('pt#16', meas(13,:), [zeros(1,12),1,-1,0,0] );
   unit_test_cmp('pt#17', pat(10).stim_pattern, [zeros(9,1);-1;1;zeros(5,1)] );

   meas= pat(10).meas_pattern;
   unit_test_cmp('pt#18', size(meas), [16 16] );
   unit_test_cmp('pt#19', meas(1,:), [1,-1,zeros(1,14)] );
   unit_test_cmp('pt#20', meas(13,:), [zeros(1,12),1,-1,0,0] );


function ok= test_adj_rotate(pat)
   %%%% test adjacent current pattern (rotate)

   unit_test_cmp('pt#21', length(pat), 16);
   unit_test_cmp('pt#22', pat(1).stimulation, 'Amp');
   unit_test_cmp('pt#23', pat(1).stim_pattern, [-1;1;zeros(14,1)]);

   meas= pat(1).meas_pattern;
   unit_test_cmp('pt#24', size(meas), [13 16] );
   unit_test_cmp('pt#25', meas(1,:), [0,0,1,-1,zeros(1,12)] );
   unit_test_cmp('pt#26', meas(13,:), [zeros(1,14),1,-1] );
   unit_test_cmp('pt#27', pat(10).stim_pattern, [zeros(9,1);-1;1;zeros(5,1)] );

   meas= pat(10).meas_pattern;
   unit_test_cmp('pt#28', size(meas), [13 16] );
   unit_test_cmp('pt#29', meas(1,:), [zeros(1,11),1,-1,zeros(1,3)] );
   unit_test_cmp('pt#30', meas(13,:), [zeros(1,7),1,-1,zeros(1,7)] );

function ok= test_adj_no_redund(pat)
   %%% test adjacent current pattern (rotate)

   unit_test_cmp('pt#31', length(pat), 14);
   unit_test_cmp('pt#32', pat(1).stimulation, 'Amp');
   unit_test_cmp('pt#33', pat(1).stim_pattern, [-1;1;zeros(14,1)]);

   meas= pat(1).meas_pattern;
   unit_test_cmp('pt#34', size(meas), [13 16] );
   unit_test_cmp('pt#35', meas(1,:), [0,0,1,-1,zeros(1,12)] );
   unit_test_cmp('pt#36', meas(13,:), [zeros(1,14),1,-1] );
   unit_test_cmp('pt#37', pat(10).stim_pattern, [zeros(9,1);-1;1;zeros(5,1)] );

   meas= pat(10).meas_pattern;
   unit_test_cmp('pt#38', size(meas), [5 16] );
   unit_test_cmp('pt#39', meas(1,:), [zeros(1,11),1,-1,zeros(1,3)] );
   unit_test_cmp('pt#40', meas(5,:), [-1,zeros(1,14),1] );

