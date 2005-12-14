function ok= mk_stim_patterns_test
% Verify mk_stim_patterns function

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: mk_stim_patterns_test.m,v 1.2 2005-12-14 16:33:05 aadler Exp $

ok= 1;

pat= mk_stim_patterns(16,1,'{ad}','{ad}');
ok= ok& test_adj(pat);

options= {'no_rotate_meas'};
pat= mk_stim_patterns(16,1,'{ad}','{ad}', options);
ok= ok& test_adj(pat);

options= {'no_rotate_meas', 'no_meas_current'};
pat= mk_stim_patterns(16,1,'{ad}','{ad}', options);
ok= ok& test_adj(pat);

options= {'no_rotate_meas', 'meas_current'};
pat= mk_stim_patterns(16,1,'{ad}','{ad}', options);
ok= ok& test_adj_full(pat);

options= {'meas_current'};
pat= mk_stim_patterns(16,1,'{ad}','{ad}', options);
ok= ok& test_adj_full(pat);

options= {'rotate_meas'};
pat= mk_stim_patterns(16,1,'{ad}','{ad}', options);
ok= ok& test_adj_rotate(pat);

options= {'rotate_meas', 'no_meas_current'};
pat= mk_stim_patterns(16,1,'{ad}','{ad}', options);
ok= ok& test_adj_rotate(pat);

options= {'rotate_meas','no_redundant', 'no_meas_current'};
pat= mk_stim_patterns(16,1,'{ad}','{ad}', options);
ok= ok& test_adj_no_redund(pat);


function ok= test_adj(pat)
   eidors_msg('test adjacent current pattern',2);

   ok=1;
   if length(pat) ~= 16
      ok=0; warning('Fail at pt#01'); return; end

   if ~strcmp( pat(1).stimulation, 'mA');
      ok=0; warning('Fail at pt#02'); return; end

   % Stim pattern # 1
   if pat(1).stim_pattern ~= [-1;1;zeros(14,1)]
      ok=0; warning('Fail at pt#03'); return; end

   meas= pat(1).meas_pattern;

   if any( size(meas)~= [13 16] )
      ok=0; warning('Fail at pt#04'); return; end

   if any( meas(1,:)~= [0,0,1,-1,zeros(1,12)] )
      ok=0; warning('Fail at pt#05'); return; end

   if any( meas(13,:)~= [zeros(1,14),1,-1] )
      ok=0; warning('Fail at pt#06'); return; end

   % Stim pattern # 10
   if any( pat(10).stim_pattern ~= [zeros(9,1);-1;1;zeros(5,1)] )
      ok=0; warning('Fail at pt#07'); return; end

   meas= pat(10).meas_pattern;

   if any( size(meas)~= [13 16] )
      ok=0; warning('Fail at pt#08'); return; end

   if any( meas(1,:)~= [1,-1,zeros(1,14)] )
      ok=0; warning('Fail at pt#09'); return; end

   if any( meas(13,:)~= [-1,zeros(1,14),1] )
      ok=0; warning('Fail at pt#10'); return; end

function ok= test_adj_full(pat)
   eidors_msg('test adjacent current pattern (full)',2);

   ok=1;
   if length(pat) ~= 16
      ok=0; warning('Fail at pt#11'); return; end

   if ~strcmp( pat(1).stimulation, 'mA');
      ok=0; warning('Fail at pt#12'); return; end

   % Stim pattern # 1
   if pat(1).stim_pattern ~= [-1;1;zeros(14,1)]
      ok=0; warning('Fail at pt#13'); return; end

   meas= pat(1).meas_pattern;

   if any( size(meas)~= [16 16] )
      ok=0; warning('Fail at pt#14'); return; end

   if any( meas(1,:)~= [1,-1,zeros(1,14)] )
      ok=0; warning('Fail at pt#15'); return; end

   if any( meas(13,:)~= [zeros(1,12),1,-1,0,0] )
      ok=0; warning('Fail at pt#16'); return; end

   % Stim pattern # 10
   if any( pat(10).stim_pattern ~= [zeros(9,1);-1;1;zeros(5,1)] )
      ok=0; warning('Fail at pt#17'); return; end

   meas= pat(10).meas_pattern;

   if any( size(meas)~= [16 16] )
      ok=0; warning('Fail at pt#18'); return; end

   if any( meas(1,:)~= [1,-1,zeros(1,14)] )
      ok=0; warning('Fail at pt#19'); return; end

   if any( meas(13,:)~= [zeros(1,12),1,-1,0,0] )
      ok=0; warning('Fail at pt#20'); return; end


function ok= test_adj_rotate(pat)
   eidors_msg('test adjacent current pattern (rotate)',2);

   ok=1;
   if length(pat) ~= 16
      ok=0; warning('Fail at pt#21'); return; end

   if ~strcmp( pat(1).stimulation, 'mA');
      ok=0; warning('Fail at pt#22'); return; end

   % Stim pattern # 1
   if pat(1).stim_pattern ~= [-1;1;zeros(14,1)]
      ok=0; warning('Fail at pt#23'); return; end

   meas= pat(1).meas_pattern;

   if any( size(meas)~= [13 16] )
      ok=0; warning('Fail at pt#24'); return; end

   if any( meas(1,:)~= [0,0,1,-1,zeros(1,12)] )
      ok=0; warning('Fail at pt#25'); return; end

   if any( meas(13,:)~= [zeros(1,14),1,-1] )
      ok=0; warning('Fail at pt#26'); return; end

   % Stim pattern # 10
   if any( pat(10).stim_pattern ~= [zeros(9,1);-1;1;zeros(5,1)] )
      ok=0; warning('Fail at pt#27'); return; end

   meas= pat(10).meas_pattern;

   if any( size(meas)~= [13 16] )
      ok=0; warning('Fail at pt#28'); return; end

   if any( meas(1,:)~= [zeros(1,11),1,-1,zeros(1,3)] )
      ok=0; warning('Fail at pt#29'); return; end

   if any( meas(13,:)~= [zeros(1,7),1,-1,zeros(1,7)] )
      ok=0; warning('Fail at pt#30'); return; end

function ok= test_adj_no_redund(pat)
   eidors_msg('test adjacent current pattern (rotate)',2);

   ok=1;
   if length(pat) ~= 14
      ok=0; warning('Fail at pt#31'); return; end

   if ~strcmp( pat(1).stimulation, 'mA');
      ok=0; warning('Fail at pt#32'); return; end

   % Stim pattern # 1
   if pat(1).stim_pattern ~= [-1;1;zeros(14,1)]
      ok=0; warning('Fail at pt#33'); return; end

   meas= pat(1).meas_pattern;

   if any( size(meas)~= [13 16] )
      ok=0; warning('Fail at pt#34'); return; end

   if any( meas(1,:)~= [0,0,1,-1,zeros(1,12)] )
      ok=0; warning('Fail at pt#35'); return; end

   if any( meas(13,:)~= [zeros(1,14),1,-1] )
      ok=0; warning('Fail at pt#36'); return; end

   % Stim pattern # 10
   if any( pat(10).stim_pattern ~= [zeros(9,1);-1;1;zeros(5,1)] )
      ok=0; warning('Fail at pt#37'); return; end

   meas= pat(10).meas_pattern;

   if any( size(meas)~= [5 16] )
      ok=0; warning('Fail at pt#38'); return; end

   if any( meas(1,:)~= [zeros(1,11),1,-1,zeros(1,3)] )
      ok=0; warning('Fail at pt#39'); return; end

   if any( meas(5,:)~= [-1,zeros(1,14),1] )
      ok=0; warning('Fail at pt#40'); return; end

