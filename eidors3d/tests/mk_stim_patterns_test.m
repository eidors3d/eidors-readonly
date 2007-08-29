function ok= mk_stim_patterns_test
% Verify mk_stim_patterns function

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: mk_stim_patterns_test.m,v 1.14 2007-08-29 09:24:12 aadler Exp $

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
      ok=0; eidors_msg('Fail at pt#01',1); return; end

   if ~strcmp( pat(1).stimulation, 'mA');
      ok=0; eidors_msg('Fail at pt#02',1); return; end

   % Stim pattern # 1
   if pat(1).stim_pattern ~= [-1;1;zeros(14,1)]
      ok=0; eidors_msg('Fail at pt#03',1); return; end

   meas= pat(1).meas_pattern;

   if any( size(meas)~= [13 16] )
      ok=0; eidors_msg('Fail at pt#04',1); return; end

   if any( meas(1,:)~= [0,0,-1,1,zeros(1,12)] )
      ok=0; eidors_msg('Fail at pt#05',1); return; end

   if any( meas(13,:)~= [zeros(1,14),-1,1] )
      ok=0; eidors_msg('Fail at pt#06',1); return; end

   % Stim pattern # 10
   if any( pat(10).stim_pattern ~= [zeros(9,1);-1;1;zeros(5,1)] )
      ok=0; eidors_msg('Fail at pt#07',1); return; end

   meas= pat(10).meas_pattern;

   if any( size(meas)~= [13 16] )
      ok=0; eidors_msg('Fail at pt#08',1); return; end

   if any( meas(1,:)~= [-1,1,zeros(1,14)] )
      ok=0; eidors_msg('Fail at pt#09',1); return; end

   if any( meas(13,:)~= [1,zeros(1,14),-1] )
      ok=0; eidors_msg('Fail at pt#10',1); return; end

function ok= test_adj_full(pat)
   eidors_msg('test adjacent current pattern (full)',2);

   ok=1;
   if length(pat) ~= 16
      ok=0; eidors_msg('Fail at pt#11',1); return; end

   if ~strcmp( pat(1).stimulation, 'mA');
      ok=0; eidors_msg('Fail at pt#12',1); return; end

   % Stim pattern # 1
   if pat(1).stim_pattern ~= [-1;1;zeros(14,1)]
      ok=0; eidors_msg('Fail at pt#13',1); return; end

   meas= pat(1).meas_pattern;

   if any( size(meas)~= [16 16] )
      ok=0; eidors_msg('Fail at pt#14',1); return; end

   if any( meas(1,:)~= [-1,1,zeros(1,14)] )
      ok=0; eidors_msg('Fail at pt#15',1); return; end

   if any( meas(13,:)~= [zeros(1,12),-1,1,0,0] )
      ok=0; eidors_msg('Fail at pt#16',1); return; end

   % Stim pattern # 10
   if any( pat(10).stim_pattern ~= [zeros(9,1);-1;1;zeros(5,1)] )
      ok=0; eidors_msg('Fail at pt#17',1); return; end

   meas= pat(10).meas_pattern;

   if any( size(meas)~= [16 16] )
      ok=0; eidors_msg('Fail at pt#18',1); return; end

   if any( meas(1,:)~= [-1,1,zeros(1,14)] )
      ok=0; eidors_msg('Fail at pt#19',1); return; end

   if any( meas(13,:)~= [zeros(1,12),-1,1,0,0] )
      ok=0; eidors_msg('Fail at pt#20',1); return; end


function ok= test_adj_rotate(pat)
   eidors_msg('test adjacent current pattern (rotate)',2);

   ok=1;
   if length(pat) ~= 16
      ok=0; eidors_msg('Fail at pt#21',1); return; end

   if ~strcmp( pat(1).stimulation, 'mA');
      ok=0; eidors_msg('Fail at pt#22',1); return; end

   % Stim pattern # 1
   if pat(1).stim_pattern ~= [-1;1;zeros(14,1)]
      ok=0; eidors_msg('Fail at pt#23',1); return; end

   meas= pat(1).meas_pattern;

   if any( size(meas)~= [13 16] )
      ok=0; eidors_msg('Fail at pt#24',1); return; end

   if any( meas(1,:)~= [0,0,-1,1,zeros(1,12)] )
      ok=0; eidors_msg('Fail at pt#25',1); return; end

   if any( meas(13,:)~= [zeros(1,14),-1,1] )
      ok=0; eidors_msg('Fail at pt#26',1); return; end

   % Stim pattern # 10
   if any( pat(10).stim_pattern ~= [zeros(9,1);-1;1;zeros(5,1)] )
      ok=0; eidors_msg('Fail at pt#27',1); return; end

   meas= pat(10).meas_pattern;

   if any( size(meas)~= [13 16] )
      ok=0; eidors_msg('Fail at pt#28',1); return; end

   if any( meas(1,:)~= [zeros(1,11),-1,1,zeros(1,3)] )
      ok=0; eidors_msg('Fail at pt#29',1); return; end

   if any( meas(13,:)~= [zeros(1,7),-1,1,zeros(1,7)] )
      ok=0; eidors_msg('Fail at pt#30',1); return; end

function ok= test_adj_no_redund(pat)
   eidors_msg('test adjacent current pattern (rotate)',2);

   ok=1;
   if length(pat) ~= 14
      ok=0; eidors_msg('Fail at pt#31',1); return; end

   if ~strcmp( pat(1).stimulation, 'mA');
      ok=0; eidors_msg('Fail at pt#32',1); return; end

   % Stim pattern # 1
   if pat(1).stim_pattern ~= [-1;1;zeros(14,1)]
      ok=0; eidors_msg('Fail at pt#33',1); return; end

   meas= pat(1).meas_pattern;

   if any( size(meas)~= [13 16] )
      ok=0; eidors_msg('Fail at pt#34',1); return; end

   if any( meas(1,:)~= [0,0,-1,1,zeros(1,12)] )
      ok=0; eidors_msg('Fail at pt#35',1); return; end

   if any( meas(13,:)~= [zeros(1,14),-1,1] )
      ok=0; eidors_msg('Fail at pt#36',1); return; end

   % Stim pattern # 10
   if any( pat(10).stim_pattern ~= [zeros(9,1);-1;1;zeros(5,1)] )
      ok=0; eidors_msg('Fail at pt#37',1); return; end

   meas= pat(10).meas_pattern;

   if any( size(meas)~= [5 16] )
      ok=0; eidors_msg('Fail at pt#38',1); return; end

   if any( meas(1,:)~= [zeros(1,11),-1,1,zeros(1,3)] )
      ok=0; eidors_msg('Fail at pt#39',1); return; end

   if any( meas(5,:)~= [1,zeros(1,14),-1] )
      ok=0; eidors_msg('Fail at pt#40',1); return; end

