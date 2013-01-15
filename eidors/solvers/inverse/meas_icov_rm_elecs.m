function meas_icov = meas_icov_rm_elecs( imdl, elec_list)
% MEAS_ICOV_RM_ELECS: remove electrodes from consideration
% meas_icov = meas_icov_rm_elecs( inv_model, elec_list )
%
% PARAMETERS:
% - elec_list = numbers of elecs to remove
%     or as imdl.meas_icov_rm_elecs.elec_list
%
% - imdl.meas_icov_rm_elecs.exponent - exponent
% - imdl.meas_icov_rm_elecs.SNR      - SNR to add (default inf)
% - imdl.meas_icov_rm_elecs.replace_value (default 0)
%      Default is to modify the current meas_icov value, if replace_value==1,
%      then a new value is calculated without reference to the current
%
% meas_icov_rm_elecs can also accept a fwd_model parameter
% 
% Reference Accounting for erroneous electrode data in EIT
% A. Adler Physiological Measurement, 25(1):227-238, 2004. 

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isstr(imdl) && strcmp(imdl,'UNIT_TEST'); do_unit_test; return; end

switch imdl.type,
  case 'inv_model'; fmdl = imdl.fwd_model;
  case 'fwd_model'; fmdl = imdl;
                    imdl.meas_icov_rm_elecs.replace_value = 1;
  otherwise;        error('meas_icov_rm_elecs: require inv- or fwd-model');
end

if nargin< 2
   elec_list = imdl.meas_icov_rm_elecs.elec_list;
end

     NSR = 0;
try; NSR = (imdl.meas_icov_rm_elecs.SNR)^(-1);
end

    exponent = 1;
try;exponent = imdl.meas_icov_rm_elecs.exponent;
end

    replace_value = 0;
try;replace_value = imdl.meas_icov_rm_elecs.replace_value;
end
 

meas_icov = [];
for stim = fmdl.stimulation(:)'
   mp = stim.meas_pattern;
   sp = stim.stim_pattern;
   icovi = ones(size(mp,1),1);
   if any(sp(elec_list) ~= 0)
      icovi = 0*icovi;
   else
      icovi = ~any( mp(:,elec_list) ~= 0, 2);
   end

   meas_icov = [meas_icov; icovi];
end

n = length(meas_icov);

meas_icov(meas_icov == 0) = NSR;
meas_icov = spdiags( meas_icov.^exponent, 0, n,n );

if replace_value == 0
   meas_icov = calc_meas_icov(imdl)*meas_icov;
end

function do_unit_test
   imdl = mk_common_model('a2c0',4);
   covar = meas_icov_rm_elecs(imdl,1);
   unit_test_cmp('ire # 1', covar, zeros(size(covar))); 

   imdl = mk_common_model('a2c0',8);
   covar = meas_icov_rm_elecs(imdl,1);
   ff =    find( diag(covar)~= 1);
   ffcmp = [1;2;3;4;5;10;11;15;16;20;21;25;26;30;31;36;37;38;39;40];
   unit_test_cmp('ire # 2', ff,ffcmp);

   covar = meas_icov_rm_elecs(imdl.fwd_model,1);
   unit_test_cmp('ire # 3', ff,ffcmp);

   imdl.meas_icov_rm_elecs.elec_list = 1;
   covar = meas_icov_rm_elecs(imdl);
   unit_test_cmp('ire # 4', ff,ffcmp);
   ff =    find( diag(covar)~= 1);

   
   covar = meas_icov_rm_elecs(imdl,[1,2]);
   ff =    find( diag(covar)==1);
   ffcmp = [12;13;14;18;19;23;24;28;29;33;34;35];
   unit_test_cmp('ire # 5', ff,ffcmp);

   covar = meas_icov_rm_elecs(imdl,[1,2,3]);
   ff =    find( diag(covar)==1);
   ffcmp = [18;19;24;29;34;35];
   unit_test_cmp('ire # 6', ff,ffcmp);

   covar = meas_icov_rm_elecs(imdl,[1,2,3,4]);
   ff =    find( diag(covar)==1);
   ffcmp = [24;35];
   unit_test_cmp('ire # 7', ff,ffcmp);

% NOW CHECK OTHER MODES
   imdl.meas_icov_rm_elecs.elec_list = 1:4;
   imdl.meas_icov_rm_elecs.exponent = 1;
   imdl.meas_icov_rm_elecs.SNR      = 0;
   covar=  diag( meas_icov_rm_elecs(imdl) );
   ffcmp = [24;35];
   unit_test_cmp('ire # 8', find(covar==1),ffcmp);

   imdl.meas_icov_rm_elecs.SNR      = 100;
   covar=  diag( meas_icov_rm_elecs(imdl) );
   unit_test_cmp('ire # 9', find(covar==1),ffcmp);
   ff =    find( covar~=1);
   unit_test_cmp('ire #10', covar(ff),1/100);

   imdl.meas_icov_rm_elecs.exponent = -1;
   covar=  diag( meas_icov_rm_elecs(imdl) );
   unit_test_cmp('ire #11', find(covar==1),ffcmp);
   ff =    find( covar~=1);
   unit_test_cmp('ire #12', covar(ff),100);

   imdl = mk_common_model('a2c0',8);
   imdl.meas_icov = spdiag((1:40)');
   covar=  diag( meas_icov_rm_elecs(imdl,[]) );
   unit_test_cmp('ire #13', covar,(1:40)');

   imdl.meas_icov_rm_elecs.replace_value = 0;
   covar=  diag( meas_icov_rm_elecs(imdl,[]) );
   unit_test_cmp('ire #14', covar,(1:40)');

   imdl.meas_icov_rm_elecs.replace_value = 1;
   covar=  diag( meas_icov_rm_elecs(imdl,[]) );
   unit_test_cmp('ire #15', covar,1);
    
