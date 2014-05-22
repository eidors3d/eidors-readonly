function params_stats
% Statistical tests as per Inez' paper chart
% $Id$

params_out; % Get data
fid= fopen(fnameout,'w');
fprintf(fid,'%% ALGS'); fn=alg_table; fprintf(fid,',%s',fn{:}); 
fprintf(fid,'\n'); fclose(fid);


zeep1_21 = 1; zeep2_21 = 5; zeep_100= 3;
peep1_21 = 2; peep2_21 = 6; peep_100= 4;
ceelv_z_p_21_1 = 12;
ceelv_z_p_21_2 = 13;
ceelv_z_p_100 = 14;
ceelv_z_z = 15; ceelv_p_p = 16;
ceelv1_z21_z100= 17;
ceelv2_z21_z100= 18;

param_Amp = 1;
param_CoV = 2;


S={'VT is reproducible             VT EIT ZEEP1 = VT EIT ZEEP2', 'VTZAZB'};
equal_test(S, data, zeep1_21, zeep2_21, param_Amp);
S={'VT is reproducible             VT EIT PEEP1 = VT EIT PEEP2', 'VTPAPB'};
equal_test(S, data, peep1_21, peep2_21, param_Amp);

S={'VT is independent of PEEP      VT EIT ZEEP = VT EIT PEEP  (21#1)','VTZPA'};
% equal_test(S, data, zeep1_21, peep1_21, param_Amp);
S={'VT is independent of PEEP      VT EIT ZEEP = VT EIT PEEP  (100)','VTZPO'};
equal_test(S, data, zeep_100, peep_100, param_Amp);
S={'VT is independent of PEEP      VT EIT ZEEP = VT EIT PEEP  (21#2)','VTZPB'};
% equal_test(S, data, zeep2_21, peep2_21, param_Amp);
S={'VT is independent of PEEP      VT EIT ZEEP = VT EIT PEEP  (21#1&2)','VTZPC'};
equal_test(S, data, [zeep1_21,zeep2_21], [peep1_21,peep2_21], param_Amp);


S={'VT is independent of FIO2      VT EIT ZEEP 21 = VT EIT ZEEP 100  (21#1)', 'VTZZA'};
% equal_test(S, data, zeep1_21, zeep_100, param_Amp);
S={'VT is independent of FIO2      VT EIT ZEEP 21 = VT EIT ZEEP 100  (21#2)', 'VTZZB'};
% equal_test(S, data, zeep2_21, zeep_100, param_Amp);
S={'VT is independent of FIO2      VT EIT ZEEP 21 = VT EIT ZEEP 100  (21#1&2)', 'VTZZC'};
equal_test(S, data, [zeep1_21,zeep2_21], [zeep_100,zeep_100], param_Amp);

S={'VT is independent of FIO2      VT EIT PEEP 21 = VT EIT PEEP 100  (21#1)', 'VTPPA'};
% equal_test(S, data, peep1_21, peep_100, param_Amp);
S={'VT is independent of FIO2      VT EIT PEEP 21 = VT EIT PEEP 100  (21#2)', 'VTPPB'};
% equal_test(S, data, peep2_21, peep_100, param_Amp);
S={'VT is independent of FIO2      VT EIT PEEP 21 = VT EIT PEEP 100  (21#1&2)', 'VTPPC'};
equal_test(S, data, [peep1_21,peep2_21], [peep_100,peep_100], param_Amp);

S={'CoG is reproducible             VT EIT ZEEP1 = VT EIT ZEEP2', 'CGZAZB'};
CoG_equal_test(S, data, zeep1_21, zeep2_21, param_CoV);
S={'CoG is reproducible             VT EIT PEEP1 = VT EIT PEEP2', 'CGPAPB'};
CoG_equal_test(S, data, peep1_21, peep2_21, param_CoV);
S={'CoG + with PEEP     CoG ZEEP > CoG PEEP (21#1)', 'CGZAPA'};
% agtb_test(S, data, zeep1_21, peep1_21, param_CoV);
S={'CoG + with PEEP     CoG ZEEP > CoG PEEP (21#2)', 'CGZBPB'};
% agtb_test(S, data, zeep2_21, peep2_21, param_CoV);
S={'CoG + with PEEP     CoG ZEEP > CoG PEEP (21#1&2)', 'CGZCPC'};
agtb_test(S, data, [zeep1_21,zeep2_21], [peep1_21,peep2_21], param_CoV);
S={'CoG + with PEEP     CoG ZEEP > CoG PEEP (100)','CGZOPO'};
agtb_test(S, data, zeep_100, peep_100, param_CoV);

S={'CoG + with FIO2     CoG ZEEP 100 > CoG ZEEP (21#1)', 'CGZOZA'};
% agtb_test(S, data, zeep_100, zeep1_21, param_CoV);
S={'CoG + with FIO2     CoG ZEEP 100 > CoG ZEEP (21#2)', 'CGZOZB'};
% agtb_test(S, data, zeep_100, zeep2_21, param_CoV);
S={'CoG + with FIO2     CoG ZEEP 100 > CoG ZEEP (21#1&2)', 'CGZOZC'};
agtb_test(S, data, [zeep_100,zeep_100], [zeep1_21,zeep2_21], param_CoV);

S={'CoG + with FIO2     CoG PEEP 100 > CoG PEEP (21#1)', 'CGPOPA'};
% agtb_test(S, data, peep_100, peep1_21, param_CoV);
S={'CoG + with FIO2     CoG PEEP 100 > CoG PEEP (21#2)', 'CGPOPB'};
% agtb_test(S, data, peep_100, peep2_21, param_CoV);
S={'CoG + with FIO2     CoG PEEP 100 > CoG PEEP (21#1&2)', 'CGPOPC'};
agtb_test(S, data, [peep_100,peep_100], [peep1_21,peep2_21], param_CoV);

% Can't do this because there is too much drift.

%S={'Global EELV is PEEP dependent               EELV EIT ZEEP < EELV EIT PEEP (#1)', 'EZAEPA'}
%agtb_test(S, data, ceelv_z_p_21_1, 0, param_Amp);

return %%%%%%%%%%%%%%%



function algs = alg_table
   algs = { 'R0', 'R3D', 'Rdif', 'Rbkg', 'Rmv', 'Rn', 'RL1', 'Rlpf','RTSVD', 'RTV',  'RGR', 'RSBP'};

function fn = fnameout; fn =  '../paper/stat_tests_out.tex';

function equal_test(str, data,  a, b, paramidx)
   stat_test(str, data, a, b, paramidx, @article_equal);

function CoG_equal_test(str, data,  a, b, paramidx)
   stat_test(str, data, a, b, paramidx, @CoG_equal);

% we accept a<b if we reject a>b
function altb_test(str, data, a, b, paramidx)
   stat_test(str, data, a, b, paramidx, @a_lt_b_t_test)

% we accept a>b if we reject a<b
function agtb_test(str, data, a, b, paramidx)
   stat_test(str, data, a, b, paramidx, @a_gt_b_t_test)

% test equality by (a-b)/mean(a,b) < 0.1
function aeqb_test(str, data, roiidx, depstr, a, b)
   stat_test(str, data, roiidx, depstr, a, b, @a_eq_b_t_test)


function stat_test(str, data, a, b, paramidx, test_ptr)
   TN = str{2}; str= str{1};
   PRINT_DEBUG_STATS = 0;

   fprintf('TEST:%s',str);

   % print to logfile
   fid= fopen(fnameout,'a');
   fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n%% "%s"\n',str);

   sel = 5:12; % only do these animals
   STATLINE = sprintf('\\def\\STATLINE%s{',TN);
   for alg= alg_table; % SPECIFY OUR ORDER
      data_ar = [data(sel,a).( alg{1} )];
      data_ar = data_ar(paramidx,:);
      if b==0;
         data_br = zeros(size(data_ar));
      else
         data_br = [data(sel,b).( alg{1} )];
         data_br = data_br(paramidx,:);
      end

      p= feval( test_ptr, data_ar, data_br );

      if p<.05; res= 'PASS (p<.05)';
      else    ; res= 'FAIL (p>.05)';
      end

      if PRINT_DEBUG_STATS
         fprintf(fid,'%%   ALG=%-6s:  p=%6.4f %s', alg{1}, p, res);
      end

      if strcmp(alg{1},'RGR') || strcmp(alg{1},'RGRb');
			       DD= 1e-1*[data_ar;data_br];
      else                     DD= 1e1*[data_ar;data_br]; end

      switch alg{1}
        case 'RGR';  fprintf('  GR=%3.6f',p);
        case 'R0';  fprintf('   R0=%3.2f',p);
        case 'RSBP'; fprintf('  BP=%3.2f\n',p);
      end


      D1 = DD(1,:); D2 = DD(2,:);
      Dd = ([1,-1]*DD);
      if PRINT_DEBUG_STATS
         fprintf(fid,' A-B=(');
         fprintf(fid,'%4d-%4d;',[D1;D2]);
         fprintf(fid,') (M=%4d;S=%4d)\n', mean(Dd), std(Dd));
      end
%       join   = ''; if ~strcmp(alg{1},'RGRb'), join='&';end
      join   = ''; if ~strcmp(alg{1},'RSBP'), join='&';end
      signif = '$^{~}$'; if p<.05, signif='$^*$';end
      STATLINE = [STATLINE, sprintf(' %1.3f%s %s',p,signif,join)];
      
   end
   STATLINE = [STATLINE, sprintf('} %% "%s" \n',str)];
   fprintf(fid,'%s',STATLINE);
   fclose(fid);

%to test for equality ($a=b$), we use the two-tailed
%$t$-test to test rejection of the null hypothesis,
% H$_0$:$|a-b|>0.1\times\frac{1}{2}(a+b)$.
%  stat_test(str, data, roiidx, depstr, a, b, @article_equal)
function p = article_equal( d1, d2)
   dd = d1 - d2;
   da = 0.1*(abs(d1+d2)/2);

% This is wrong, we need to compare the abs difference againt the chosen threshold
%  p = paired_t_test( dd, da );

   p = a_gt_b_t_test( da, abs(dd) );

function p = CoG_equal( d1, d2)
   dd = d1 - d2;
   da = 0.03 *ones(size(dd));  % FIXME

   p = a_gt_b_t_test( da, abs(dd) );

% Output is the p( d1 = d2 can be explained by the data )
function p = paired_t_test( d1, d2 );
   dd = d1 - d2;   
   m_diff= mean(dd,2);
   s_diff= std( dd,[],2);
   N = size(dd,2);
   DF= N - 1;
   T = m_diff / ( s_diff / sqrt(N) );
   p= two_sided(T,DF);

% Output is the p that d1>=d2 cannot be explained by the data
function p = a_gt_b_t_test( d1, d2 );
   dd = d1 - d2;   
   m_diff= mean(dd,2);
   s_diff= std( dd,[],2);
   N = size(dd,2);
   DF= N - 1;
   T = m_diff / ( s_diff / sqrt(N) );
   p= one_sided(T,DF);

function p = a_lt_b_t_test( d1, d2 );
  p = a_gt_b_t_test( d2, d1 );


% test equality by (a-b)/mean(a,b) < 0.1
function p = a_eq_b_t_test( d1, d2 );
   LIM = 0.10;

   dd = (d1 - d2) ./ (0.5*d1+0.5*d2);
   N = size(dd,2);
   DF= N - 1;

   m_diff= mean(dd,2);
   s_err = std( dd,[],2) / sqrt(N);
   

   % p is prob dd < LIM
   T1 = ( m_diff - LIM) / s_err; p1= one_sided(T1, DF);
   T2 = (-m_diff - LIM) / s_err; p2= one_sided(T2, DF);
   % p is prob data can be explained by > LIM
   p= 1 - min([p1,p2]);

function p= two_sided(T,DF)
   p= 1.0*betainc (1 / (1 + T^2/DF), DF/2, 0.5);

function p= one_sided(T,DF)
   p= 0.5*betainc (1 / (1 + T^2/DF), DF/2, 0.5);
   if (T<0); p = 1 - p; end

function old_tests

S={'Regional VT is PEEP dependent         nondep VT EIT ZEEP > nondep VT EIT PEEP (#1)', 'nVZnPA'};
agtb_test(S, data, roi, 'nondep', zeep1_21, peep1_21)
S={'Regional VT is PEEP dependent         nondep VT EIT ZEEP > nondep VT EIT PEEP (#2)', 'nVZnPO'};
agtb_test(S, data, roi, 'nondep', zeep_100, peep_100)
S={'Regional VT is PEEP dependent         nondep VT EIT ZEEP > nondep VT EIT PEEP (#3)', 'nVZnPB'};
agtb_test(S, data, roi, 'nondep', zeep2_21, peep2_21)

S={'Regional VT is PEEP dependent            dep VT EIT ZEEP <    dep VT EIT PEEP (#1)', 'dVZdPA'};
altb_test(S, data, roi, 'dep', zeep1_21, peep1_21)
S={'Regional VT is PEEP dependent            dep VT EIT ZEEP <    dep VT EIT PEEP (#2)', 'dVZdPO'};
altb_test(S, data, roi, 'dep', zeep_100, peep_100)
S={'Regional VT is PEEP dependent            dep VT EIT ZEEP <    dep VT EIT PEEP (#3)', 'dVZdPB'};
altb_test(S, data, roi, 'dep', zeep2_21, peep2_21)

S={'Regional VT is FIO2 dependent         nondep VT EIT ZEEP 21 < nondep VT EIT ZEEP 100 (#1)', 'nVZAnVZO'};
altb_test(S, data, roi, 'nondep', zeep1_21, zeep_100)
S={'Regional VT is FIO2 dependent         nondep VT EIT ZEEP 21 < nondep VT EIT ZEEP 100 (#2)', 'nVZBnVZO'};
altb_test(S, data, roi, 'nondep', zeep2_21, zeep_100)

S={'Regional VT is FIO2 dependent         nondep VT EIT PEEP 21 < nondep VT EIT PEEP 100 (#1)', 'nVPAnVPO'};
altb_test(S, data, roi, 'nondep', peep1_21, peep_100)
S={'Regional VT is FIO2 dependent         nondep VT EIT PEEP 21 < nondep VT EIT PEEP 100 (#2)', 'nVPAnVPO'};
altb_test(S, data, roi, 'nondep', peep2_21, peep_100)

S={'Regional VT is FIO2 dependent            dep VT EIT ZEEP 21 > dep VT EIT ZEEP 100 (#1)', 'dVZAdVZO'};
agtb_test(S, data, roi, 'dep', zeep1_21, zeep_100)
S={'Regional VT is FIO2 dependent            dep VT EIT ZEEP 21 > dep VT EIT ZEEP 100 (#2)', 'dVZBdVZO'};
agtb_test(S, data, roi, 'dep', zeep2_21, zeep_100)

S={'Regional VT is FIO2 dependent            dep VT EIT PEEP 21 > dep VT EIT PEEP 100 (#1)', 'dVPAdVPO'};
agtb_test(S, data, roi, 'dep', peep1_21, peep_100)
S={'Regional VT is FIO2 dependent            dep VT EIT PEEP 21 > dep VT EIT PEEP 100 (#2)', 'dVPBdVPO'};
agtb_test(S, data, roi, 'dep', peep2_21, peep_100)


S={'Regional VT is reproducible           nondep VT EIT ZEEP1 = nondep VT EIT ZEEP2 (#1)', 'nVZAnDZB'};
equal_test(S, data, roi, 'nondep', zeep1_21, zeep2_21)

S={'Regional VT is reproducible           nondep VT EIT PEEP1 = nondep VT EIT PEEP2 (#1)', 'nVPAdVPB'};
equal_test(S, data, roi, 'nondep', peep1_21, peep2_21)

S={'Regional VT is reproducible              dep VT EIT ZEEP1 = dep VT EIT ZEEP2 (#1)', 'dVZAdVZB'};
equal_test(S, data, roi, 'dep', zeep1_21, zeep2_21)
S={'Regional VT is reproducible              dep VT EIT PEEP1 = dep VT EIT PEEP2 (#1)', 'dVPAdVPB'};
equal_test(S, data, roi, 'dep', peep1_21, peep2_21)


S={'Global EELV is PEEP dependent               EELV EIT ZEEP < EELV EIT PEEP (#1)', 'EZAEPA'}
altb_test(S, data, roi, 'global', ceelv_z_p_21_1, 0)
S={'Global EELV is PEEP dependent               EELV EIT ZEEP < EELV EIT PEEP (100)', 'EZOEPO'}
altb_test(S, data, roi, 'global', ceelv_z_p_100, 0)
S={'Global EELV is PEEP dependent               EELV EIT ZEEP < EELV EIT PEEP (#2)', 'EZBEPB'}
altb_test(S, data, roi, 'global', ceelv_z_p_21_2, 0)

S={'Global EELV is reproducible                EELV EIT ZEEP1 = EELV EIT ZEEP2 (#1)', 'EZAEZB'};
equal_test(S, data, roi, 'global', ceelv_z_z, 0)
S={'Global EELV is reproducible                EELV EIT PEEP1 = EELV EIT PEEP2 (#1)', 'EPAEPB'};
equal_test(S, data, roi, 'global', ceelv_p_p, 0)

S={'Regional EELV is FIO2 dependent       nondep EELV EIT ZEEP 21 < nondep EELV EIT ZEEP 100 (#1)', 'nEZAnEZO'};
altb_test(S, data, roi, 'nondep', ceelv1_z21_z100, 0);
S={'Regional EELV is FIO2 dependent       nondep EELV EIT ZEEP 21 < nondep EELV EIT ZEEP 100 (#2)', 'nEZBnEZO'};
altb_test(S, data, roi, 'nondep', ceelv2_z21_z100, 0);
S={'Regional EELV is FIO2 dependent          dep EELV EIT ZEEP 21 > dep EELV EIT ZEEP 100 (#1)', 'dEZAdEZO'};
agtb_test(S, data, roi, 'dep', ceelv2_z21_z100, 0);
S={'Regional EELV is FIO2 dependent          dep EELV EIT ZEEP 21 > dep EELV EIT ZEEP 100 (#2)', 'dEZBdEZO'};
agtb_test(S, data, roi, 'dep', ceelv2_z21_z100, 0);


S={'Regional EELV is reproducible         nondep EELV EIT ZEEP1 = nondep EELV EIT ZEEP2', 'nEZAnEZB'};
equal_test(S, data, roi, 'nondep', zeep1_21, zeep2_21)
S={'Regional EELV is reproducible            dep EELV EIT ZEEP1 =    dep EELV EIT ZEEP2', 'dEZAdEZB'};
equal_test(S, data, roi, 'dep', zeep1_21, zeep2_21)
S={'Regional EELV is reproducible         nondep EELV EIT PEEP1 = nondep EELV EIT PEEP2', 'nEPAnEPB'};
equal_test(S, data, roi, 'nondep', peep1_21, peep2_21)
S={'Regional EELV is reproducible            dep EELV EIT PEEP1 =    dep EELV EIT PEEP2', 'dEPAdEPB'};
equal_test(S, data, roi, 'dep', peep1_21, peep2_21)
