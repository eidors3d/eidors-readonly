function ok= demo_real_test
% Perform tests based on the demo_real function

% (C) 2005 Andy Adler + Nick Polydorides. License: GPL version 2 or version 3
% $Id: demo_real_test.m,v 1.22 2007-08-29 09:26:40 aadler Exp $

isOctave= exist('OCTAVE_VERSION');

datareal= 'datareal.mat';
datacom=  'datacom.mat';
drt=      'demo_real_test.mat';
if isOctave
    datareal= file_in_loadpath(datareal);
    datacom=  file_in_loadpath(datacom);
    drt    =  file_in_loadpath(drt);
    page_screen_output= 0;
end

load(datareal);
perm_sym= '{y}';

[I,Ib] = set_3d_currents(protocol,elec,vtx,gnd_ind,no_pl);

mat_ref = 1*ones(828,1);

%Set the tolerance for the forward solver
tol = 1e-5;

[Eref,D,Ela,ppr] = fem_master_full(vtx,simp,mat_ref,gnd_ind,elec,zc,perm_sym);
[Vref] = forward_solver(Eref,I,tol,ppr);
[refH,refV,indH,indV,dfr]=get_3d_meas(elec,vtx,Vref,Ib,no_pl);
dfr = dfr(1:2:end); %Taking just the horrizontal measurements

mat=mat_ref;
load(datacom,'A','B'); %Indices of the elements to represent the inhomogeneity
sA = mat_ref(A(1))+0.15;
sB = mat_ref(B(1))-0.20;
mat(A) = sA;
mat(B) = sB;

[En,D,Ela,ppn] = fem_master_full(vtx,simp,mat,gnd_ind,elec,zc,perm_sym);
[Vn] = forward_solver(En,I,tol,ppn,Vref);
[voltageH,voltageV,indH,indV,dfv]=get_3d_meas(elec,vtx,Vn,Ib,no_pl);
dfv = dfv(1:2:end);

if size(dfr)~= size(dfv)
   error('Mismatched measurements')
end

[v_f] = m_3d_fields(vtx,32,indH,Eref,tol,gnd_ind);

[J] = jacobian_3d(I,elec,vtx,simp,gnd_ind,mat_ref,zc,v_f,dfr,tol,perm_sym);

[Reg] = iso_f_smooth(simp,vtx,3,1);

tfac = 1e-8;

sol = (J'*J +  tfac*Reg'*Reg)\J' * (voltageH - refH);

Diag_Reg_012= [diag(Reg,0),[diag(Reg,1);0],[diag(Reg,2);0;0]];
Jcolsby100=J(:,1:100:size(J,2));

%% verifications

load(drt);

% Need to divide by 2 since code bugs are fixed
compare_tol( drt.voltageH, voltageH, 'voltageH' )
compare_tol( drt.voltageV, voltageV, 'voltageV' )
compare_tol( drt.sol, sol, 'sol' )
compare_tol( drt.Jcolsby100, Jcolsby100, 'Jcolsby100' )
compare_tol( drt.Diag_Reg_012, Diag_Reg_012, 'Diag_Reg_012' )

ok=1;


function compare_tol( cmp1, cmp2, errtext )
% compare matrices and give error if not equal
fprintf(2,'testing parameter: %s ...\n',errtext);

tol= 1e-4;

vd= mean(mean( abs(cmp1 - cmp2) ));
vs= mean(mean( abs(cmp1) + abs(cmp2) ));
if vd/vs > tol
   eidors_msg( ...
     'parameter %s exceeds tolerance %g (=%g)', errtext, tol, vd/vs,1 );
end

