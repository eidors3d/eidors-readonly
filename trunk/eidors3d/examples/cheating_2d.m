% code to simulate inverse crimes in EIT
% $Id: cheating_2d.m,v 1.2 2005-10-18 15:25:23 aadler Exp $

function out=cheating_2d

   [vis,vhs,s_mdl]= small_2d_mdl;
    i_mdl= large_inv_model;

   im= inv_solve( i_mdl, vis, vhs );
   show_fem(im);
   reconst_with_noise(i_mdl, vis, vhs, 4)


return;

% simulate 'sad' data for small model
function [vis,vhs,mdl]= small_2d_mdl
   params= mk_circ_tank(8, [], 16 ); 
   params.stimulation= mk_stim_patterns(16, 1, '{ad}','{ad}', ...
                         {'no_meas_current','no_rotate_meas'}, 1);
   params.solve=      'aa_fwd_solve';
   params.system_mat= 'aa_calc_system_mat';
   mdl= eidors_obj('fwd_model', params);

   mat= ones( size(mdl.elems,1), 1);

   % homogeneous data
   vhs= fwd_solve( eidors_obj('image','name',  ...
                     'elem_data', mat, 'fwd_model', mdl ));

   % inhomogeneous data for sad face model
   leye= [78,97,98,117,118,141];
   reye= [66,82,83,102,103,123];
   mouth= [28,31,40,43,58,57,53,54,74,73,92,93,113, ...
            112,135,69,87,70,107,88,108,129];
   mat(leye)= 2;
   mat(reye)= 2;
   mat(mouth)=1.5;

   vis= fwd_solve( eidors_obj('image','name',  ...
                     'elem_data', mat, 'fwd_model', mdl ));

function i_mdl= large_inv_model;
   params= mk_circ_tank(12, [], 16 ); 
   params.stimulation= mk_stim_patterns(16, 1, '{ad}','{ad}', ...
                         {'no_meas_current','no_rotate_meas'}, 1);
   params.solve=      'aa_fwd_solve';
   params.system_mat= 'aa_calc_system_mat';
   params.jacobian  = 'aa_calc_jacobian';
   l_mdl= eidors_obj('fwd_model', params);

% create inverse model
  %hparam.value = 1e-8;
   hparam.func = 'aa_calc_noise_figure';
   hparam.noise_figure= 4;
   hparam.tgt_elems= 1:4;

   img_prior.func = 'tikhonov_image_prior';
  %img_prior.func = 'aa_calc_image_prior';

   i_mdl= eidors_obj( ...
          'inv_model', 'large 2D inverse', ...
          'hyperparameter', hparam, 'image_prior', img_prior, ...
          'reconst_type', 'difference', ...
          'fwd_model', l_mdl, 'solve', 'aa_inv_solve' );


function reconst_with_noise(i_mdl, vis, vhs, num_tries)
    noise = .0002*randn(size(vis.meas,1),num_tries);

    vi_n(1:num_tries) = vis;
    for i= 1:num_tries % stupid matlab doesn't allow easy vectorization
       vi_n(i).meas = vi_n(i).meas + noise(:,i);
    end
    show_slices( inv_solve( i_mdl, vi_n, vhs ));

function boo;

global vhs vis;

if ~exist('vhs') || isempty(vhs)
    resetup('d0')

    leye= [78,97,98,117,118,141];
    reye= [66,82,83,102,103,123];
    mouth= [28,31,40,43,58,57,53,54,74,73,92,93,113, ...
            112,135,69,87,70,107,88,108,129];
    sad= zeros(1,256);
    vhs= prob_dir( sad );
    sad(leye)=1;
    sad(reye)=1;
    sad(mouth)=.5;
    vis= prob_dir( sad );
end

function cleancolourmap
  s= 1-[flipud(hot(64)) ;fliplr(hot(64))];
  s= [flipud(fliplr(hot(64))); hot(64)];
  s= [s([1:2:63 64:2:128],:)*.7+.2;[1 1 1]];
  colormap(s);
endfunction



resetup('e0'); cleancolourmap;



function noise_test(vhs,vis)
    % Noise Test
    nn= 400;
    nvec= .0002*randn(208,nn);
    vn = vis*ones(1,nn) + nvec;
    image(imgr(irec(vn,vhs,0,4),50),1);
endfunction

noise1= 1e-5*[-27.76546; -7.59427; -1.93711; -27.05896; -12.85904; -11.46233;
    1.71458; -15.84886; 9.49042; 26.54607; 26.12679; 10.28903; 11.27190;
    5.94225; 6.44672; -19.99604; -22.04821; -8.67624; 11.81035; 4.80488;
   29.03764; -22.83378; 32.48301; 30.09973; 33.70337; -0.47085; 3.85125;
    9.05904; 12.30617; 14.14135; -10.76137; 11.01337; 17.39058; 37.05449;
   28.25643; 18.98459; 26.58454; 11.32970; -35.71133; -1.36525; -5.43171;
    9.98939; -3.04222; 21.29392; -15.11127; -8.27054; 13.01116; 5.32276;
   19.24881; -16.71019; -15.73683; -4.06076; 2.37260; 6.83274; 5.26589;
    8.83536; -6.29922; 0.63870; -20.85604; -5.73974; 11.16503; -9.94519;
    8.64528; 21.63097; -4.20898; 8.17845; -17.65653; 10.60018; 8.51757;
   15.55584; 30.22712; 0.71698; 18.94148; -15.84283; 26.82607; -12.08514;
  -21.82176; -10.15707; -15.43742; -1.46121; -31.93084; -35.00449; 43.65366;
  -28.86313; -8.89273; -4.88046; 12.11069; 13.11200; 24.02962; 25.25127;
   -6.39429; -18.28042; -22.45900; 11.12896; -16.43752; 11.79022; 26.50192;
    6.81792; 14.81372; -15.59964; -18.31752; 2.16775; 25.24500; 20.77933;
   -2.09474; 13.79688; -31.40512; 9.63977; 13.49541; -35.61789; -9.74602;
   15.56013; 12.24740; -2.10360; 11.13246; -5.85497; 0.79389; 36.91269;
   35.25533; 43.20230; -1.32461; -15.47513; 29.92197; 6.57967; -0.61718;
   12.68731; 39.90042; 10.47227; -38.75338; -35.50139; -35.01525; -14.01510;
   -6.43549; 11.36402; 11.05323; 18.26093; 6.74087; -25.88768; -36.80022;
   15.95831; 6.77176; 20.80731; 21.31454; 27.03949; -36.53713; -0.63354;
   37.45211; 21.17920; 1.71484; 10.29992; -37.99016; 6.85604; -16.11494;
   -7.55474; 20.94234; 12.73071; 16.33950; -2.41776; -16.23419; 25.60348;
  -13.22026; -0.41706; -6.05165; -13.65678; -40.08519; -9.26705; 4.14438;
    9.65111; 21.48499; 9.51472; -8.41532; -6.05089; -2.66941; -13.33649;
   -6.67147; 22.97640; -0.21244; 1.87677; 37.91273; -21.12891; 11.22574;
   13.86916; 11.60222; -14.87688; -16.05350; -19.14326; -11.15221; 3.82456;
   10.36910; 11.75312; -30.31466; -28.21045; 17.44571; -38.11927; 1.55531;
   11.29632; -13.01711; 23.26577; 28.83256; -44.26250; -2.35528; -14.15714;
  -28.49528; 34.63137; -19.37650; -14.00479; 11.41968; -8.54238]; 
noise3=1e-5*[ 5.89382; -11.76715; -20.38383; 25.08106; -12.05456;
   42.44387; 30.43318; 6.56803; 16.58853; 1.97596; 24.88580;
    0.83537; -1.81821; -21.13378; 17.27234; 7.04205; -35.14670;
   15.98648; 19.50310; -7.40796; -11.66726; -23.45204; 44.69655;
    0.17090; 31.43564; -27.55132; 2.21016; -21.08567; 16.25432;
  -11.63719; 12.41088; 14.11388; 7.72355; -5.39017; -0.12547;
  -38.58446; -20.10008; -15.15861; -15.86851; -5.84053; -30.07290;
  -37.36250; -2.53053; 2.54702; -43.21144; 35.84945; 20.05358;
  -11.16362; -2.38681; 20.40710; 18.91450; 30.85527; 23.66431;
   -9.11852; 3.31994; -10.97963; 29.01126; -40.09789; -2.74479;
    9.48409; -17.94402; 36.94314; 31.88554; -25.40466; -13.15508;
  -47.24235; -17.08644; -46.53377; -10.74565; 11.32868; 18.86069;
   22.13141; 5.42699; 1.41991; -36.29532; -20.18932; -30.89348;
   20.54996; 16.83219; -33.43325; 20.56161; -19.14693; -6.48806;
    5.46298; -9.15376; -3.86233; 9.63450; 11.44795; 37.29246;
   43.19218; 22.32386; -28.51360; 9.94637; -6.16323; -48.93215;
  -12.29428; 17.79531; -8.09356; 17.83567; 28.58079; 13.56182;
   -7.04320; -28.85897; -13.93723; 15.28342; -13.16094; 40.85824;
    9.99891; -35.85101; 8.10981; 0.66691; 34.24551; 6.44166;
   15.34364; -2.77896; 42.85408; 6.59410; 11.76109; 11.95449;
   15.50520; -4.78434; 15.30083; -46.30672; 39.56400; -3.20849;
   -9.92821; -9.78347; 34.56455; 7.89952; 66.45235; -32.17386;
    1.68791; 39.34762; -17.82598; 32.91579; 15.01798; -10.67046;
    6.50035; -19.61355; 5.31837; 34.16613; 8.41203; -9.61931;
   -7.23913; 27.83542; -9.88770; 33.82686; -23.26859; 5.41904;
    9.06312; 16.01742; 8.79307; -6.71841; 8.02480; 1.07434;
    0.70523; -6.28717; 4.78425; -57.24199; -29.81538; -19.67294;
   27.44577; 18.46645; -1.85555; -16.30453; -24.93155; 10.69559;
  -29.67271; 6.18243; 14.99266; 22.73511; -26.77038; -13.64186;
   -3.83900; 1.00810; 63.84263; -7.44160; -19.39190; 2.35079;
    8.04457; 24.87581; 19.30803; -0.26582; 6.46934; -15.08571;
    2.17275; 15.56750; 18.75319; 2.68343; 18.73000; 6.30140;
  -26.80381; -34.22573; 14.17559; -8.28483; 5.89168; -2.67605;
   20.63114; 2.23133; -45.23000; 0.21635; 4.08253; -29.12417;
  -11.76195; -24.50549; 26.35054; 4.04918; -13.95879]; 

%
% STEP 1: Create image of 'random' noise which happens
%    to result in a happytransform
%
ii= imgr(irec([vis+noise1,vis+noise3],vhs,0,4));
image(ii)
imwrite('happynoise.png',ii,colormap);


% Create Tikhonov reconstruction matrix
%    p_prior        ( 0 -> I         = d^0
%                     1 -> sqrt(S'S) = d^.5
%                     2 -> S'S       = d 
%    p_noise        ( 0 -> Hi = 1 
%                     1 -> Hi = 1/homg
%                     2 -> Hi = 1/homg^2 
%
%    Prior Weighting  p_ww (default = ones)
function Z= mkreconst( lamda, p_prior, p_noise, p_ww );
    global DVV;
    [m,e]= size(DVV);

    if ~exist('p_ww'); p_ww = ones(1,e); end

    n_var= 1./prob_dir( zeros(1,e) );
    W= sparse(1:m,1:m, n_var.^(-p_noise) );
    pond=sparse(1:e,1:e, ...
          ((n_var.^(-p_noise) )'*(DVV.^2)).^(p_prior/2) .* p_ww );
    dd= DVV'*W*DVV;
    dx= DVV'*W;
    Z= (dd+ lamda*pond)\dx;
endfunction

function nf= calc_nf( Z )
    global AIRE;
    global DVV; [m,e]= size(DVV);

    nfg= sqrt(e/m);         % noise figure gain
    n_var= 1./prob_dir( zeros(1,e) );
    sig= mean(DVV(:,1:4)')';
    nf= sum(sig)*sqrt(sum(sum( ((AIRE*n_var').*Z).^2 ))) / ...
        sum(AIRE'*Z*sig)/sqrt(n_var'*n_var)*nfg;
   
endfunction

function ii=bad_reg_icrime
    global vis vhs
    resetup('d0'); cleancolourmap;

    ll= .01;
    zg1=  mkreconst(ll,1,1);
    calc_nf(zg1)
    rg1= irec(vis,vhs,zg1);

    leye= [78,97,98,117,118,141];
    reye= [66,82,83,102,103,123];
    rsmile= [40:41, 53:55, 69:71, 86:88];
    lsmile= [43:44, 57:59, 73:75, 91:93];
    lsad = [31,43,58,57,74,73,92,93,113,112,135];
    sad = [28,31,40,43,58,57,53,54,74,73,92,93,113, ...
            112,135,69,87,70,107,88,108,129];

    p_ww= ones(1,256);
    p_ww([leye,reye,sad])= .5;
    zg3=  mkreconst(ll,1,1,p_ww);
    calc_nf(zg3)
    rg3= irec(vis,vhs,zg3);

    p_ww= ones(1,256);
    p_ww([leye,reye,rsmile,lsad])= .5;
    zb1=  mkreconst(ll,1,1,p_ww);
    calc_nf(zb1)
    rb1= irec(vis,vhs,zb1);

    p_ww= ones(1,256);
    p_ww([leye,reye,rsmile,lsmile])= .5;
    zb3=  mkreconst(ll,1,1,p_ww);
    calc_nf(zb3)
    rb3= irec(vis,vhs,zb3);

    ii=imgr([rg1,  rg3, rb1,  rb3]);
endfunction


%
% STEP 2: Tikhonov regularization happytransform WITH i_crime
%

ii=bad_reg_icrime;
imwrite('tikhonov-icrime.png',ii,colormap);
image(ii);

%
% STEP 3: Tikhonov regularization happytransform WITHOUT i_crime
%
sad_e0= [ 53; 57; 69; 73; 86; 87; 91; 92; 106; 107; 111; 112; 127; 128;
         129; 133; 134; 135; 147; 151; 152; 153; 157; 158; 159; 165;
         171; 172; 177; 178; 179; 184; 185; 186; 192; 193; 199; 200;
         205; 206; 207; 212; 213; 214; 220; 221; 227; 228; 229; 235;
         236; 243; 244; 251; 252; 253; 259; 260; 261; 267; 268; 275;
         276; 283; 284; 285; 292; 293; 301; 310; 319; 320]; 

happy_e0=[ 69; 70; 71; 73; 74; 75; 86; 87; 88; 89; 91; 92; 93; 94; 106;
       107; 108; 109; 111; 112; 113; 114; 127; 128; 129; 130; 131; 133;
       134; 135; 136; 137; 147; 151; 152; 153; 154; 155; 157; 158; 159;
       160; 161; 165; 171; 172; 176; 177; 178; 179; 180; 183; 184; 185;
       186; 187; 192; 193; 199; 200; 220; 221; 227; 228; 229; 251; 252;
       253; 259; 260; 261; 283; 284; 285; 292; 293; 319; 320]; 

halfy_e0=[ 57; 69; 70; 71; 73; 86; 87; 88; 89; 91; 92; 106; 107; 108; 109;
          111; 112; 127; 128; 129; 130; 131; 133; 134; 135; 147; 151; 152;
          153; 154; 155; 157; 158; 159; 165; 171; 172; 176; 177; 178; 179;
          180; 184; 185; 186; 192; 193; 199; 200; 212; 213; 214; 220; 221;
          227; 228; 229; 243; 244; 251; 252; 253; 259; 260; 261; 275; 276;
          283; 284; 285; 292; 293; 310; 319; 320];


resetup('e0'); cleancolourmap;
rr=zeros(1,576); rr(happy_e0) = 1; imwrite('happy.png',imgr(rr),colormap);
rr=zeros(1,576); rr(sad_e0) = 1; imwrite('sad.png',imgr(rr),colormap);

ll= .0025;
z001=  mkreconst(ll,1,1);
calc_nf(z001)
rg1= irec(vis,vhs,z001);

p_ww= ones(1,576); p_ww(sad_e0)= .65;
b001=  mkreconst(ll,1,1,p_ww);
calc_nf(b001)
rg3= irec(vis,vhs,b001);

p_ww= ones(1,576); p_ww(halfy_e0)= .65;
b001=  mkreconst(ll,1,1,p_ww);
calc_nf(b001)
rb1= irec(vis,vhs,b001);


p_ww= ones(1,576); p_ww(happy_e0)= .65;
b001=  mkreconst(ll,1,1,p_ww);
calc_nf(b001)
rb3= irec(vis,vhs,b001);

ii=imgr([rg1,  rg3, rb1,  rb3]);
imwrite('tikhonov.png',ii,colormap);
image(ii);

function yerr= dataprior( ii, vis, vhs)
   zz= (vis - vhs) ./ (vis + vhs) /2;
   vh= prob_dir( zeros(size(ii)));
   vi= prob_dir( ii );
   zi= (vi - vh) ./ (vi + vh) *2;
   yerr= sumsq(zz-zi) / sumsq(zz);
endfunction

printf('dataprior rg1 = %f\n', dataprior(rg1, vis,vhs));
printf('dataprior rg3 = %f\n', dataprior(rg3, vis,vhs));
printf('dataprior rb1 = %f\n', dataprior(rb1, vis,vhs));
printf('dataprior rb3 = %f\n', dataprior(rb3, vis,vhs));
keyboard

% deform a model by Delta
% A1 = A0 + k1*sin(A0) + k2*cos(A0) + k3*sin(2*A0) ... 
function node1= angl_deform(node0, kk)
   A0 = atan2( node0(2,:), node0(1,:) );
   R0 = sqrt( sum( node0.^2 ));
   A1= A0;
   for k= 2:2:length(kk)
      m = k/2;
      A1= A1 + kk(k-1)*cos(m*A0) + kk(k)*sin(m*A0);
   end
   node1 = [R0.*cos(A1); R0.*sin(A1) ];
endfunction

function [vi,vh,dist] = make_deform( deform );
    global ChoiX;
    if ~strcmp(ChoiX, 'd001') 
       resetup('d0'); cleancolourmap;
    end
    global NODE; node0 = NODE;

    leye= [78,97,98,117,118,141];
    reye= [66,82,83,102,103,123];
    mouth= [28,31,40,43,58,57,53,54,74,73,92,93,113, ...
            112,135,69,87,70,107,88,108,129];

    node1= angl_deform(node0, deform);
    sad= zeros(1,256);
    vh= prob_dir( sad, node1 );
    sad(leye)=1;
    sad(reye)=1;
    sad(mouth)=.5;
    vi= prob_dir( sad, node1 );
% Find deformation of nodes at boundary
    bdy = 114:145;
    dist= mean(sumsq( node1(:,bdy) - node0(:,bdy) ));
endfunction

%    deform = .03*(ones(400,1)*(8:-1:1)).*randn(400,8);
     deform = [ ...
-0.3269493,-0.3918053,-0.0937888,-0.0018707,-0.0673451, 0.0110516,0.0718761   0.0170832;
-0.3455931, 0.0084441, 0.1094172,-0.1296894,-0.2014787, 0.0140448,-0.0650590   0.0236333;
-0.0898745,-0.1456434,-0.0411748,-0.0883531,-0.2377322,-0.0182819,0.1011687   0.0258459;
 0.2951051,-0.0881201,-0.2394838, 0.0369450,-0.0692927,-0.0078547,-0.0438427   0.0123168];

%
% STEP 4: Create image of 'random' deformations which happens
%    to result in a happytransform
%
    for i=1:size(deform,1)
        [vi,vh,dist] = make_deform( deform(i,:));
        vhd(:,i)= vh;
        vid(:,i)= vi;
        dists(i)= dist;
    end

  

    resetup('e0'); cleancolourmap;
ii=imgr(irec(vid,vhd,0,4));
image(ii,1);
imwrite('deform.png',ii,colormap);


% find numbers of elements which are adjoining given element
%     find_adjoin(3) => [7, 2, 4]
function elems= find_adjoin(ee)
   global ELEM;
   nn= ELEM(:,ee);
   [d,e]= size(ELEM);
   ss= zeros(1,e);
   for i=1:d
     ss= ss+ any(ELEM==nn(i));
   end
   elems= find(ss==d-1);
endfunction

% make pseudo laplacian filter
%   roielems = region to exclude boundary
%   if both ii and jj are in or out of roielems
%   then include the term
function D= mk_plaplace( roielems, ff )
   global ELEM;
   [d,e]= size(ELEM);
   ROI = zeros(1,e); ROI(roielems) = 1;
   Iidx= [];
   Jidx= [];
   Vidx= [];
   for ii=1:e
     el_adj = find_adjoin( ii );
%        ll= ones(size(el_adj));
%        Iidx= [Iidx, i, i*ll];
%        Jidx= [Jidx, i, el_adj];
%        Vidx= [Vidx, 3, -ll];
     for jj=el_adj(:)'
         if (ROI(ii) + ROI(jj)) == 1 %one or both
            fac= ff*.5;
         else 
            fac = .5;
         end
         Iidx= [Iidx,      ii, ii, jj, jj];
         Jidx= [Jidx,      ii, jj, ii, jj];
         Vidx= [Vidx, fac*([1, -1, -1,  1]) ];
     end
   end
   D= sparse(Iidx,Jidx, Vidx,e,e,'sum');
endfunction
   
% Create Reconstruction matrix
%    p_noise        ( 0 -> Hi = 1 
%                     1 -> Hi = 1/homg
%                     2 -> Hi = 1/homg^2 
%
%    Prior Weighting  p_ww (default = ones)
function Z= mkreconst_R( lamda, R, p_noise)
    global DVV;
    [m,e]= size(DVV);

    n_var= 1./prob_dir( zeros(1,e) );
    W= sparse(1:m,1:m, n_var.^(-p_noise) );
    dd= DVV'*W*DVV;
    dx= DVV'*W;
    Z= (dd+ lamda*R)\dx;
endfunction

% calc functions with laplacian regulies
% ff is amount of Non-penalty for elements
function ii=bad_reg_icrime_laplace(ff, ll);
    global vis vhs
    resetup('d0'); cleancolourmap;

    if ~exist('ll'); ll= .00025; end
    D= mk_plaplace([],ff);
    zg1=  mkreconst_R(ll,D,1);
    calc_nf(zg1)
    rg1= irec(vis,vhs,zg1);

    leye= [78,97,98,117,118,141];
    reye= [66,82,83,102,103,123];
    rsmile= [40:41, 53:55, 69:71, 86:88];
    lsmile= [43:44, 57:59, 73:75, 91:93];
    lsad = [31,43,58,57,74,73,92,93,113,112,135];
    sad = [28,31,40,43,58,57,53,54,74,73,92,93,113, ...
            112,135,69,87,70,107,88,108,129];

    DD= mk_plaplace([leye,reye,sad],ff);
    zg3=  mkreconst_R(ll,DD,1);
    calc_nf(zg3)
    rg3= irec(vis,vhs,zg3);

    DD= mk_plaplace([leye,reye,rsmile,lsad],ff);
    zb1=  mkreconst_R(ll,DD,1);
    calc_nf(zb1)
    rb1= irec(vis,vhs,zb1);

    DD= mk_plaplace([leye,reye,rsmile,lsmile],ff);
    zb3=  mkreconst_R(ll,DD,1);
    calc_nf(zb3)
    rb3= irec(vis,vhs,zb3);

    ii=imgr([rg1,  rg3, rb1,  rb3]);
endfunction

%
% STEP 5: Laplace regularization happytransform WITH i_crime
%
    ii=bad_reg_icrime_laplace(0,.00025);
    imwrite('laplace0-icrime.png',ii,colormap);
    image(ii);
    ii=bad_reg_icrime_laplace(.3,.00025);
    imwrite('laplace3-icrime.png',ii,colormap);
    image(ii);

%
% STEP 3: Tikhonov regularization happytransform WITHOUT i_crime
%
resetup('e0'); cleancolourmap;

ff= .3;
ll= .00025;
DD= mk_plaplace([],ff);
zg1=  mkreconst_R(ll,DD,1);
calc_nf(zg1)
rg1= irec(vis,vhs,zg1);

DD= mk_plaplace([sad_e0],ff);
zg3=  mkreconst_R(ll,DD,1);
calc_nf(zg3)
rg3= irec(vis,vhs,zg3);

DD= mk_plaplace([halfy_e0],ff);
zb1=  mkreconst_R(ll,DD,1);
calc_nf(zb1)
rb1= irec(vis,vhs,zb1);


DD= mk_plaplace([happy_e0],ff);
zb3=  mkreconst_R(ll,DD,1);
calc_nf(zb3)
rb3= irec(vis,vhs,zb3);

ii=imgr([rg1,  rg3, rb1,  rb3]);
imwrite('laplace3.png',ii,colormap);
image(ii);
