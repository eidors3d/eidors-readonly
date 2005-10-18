% code to simulate inverse crimes in EIT
% $Id: cheating_2d.m,v 1.4 2005-10-18 16:49:37 aadler Exp $

%TODO: calculate how well data matches priors
function out=cheating_2d

   [vis,vhs,s_mdl]= small_2d_mdl;

   il_g = make_inv_model( 12 ); %large model
   disp('Approach #1: reconstruct with noise');
%  reconst_with_noise(il_g, vis, vhs, 10); pause

   disp('Approach #2: reconstruct with tikhonov cheat - with inv crime');
   pp= small_face;
   % homog (normal) model
   image_prior.func= @cheat_tikhonov;
   image_prior.cheat_elements= [];
   image_prior.cheat_weight = 0.5;
   is_n = make_inv_model( 8 , image_prior ); 
   % sad model
   image_prior.cheat_elements= [pp.eyes, pp.sad];
   is_s = make_inv_model( 8 , image_prior ); 
   % happy model
   image_prior.cheat_elements= [pp.eyes, pp.smile];
   is_h = make_inv_model( 8 , image_prior ); 
   % happy/sad (medium) model
   image_prior.cheat_elements= [pp.eyes, pp.rsmile, pp.lsad];
   is_m = make_inv_model( 8 , image_prior ); 

   show_slices( [ inv_solve( is_n, vis, vhs ), ... 
                  inv_solve( is_s, vis, vhs ), ... 
                  inv_solve( is_h, vis, vhs ), ... 
                  inv_solve( is_m, vis, vhs ) ] ); pause



   disp('Approach #3: reconstruct with tikhonov cheat - without inv crime');
   pp= large_face;
   % homog (normal) model
   image_prior.func= @cheat_tikhonov;
   image_prior.cheat_elements= [];
   image_prior.cheat_weight = 0.5;
   il_n = make_inv_model(12 , image_prior ); 
   % sad model
   image_prior.cheat_elements= pp.sad;
   il_s = make_inv_model(12 , image_prior ); 
   % happy model
   image_prior.cheat_elements= pp.happy;
   il_h = make_inv_model(12 , image_prior ); 
   % happy/sad (medium) model
   image_prior.cheat_elements= pp.halfy;
   il_m = make_inv_model(12 , image_prior ); 

   show_slices( [ inv_solve( il_n, vis, vhs ), ... 
                  inv_solve( il_s, vis, vhs ), ... 
                  inv_solve( il_h, vis, vhs ), ... 
                  inv_solve( il_m, vis, vhs ) ] ); pause

return;

function p= small_face;
   p.leye=   [78,97,98,117,118,141];
   p.reye=   [66,82,83,102,103,123];
   p.rsmile= [40:41, 53:55, 69:71, 86:88];
   p.lsmile= [43:44, 57:59, 73:75, 91:93];
   p.lsad =  [31,43,58,57,74,73,92,93,113,112,135];
   p.sad =   [28,31,40,43,58,57,53,54,74,73,92,93,113, ...
              112,135,69,87,70,107,88,108,129];
   p.eyes= [p.leye,p.reye];
   p.smile= [p.rsmile, p.lsmile];

function p=large_face
p.sad=  [ 53; 57; 69; 73; 86; 87; 91; 92; 106; 107; 111; 112; 127; 128;
         129; 133; 134; 135; 147; 151; 152; 153; 157; 158; 159; 165;
         171; 172; 177; 178; 179; 184; 185; 186; 192; 193; 199; 200;
         205; 206; 207; 212; 213; 214; 220; 221; 227; 228; 229; 235;
         236; 243; 244; 251; 252; 253; 259; 260; 261; 267; 268; 275;
         276; 283; 284; 285; 292; 293; 301; 310; 319; 320]; 

p.happy= [ ...
        69; 70; 71; 73; 74; 75; 86; 87; 88; 89; 91; 92; 93; 94; 106;
       107; 108; 109; 111; 112; 113; 114; 127; 128; 129; 130; 131; 133;
       134; 135; 136; 137; 147; 151; 152; 153; 154; 155; 157; 158; 159;
       160; 161; 165; 171; 172; 176; 177; 178; 179; 180; 183; 184; 185;
       186; 187; 192; 193; 199; 200; 220; 221; 227; 228; 229; 251; 252;
       253; 259; 260; 261; 283; 284; 285; 292; 293; 319; 320]; 

p.halfy= [ ...
        57; 69; 70; 71; 73; 86; 87; 88; 89; 91; 92; 106; 107; 108; 109;
       111; 112; 127; 128; 129; 130; 131; 133; 134; 135; 147; 151; 152;
       153; 154; 155; 157; 158; 159; 165; 171; 172; 176; 177; 178; 179;
       180; 184; 185; 186; 192; 193; 199; 200; 212; 213; 214; 220; 221;
       227; 228; 229; 243; 244; 251; 252; 253; 259; 260; 261; 275; 276;
       283; 284; 285; 292; 293; 310; 319; 320];

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
   pp= small_face;
   mat(pp.eyes)= 2;
   mat(pp.sad)=1.5;

   vis= fwd_solve( eidors_obj('image','name',  ...
                     'elem_data', mat, 'fwd_model', mdl ));

% create inv_model, and specify img_prior if required
function i_mdl= make_inv_model( n_rings, img_prior );
   params= mk_circ_tank(12, [], 16 ); 
   params.stimulation= mk_stim_patterns(16, 1, '{ad}','{ad}', ...
                         {'no_meas_current','no_rotate_meas'}, 1);
   params.solve=      'aa_fwd_solve';
   params.system_mat= 'aa_calc_system_mat';
   params.jacobian  = 'aa_calc_jacobian';
%  params.normalize_measurements  = 1; TODO: we have a bug here
   l_mdl= eidors_obj('fwd_model', params);

% create inverse model
   hparam.value = 1e-8;
  %hparam.func = 'aa_calc_noise_figure';
  %hparam.noise_figure= 1;
  %hparam.tgt_elems= 1:4;

   if nargin < 2
      img_prior.func = 'tikhonov_image_prior';
     %img_prior.func = 'aa_calc_image_prior';
   end

   i_mdl= eidors_obj( ...
          'inv_model', '2D inverse', ...
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

% Image prior with Tikhonov + cheating elements
function Reg= cheat_tikhonov( inv_model, elems, weight );
% Reg= cheat_tikhonov( inv_model )
% Reg        => output regularization term
% Parameters:
%   elems    = inv_model.image_prior.cheat_elements;
%            => elements weights to modify
%   weight   = inv_model.image_prior.cheat_weight;
%            => new weight to set elements to


pp= aa_fwd_parameters( inv_model.fwd_model );
idx= 1:pp.n_elem;
weight= ones(1,pp.n_elem);
weight( inv_model.image_prior.cheat_elements ) = ...
        inv_model.image_prior.cheat_weight;

Reg = sparse( idx, idx, weight );

function boo;


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
