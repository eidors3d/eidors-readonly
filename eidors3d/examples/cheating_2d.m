% code to simulate inverse crimes in EIT

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: cheating_2d.m,v 1.8 2005-10-31 01:52:11 aadler Exp $

%TODO: calculate how well data matches priors
function out=cheating_2d

   [vis,vhs,s_mdl]= small_2d_mdl;

   il_g = make_inv_model( 12 ); %large model

%
% APPROACH 1
%
   disp('Approach #1: reconstruct with noise');

   num_tries=12;

   noise = .0002*randn(size(vis.meas,1),num_tries);

   vi_n(1:num_tries) = vis;
   for i= 1:num_tries % stupid matlab doesn't allow easy vectorization
      vi_n(i).meas = vi_n(i).meas + noise(:,i);
   end
   show_slices( inv_solve( il_g, vhs, vi_n ));

   levels= [0,0,0,1,1];
%
% APPROACH 2
%
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

   show_slices( [ inv_solve( is_n, vhs, vis ), ... 
                  inv_solve( is_s, vhs, vis ), ... 
                  inv_solve( is_h, vhs, vis ), ... 
                  inv_solve( is_m, vhs, vis ) ], levels ); pause



%
% APPROACH 2B
%
   disp('Approach #2B: reconstruct with tikhonov cheat - without inv crime');
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

   show_slices( [ inv_solve( il_n, vhs, vis ), ... 
                  inv_solve( il_s, vhs, vis ), ... 
                  inv_solve( il_h, vhs, vis ), ... 
                  inv_solve( il_m, vhs, vis ) ] ); pause

%
% APPROACH 3
%
   disp('Approach #3: reconstruct with Laplace filter cheat - with inv crime');
   pp= small_face;
   % homog (normal) model
   image_prior.func= @cheat_laplace;
   image_prior.cheat_elements= [];
   image_prior.cheat_weight = 0.2;
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

   show_slices( [ inv_solve( is_n, vhs, vis ), ... 
                  inv_solve( is_s, vhs, vis ), ... 
                  inv_solve( is_h, vhs, vis ), ... 
                  inv_solve( is_m, vhs, vis ) ] ); pause


%
% APPROACH 3B
%
   disp('Approach #3B: reconstruct with Laplace cheat - without inv crime');
   pp= large_face;
   % homog (normal) model
   image_prior.func= @cheat_laplace;
   image_prior.cheat_elements= [];
   image_prior.cheat_weight = 0.2;
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

   show_slices( [ inv_solve( il_n, vhs, vis ), ... 
                  inv_solve( il_s, vhs, vis ), ... 
                  inv_solve( il_h, vhs, vis ), ... 
                  inv_solve( il_m, vhs, vis ) ] ); pause


%
% APPROACH 4
%
   disp('Approach #4: deform the model');

   params= mk_circ_tank(8, [], 16 ); 
   params.stimulation= mk_stim_patterns(16, 1, '{ad}','{ad}', ...
                         {'no_meas_current','no_rotate_meas'}, 1);
   params.solve=      'aa_fwd_solve';
   params.system_mat= 'aa_calc_system_mat';

   mat= ones( size(params.elems,1), 1);
   pp= small_face;
   mat(pp.eyes)= 2;
   mat(pp.sad)=1.5;

   for i= 1:num_tries % stupid matlab doesn't allow easy vectorization
      def_mdl = eidors_obj('fwd_model', angl_deform(params ) );
      vi_m(i)= fwd_solve( eidors_obj('image','name',  ...
                     'elem_data', mat, 'fwd_model', def_mdl ));
   end
   show_slices( inv_solve( il_g, vhs, vi_m ));



   
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
   params= mk_circ_tank(n_rings, [], 16 ); 
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
function Reg= cheat_tikhonov( inv_model )
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

% find elems which are connected to elems ee
function elems= find_adjoin(ee, ELEM)
   nn= ELEM(:,ee);
   [d,e]= size(ELEM);
   ss= zeros(1,e);
   for i=1:d
     ss= ss+ any(ELEM==nn(i));
   end
   elems= find(ss==d-1);

% Image prior with FEM based Laplacian + cheating
function Reg= cheat_laplace( inv_model )
% Reg= cheat_laplace( inv_model )
% Reg        => output regularization term
% Parameters:
%   elems    = inv_model.image_prior.cheat_elements;
%            => elements weights to modify
%   weight   = inv_model.image_prior.cheat_weight;
%            => new weight to set elements to

% make pseudo laplacian filter
%   roielems = region to exclude boundary
%   if both ii and jj are in or out of roielems
%   then include the term
   pp= aa_fwd_parameters( inv_model.fwd_model );

   ROI = zeros(1,pp.n_elem);
   ROI( inv_model.image_prior.cheat_elements ) = 1;

   Iidx= [];
   Jidx= [];
   Vidx= [];
   for ii=1:pp.n_elem
     el_adj = find_adjoin( ii, pp.ELEM );
     for jj=el_adj(:)'
         if (ROI(ii) + ROI(jj)) == 1 %one only
            fac= inv_model.image_prior.cheat_weight *.5;
         else 
            fac = .5;
         end
         Iidx= [Iidx,      ii, ii, jj, jj];
         Jidx= [Jidx,      ii, jj, ii, jj];
         Vidx= [Vidx, fac*([1, -1, -1,  1]) ];
     end
   end
   Reg = sparse(Iidx,Jidx, Vidx, pp.n_elem, pp.n_elem );

% deform a model by Delta
% A1 = A0 + k1*sin(A0) + k2*cos(A0) + k3*sin(2*A0) ... 
function mdl1 = angl_deform(mdl0 );

   node0= mdl0.nodes';
   deform = .0001*(ones(400,1)*(8:-1:1)).*randn(400,8);

   A0 = atan2( node0(2,:), node0(1,:) );
   R0 = sqrt( sum( node0.^2 ));
   A1= A0;
   for k= 2:2:length(deform)
      m = k/2;
      A1= A1 + deform(k-1)*cos(m*A0) + deform(k)*sin(m*A0);
   end
   node1 = [R0.*cos(A1); R0.*sin(A1) ];

   mdl1 = mdl0;
   mdl1.nodes = node1';

function boo;


   
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
