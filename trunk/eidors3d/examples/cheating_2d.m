% code to simulate inverse crimes in EIT

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: cheating_2d.m,v 1.29 2007-08-29 09:21:13 aadler Exp $

%TODO: calculate how well data matches priors
function out=cheating_2d( figno, rand_seed )
   [vis,vhs,s_mdl]= small_2d_mdl;

   if nargin<2; rand_seed= []; end

   il_g = make_inv_model( 12 ); %large model

   if nargin==0; figno= {'1','2a','2b','3a','3b','4'}; end 
   if isstr(figno); figno= {figno}; end

   for idx= 1:length(figno)
       if idx~=1; pause; end

       if     strcmp( figno{idx}, '1' )
           approach1(vis, vhs, s_mdl, il_g, rand_seed)
       elseif strcmp( figno{idx}, '2a' )
           approach2a(vis, vhs, s_mdl, il_g)
       elseif strcmp( figno{idx}, '2b' )
           approach2b(vis, vhs, s_mdl, il_g)
       elseif strcmp( figno{idx}, '3a' )
           approach3a(vis, vhs, s_mdl, il_g)
       elseif strcmp( figno{idx}, '3b' )
           approach3b(vis, vhs, s_mdl, il_g)
       elseif strcmp( figno{idx}, '4' )
           approach4(vis, vhs, s_mdl, il_g, rand_seed)
       elseif strcmp( figno{idx}, '4b' )
	   mdl= mk_common_model('b2c');
           if ~isempty(rand_seed); randn('state', rand_seed(i)); end
           def_mdl = eidors_obj('fwd_model', angl_deform(mdl.fwd_model, .02));
	   show_fem(def_mdl);
       end
   end
%
% APPROACH 1
%
function approach1(vis, vhs, s_mdl, il_g, rand_seed)
   disp('Approach #1: reconstruct with noise');

   num_tries=6;
   levels= [];
   if ~isempty(rand_seed);
      num_tries= length(rand_seed);
      levels= [0,0,0,1,1];
   end

   vi_n(1:num_tries) = vis;
   for i= 1:num_tries % stupid matlab doesn't allow easy vectorization
      if ~isempty(rand_seed); randn('state', rand_seed(i)); end
      noise = .0002*randn(size(vis.meas));
      vi_n(i).meas = vi_n(i).meas + noise;
   end
   show_slices( inv_solve( il_g, vhs, vi_n ), levels);


%
% APPROACH 2A
%
function approach2a(vis, vhs, s_mdl, il_g)
   levels= [0,0,0,1,1];
   disp('Approach #2: reconstruct with tikhonov cheat - with inv crime');
   pp= small_face;
   % homog (normal) model
   RtR_prior= @cheat_tikhonov;
   param_vals.cheat_elements= [];
   param_vals.cheat_weight = 0.5;
   is_n = make_inv_model( 8 , RtR_prior, 'cheat_tikhonov', param_vals ); 
   % sad model
   param_vals.cheat_elements= [pp.eyes, pp.sad];
   is_s = make_inv_model( 8 , RtR_prior, 'cheat_tikhonov', param_vals ); 
   % happy model
   param_vals.cheat_elements= [pp.eyes, pp.smile];
   is_h = make_inv_model( 8 , RtR_prior, 'cheat_tikhonov', param_vals ); 
   % happy/sad (medium) model
   param_vals.cheat_elements= [pp.eyes, pp.rsmile, pp.lsad];
   is_m = make_inv_model( 8 , RtR_prior, 'cheat_tikhonov', param_vals ); 

   show_slices( [ inv_solve( is_n, vhs, vis ), ... 
                  inv_solve( is_s, vhs, vis ), ... 
                  inv_solve( is_h, vhs, vis ), ... 
                  inv_solve( is_m, vhs, vis ) ], levels );



%
% APPROACH 2B
%
function approach2b(vis, vhs, s_mdl, il_g)
   disp('Approach #2B: reconstruct with tikhonov cheat - without inv crime');
   levels= [0,0,0,1,1];
   pp= large_face;
   % homog (normal) model
   RtR_prior= @cheat_tikhonov;
   param_vals.cheat_elements= [];
   param_vals.cheat_weight = 0.5;
   il_n = make_inv_model(12 , RtR_prior, 'cheat_tikhonov', param_vals ); 
   % sad model
   param_vals.cheat_elements= pp.sad;
   il_s = make_inv_model(12 , RtR_prior, 'cheat_tikhonov', param_vals ); 
   % happy model
   param_vals.cheat_elements= pp.happy;
   il_h = make_inv_model(12 , RtR_prior, 'cheat_tikhonov', param_vals ); 
   % happy/sad (medium) model
   param_vals.cheat_elements= pp.halfy;
   il_m = make_inv_model(12 , RtR_prior, 'cheat_tikhonov', param_vals ); 

   show_slices( [ inv_solve( il_n, vhs, vis ), ... 
                  inv_solve( il_s, vhs, vis ), ... 
                  inv_solve( il_h, vhs, vis ), ... 
                  inv_solve( il_m, vhs, vis ) ], levels );

%
% APPROACH 3A
%
function approach3a(vis, vhs, s_mdl, il_g)
   disp('Approach #3: reconstruct with Laplace filter cheat - with inv crime');
   levels= [0,0,0,1,1];
   pp= small_face;
   % homog (normal) model
   RtR_prior= @cheat_laplace;
   param_vals.cheat_elements= [];
   param_vals.cheat_weight = 0.2;
   is_n = make_inv_model( 8 , RtR_prior, 'cheat_laplace', param_vals ); 
   % sad model
   param_vals.cheat_elements= [pp.eyes, pp.sad];
   is_s = make_inv_model( 8 , RtR_prior, 'cheat_laplace', param_vals ); 
   % happy model
   param_vals.cheat_elements= [pp.eyes, pp.smile];
   is_h = make_inv_model( 8 , RtR_prior, 'cheat_laplace', param_vals ); 
   % happy/sad (medium) model
   param_vals.cheat_elements= [pp.eyes, pp.rsmile, pp.lsad];
   is_m = make_inv_model( 8 , RtR_prior, 'cheat_laplace', param_vals ); 

   show_slices( [ inv_solve( is_n, vhs, vis ), ... 
                  inv_solve( is_s, vhs, vis ), ... 
                  inv_solve( is_h, vhs, vis ), ... 
                  inv_solve( is_m, vhs, vis ) ], levels );


%
% APPROACH 3B
%
function approach3b(vis, vhs, s_mdl, il_g)
   disp('Approach #3B: reconstruct with Laplace cheat - without inv crime');
   levels= [0,0,0,1,1];
   pp= large_face;
   % homog (normal) model
   RtR_prior= @cheat_laplace;
   param_vals.cheat_elements= [];
   param_vals.cheat_weight = 0.2;
   il_n = make_inv_model(12 , RtR_prior, 'cheat_laplace', param_vals ); 
   % sad model
   param_vals.cheat_elements= pp.sad;
   il_s = make_inv_model(12 , RtR_prior, 'cheat_laplace', param_vals ); 
   % happy model
   param_vals.cheat_elements= pp.happy;
   il_h = make_inv_model(12 , RtR_prior, 'cheat_laplace', param_vals ); 
   % happy/sad (medium) model
   param_vals.cheat_elements= pp.halfy;
   il_m = make_inv_model(12 , RtR_prior, 'cheat_laplace', param_vals ); 

   show_slices( [ inv_solve( il_n, vhs, vis ), ... 
                  inv_solve( il_s, vhs, vis ), ... 
                  inv_solve( il_h, vhs, vis ), ... 
                  inv_solve( il_m, vhs, vis ) ], levels );


%
% APPROACH 4
%
function approach4(vis, vhs, s_mdl, il_g, rand_seed)
   disp('Approach #4: deform the model');
   num_tries=6;
   levels= [];
   if ~isempty(rand_seed);
      num_tries= length(rand_seed);
      levels= [0,0,0,1,1];
   end

   params= mk_circ_tank(8, [], 16 ); 
   params.stimulation= mk_stim_patterns(16, 1, '{ad}','{ad}', ...
                         {'no_meas_current','no_rotate_meas'}, 1);
   params.solve=      'aa_fwd_solve';
   params.system_mat= 'aa_calc_system_mat';

   mat= ones( size(params.elems,1), 1);
   pp= small_face;
   mat(pp.eyes)= 2;
   mat(pp.sad)=1.5;

   def_amount= .0035;

   for i= 1:num_tries % stupid matlab doesn't allow easy vectorization
      if ~isempty(rand_seed); randn('state', rand_seed(i)); end

      def_mdl = eidors_obj('fwd_model', angl_deform(params, def_amount ) );
      vi_m(i)= fwd_solve( eidors_obj('image','name',  ...
                     'elem_data', mat, 'fwd_model', def_mdl ));
   end
   show_slices( inv_solve( il_g, vhs, vi_m ), levels);



   
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
function i_mdl= make_inv_model( n_rings, img_prior, param_name, param_vals );
   params= mk_circ_tank(n_rings, [], 16 ); 
   params.stimulation= mk_stim_patterns(16, 1, '{ad}','{ad}', ...
                         {'no_meas_current','no_rotate_meas'}, 1);
   params.solve=      'aa_fwd_solve';
   params.system_mat= 'aa_calc_system_mat';
   params.jacobian  = 'aa_calc_jacobian';
%  params.normalize_measurements  = 1; TODO: we have a bug here
   l_mdl= eidors_obj('fwd_model', params);

% create inverse model
   hparam.value = 3e-4;
  %hparam.func = 'select_noise_figure';
  %hparam.noise_figure= 1;
  %hparam.tgt_elems= 1:4;

   if nargin < 2
      img_prior = 'tikhonov_image_prior';
      param_name= 'jnk___'; param_vals= 'jnk___';
   end

   i_mdl= eidors_obj( ...
          'inv_model', '2D inverse', ...
          'hyperparameter', hparam, 'RtR_prior', img_prior, ...
          param_name, param_vals, ...
          'reconst_type', 'difference', ...
          'fwd_model', l_mdl, 'solve', 'aa_inv_solve' );
   i_mdl.jacobian_bkgnd.value= 1;


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
%   elems    = inv_model.RtR_prior.cheat_elements;
%            => elements weights to modify
%   weight   = inv_model.RtR_prior.cheat_weight;
%            => new weight to set elements to


pp= aa_fwd_parameters( inv_model.fwd_model );
idx= 1:pp.n_elem;
weight= ones(1,pp.n_elem);
weight( inv_model.cheat_tikhonov.cheat_elements ) = ...
        inv_model.cheat_tikhonov.cheat_weight;

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
%   elems    = inv_model.RtR_prior.cheat_elements;
%            => elements weights to modify
%   weight   = inv_model.RtR_prior.cheat_weight;
%            => new weight to set elements to

% make pseudo laplacian filter
%   roielems = region to exclude boundary
%   if both ii and jj are in or out of roielems
%   then include the term
   pp= aa_fwd_parameters( inv_model.fwd_model );

   ROI = zeros(1,pp.n_elem);
   ROI( inv_model.cheat_laplace.cheat_elements ) = 1;

   Iidx= [];
   Jidx= [];
   Vidx= [];
   for ii=1:pp.n_elem
     el_adj = find_adjoin( ii, pp.ELEM );
     for jj=el_adj(:)'
         if (ROI(ii) + ROI(jj)) == 1 %one only
            fac= inv_model.cheat_laplace.cheat_weight *.5;
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
function mdl1 = angl_deform(mdl0, def_amount );

   node0= mdl0.nodes';
   deform = def_amount*(8:-1:1).*randn(1,8);

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

