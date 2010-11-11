function inv_mdl= select_imdl( mdl, options )
% SELECT_IMDL: select pre-packaged inverse model features
% inv_mdl= select_imdl( mdl, options )
%
%   mdl - inv_model structure - parameters are replaced as specified
%    OR
%   mdl - fwd_model - a basic GN difference model is created, and parameters replaced
%
% OPTIONS => {'opt1','opt2'} options are processed in the order specified
%
% Available options are:
%
% 'Basic GN dif';   Basic GN one step difference solver with Laplace prior
% 'Basic GN abs';   Basic Gauss-Newton absolute solver with Laplace prior
% 'NOSER dif';      Basic GN one step difference solver with NOSER prior 
% 'Nodal GN dif';   Basic GN solver, solves onto nodes
% 'TV solve dif';   Total Variation PDIPM difference solver 
% 'Elec Move GN';   One step GN difference solver with compensation for electrode movement
% 'Choose NF=1.0';  Choose hyperparameter value appropriate for specified noise figure (NF)

% (C) 2010 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isstr(mdl) && strcmp(mdl,'UNIT_TEST'); do_unit_test; return; end

if nargin == 1; options = {}; end

switch mdl.type
  case 'inv_model'; inv_mdl = mdl;
  case 'fwd_model'; inv_mdl = basic_imdl( mdl );
  otherwise;        error('select_imdl: expects inv_model or fwd_model input');
end

for i=1:length(options);
  % split the option on an equals sign
  [s,f,tok]= regexp(options{i},'(.[^=]*)=?(.*)');
  tok = tok{1};
  for t=1:size(tok,1);
     opt{t} = options{i}(tok(t,1):tok(t,2));
  end
% This is the Mat7 only code:
%   opt= regexp(options{i},'(.[^=]*)=?(.*)','tokens');
%   opt= opt{1};
  switch opt{1}
    case 'NOSER dif';       inv_mdl = NOSER_dif( inv_mdl );
    case 'Basic GN dif';    inv_mdl = Basic_GN_Dif( inv_mdl );
    case 'Basic GN abs';    inv_mdl = Basic_GN_Abs( inv_mdl );
    case 'Nodal GN dif';    inv_mdl = Nodal_GN_Dif( inv_mdl );
    case 'TV solve dif';    inv_mdl = TV_solve_Dif( inv_mdl );
    case 'Elec Move GN';    inv_mdl = Elec_Move_GN( inv_mdl );
    case 'Choose NF';       inv_mdl = Choose_NF( inv_mdl, str2num(opt{2}) );
    
    otherwise; error('option {%s} not understood', options{i});
  end
end

function imdl = basic_imdl( fmdl );
   imdl.name= 'Basic imdl from select_imdl';
   imdl.type= 'inv_model';

   imdl.solve= @aa_inv_solve;
   imdl.hyperparameter.value = .01;
   imdl.RtR_prior = @laplace_image_prior;
   imdl.jacobian_bkgnd.value = 1;
   imdl.reconst_type= 'difference';
   imdl.fwd_model = fmdl;

function imdl = NOSER_dif( imdl );
   imdl.RtR_prior = @noser_image_prior;
   try; imdl = rmfield(imdl,'R_prior'); end
   imdl.hyperparameter.value = .03;
   imdl.solve= @aa_inv_solve;
   imdl.reconst_type= 'difference';

function imdl = Basic_GN_Dif( imdl );
   imdl.RtR_prior = @laplace_image_prior;
   try; imdl = rmfield(imdl,'R_prior'); end
   imdl.solve= @aa_inv_solve;
   imdl.reconst_type= 'difference';

function imdl = Basic_GN_Abs( imdl );
   imdl.RtR_prior = @laplace_image_prior;
   try; imdl = rmfield(imdl,'R_prior'); end
   imdl.solve= @GN_abs_solve;
   imdl.parameters.max_iterations= 10;
   imdl.reconst_type= 'absolute';

function imdl = TV_solve_Dif( imdl );
   imdl.R_prior = @ab_calc_tv_prior;
   try; imdl = rmfield(imdl,'RtR_prior'); end
   imdl.solve= @ab_tv_diff_solve;
   imdl.reconst_type= 'difference';
   imdl.parameters.max_iterations= 15;
   imdl.hyperparameter.value = 1e-1;
   imdl.parameters.term_tolerance = 1e-3;

function imdl = Elec_Move_GN( imdl );
   % keep previous model as conductivity jacobian, so it should be ok
   imdl.fwd_model.conductivity_jacobian = imdl.fwd_model.jacobian; 
   imdl.fwd_model.jacobian = @aa_e_move_jacobian;
   imdl.RtR_prior =          @aa_e_move_image_prior;
   imdl.solve= @aa_inv_solve;

   MV_prior = 1./mean( std( imdl.fwd_model.nodes ));
   imdl.aa_e_move_image_prior.parameters = MV_prior;
   imdl.aa_e_move_image_prior.RegC.func = @gaussian_HPF_prior;
   imdl.hyperparameter.value = .03;

   n_elems = size(imdl.fwd_model.elems,1);
   imdl.inv_solve.select_parameters = 1:n_elems;

   imdl.prior_use_fwd_not_rec = 1; % for c2f mapping

function imdl = Nodal_GN_Dif( imdl );
   imdl.solve = @nodal_solve;


function imdl = Choose_NF( imdl, NF_req );
   if ~strcmp(imdl.reconst_type, 'difference');
      error('Choose NF only works for difference solvers right now');
   end

% Find 4 elems in mesh ctr to be the NF target elems
   xyz_elems = interp_mesh( imdl.fwd_model );
   ctr_elems = mean(xyz_elems, 1);
   xyz_elems = xyz_elems - ones(size(xyz_elems,1),1)*ctr_elems;
   d_elems   = sqrt( sum( xyz_elems.^2, 2 ));
   [jnk, e_idx] = sort(d_elems);

   imdl.hyperparameter.tgt_elems = e_idx(1:4);
   imdl.hyperparameter.noise_figure = NF_req;

   sv_log = eidors_msg('log_level'); eidors_msg('log_level', 1);
   HP = choose_noise_figure( imdl );
   eidors_msg('log_level', sv_log);

   imdl.hyperparameter.value = HP;


function do_unit_test
% Test difference solvers on lung images
   load montreal_data_1995;
   imdl = mk_common_model('b2t3',16); 

   for i=1:100; % Break when finished
      vh = zc_resp(:,1); vi= zc_resp(:,23);
      switch i
         case 01;
            imdl0 = select_imdl( imdl );
         case 02;
            imdl0 = select_imdl( imdl.fwd_model );
         case 03;
            imdl0 = select_imdl( imdl, {'NOSER dif'} );
         case 04;
            imdl0 = select_imdl( imdl, {'NOSER dif','Choose NF=1.1'});
         case 05;
            imdl0 = select_imdl( imdl, {'Basic GN dif'} );
         case 06;
            imdl0 = select_imdl( imdl, {'TV solve dif'} );
         case 07;
            imdl0 = select_imdl( imdl.fwd_model, {'Basic GN dif','Elec Move GN','Choose NF=0.5'} );
         case 08;
            imdl0 = select_imdl( imdl, {'Elec Move GN'} );
         case 09;
            imdl0 = mk_common_model('b2C2',16); 
            imdl0 = select_imdl( imdl0, {'Elec Move GN'} );
         case 10;
            imdl0 = select_imdl( imdl, {'Nodal GN dif'} );
         case 11;
            imdl0 = select_imdl( imdl, {'Nodal GN dif', 'Choose NF=0.50'} );
         case 12;
            imdl0 = mk_common_model('b2C2',16); 
            imdl0 = select_imdl( imdl0, {'Basic GN dif', 'TV solve dif'} );
            imdl0.parameters.max_iterations= 2;
            imdl0 = select_imdl( imdl0, {'Choose NF=0.8'} );
            [vh,vi] = simulate_movement(mk_image(imdl0), [0;0.5;0.35]);
         case 13;
            imdl0 = mk_common_model('b2C2',16); 
            imdl0 = select_imdl( imdl0, {'Nodal GN dif', 'Choose NF=0.50'} );
            [vh,vi] = simulate_movement(mk_image(imdl0), [0;0.5;0.05]);
         case 14;
            imdl0 = mk_common_model('b2C2',16); 
            imdl0 = select_imdl( imdl0, {'Basic GN abs'} );
            [vh,vi] = simulate_movement(mk_image(imdl0), [0;0.5;0.05]);
         case 15; break
      end;

%     disp([i,imdl0.hyperparameter.value]);
      if strcmp( imdl0.reconst_type, 'absolute')
         imgr = inv_solve( imdl0, vi);
      else
         imgr = inv_solve( imdl0, vh, vi);
      end
      subplot(4,4,i); show_slices( imgr );
   end
