function make_figures(figno)

   switch figno
     case 1;
       imdl = mk_common_model('c2c2',16);
       imdl.fwd_model.normalize_measurements= 0;
       imdl.hyperparameter.value = 0.0439;

     case 2;
       imdl = mk_common_model('c2c2',16);
       imdl.fwd_model.normalize_measurements= 1;
       imdl.hyperparameter.value = 0.0189;

     case 3;
       imdl = mk_common_model('c2c2',16);
       [f_mdl,c2f]= fmdl3d( imdl.fwd_model);
       imdl.rec_model = imdl.fwd_model;
       imdl.fwd_model = f_mdl;
       imdl.fwd_model.coarse2fine = c2f;
       imdl.fwd_model.normalize_measurements= 0;
       imdl.hyperparameter.value = 0.00449;

     case 4;
       imdl = mk_common_model('c2c2',16);
       [f_mdl,c2f]= fmdl3d( imdl.fwd_model);
       imdl.rec_model = imdl.fwd_model;
       imdl.fwd_model = f_mdl;
       imdl.fwd_model.coarse2fine = c2f;
       imdl.fwd_model.normalize_measurements= 1;
       imdl.hyperparameter.value = 0.00190;


     otherwise
       error('huh?')
   end
   calc_noise_figure( set_tgts( imdl ));
%  fprintf('hp=%1.8g\n', choose_noise_figure( set_tgts(imdl, 0.79)) );

   [eelv,eilv] = pig_data;
    vh = eelv.p20; 
    vi = eilv.p20; 


   img= inv_solve(imdl, vh, vi);
   show_fem(img);

   
function imdl = set_tgts( imdl, nf )
   imdl00 = mk_common_model('a2c0',16);
   meas_sel = imdl00.fwd_model.meas_select;   
   [eelv,eilv] = pig_data;
    vh = eelv.p20(meas_sel); 
    vi = eilv.p20(meas_sel); 
   imdl.hyperparameter.tgt_data.meas_t1= vh;
   imdl.hyperparameter.tgt_data.meas_t2= vi;

   if nargin==2;
      imdl.hyperparameter.noise_figure = nf;
   end
 


function [f_mdl, c2f] = fmdl3d( c_mdl );
   if ~exist('ng_mdl_16x1_coarse.mat','file')
      !wget http://eidors3d.sf.net/data_contrib/netgen_moving_ball/ng_mdl_16x1_coarse.7z
      !C:\progra~1\7-zip\7z e ng_mdl_16x1_coarse.7z ng_mdl_16x1_coarse.mat
   end

   load ng_mdl_16x1_coarse;
   f_mdl = ng_mdl_16x1_coarse;

   f_mdl.stimulation = c_mdl.stimulation;
   f_mdl.meas_select = c_mdl.meas_select;
   f_mdl.solve = @aa_fwd_solve;
   f_mdl.jacobian = @aa_calc_jacobian;
   f_mdl.system_mat = @aa_calc_system_mat;

   scl= 15;
   c_mdl.mk_coarse_fine_mapping.f2c_offset = [0,0,scl];
   c_mdl.mk_coarse_fine_mapping.f2c_project = (1/scl)*speye(3);
   c_mdl.mk_coarse_fine_mapping.z_depth = inf;
   c2f= mk_coarse_fine_mapping( f_mdl, c_mdl);

function [eelv,eilv] = pig_data
   if ~exist('p1130107.get','file');
      !wget http://eidors3d.sf.net/data_contrib/if-peep-acute-lung-injury/if_data_2003.zip
      !unzip if_data_2003.zip p1130107.get
   end
   vv= eidors_readdata('p1130107.get');
   eelv.p0  = vv(:,72);  eilv.p0  = vv(:,80);
   eelv.p20 = vv(:,554); eilv.p20 = vv(:,561);
   eelv.p30 = vv(:,795); eilv.p30 = vv(:,803);

