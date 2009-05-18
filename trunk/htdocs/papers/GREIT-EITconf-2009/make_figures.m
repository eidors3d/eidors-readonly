function make_figures( fig_list)
% Make figures for EIT 2009 paper 
   if nargin==0; fig_list = 0:4; end
   
   clf;
   for i= fig_list
      axes('position',[i*0.2,0.42,0.2,0.58]);
      make_fig(i,0) 
      axis normal
      axes('position',[i*0.2,-0.06,0.2,0.58]);
      make_fig(i,1) 
      axis normal
   end
   ax= get(gcf,'paperposition')
   set(gcf,'paperposition',[ax(1:3),ax(3)/2.5])
   print -depsc2 fig2.eps
   set(gcf,'paperposition',[ax])
   


function make_fig(figno,imgno)
   switch figno
     case 0;
       imdl = mk_common_gridmdl('backproj');

     case {1,2};
       imdl = mk_common_model('c2c2',16);
       load ng_mdl_16x1_coarse;  ngm = ng_mdl_16x1_coarse;
%      load ng_mdl_16x1_fine;    ngm = ng_mdl_16x1_fine;
%      load ng_mdl_16x1_vfine;   ngm = ng_mdl_16x1_vfine;
       imdl.fwd_model.elems      = ngm.elems;
       imdl.fwd_model.nodes      = ngm.nodes(:,[2,1,3])/15; % Flip orientation
       imdl.fwd_model.nodes(:,3) =(ngm.nodes(:,3) - mean(ngm.nodes(:,3)))/15;
       imdl.fwd_model.electrode  = ngm.electrode;
       imdl.fwd_model.boundary   = ngm.boundary;
       weight = 25;


     case {3,4};
       imdl = eidors_obj('get-cache', 'paper09', 'd2d4c_model');
       if ~isempty(imdl)
          imdl = mk_common_model('d2d4c',16);
          eidors_obj('set-cache', 'paper09', 'd2d4c_model', imdl);
       end
       weight = 25;

     otherwise
       error('huh?')
   end

   switch figno;
     case 1; weight = 16; % NF 0.71
     case 2; weight = 23; % NF 0.78
     case 3; weight = 30; % NF 0.74
     case 4; weight = 40; % NF 0.75 
   end

   switch figno;
      case {0,1,3}; imdl.fwd_model.normalize_measurements= 1;
      otherwise;    imdl.fwd_model.normalize_measurements= 0;
   end

   if figno>0
      if imgno == 0;
         imdl.jacobian_bkgnd.value = 1;
      else
         mdl_pts = interp_mesh( imdl.fwd_model, 5);
         out_circ= mdl_pts(:,1,:).^2 + mdl_pts(:,2,:).^2 > 0.75^2;
         out_circ= mean(out_circ, 3);
         sigma   = 2.2 - 1.2*out_circ;
         imdl.jacobian_bkgnd.value = sigma;
      end

      img = calc_jacobian_bkgnd( imdl );
      imdl= GREIT_mdl( img, weight, imdl.fwd_model.normalize_measurements );
if 1
   calc_noise_figure( set_tgts( imdl ));
%  fprintf('hp=%1.8g\n', choose_noise_figure( set_tgts(imdl, 0.79)) );
end
   end

      

   if figno>0
%  calc_noise_figure( set_tgts( imdl ));
%  fprintf('hp=%1.8g\n', choose_noise_figure( set_tgts(imdl, 0.79)) );
   end

   [eelv,eilv] = pig_data; vh = eelv.p0;  vi = eilv.p0; 
%  [eelv,eilv] = pig_data; vh = eelv.p20;  vi = eilv.p20; 
if 0
   switch imgno;
       case 0; vh = eelv.p0;  vi = eilv.p0; 
       case 1; vh = eelv.p20; vi = eilv.p20; 
       otherwise; error('huh?')
    end
end

   img= inv_solve(imdl, vh, vi);
   show_slices(img);

function imdl = GREIT_mdl( img, noiseampl, normalize );
   switch size(img.fwd_model.nodes,2)
      case 3; % 3D Model
         [xc,yc,zc]= ndgrid( linspace(-1,1,25), linspace(-1,1,25), [-.05,0,05] );
         elim = xc.^2 + yc.^2 > 0.9.^2; xc(elim)= []; yc(elim)= []; zc(elim)= []; 
         xyzr_pt = [xc(:),yc(:),zc(:), 0.05 + 0*zc(:)]';
      case 2; % 3D Model
         [xc,yc]= ndgrid( linspace(-1,1,75), linspace(-1,1,75) );
         elim = xc.^2 + yc.^2 > 0.9.^2; xc(elim)= []; yc(elim)= [];
         xyzr_pt = [xc(:),yc(:), 0.05 + 0*xc(:)]';
     otherwise
       error('huh?')
   end

   [vh,vi] = simulate_movement(img, xyzr_pt); 
   RM= calc_GREIT_RM(vh,vi,xyzr_pt,0.25, noiseampl, normalize);
   imdl= mk_common_gridmdl('b2c',RM);
   imdl.fwd_model.normalize_measurements = normalize;

 

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
      !7z e ng_mdl_16x1_coarse.7z ng_mdl_16x1_coarse.mat
   end

   load ng_mdl_16x1_coarse;
   f_mdl = ng_mdl_16x1_coarse;
   f_mdl.nodes = f_mdl.nodes(:,[2,1,3]); % flip x,y

   f_mdl.stimulation = c_mdl.stimulation;
   f_mdl.meas_select = c_mdl.meas_select;
   f_mdl.solve = @aa_fwd_solve;
   f_mdl.jacobian = @aa_calc_jacobian;
   f_mdl.system_mat = @aa_calc_system_mat;

   scl= 15;
   c_mdl.mk_coarse_fine_mapping.f2c_offset = [0,0,scl];
   c_mdl.mk_coarse_fine_mapping.f2c_project = (1/scl)*speye(3);
   c_mdl.mk_coarse_fine_mapping.z_depth = scl/5;
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

