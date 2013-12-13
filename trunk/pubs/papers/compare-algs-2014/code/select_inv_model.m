function [imdl, ROI, mdlname] = select_inv_model(modelno,dir,W );
% call as:  select_inv_model( #, 'S12', W)
%           select_inv_model('paper-2D' | 'paper-3D');

   if isstr(modelno); switch modelno
      case 'paper-2D';       imdl = mk_paper_fem(2);
      case 'paper-2D-lungs'; imdl = mk_paper_fem(2.1);
      case 'paper-3D';       imdl = mk_paper_fem(3);
      case 'paper-3D-lungs'; imdl = mk_paper_fem(3.1);
      case 'test-NF';        test_NF(dir);
      otherwise;       error('huh?');
   end; return; end

   [imdl, ROI, mdlname] = get_model(modelno,dir,W );
   imdl = modify_if_necessary(imdl,dir);
   imdl = set_hyperparams(imdl, modelno);

function test_NF(mdl_list)
warning off EIDORS:Deprecated
warning off EIDORS:DeprecatedInterface
warning off EIDORS:calc_jacobian_input_params
   fmdl = mk_paper_fem(4); % dense cyl
   fmdl.stimulation = SBP_stim;
   img = mk_image(fmdl,1); % 3D case
if 1
   [vh,vi] = simulate_movement(img ,[0.3536;0.3536;0.5;0.05]);
else
         dd= eidors_readdata('S01/S01-001.get');
         [jnk, sel] = mk_stim_patterns(16,1,[0,1],[0,1],{},1);
         vh= double(dd(sel,180));
         vi= double(dd(sel,152));
end
   for mdl = mdl_list;
       if mdl == 1003
           bvh = vh; bvi = vi;
           fmdl = mk_paper_fem(4.1); % bkgnd
           fmdl.stimulation = SBP_stim;
           img = mk_image(fmdl,1);
           img.elem_data(fmdl.mat_idx{2}) = 0.3;
           [vh,vi] = simulate_movement(img ,[0.3536;0.3536;0.5;0.05]);
       end
      [imdl,jnk,name] = select_inv_model(mdl,'S12',[]);
      imdl.hyperparameter.tgt_data.meas_t1 = vh;
      imdl.hyperparameter.tgt_data.meas_t2 = vi;
      eidors_msg('log_level',0);
      if mdl == 1010
          NF = calc_noise_figure(imdl,NaN,1000);
      else
          NF = calc_noise_figure(imdl,imdl.hyperparameter.value,100);
      end
      fprintf('%s:\tHP=%1.8f\tNF=%f\n',name,imdl.hyperparameter.value,NF);
      subplot(4,4,mdl-999) 
      rimg = inv_solve(imdl,vh,vi);
      rimg.calc_colours.ref_level = 0;
      if size(rimg.fwd_model.nodes,2) == 3
          show_slices(rimg, [inf inf 0.5]);
      else
          show_slices(rimg)
      end
      title(sprintf('%s NF = %f',name,NF));
      eidors_msg('log_level',2);
      
      if mdl == 1003
          vi = bvi; vh = bvh;
      end
   end
warning on EIDORS:Deprecated
warning on EIDORS:DeprecatedInterface
warning on EIDORS:calc_jacobian_input_params



% returns the fmdl as imdl
function  str = lung_geometry(name, ctr, th1,l1,th2,l2);
   ll= 'solid %s= ellipticcylinder( %f,%f,0; %f,%f,0; %f,%f,0 );';
   p1 = [cos(pi*th1/180); sin(pi*th1/180)]*l1;
   p2 = [cos(pi*th2/180); sin(pi*th2/180)]*l2;
   str=sprintf(ll, name, ctr, p1, p2);
  

function  fmdl = mk_paper_fem(Dims);
%  lungs=['solid rL= ellipticcylinder( 0.3,0.2,0;0.4,0,0;0,0.5,0);' ...
%         'solid lL= ellipticcylinder(-0.3,0.2,0;0.4,0,0;0,0.5,0);' ...
%         'solid lungs = (rL or lL) and ob;'];
   lungs=[lung_geometry('rL',[ 0.3,0.1], 20,0.5, 110,0.6), ...
          lung_geometry('lL',[-0.3,0.1],-20,0.5,-110,0.6), ...
          'solid lungs = (rL or lL) and ob;'];

switch Dims; 
   case 2.0;
     [fmdl,mat_idx]= ng_mk_cyl_models([0,1,.10],[16],[0.05]);
   case 2.1;
     extra={'lungs', ...
           ['solid ob = orthobrick(-2,-2,0;2,2,.05);',lungs]};
     [fmdl,mat_idx]= ng_mk_cyl_models([0,1,.10],[16],[0.05], extra); 
   case 3.0;
     [fmdl,mat_idx]= ng_mk_cyl_models([1,1,.10],[16,.5],[0.05]); 
   case 3.1;
     extra={'lungs', ...
           ['solid ob = orthobrick(-2,-2,0;2,2,1);',lungs]};
     [fmdl,mat_idx]= ng_mk_cyl_models([1,1,.10],[16,.5],[0.05],extra); 
   case 4
     [fmdl mat_idx] = ng_mk_cyl_models([1,1,.05],[16,.5],[0.05 0 0.01]); 
   case 4.1;
     extra={'lungs', ...
           ['solid ob = orthobrick(-2,-2,0;2,2,1);',lungs]};
     [fmdl,mat_idx]= ng_mk_cyl_models([1,1,.05],[16,.5],[0.05 0 0.01],extra); 
   otherwise;       error('huh?');
end
fmdl.mat_idx = mat_idx;

%fmdl= mk_image(fmdl,1);
%fmdl.elem_data(mat_idx{2}) = 0.9;
%fmdl.calc_colours.ref_level = 1;
  
function hp = get_hp(modelno)
%---------------------- NF = 1.0 -----------------------------------------%
if 1
switch modelno
    case 1000
        hp = 0.0064;
    case 1001
%         hp = 0.00282; % HPF
        hp = 0.0066;
    case 1002
        hp = 0.0555;
    case 1003
%         hp = 0.0078;%0.0073; %HPF
        hp = 0.0127;
    case 1004
        hp = 0.310; % for m_hp = 1e5
    case 1005
        hp = 0.32;
    case 1006
        hp = 0.0000493;
    case 1007
%         hp = 0.00268;%diam_frac = 0.1
        hp = 0.00188;%diam_frac = 0.2
    case 1008
        hp = 0.8295;
    case 1009
        hp = 3.78;%4;
    case 1011
        hp = 5.6;
    case 1013
        hp = 0.00185;
    case 1014
        hp = 0.0525;
    otherwise
        error('No hp for %d',modelno);
end
end
%---------------------- NF = 0.5 -----------------------------------------%
if 0
switch modelno
    case 1000
        hp = 0.01615;
    case 1001
        hp = 0.0075;
    case 1002
        hp = 0.155;
    case 1003
        hp = 0.01285;
    case 1004
        hp = 0.86;
    case 1005
        hp = 0.9;
    case 1006
        hp = 0.0000268;
    case 1007
        hp = 0.00715; 
    case 1008
        hp = 1.1;
    case 1009
        hp = 12.85;
    case 1011
        hp = 11.8;
    otherwise
        error('No hp for %d',modelno);
end
end

end
% Sheffield stim patterns used
function [stim,meas_sel] = SBP_stim;
   [stim,meas_sel] =  ...
       mk_stim_patterns(16,1,0:1,0:1,{'no_meas_current'},1);

function [imdl, ROI, mdlname] = get_model(modelno,fname,W );
   fmdl = mk_paper_fem(2);
   fmdl.stimulation = SBP_stim;
   fmdl.normalize_measurements = 1;
   iR0 = select_imdl(fmdl,{'Basic GN dif'});
   iR0.jacobian_bkgnd.value = 1;
   iR0.RtR_prior = @noser_image_prior;
   ROI = (interp_mesh(fmdl)*[0;1]) <0;
   switch modelno
% 1010  RSBP = Sheffield backprojection
      case 1010; mdlname= 'RSBP';
         imdl= mk_common_gridmdl('b2c','backproj');
         imdl.fwd_model.stimulation = SBP_stim;
         ROI= [zeros(1,856/2),ones(1,856/2)]; % lower half
         imdl.fwd_model.nodes = imdl.fwd_model.nodes;%*[1,0;0,-1];
         imdl.hyperparameter.value = NaN;
% 1000  R0   = 2D, norm, homog, nomove, noise=1, L2, NOSER
      case 1000; mdlname= 'R0';
         imdl = iR0;
         imdl.hyperparameter.value = get_hp(modelno);

% 1002  Rdif =     nonorm
      case 1002; mdlname = 'Rdif';
         imdl = iR0;
         imdl.fwd_model.normalize_measurements = 0;
         imdl.hyperparameter.value = get_hp(modelno);

% 1007  RLPF                                         Gaussian
      case 1007; mdlname = 'Rlpf';
         imdl = iR0;
         imdl.RtR_prior = @prior_gaussian_HPF;
         imdl.fwd_model.prior_gaussian_HPF.diam_frac = 0.2;
         imdl.hyperparameter.value = get_hp(modelno);
         
% 1013  RLP                                          Laplace
      case 1012; mdlname = 'Rlp';
         imdl = iR0;
         imdl.RtR_prior = @prior_laplace;
         imdl.hyperparameter.value = get_hp(modelno);

% 1004  Rmv                      move
      case 1004; mdlname = 'Rmv';
         m_hp = 10000; hp = get_hp(modelno);
m_hp = 1e+5;
         imdl = iR0; img = mk_image(imdl);
         Jc = calc_jacobian( img );
         iRtR = inv(noser_image_prior( imdl ));
         iRN = hp^2 * speye(size(Jc,1));
         Jm = calc_Jm( img );
         RM = iRtR*Jc'/(Jc*iRtR*Jc' + m_hp*Jm*Jm' + iRN);
         imdl.solve = @solve_use_matrix; 
         imdl.solve_use_matrix.RM  = RM;
         imdl.hyperparameter.value = hp; % just for info


% 1005  Rn        reciproc 
      case 1005; mdlname= 'Rn';
         imdl = iR0;
         imdl.hyperparameter.value = get_hp(modelno);
         % get data 
         [jnk,msel] = SBP_stim;
         ff=dir([fname,'*.mat']); vv=[];
         for fi = 1:length(ff); 
           if ff(fi).name(9) == 'R', continue, end % Don't look at saved image files
           % TODO continue if we don't have dsv
           dd= load(ff(fi).name); vv=[vv,dd.dsv(msel,:)];
         end
         imdl.calc_reciproc_error.tau = 1e-5;
         % for this value for S12, mean(diag(imdl.meas_icov))=0.6604
         imdl.meas_icov = calc_reciproc_error( imdl, vv );
         fprintf('Mean meas_icov = %f\n', ...
             full(mean(diag(imdl.meas_icov))) );

% 1006  RTV       TV
      case 1006; mdlname= 'RTV';
         imdl= rmfield(iR0,'RtR_prior');
         imdl.R_prior= @ab_calc_tv_prior;
         imdl.fwd_model.normalize_measurements= 1;
         imdl.solve = @pdipm_diff;
         imdl.hyperparameter.value = get_hp(modelno);
         imdl.pdipm_diff.beta = 1e-10;
         imdl.pdipm_diff.norm_data =2;
         imdl.pdipm_diff.norm_image=1;
%        imdl.pdipm_diff.beta     (default 1e-6)
         imdl.parameters.max_iterations = 15;

% 1008  RL1       L1n
      case 1008; mdlname= 'RL1';
         imdl= iR0;
         imdl.fwd_model.normalize_measurements= 1;
         imdl.solve = @pdipm_diff;
         imdl.hyperparameter.value = get_hp(modelno);
         imdl.pdipm_diff.beta = 1e-10;
         imdl.pdipm_diff.norm_data =1;
         imdl.pdipm_diff.norm_image=2;
         imdl.parameters.max_iterations = 10;

% 1001  R3D  = 3D
      case 1001; mdlname = 'R3D';
         fmdl = mk_paper_fem(3);
%        fmdl.nodes(:,1) = -fmdl.nodes(:,1); % reverse x-ax
%  no longer needed
         fmdl.stimulation = SBP_stim;
         fmdl.normalize_measurements = 1;
         imdl = select_imdl(fmdl,{'Basic GN dif'});
         imdl.rec_model = mk_paper_fem(2);
         c2f= mk_coarse_fine_mapping( fmdl, imdl.rec_model );
         imdl.fwd_model.coarse2fine = c2f;
         imdl.RtR_prior = @noser_image_prior; 
         imdl.prior_use_fwd_not_rec = 1;
       %-- singular
%          imdl.RtR_prior = @gaussian_HPF_prior;
         imdl.hyperparameter.value = get_hp(modelno);
         ROI = (interp_mesh(imdl.rec_model)*[0;1]) <0;
% 1014 R3Ddiff         
      case 1014; mdlname = 'R3Ddiff';
         [imdl, ROI] = get_model(1001,fname);
         imdl.fwd_model.normalize_measurements = 0;
         imdl.hyperparameter.value = get_hp(modelno);

% 1003  Rbkg            nohomg
      case 1003; mdlname = 'Rbkg';
         fmdl = mk_paper_fem(3.1);
%        fmdl.nodes(:,1) = -fmdl.nodes(:,1); %reverse x-ax
         fmdl.stimulation = SBP_stim;
         fmdl.normalize_measurements = 1;
         imdl = select_imdl(fmdl,{'Basic GN dif'});
         imdl.jacobian_bkgnd.value = ones(num_elems(fmdl),1);
         imdl.jacobian_bkgnd.value(fmdl.mat_idx{2}) = 0.3;
         imdl.rec_model = mk_paper_fem(2);
         c2f= mk_coarse_fine_mapping( fmdl, imdl.rec_model );
         imdl.fwd_model.coarse2fine = c2f;
         imdl.RtR_prior = @noser_image_prior; %-- singular
         imdl.prior_use_fwd_not_rec = 1;
%          imdl.RtR_prior = @gaussian_HPF_prior;
 % Not real, but NF makes less sense here;
         imdl.hyperparameter.value = get_hp(modelno);
         ROI = (interp_mesh(imdl.rec_model)*[0;1]) <0;

% 1009  RGREIT
      case 1009; mdlname = 'RGR';
         fmdl = mk_paper_fem(3);
         fmdl.stimulation = SBP_stim;
         fmdl.normalize_measurements = 1;
         img = mk_image(fmdl,1);
         opt.imgsz = [32 32];
         opt.distr = 3; % non-random, uniform
         opt.Nsim = 1000;
         opt.target_size = 0.05; % Target size (frac of medium)
%          opt.noise_figure = 0.75; % Recommended NF=0.5;
         imdl = mk_GREIT_model(img, 0.25, get_hp(modelno), opt);
         ROI = (interp_mesh(imdl.rec_model)*[0;1]) <0;
         imdl.hyperparameter.value = get_hp(modelno); %only for info
% 1012 RGREIT Background         
      case 1012; mdlname = 'RGRb';
         fmdl = mk_paper_fem(3.1);
         fmdl.stimulation = SBP_stim;
         fmdl.normalize_measurements = 1;
         img = mk_image(fmdl,1);
         img.elem_data(fmdl.mat_idx{2}) = 0.3;
         opt.imgsz = [32 32];
         opt.distr = 3; % non-random, uniform
         opt.Nsim = 1000;
         opt.target_size = 0.05; % Target size (frac of medium)
         opt.noise_figure = 0.79; % Recommended NF=0.5;
         imdl = mk_GREIT_model(img, 0.25, [], opt);
         ROI = (interp_mesh(imdl.rec_model)*[0;1]) <0;
% 1011 TSVD
       case 1011; mdlname = 'RTSVD';
          fmdl = mk_paper_fem(2);
          fmdl.stimulation = SBP_stim;
          fmdl.normalize_measurements = 1;
          imdl = select_imdl(fmdl,{'Basic GN dif'});
          imdl.solve = @inv_solve_TSVD;
          imdl.hyperparameter.value = get_hp(modelno);


      otherwise 
         error(['dont know what to do with model #',modelno]);
   end

function imdl = modify_if_necessary(imdl,dir)
   switch dir  % special handling of any directories
% S07-S12 Electrodes attached in a reversed order (i.e. electrode 5 on the right chest, electrode 13 on the left chest):  
      case {'SS1', 'SS2'}
         % Simulation - do nothing

      case {'S07','S08','S09','S10','S11','S12'}
          if isfield(imdl,'rec_model');
              imdl.rec_model.nodes(:,1)= -imdl.rec_model.nodes(:,1);
          else
              imdl.fwd_model.nodes(:,1)= -imdl.fwd_model.nodes(:,1);
          end
      case {'S05','S06','S13','S14','S15','S16','S17','S18', ...
            'S19','S20','S21','S22'}
         % model is ok. Do nothing

      otherwise,
         error('shouldn`t get here');
   end

function ROI= c2c2ROI(imdl);
   xyel = interp_mesh(imdl.fwd_model);
   r_el = sqrt(sum(xyel.^2,2));
   ROI   = (xyel(:,2)<0) & (r_el < max(r_el)*.8);
   imgroi= eidors_obj('image','','elem_data',ROI,'fwd_model',imdl.fwd_model);
   pp= aa_fwd_parameters(imdl.fwd_model);
   ROI= (ROI.*pp.VOLUME)';

function img= inv_solve_sim_norm( imdl, data1, data2)
   img= eidors_obj('image','','fwd_model',imdl.fwd_model);
   vh= fwd_solve(calc_jacobian_bkgnd(imdl));
   d2 = data2./data1.*vh.meas; 
   d1 = vh.meas;
   img= aa_inv_solve( imdl, d1, d2);

function Jm = calc_Jm( img);
   cache_obj = img; cache_obj.type='cache-image';
   Jm = eidors_obj('get-cache',cache_obj,'calc_Jm');
   if ~isempty(Jm); return; end
eidors_msg('recalc Jm',1);

   Jm= []; delta= 1e-6;

   d0= fwd_solve( img );
   node0 = img.fwd_model.nodes;
   for d= 1:size(node0,2);
      for i= 1:size(node0,1);
         img.fwd_model.nodes( i, d)= node0(i,d) + delta;
         di= fwd_solve( img );
         img.fwd_model.nodes( i, d)= node0(i,d);

         Jm= [Jm, (1/delta) * (di.meas - d0.meas)];
      end
   end
   eidors_obj('set-cache', cache_obj, 'calc_Jm', Jm);