function img= inv_solve_abs_GN_diff_tangential( inv_model, data0);
%INV_SOLVE_ABS_GN_TANGENTIAL Absolute solver using GN
% img= inv_solve_abs_GN( inv_model, data0)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data0      => EIT data

%Pass options on inverse model
opt = parse_options(inv_model);

%Passed in Jacobian Adjoint with a c2f on inv_model
img = initial_estimate( inv_model, data0 );

%Calculate hyperparameter and RtR (movement prior!!!)
hp  = calc_hyperparameter( inv_model );
RtR = calc_RtR_prior( inv_model ); hp2RtR= hp^2*RtR;

% Physics and initial estimates
img0 = physics_data_mapper(img);
img0.elem_data=inv_model.prior_c;
img.elem_data=inv_model.prior_c;

%get the number of electrodes and coarse elements
n_elec=length(inv_model.img_fine.fwd_model.electrode);
n_elemsc=inv_model.fwd_model.n_elemsc;

%We want to get the prior conductivity and tangential change 
prior_ce=img0.elem_data; 
prior_ce(n_elemsc+1:n_elemsc+2*n_elec)=inv_model.prior_e;
cur_ce=prior_ce; 
cur_move=prior_ce(n_elemsc+1:n_elemsc+2*n_elec);

%Get the fine image and coarse model 
img_h=inv_model.img_fine;
mdl_c=inv_model.mdl_coarse;

%Calculate the electrode components on fine model
elec_comp_h=calc_electrode_components(img_h.fwd_model);
for jj=1:length(elec_comp_h)
    elec_posHS(jj,:)=elec_comp_h{jj}.com; 
end
elec_posH = elec_posHS;

%Compute the dat and regularisaiton residuals
vsim=fwd_solve(img_h);
data_res=zeros(opt.max_iter+1,1); 
data_res(1) = 0.5*(vsim.meas-data0)'*(vsim.meas-data0);
regu_res=zeros(opt.max_iter+1,1); 
regu_res(1) = 0;

%Movement and conductivity Jacobian
Jm=jac_move_eidors_elec_tang_only(inv_model.fwd_model,img_h); 
J = calc_jacobian( img ); J=[J,Jm]; resolve_jacobian=0;

%Loop over iterations
for i = 1:opt.max_iter
  %Simulate measured data on FINE image and calc diff data
  vsim = fwd_solve( img_h );
  vsim.meas=vsim.meas+Jm*cur_move;
  dv = calc_difference_data( vsim , data0, img.fwd_model);  

  %Jacobian here  
  if(resolve_jacobian==1)
    %Jacobian adjoint and coarse2fine field to get movement
    J = calc_jacobian( img ); J=[J,Jm];
  end
 
  %Compute the search direction -inv(H)*J'
  RDx = hp2RtR*(prior_ce - cur_ce);
  dx = (J'*J + hp2RtR)\(J'*dv + RDx); 
  
  %Perform a linesearch on this using polynomial search
  opt.line_optimize.hp2RtR = hp2RtR;
  img.elem_data=cur_ce; %Elem_data with conductivity AND move
  img0.elem_data=prior_ce;
  img = feval(opt.line_optimize_func,img, dx, data0, ...
      img_h, inv_model.c2f2,n_elemsc,Jm,img0,opt.line_optimize);

  %Reassign tangential and conductivity (fine already sorted out)
  cur_ce = img.elem_data;
  cur_move=cur_ce(n_elemsc+1:n_elemsc+2*n_elec);
  img.elem_data = cur_ce(1:n_elemsc);

  %Add coarse data convert to fine differential as a test
  fine_from_c = inv_model.c2f2*img.elem_data;
  img_h=mk_image(inv_model.fwd_model,fine_from_c); 
  
  %Calculate the electrode components
  for jj=1:length(elec_comp_h)
      %We can then decompose this to update the direction
      a_i_elec_ii = cur_move(jj)*elec_comp_h{jj}.tangent(:,1) + ...
          cur_move(jj+n_elec)*elec_comp_h{jj}.tangent(:,2);    
      %Calculate the new coordinates
      elec_pos_NEW(jj,1:3) = elec_posHS(jj,1:3) + a_i_elec_ii';   
  end
  %Get the normals
  elec_pos_NEW(:,4:6)=elec_pos_NEW(:,1:3);

  %Compute the data residual
  vsim=fwd_solve(img_h); %This is conductivity change
  vsim.meas=vsim.meas+Jm*cur_move; 
  
  %The new residuals and message
  data_res(i+1) = 0.5*(vsim.meas-data0)'*(vsim.meas-data0);
  regu_res(i+1) = 0.5*(cur_ce-prior_ce)'*hp2RtR*(cur_ce-prior_ce);   
  eidors_msg('#%02d squared summed residual=%.3g', i, ...
      data_res(i+1)+regu_res(i+1), 1);    
end

%Attach the movement
img.movement_data = cur_move; %this is on tangents

function val = GN_objective_function(data0, data, img_k,img0,opt)
   dv = calc_difference_data(data, data0, img0.fwd_model);
   de = img_k.elem_data - img0.elem_data;
   val = 0.5*( dv'*dv + de' * opt.hp2RtR * de);   


function img = initial_estimate( imdl, data )
   img = calc_jacobian_bkgnd( imdl );
   vs = fwd_solve(img);
   pf = polyfit(data,vs.meas,1);
   img = physics_data_mapper(img);
   if isfield(img.fwd_model,'coarse2fine');
      nc = size(img.fwd_model.coarse2fine,2);
      img.elem_data = mean(img.elem_data)*ones(nc,1)*pf(1);
   else
      img.elem_data = img.elem_data*pf(1);
   end
   % remove elem_data
   img = physics_data_mapper(img,1);  
   
function opt = parse_options(imdl)
   try
       % for any general options
       opt = imdl.parameters;
   end
   opt.max_iter = 1;
   try
      opt.max_iter = imdl.parameters.max_iterations;
   end
   
   if isfield(imdl, 'inv_solve_abs_GN');
      fnames = fieldnames(imdl.inv_solve_abs_GN);
      for i = 1:length(fnames)
          opt.(fnames{i}) = imdl.inv_solve_abs_GN.(fnames{i});
      end
   end
   
   if ~isfield(opt,'line_optimize_func')
      opt.line_optimize_func = @line_optimize_tangential;
   end
   
   if ~isfield(opt,'line_optimize')
      opt.line_optimize = [];
   end
   
   if ~isfield(opt, 'line_optimize') || ...
      ~isfield(opt.line_optimize, 'objective_func')
    % not sure this should be allowed to change
      opt.line_optimize.objective_func = @GN_objective_function;       
   end
   
   if ~isfield(opt,'do_starting_estimate')
       opt.do_starting_estimate = 1;
   end