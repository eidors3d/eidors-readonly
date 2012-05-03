function img= mc_e_move_inv_solve( inv_model, data1, data2)
% MC_INV_SOLVE inverse solver 
% img= mc_e_move_inv_solve( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%
% both data1 and data2 may be matrices (MxT) each of
%  M measurements at T times
% if either data1 or data2 is a vector, then it is expanded
%  to be the same size matrix
%
%NB : Jacobian abd prior calculation for the movement are taken from Andy's
%code (aa_e_move_jacobian, aa_e_move_image_prior). I have just put the code
%here as a wrapper for soem scripts I am running

%Calculate the voltage difference data
diff_volt = calc_difference_data( data1, data2, inv_model.fwd_model);

%Get the change in conductivity from generalised tikhonov regularisation
sol = tikhonov_one_step_inv(inv_model) * diff_volt;

%Create a data structure to return (difference conductivity!!!!)
img.name= 'solved by mc_GN_solve';
img.elem_data = sol;
img.fwd_model= inv_model.fwd_model;
img.type='image';

function tikhonov_inv = tikhonov_one_step_inv(inv_model)
% The one_step reconstruction matrix is cached
   tikhonov_inv = eidors_obj('get-cache', inv_model, 'mc_e_move_one_step_inv');
   if ~isempty(tikhonov_inv)
       eidors_msg('mc_e_move_inv_solve: using cached value', 3);
   else     
       %Calculate the background for the Jacobian
       img_bkgnd= calc_jacobian_bkgnd( inv_model );
       
       %Calculate conductivity and movement Jacobian separately
       Jc = calc_jacobian( inv_model.fwd_model, img_bkgnd);
       Jm = perturb_movement_jacobian(img_bkgnd); 
       J=[Jc,Jm]; %Concatenate
       
       %Calculate conductivity and movement image prior separately
       Rc = calc_R_prior( inv_model );
       if(strcmp(inv_model.movement_prior,'@tikhonov_movement_prior'))
           Rm = tikhonov_movement_prior(img_bkgnd);
       elseif(strcmp(inv_model.movement_prior,'@laplace_movement_prior'))
           Rm = laplace_movement_prior(img_bkgnd);
       end
       
       %Concatenate into matrix (hp_move^2=ratio of cond/elec variance)  
       n_cond=size(Rc,1); n_extra=size(Rm,1); R0=zeros(n_cond,n_extra);
       hp_move=inv_model.cond_movement_param;
       R= [Rc, R0; R0', hp_move^2*Rm ]; 
           
       %Calculate total hyperparameter and compute inverse
       hp= calc_hyperparameter( inv_model );
       tikhonov_inv= (J'*J +  hp^2*(R'*R))\J';
       
       %Create an eidors object for iterations
       eidors_obj('set-cache', inv_model, 'mc_e_move_one_step_inv', tikhonov_inv);
       eidors_msg('mc_e_move_inv_solve: setting cached value', 3);
   end
end
  
end

function J= perturb_movement_jacobian(img )
delta=1e-6;
%Find electrode stucture and no.of electrodes 
elecstruc=img.fwd_model.electrode; nelecs=size(elecstruc,2);
stimstruc=img.fwd_model.stimulation; nstims=size(stimstruc,2); 
nodestruc=img.fwd_model.nodes; ndims=size(nodestruc,2); 

%Find total number of measurements
nmeass=0;
for k=1:nstims
    stimkmeasmatrix = stimstruc(k).meas_pattern;
    nmeass=nmeass+size(stimkmeasmatrix,1);
end

J = zeros( nmeass, nelecs*ndims );

d0= fwd_solve( img );
for d= 1:ndims
   for i= 1:nelecs
      idx= img.fwd_model.electrode(i).nodes;

      img.fwd_model.nodes( idx, d)= nodestruc(idx,d) + delta;
      di= fwd_solve( img );
      img.fwd_model.nodes( idx, d)= nodestruc(idx,d);

      J_idx = nelecs*(d-1) + i;
      J(:,J_idx) = (1/delta) * (di.meas - d0.meas);
   end
end

end

function RegM = tikhonov_movement_prior(img )
elecstruc=img.fwd_model.electrode; nelecs=size(elecstruc,2);
nodestruc=img.fwd_model.nodes; ndims=size(nodestruc,2); 
RegM=eye(nelecs*ndims);
end

function RegM = laplace_movement_prior(img )
elecstruc=img.fwd_model.electrode; nelecs=size(elecstruc,2);
node0=img.fwd_model.nodes; ndims=size(node0,2); 
% movement constraint in each dimention
idx =(0:nelecs-1)';
im1= rem(idx-1+nelecs,nelecs);
ip1= rem(idx+1,nelecs); 
mv= sparse([idx,idx,idx]+1,[im1,idx,ip1]+1,ones(nelecs,1)*[-1,2.1,-1]);

RegM= spalloc(ndims*nelecs,ndims*nelecs, 3*ndims*nelecs);

for i=0:ndims-1;
    m_idx= idx + i*nelecs + 1;
    RegM( m_idx, m_idx ) = mv;
end

end