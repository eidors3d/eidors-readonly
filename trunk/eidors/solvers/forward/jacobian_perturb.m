function J= jacobian_perturb( fwd_model, img)
% JACOBIAN_PERTURB: J= jacobian_perturb( fwd_model, img)
% Calculate Jacobian Matrix, based on small perturbations
%   in the forward model. This will tend to be slow, but
%   should be best used to 'sanity check' other code
%
% J         = Jacobian matrix
% fwd_model = forward model
% fwd_model.jacobian_perturb.delta   - delta perturbation to use
% img = image background for jacobian calc

% (C) 2006 Andy Adler. License: GPL version 2 or version 3
% $Id$

if ischar(fwd_model) && strcmp(fwd_model,'UNIT_TEST'); do_unit_test; return;end

if isfield(fwd_model,'jacobian_perturb')
   delta = fwd_model.jacobian_perturb.delta;
else
   delta= 1e-6; % tests indicate this is a good value
end

n_elem = size(fwd_model.elems,1);
% force image to use provided fwd_model
img.fwd_model= fwd_model;

jnk = data_mapper(img);
curprms = jnk.current_params;

if strcmp(curprms, 'conductivity') && ~isfield(img,'conductivity')
   curprms = '';
end
% solve one time to get the size
d0= fwd_solve( img );

if isfield(img.fwd_model,'coarse2fine');
   Jcol= perturb_c2f(img, 1, delta, d0, curprms);
   Jrows= size(img.fwd_model.coarse2fine,2);
   J= zeros(length(Jcol), Jrows );
   J(:,1)= Jcol;
   for i=2:size(img.fwd_model.coarse2fine,2);
     J(:,i)= perturb_c2f(img, i, delta, d0, curprms);
     if rem(i,50)==0; fprintf('+'); end
   end
else
   Jcol= perturb(img, 1, delta, d0, curprms);

   J= zeros(length(Jcol), n_elem);
   J(:,1)= Jcol;

   for i=2:n_elem
     J(:,i)= perturb(img, i, delta, d0, curprms);
   end
end


function Jcol= perturb( img, i, delta, d0, curprms)
   if ~isempty(curprms)
      img.(curprms).elem_data(i)= img.(curprms).elem_data(i) + delta;
   else
      img.elem_data(i)= img.elem_data(i) + delta;
   end
   di= fwd_solve( img );
   Jcol = (1/delta) * (di.meas - d0.meas);

function Jcol= perturb_c2f( img, i, delta, d0, curprms)
   if ~isempty(curprms)
      img.(curprms).elem_data= img.(curprms).elem_data + delta*img.fwd_model.coarse2fine(:,i);
   else
       img.elem_data(i)= img.elem_data(i) + delta; %*img.fwd_model.coarse2fine(:,i);
%       img.elem_data= img.fwd_model.coarse2fine*img.elem_data + delta*img.fwd_model.coarse2fine(:,i);
   end
   di= fwd_solve( img );
   Jcol = (1/delta) * (di.meas - d0.meas);


function do_unit_test
% Perturbation Jacobians
% $Id$

  models ={'c2c2','a3cr','n3r2'};
  for i = 1:length(models)
     imdl= mk_common_model(models{i},16);
     imdl.fwd_model.nodes = imdl.fwd_model.nodes*.25;
     img= calc_jacobian_bkgnd(imdl);

     img.fwd_model = mdl_normalize(img.fwd_model, 0);

     img.fwd_model.jacobian=   @jacobian_adjoint;
     img.fwd_model.system_mat= @system_mat_1st_order;
     img.fwd_model.solve=      @fwd_solve_1st_order;
     J_aa= calc_jacobian( img ); % 2 for bug in my code

     img.fwd_model.jacobian=   @jacobian_perturb;
     J_aa_p= calc_jacobian( img ); % 2 for bug in my code

     unit_test_cmp('J_ - J_p', norm(J_aa - J_aa_p,'fro')/norm(J_aa,'fro'), 0, 1e-3);
   end


function display_jacobian_deltas
% Calculate the perturbation Jacobian
  vh= fwd_solve(img);

  for el= 1:50:size(img.fwd_model.elems,1);
     fprintf('\nelem#%03d: ',el);
     for delta= [1e-4,1e-6,1e-8] % delta is perturb
        img_i = img;
        img_i.elem_data(el)= img_i.elem_data(el) + delta;
        vi= fwd_solve(img_i);
        J_p = (vi.meas - vh.meas)/delta; % perturb Jacobian
        fprintf(' %8.6f', norm(J_p - J_np(:,el))/norm(J_np(:,el)));
     end
  end

