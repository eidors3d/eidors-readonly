function [e res] = calc_solution_error(imgc, imdl, vh, vi)
%CALC_SOLUTION_ERROR Calculate residuals for a solution
% E = CALC_SOLUTION_ERROR(SOL, IMDL, DATA) calculates the residual error
% where 
%     SOL   -  EIDORS 'image' structure containing the solution
%     IMDL  -  EIDORS 'inv_model' structure used to calculate the solution
%     DATA  -  the data to be fitted (either a vector, or EIDORS 'data'
%              struct)
% and 
%     E = norm(DATA - SOL_MEAS) / norm( DATA ),
% where SOL_MEAS is the simulated data obtained from SOL.
%
% E = CALC_SOLUTION_ERROR(SOL, IMDL, VH, VI) allows specifying two data 
% inputs for difference solvers, such that DATA = VI - VH
%
% [E RES] = CALC_SOLUTION_ERROR(...) also returns the vector of residuals RES
%
% See also INV_SOLVE, CALC_DIFFERENCE_DATA

% (C) 2013 Bartlomiej Grychtol, License: GPL version 2 or 3.
% $Id$

% TODO: Generalize coarse2fine


switch imdl.reconst_type
   case 'difference'
      if nargin==4
         data = calc_difference_data(vh,vi,imdl.fwd_model);
      else
         if isstruct(vh)
            data = vh.meas;   % eidors object
         else
            data = vh;        % vector of measurements
         end
      end
      res = calc_diff_residual(imgc,imdl,data);
   case {'absolute', 'static'}
      if isstruct(vh)
         data = vh.meas;
      else
         data = vh;
      end
      res = calc_abs_residual(imgc,imdl,data);
   otherwise
      error('reconst_type not recognized');
end


% avarage error
e = norm(res)/norm(data);


function res = calc_abs_residual(imgc,imdl,data)
fmdl = imdl.fwd_model;

% make sure to have elem_data irrespective of physics
img = calc_jacobian_bkgnd(imdl); 
img = physics_data_mapper(img);
imgc = physics_data_mapper(imgc);

% account for coarse2fine
if size(img.elem_data,1) == size(imgc.elem_data,1)
   img.elem_data = imgc.elem_data;
else
   img.elem_data = fmdl.coarse2fine*imgc.elem_data;
end
img = physics_data_mapper(img,1);

% simualate data from solution
sim = fwd_solve(img);

% residuals
res = sim.meas - data;



function res = calc_diff_residual(sol,imdl,data)

% protect agains legacy movement image
sol = convert_legacy_movement_img(sol);

fmdl = imdl.fwd_model;
img = calc_jacobian_bkgnd(imdl);
if isfield(sol, 'movement') && ~isfield(img,'movement')
   try
      n_c2f = size(img.fwd_model.coarse2fine,2);
   catch
      n_c2f = 0;
   end
   n_elem = size(img.fwd_model.elems,1);
   n_elem_data = length(img.elem_data);
   if n_elem_data == n_elem || n_elem_data == n_c2f
      % jacobian background misses movement
      n_elec = length(img.fwd_model.electrode);
      n_dim  = mdl_dim(img.fwd_model);
      img.elem_data(end + (1:n_elec*n_dim),:) = 0;
   end
   img = convert_legacy_movement_img(img);  
end
simh = fwd_solve(img);

% sort out the size of any physics that becomes elem_data
img = physics_data_mapper(img);
sol = physics_data_mapper(sol);
if isfield(img, 'elem_data');
   try
      n_elems = size(sol.fwd_model.coarse2fine,2);
   catch
      n_elems = num_elems(sol.fwd_model);
   end
   if size(img.elem_data,1) ~= size(sol.elem_data,1);
      try
         sol.elem_data = fmdl.coarse2fine*sol.elem_data(1:n_elems,:);
      catch e
         eidors_msg(['dimensions of solution and jacobian background don''t agree' ...
            ' and using coarse2fine failed.']);
         rethrow(e);
      end
   end
end
sol = physics_data_mapper(sol,1);
img  = physics_data_mapper(img,1);

% add solution to jacobian background in param space
img  = physics_param_mapper(img);
sol = physics_param_mapper(sol);
initial = repmat(img.params,1,size(sol.params,2));

% add the solution params to the initial params from jacobian bkgnd
img.params = initial + sol.params;

% deal with arrays
jnk = img;
S = size(sol.params,2);
if S > 10
   fprintf(1, '   Error calculation progress:     ');
   step = ceil(S/100);
   for i = fliplr(1:S);
      if ~mod(i,step)
         fprintf(1, '\b\b\b%02d%%',100*(S-i)/S);
      end
      jnk.params = img.params(:,i);
      %    jnk = physics_param_mapper(jnk,1);
      tmp = fwd_solve(physics_param_mapper(jnk,1));
      simi(:,i) = tmp.meas;
   end
   fprintf(1, '\b\b\b100%%\n',100*(S-i/S));
end
sim = calc_difference_data(simh,simi,fmdl);

% residuals
res = data - sim;
