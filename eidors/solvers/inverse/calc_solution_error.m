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
img = physics_data_mapper(img);

% simualate data from solution
sim = fwd_solve(img);

% residuals
res = sim.meas - data;



function res = calc_diff_residual(imgc,imdl,data)

fmdl = imdl.fwd_model;
img = calc_jacobian_bkgnd(imdl);
simh = fwd_solve(img);

% map physics to elem_data
img = physics_data_mapper(img);
imgc = physics_data_mapper(imgc);

if ~isfield(imgc,'elem_data') && isfield(imgc,'node_data')
   eidors_msg('Solution error calculation for nodal solvers not supported (yet).',2);
   res = NaN;
   return
end
   
% add solution to jacobian background
e_data = repmat(img.elem_data,1,size(imgc.elem_data,2));
has_c2f = isfield(imgc.fwd_model,'coarse2fine');
if ~has_c2f
   n_elems = num_elems(imgc.fwd_model);
   img.elem_data = e_data + imgc.elem_data(1:n_elems,:);
else
   n_elems = size(imgc.fwd_model.coarse2fine,2);
   img.elem_data = e_data + fmdl.coarse2fine*imgc.elem_data(1:n_elems,:);
end
img = physics_data_mapper(img,1);

% simulate data from solution
jnk = img;
for i = fliplr(1:size(imgc.elem_data,2));
   jnk.elem_data = img.elem_data(:,i);
   tmp = fwd_solve(jnk);
   simi(:,i) = tmp.meas;
end
sim = calc_difference_data(simh,simi,fmdl);

% residuals
res = data - sim;
