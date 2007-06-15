function [data,img] = dg_fwd_solve(img)
%% Solve the forward problem with a mapping function
%
% This function first applies a mapping function which gives the
% conductivities in each fine elements of the FEM structure from a set of
% parameters and auxiliary data:
%
% img.params_mapping.function = handle pointing onto the mapping function
% img.params_mapping.params = the parameters from which the conductivites
%                             are obtained
% imp.params_mapping.data = auxiliary data needed by the mapping function
%
% Dominique Gibert, April 2007
%
%%
mapping_function= img.params_mapping.function;
if ~isempty(mapping_function)
    img= feval(mapping_function,img);
end
fwd_model= img.fwd_model;
pp= dg_fwd_parameters(fwd_model);
s_mat= dg_calc_system_mat(fwd_model,img);
idx= 1:pp.n_node;
idx(fwd_model.gnd_node) = [];
v= zeros(pp.n_node,pp.n_stim);
v(idx,:)= s_mat(idx,idx)\pp.QQ(idx,:);
%vv= v(MES(1,:),:)- v(MES(2,:),:);
%vv= vv(ELS);

% calc voltage on electrodes
v_els= pp.N2E * v;

% measured voltages from v
vv = zeros(pp.n_meas,1);
idx=0;
for i=1:pp.n_stim
   meas_pat= fwd_model.stimulation(i).meas_pattern;
   n_meas  = size(meas_pat,1);
   vv(idx+(1:n_meas)) = meas_pat*v_els(:,i);
   idx= idx+ n_meas;
end

% create a data structure to return
data.meas= vv;
data.time= -1; % unknown
data.name= 'solved by dg_fwd_solve';
% TODO: figure out how to describe measurement pattern
data.configuration='unknown';
