function [m2Dto3D,nodes2D,elems2D]= mk_2Dto3D_mapping(fwd_model)
%% Makes the correspondance between the elements of a 2D FEM structure and
%  its extruded 3D FEM structure (initially made for FEM gallery model)
%
% m2Dto3D = matrix such that elems_3D = m2Dto3D*elems_2D
%
% Dominique Gibert, April 2007
%
%%
nodes2D= fwd_model.nodes(fwd_model.nodes(:,3)==0,1:2);
n_nodes2D= size(nodes2D,1);
elems2D= fwd_model.elems(:,2:4);
elems2D= elems2D(elems2D(:,1)<=n_nodes2D,:);
elems2D= elems2D(elems2D(:,2)<=n_nodes2D,:);
elems2D= elems2D(elems2D(:,3)<=n_nodes2D,:);
n_elems2D= size(elems2D,1);
n_elems3D= size(fwd_model.elems,1);
m2Dto3D= zeros(n_elems3D,n_elems2D);
for k= 1:n_elems2D
    m2Dto3D(fwd_model.misc.map_2Dto3D(k,:),k)= 1.;
end
m2Dto3D= sparse(m2Dto3D);
