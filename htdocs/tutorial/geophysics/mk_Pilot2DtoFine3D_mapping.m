function img= mk_Pilot2DtoFine3D_mapping(img,sparsify)
%% Maps the fine elements of a 3D FEM structure onto a set of sparse points
% of the corresponding 2D FEM structure (initially made for gallery structure)
%
% The fine-elements vector is obtained by applying a matrix onto the vector
% of coarse elements
%
% First step: computes the matrix mapping fine 3D onto fine 2D elements
%
% Second step: computes the barycenters of each 2D element and select a
% subset of them (1 every sparsify point)
%
% Third step: construct a matrix to map the 2D elements onto the 3D elements
% of the FEM structure
%
% Fourth step: construct a matrix mapping the selected (sparse) 2D barycenters
% onto the dense barycenters (hence in the fine 2D elements)
%
% Dominique Gibert, April 2007
%
%%

n_fine_elems3D= size(img.fwd_model.elems,1);
[m2Dto3D,nodes2D,elems2D]= mk_2Dto3D_mapping(img.fwd_model);
n_fine_elems2D= size(elems2D,1);
for k= 1:n_fine_elems2D
    fine_bary2D(k,:)= mean(nodes2D(elems2D(k,:),:));
end
sparse_bary2D= fine_bary2D(1:sparsify:n_fine_elems2D,:);
n_sparse_bary2D = size(sparse_bary2D,1);
ms2f= sparse2Dtodense2D_mapping(sparse_bary2D,fine_bary2D);
c2f= m2Dto3D*ms2f;
if isempty(img.elem_data)
    img.elem_data= ones(n_fine_elems3D,1);
end
img.params_mapping.params= log(img.elem_data(1:sparsify:n_fine_elems2D));
img.params_mapping.function = @dg_sparse2Dto3D_mapping;
img.params_mapping.data = c2f;
img.params_mapping.perturb= log(2)*ones(size(img.params_mapping.params,1),1);
disp(['3D FEM model of gallery created | number of elements = ' num2str(n_fine_elems3D)]);
disp(['                 Number of parameters to be inverted = ' num2str(n_sparse_bary2D)]);

function ms2d= sparse2Dtodense2D_mapping(sparseXY,denseXY)
%% merge = n -> coarse elements are obtained by merging n contiguous fine elements
% merge = [m1 m2 ... mN] -> coarse element number 1 is obtained by merging
% the m1 first fine elements, etc.. If sum(merge) < total number of fine
% elements, the remaining fine elements are merged to produce a N+1th coarse
% element.
n_keep_sparse=7;
n_sparse= size(sparseXY,1);
n_dense= size(denseXY,1);
sX= sparseXY(:,1);
sY= sparseXY(:,2);
dX= denseXY(:,1);
dY= denseXY(:,2);
ms2d=zeros(n_dense,n_sparse);
n_neighbor= min(n_keep_sparse,n_sparse);
for k= 1:n_dense
    dist= sqrt((sX(:)-dX(k)).^2 + (sY(:)-dY(k)).^2);
    [d,kdist]= sort(dist);
    if dist(kdist(1)) == 0
        dist(kdist(1))= dist(kdist(2))/2;
    end
    dist = 1./dist;
    dist_kdist = dist(kdist(1:n_neighbor));
    ms2d(k,kdist(1:n_neighbor)) = dist_kdist(:)'/sum(dist_kdist);
end
ms2d= sparse(ms2d);

function img = dg_sparse2Dto3D_mapping(img)
%% Function to be called to perform the mapping in the forward problem
mapping_matrix= img.params_mapping.data;
coarse_elems= img.params_mapping.params;
img.elem_data= exp(mapping_matrix*coarse_elems); % take exp because parameters are in log domain
