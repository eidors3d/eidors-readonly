function img= mk_Coarse2DtoFine3D_mapping(img,merge_elems)
%% Maps the fine elements of a 3D FEM structure onto coarse elements of the
%% corresponding 2D FEM structure (initially made for gallery structure)
%
% The fine elements vector is obtained by applying a matrix onto the vector
% of coarse elements
%
% First step: computes the matrix mapping fine 3D onto fine 2D elements
%
% Second step: computes the matrix mapping fine 2D onto coarse 2D elements
%
% For this last step:
%
% merge = n -> coarse elements are obtained by merging n contiguous fine elements
% merge = [m1 m2 ... mN] -> coarse element number 1 is obtained by merging
% the m1 first fine elements, etc..
% If sum(merge) < total number of fine elements, the remaining fine elements
% are merged to produce a N+1th coarse element.
%
% Dominique Gibert, April 2007
%
%%
c2f_3D= mk_2Dto3D_mapping(img.fwd_model);
c2f_2D= mk_gallery_2Dcto2Df_mapping(size(c2f_3D,2),merge_elems);
c2f= c2f_3D*c2f_2D;
n_coarse_elems = size(c2f,2);
n_fine_elems= size(img.fwd_model.elems,1);
img.elem_data= zeros(n_fine_elems,1);
img.params_mapping.function = @dg_c2f_matrix_mapping;
img.params_mapping.data = c2f;
img.params_mapping.params= zeros(n_coarse_elems,1);
disp(['3D FEM model of gallery created | number of elements = ' num2str(n_fine_elems)]);
disp(['                 Number of parameters to be inverted = ' num2str(n_coarse_elems)]);
end

function m2Dcto2Df= mk_gallery_2Dcto2Df_mapping(n_elems_2Dfine,merge)
%% merge = n -> coarse elements are obtained by merging n contiguous fine elements
% merge = [m1 m2 ... mN] -> coarse element number 1 is obtained by merging
% the m1 first fine elements, etc.. If sum(merge) < total number of fine
% elements, the remaining fine elements are merged to produce a N+1th coarse
% element.
if isscalar(merge)
    n_elems_2Dcoarse= n_elems_2Dfine/merge;
    m2Dcto2Df= zeros(n_elems_2Dfine,n_elems_2Dcoarse);
    for k= 1:n_elems_2Dcoarse
        m2Dcto2Df((1:merge)+merge*(k-1),k)= 1.;
    end
else
    n_elems_2Dcoarse= size(merge,2);
    if sum(merge) == n_elems_2Dfine
        csm= [0 cumsum(merge)];
    elseif sum(merge) < n_elems_2Dfine
        n_elems_2Dcoarse= n_elems_2Dcoarse+1;
        csm= [0 cumsum(merge) n_elems_2Dfine];
    elseif sum(merge) > n_elems_2Dfine
        csm= [0 cumsum(merge)];
        csm= [csm(find(csm < n_elems_2Dfine)) n_elems_2Dfine];
        n_elems_2Dcoarse= size(csm,2)-1;
    end
    m2Dcto2Df= zeros(n_elems_2Dfine,n_elems_2Dcoarse);
    for k= 1:n_elems_2Dcoarse
        m2Dcto2Df(csm(k)+1:csm(k+1),k)= 1.;
    end  
end
m2Dcto2Df= sparse(m2Dcto2Df);

end

function img = dg_c2f_matrix_mapping(img)
%% Function to be called to perform the mapping in the forward problem
mapping_matrix= img.params_mapping.data;
coarse_elems= img.params_mapping.params;
img.elem_data= mapping_matrix*coarse_elems;
end