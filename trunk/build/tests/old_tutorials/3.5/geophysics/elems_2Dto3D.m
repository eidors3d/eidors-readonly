function map_2Dto3D = elems_2Dto3D(n_elems_2D,n_levels)
a= (1:n_elems_2D:3*n_elems_2D*(n_levels-1));
for k= 1:n_elems_2D
    map_2Dto3D(k,:)= (k-1)+a;
end
