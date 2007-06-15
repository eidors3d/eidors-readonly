function cost = dg_gallery_penalty_simplex(elems_coarse,image,real_data)
%
maximum_resistivity= 50.;
minimum_resistivity= 5.;
lowerb= -log(maximum_resistivity);
upperb= -log(minimum_resistivity);
if isempty(find(elems_coarse < lowerb | elems_coarse > upperb))
    image.params_mapping.params= elems_coarse;
    [model_data,image]= dg_fwd_solve(image);
    cost= norm(real_data.meas-model_data.meas,2);
else
    cost= 1000000.;
end

return