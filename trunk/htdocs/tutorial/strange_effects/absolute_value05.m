% Modify Jacobian 
imdl.fwd_model.jacobian = @jacobian_absolute;
show_fem( inv_solve( imdl, vha, via));
print_convert absolute_value05a.png
