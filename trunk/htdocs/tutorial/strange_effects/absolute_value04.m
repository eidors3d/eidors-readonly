% simulate homogeneous voltages
vhomg = fwd_solve( mk_image(imdl, 1));
flip  = sign(vhomg.meas);

% recalculate the correct values
vhr = vha.*flip;
vir = via.*flip;

show_fem( inv_solve( imdl, vhr, vir));
print_convert absolute_value04a.png
