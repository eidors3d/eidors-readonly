% Simulate obj $Id$

% Simulate positions that are 0.25 above the plane
r =  linspace(0,0.9,100);
xyzr = [r;
        zeros(1,100);
     1.25*ones(1,100);
     0.05*ones(1,100)];

[vh,vi] = simulate_movement(imgs, xyzr);

% Show GREIT images
i_gr = mk_common_gridmdl('GREITc1');
imgr = inv_solve(i_gr, vh, vi(:,1:5:100));
imgr.show_slices.img_cols = 5;
show_slices(imgr);

print_convert('GREIT_test_params05a.png','-density 60');
