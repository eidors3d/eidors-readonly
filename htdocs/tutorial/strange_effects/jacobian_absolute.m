function J = jacobian_absolute( fwd_model, img);
    vh = fwd_solve(img);
    flip = sign(vh.meas);
    flip = spdiags(flip, 0, length(flip), length(flip));
    J = flip*jacobian_adjoint(img);
