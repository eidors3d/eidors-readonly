function imv = solve_volt_image(img,msm,ctr);
% SOLVE_VOLT_IMAGE: solve and create images of voltage
%  Parameters:
%     img - Image to solve to
%     Optional parameters
%     msm - model_slice_mapper field
%     ctr - contour filed

% $Id$
% (C) 2021 Andy Adler. Licence GPL v2 or v3

    img.fwd_solve.get_all_meas = true;
    vv = fwd_solve(img);
    imv = rmfield(img, 'elem_data');
    imv.node_data = vv.volt;
    if nargin==1; return; end

    imv.fwd_model.mdl_slice_mapper = msm;
    imv.show_slices.axes_msm = true;
    if nargin==2; return; end

    imv.show_slices.contour_levels = ctr;

