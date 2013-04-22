if ~exist('ng_cyl_mdl.mat')
    make_ng_cyl_mdl
end

if ~exist('sim_targets.mat')
    sim_targets('sim_targets.mat');
end

if ~exist('jacobian_cyl.mat')
    calc_jacobian_mdl;
end

if ~exist('sim_testdata.mat')
    make_sim_testdata('sim_testdata.mat');
end

%ALGS = Recon_GREIT_ndiff.m \
%       Recon_NOSER_diff.mat \
%       Recon_NOSER_ndiff.mat \
%       Recon_Sheffield_backproj.mat

mkdir('ex_img')
example_images('ex_img/ex','Recon_GREIT_ndiff')