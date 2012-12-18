function test_physics(N)
% Test the new physics interface

if nargin==0
    N = 1:4;
end
for i = N
    run_test(i);
end

function run_test(N)
eidors_msg('log_level',0);
eidors_cache off
imdl = prepare_model;
fmdl = imdl.fwd_model;
switch N
    case 1
        imgh = mk_image(fmdl,1,'Old code default');
        imgi = imgh;
        imgi.elem_data(32) = 2;
    case 2 
        imgh = mk_image(fmdl,1,'conductivity','New code default');
        imgi = imgh;
        imgi.conductivity.elem_data(32) = 2;
    case 3
        imgh = mk_image(fmdl,1,'conductivity','Conductivity');
        imgi = imgh;
        imgi.conductivity.elem_data(32) = 2;
        imdl = rmfield(imdl,'jacobian_bkgnd');
        imdl.jacobian_bkgnd.conductivity = imgh.conductivity;
    case 4
        imgh = mk_image(fmdl,1,'resistivity','Resistivity');
        imgi = imgh;
        imgi.resistivity.elem_data(32) = 1/2;
        imdl = rmfield(imdl,'jacobian_bkgnd');
        imdl.jacobian_bkgnd.resistivity = imgh.resistivity;
    otherwise
        error('No test %d',N);
end
fprintf([imgh.name '...\t']);
try
%     S = calc_system_mat(imgh); %???
    J = calc_jacobian(imgh);
    vh = fwd_solve(imgh);
    vi = fwd_solve(imgi);
    imgr = inv_solve(imdl,vh, vi);
    fprintf('OK\n');
catch err
    rethrow(err)
end
subplot(3,3,N)
show_slices(imgr);
eidors_msg('log_level',2);
eidors_cache on

function imdl = prepare_model
imdl = mk_common_model('b2C',16);



