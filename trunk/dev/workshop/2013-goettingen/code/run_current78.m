fmdl2= ng_mk_cyl_models([0,1,.1],16,0.05);
fmdl3= ng_mk_cyl_models([1,1,.1],[16,.5],0.05);
stim = mk_stim_patterns(16,1,[0,1],[0,1],{},1);
fmdl2.stimulation = stim;
fmdl3.stimulation = stim;

subplot(211);
img2= mk_image(fmdl2, 1);
J= calc_jacobian(img2);
V= get_elem_volume(fmdl2);
sens = J(5,:)'./V; 
sens = (abs( sens ) ).^.5.*sign(sens);
show_slices(mk_image(img2,sens));


subplot(212);
img3= mk_image(fmdl3, 1);
J= calc_jacobian(img3);
V= get_elem_volume(fmdl3);
sens = J(5,:)'./V; 
sens = (abs( sens ) ).^.5.*sign(sens);
show_slices(mk_image(img3,sens),1);


if 0
sl=shape_library('get','adult_male','boundary');
%fmdl = ng_mk_2d_model({flipud(sl),0.1},[16,-.17],[0.1,20]);
fmdl = ng_mk_extruded_model({1,sl,1},[16,1,0.5],[0.05 0 .02 0 60]);

img = mk_image(fmdl, 1);
img.fwd_solve.get_all_meas = 1;
    
stim =[1,2,3,4;
       1,9,3,4;
       1,9,4,5];
img.fwd_model.stimulation = ...
    stim_meas_list(stim,16);
v1 = fwd_solve(img);
    
subplot(221);
imgv = mk_image(fmdl,v1.volt(:,1));
show_slices(imgv,1);

subplot(223);
imgv = mk_image(fmdl,v1.volt(:,2));
show_slices(imgv,1);


J= calc_jacobian(img);
V= get_elem_volume(fmdl);
subplot(222);
sens = J(1,:)'./V; 
sens = (abs( sens ) ).^.2.*sign(sens);
show_slices(mk_image(img,sens),1);

subplot(224);
sens = J(2,:)'./V; 
sens = (abs( sens ) ).^.2.*sign(sens);
show_slices(mk_image(img,sens),1);


    fmdl= ng_mk_cyl_models([0,0.5,.05],4,[0.1]);
    stim.stim_pattern = [2;-1;0;0];
    stim.meas_pattern = [0,0,1,-1];
    fmdl.stimulation(1)= stim;
    img= mk_image(fmdl,1);

    img.fwd_solve.get_all_meas = 1;
    show_fem(img)
    volt = fwd_solve(img);

    imgv = mk_image(fmdl,volt.volt);
    show_fem(imgv);
end