% sl=shape_library('get','adult_male','boundary');
% fmdl = ng_mk_2d_model({flipud(sl),0.1},[16,-.17],[0.1,20]);
fmdl= ng_mk_cyl_models([0,1,.1],16,[0.1]);
fmdl = mdl_normalize(fmdl, 0);
stim = mk_stim_patterns(16,1,[0,5],[0,5],{},1);
fmdl.stimulation = stim;
fmdl.electrode= fmdl.electrode([1, 16:-1:2]);


func = inline('(x-0.5).^2+(y).^2<0.25^2','x','y','z');
left = elem_select(fmdl, func);
func = inline('(x+0.5).^2+(y).^2<0.25^2','x','y','z');
right = elem_select(fmdl, func);

imgh = mk_image(fmdl,1);
imgh.calc_colours.ref_level = 1;
imgh.calc_colours.clim = 1;
imgh.elem_data = 1 - 0.1*left + 0.1*right;
subplot(3,3,1)
show_fem(imgh);
vh = fwd_solve(imgh);
imgi = imgh;
imgi.elem_data = imgh.elem_data - 0.1;
subplot(3,3,2)
show_fem(imgh);
vi = fwd_solve(imgi);


for i = 1:2
    imdl = select_imdl(fmdl,{'Basic GN dif'});
    if i == 2
        imdl.jacobian_bkgnd = imgh;
    end
    imdl = select_imdl(imdl,{'Choose NF=1'});
    
    
    subplot(3,3,3*(i-1) + 4)
    img = mk_image(imdl);
    show_fem(img)
    
    subplot(3,3,3*(i-1) + 5)
    J= calc_jacobian(img);
    V= get_elem_volume(fmdl);
    sens = sum(J.^2,1)'./V;
    sens = (abs( sens ) ).^.5.*sign(sens);
    imgs = mk_image(img,sens);
    imgs.calc_colours.ref_level = median(imgs.elem_data);
    imgs.calc_colours.clim = .2;
    hh = show_fem(imgs);
%     eidors_colourbar(imgs);
    set(hh,'EdgeColor','none');
    
    subplot(3,3, 3*(i-1) + 6)
    imgr = inv_solve(imdl,vh, vi);
    imgr.calc_colours.ref_level = 0;
    imgr.calc_colours.clim      = 0.3;
    show_fem(imgr);
end

