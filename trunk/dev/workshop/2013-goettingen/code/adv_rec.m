fmdl= ng_mk_cyl_models([2,1],[16,1.3,0.7],[0.15]);
[stim,mpat] = mk_stim_patterns(16,2,[0,1],[0,1],{},1);
fmdl.stimulation = stim;
fmdl.meas_select = mpat;

imdl = select_imdl(fmdl, {});
[rmdl, fmdl] = mk_pixel_slice(fmdl);
imdl.rec_model = rmdl;
imdl.fwd_model = fmdl;
imdl = select_imdl(imdl, {'Choose NF=0.5'});








fmdl= ng_mk_cyl_models([2,1],[16,1.3,0.7],[0.15]);

[stim,mpat] = mk_stim_patterns(16,2,[0,1],[0,1],{},1);
fmdl.stimulation = stim;
fmdl.meas_select = mpat;
 
 func = inline('(x-0.5).^2+(y-0.2).^2 +(z-1).^2<0.25^2','x','y','z');
 ball = elem_select(fmdl, func);
 img = mk_image(fmdl,1);
 vh = fwd_solve(img);
 img.elem_data = 1 + ball;
 vi = fwd_solve(img);
%  show_fem(img)
 
 
 opt.z_depth = 0.2;
 [recmdl fmdl]= mk_pixel_slice(fmdl,1,opt);
 
 show_fem(fmdl);
 hold on
 hh = show_fem(recmdl);
 set(hh,'EdgeColor','b');
 hold off
 return
 imdl = select_imdl(fmdl,{'Basic GN dif'});
% recmdl = rmfield(recmdl,'electrode');
 imdl.rec_model = recmdl;
%  imdl.prior_use_fwd_not_rec = 1;
% keyboard
  rimg = inv_solve(imdl, vh, vi);
  show_slices(rimg,[1 inf inf]);
return 
 imdl = select_imdl(fmdl,{'Basic GN dif'});
%  show_3d_slices(rimg);