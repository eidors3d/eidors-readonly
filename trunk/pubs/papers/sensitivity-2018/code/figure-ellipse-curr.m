% contrasts01.m
sz = 5;
ellipse_x = 2;

shape_str = [ ...
 sprintf('solid ec = ellipticcylinder(0,0,0;%f,0,0;0,%f,0);\n', ...
     0.5*ellipse_x, 0.5/ellipse_x ) ...
 sprintf('solid left  = plane(-%f,0,0;-1,0,0);\n',sz) ...
 sprintf('solid right = plane( %f,0,0; 1,0,0);\n',sz) ...
 sprintf('solid brick = orthobrick(-%f,-2,-0.2;%f,2,0);\n',sz+1,sz+1) ...
 'solid cyl = ec and brick; tlo cyl;\n' ...
 'solid mainobj= left and right and (not cyl) and brick;\n'];
elec_pos = [ sz,  0,  0,   0,  0,  1;
            -sz,  0,  0,   0,  0, -1];
elec_shape=[2.0];
elec_obj = {'left','right'};
fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);
fmdl = mdl2d_from3d(fmdl);

fmdl.stimulation = stim_meas_list([1,2,1,2]);
img = mk_image(fmdl,1);
img.elem_data( fmdl.mat_idx{1} ) = 2;

img.calc_colours.ref_level = 1;

show_fem( img );
%print_convert contrasts_01a.png

% contrasts02.m
img.fwd_model.solve =      @fwd_solve_1st_order;
img.fwd_model.system_mat = @system_mat_1st_order;
[img.fwd_model.electrode(:).z_contact] = deal(1000); % Large

img.fwd_solve.get_all_meas = 1;
vv = fwd_solve(img);
imgv= rmfield(img,'elem_data');
imgv.node_data = vv.volt;
imgv.calc_colours.ref_level = mean(vv.volt);
show_fem(imgv);

%print_convert contrasts_02a.png


%contrasts03.m
sz = 5; img_idx = 'b';
for ellipse_x = [0.5,1,2];
   img = contrasts_04_modeller( sz, ellipse_x); 
   targ = img.fwd_model.mat_idx{1};
   for contrast = linspace( -3,3,7);
      img.elem_data( targ ) = exp(contrast);
      vv = fwd_solve(img);
      img_name = sprintf('04%c',img_idx); img_idx= img_idx+1;
   %contrasts03.m
      imgc = img;
      imgc.fwd_model.mdl_slice_mapper.npx = 128;
      imgc.fwd_model.mdl_slice_mapper.npy = 200;
      imgc.fwd_model.mdl_slice_mapper.level = [inf,inf,0];
      imgc.calc_colours.ref_level = 1;
      q = show_current(imgc,vv.volt);
      hh=show_fem(imgc);
      set(hh,'EdgeColor',[1,1,1]*.75);
      hold on;

      sy = linspace(-2,2,20); sx= 0*sy - sz;
      hh=streamline(q.xp,q.yp, q.xc, q.yc,-sx,sy); set(hh,'Linewidth',2);

      hold off;

      if ~exist('img_name'); img_name = '03a'; end
      print_convert(sprintf('contrasts_%s.png',img_name));
      print('-dpdf',sprintf('contrasts_%s.pdf',img_name));
   end
end

