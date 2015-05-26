function tank_recon(path)
imdl_2x16 = get_base_model;
if nargin == 0
   path = cd;
   SKIP = 4;
else
   str = strtok(path,'_');
   SKIP = str2num(str(5:end));
end

curpath = cd;

[stim, m_sel] = mk_stim_patterns(32,1,[0 SKIP+1], [0 SKIP+1], {'no_rotate_meas'}, 1);
imdl_2x16.fwd_model.stimulation = stim;
imdl_2x16.fwd_model.meas_select = m_sel;
vname = sprintf('imdl_2x16_s%d',SKIP);
if exist([vname '.mat'], 'file');
   load(vname);
else
   xvec = linspace(-.145,.145,33); xvec(end) = []; xvec = xvec + (xvec(2)-xvec(1))/2;
   [x, y, z] = ndgrid(xvec,xvec,.15:.02:.31);
   idx = ( x(:).^2 + y(:).^2 ) < (.9*.145)^2;
   gopt.distr = [x(idx) y(idx) z(idx)]';
   gopt.noise_figure = 1.0;
   eval([vname '= mk_GREIT_model(imdl_2x16_s4,0.2,[],gopt);']);
   save(vname,vname);
end
if nargin > 0
   cd(['../data/' path]);
end
ref_files = dir('*ref*.eit');
vh = [];
for rf = {ref_files.name}
   vh = [ vh mean(real(eidors_readdata(rf{1},'LQ4')),2)];
end
vh = mean(vh,2);


for L = [0 35 75]
   for R = [0 50 75 100]
      fname = sprintf('pos_L%d_R%d.eit',L,R);
      vi = mean(real(eidors_readdata(fname, 'LQ4')),2);
      if L==0 && R==0
         rimg = inv_solve(imdl_2x16_s4,vh,vi);
      else
         tmp =  inv_solve(imdl_2x16_s4,vh,vi);
         rimg.elem_data = [rimg.elem_data tmp.elem_data];
      end
   end
end

cd(curpath)

   rimg.calc_colours.ref_level = 0;
%    rimg.calc_colours.clim = max(abs(rimg.elem_data(:)))/2;
   rimg.show_slices.img_cols = size(rimg.elem_data,2);
   %%
   clear lvl
   lvl(:,3) = fliplr(.15:.02:.31);
   lvl(:,1:2) = Inf;
   clf
   show_slices(rimg,lvl);
   
   popt.supersampling_factor = 1;
   popt.resolution = 300;
   if nargin > 0
      fname = sprintf('tank_recon_%s.png',path);
   else
      fname = 'tank_recon.png';
   end
   print_convert(fname,popt);
   %%
   if nargin ~= 0 
      return
   end
   clf
   img = mk_image(imdl_2x16_s4.fwd_model,1);   
   h = show_fem(img);
   set(h, 'linewidth',.1);
   hold on
   [X, Y, Z] = sphere(40);
   X = 0.023*X; Y = 0.023*Y; Z = 0.023*Z;
   for L = [0 35 75]
      for R = [0 50 75 100]
         [x,y] = pol2cart(2*pi/32,R/100 * (0.14-0.024));
         surf(x+X,y+Y,Z+0.226+L/1000,'facecolor','b','edgecolor','none','facelighting','gouraud')
      end
   end
   l = linspace(0,2*pi,101);
   X = 0.14*sin(l);
   Y = 0.14*cos(l);
   for H = .14:.02:.32
      plot3(X,Y,H*ones(size(X)),'k','linewidth',2);
   end

   popt.supersampling_factor = 2;
   popt.resolution = 150;
   print_convert('tank.png',popt);