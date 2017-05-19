function  m2_3D_sensitivity;

   DO_LUNGS_MODEL= 0;

   opt.pagesize = [2,2.2];
   skips = [0,1,3,7];
   skips = [5];
   for s = skips+1;
      [stim,msel] = mk_stim_patterns(N_ELECS,1,[0,s],[0,s],{'no_meas_current'},1);

      imgm= do_simulation(stim,msel,1.0,1);
      print_it(sprintf('../figs/fig09-offplane-sens-1plane_%c',64+s),opt)

      imgm= do_simulation(stim,msel,1.0,2);
      print_it(sprintf('../figs/fig09-offplane-sens-2plane_%c',64+s),opt)

      imgm= do_simulation(stim,msel,1.0,2.1);
      print_it(sprintf('../figs/fig09-offplane-sens-2pl_sq%c',64+s),opt)

      if DO_LUNGS_MODEL
      imgm= do_simulation(stim,msel,0.3,1);
      print_it(sprintf('m2_3D_sensit_lung_%c',64+s),opt)
      end
   end
   if DO_LUNGS_MODEL
      imgm= do_simulation(stim,msel,1.0,1);
      show_fem_enhanced(imgm); view(0,22);
      print_it('m2_lung_model_10',opt)
      imgm= do_simulation(stim,msel,0.3,1);
      show_fem_enhanced(imgm); view(0,22);
      print_it('m2_lung_model_03',opt)
   end

function nel = N_ELECS; nel=32;

function fmdl= thorax_fmdl(n_elecs, planes)
   arms = [ ...
      'solid arm1 = cylinder( 1,0,2.3; 2,0,1.3;0.3);' ...
      'solid arm2 = cylinder(-1,0,2.3;-2,0,1.3;0.3);' ...
      'solid ob   = orthobrick(-1.2,-1.2,0;1.2,1.2,2.6);' ...
      'solid arms = ( arm1 or arm2 ) and ob;'];
   lungs= ['solid lungs = ' ...
           'ellipsoid( 0.4,0.0,1.30;0.44,0,0;0,0.5,0;0,0,1.0) or ' ...
           'ellipsoid(-0.4,0.0,1.30;0.44,0,0;0,0.5,0;0,0,1.0);'];
   extra = {'arms','lungs',[arms,lungs]};
   e_plane = 1.3;
   switch floor(planes)
      case 1; planes_line = [n_elecs, e_plane];
      case 2; planes_line = [n_elecs/2, e_plane+[0.2,-0.2]];
      otherwise; error('huh?');
   end
   fmdl= ng_mk_ellip_models([2.6,1.0,0.8,0.05],planes_line,[0.02],extra); 
   if planes==2.1 %square pattern
      curr = reshape(1:n_elecs,[],2)';
      desi = reshape(1:n_elecs,2,[]);
      desi(:,2:2:end) = flipud( desi(:,2:2:end));
      fmdl.electrode(desi) = fmdl.electrode(curr);
   end
   fmdl.nodes(:,3) = fmdl.nodes(:,3) - e_plane;

function imgm= do_simulation(stim,msel, lung_conductivity, planes);
   fmdl= thorax_fmdl(N_ELECS, planes);
   fmdl.stimulation= stim;
   fmdl.meas_select= msel;
   scl = 15;
   fmdl.nodes      = fmdl.nodes*scl; % 30cm lateral diameter
show_fem(fmdl);


   imgm = mk_image(fmdl,1); %NO Lungs
   imgm.elem_data(fmdl.mat_idx{3}) = lung_conductivity;

   img = imgm;

   J = calc_jacobian( img ) * spdiag(1./get_elem_volume(fmdl));
   img.elem_data = sqrt(sum(J.^2,1))';
   img.fwd_model.mdl_slice_mapper.npy =120*2+1;
   img.fwd_model.mdl_slice_mapper.npx =230  +1;
   filter = ones(15); filter = conv2(filter,filter);
   img.calc_slices.filter = filter;

   % remove horiz bands
   sens = calc_slices(img,[inf,0,inf]);
   sens(all(isnan(sens),2),:) = [];
   sens(:,all(isnan(sens),1)) = [];

   % scale vertical
   horz = find(~isnan(sens(end,:))); h1= horz(1); he= horz(end);   
   horz = he*ones(1,size(sens,2));
   horz(1:h1) = h1;
   horz(h1:he) = h1:he;

   maxsens = max( sens(:,horz), [], 1);
   sens = sens./( ones(size(sens,1),1)*maxsens );

   sens(end-(0:13),:)= NaN;
   
   img.calc_colours.ref_level = 0;
   img.calc_colours.backgnd = [1,1,1];
   if 0
      img.calc_colours.greylev = 0.001;
   else
      img.calc_colours.greylev =-0.201;
   end
   xax = linspace(-1.2*scl, 1.2*scl,size(sens,2));
   yax = linspace( 1.3*scl,-1.3*scl,size(sens,1));
   image(xax,yax,calc_colours(sens,img));
   hold on
   contour(xax, yax ,sens, [0.25, 0.5,0.75,0.9,0.95], ...
       'linewidth',2, 'color',[0,0,0]);
   hh=line([-15,15],[0,0]); set(hh,'linewidth',2, 'color',[0,0,0]);
   hold off;
   set(gca,'YDIR','normal');
%  axis equal
   set(gca,'box','off')
   ylim([-18,max(yax)])
   xlim([min(xax)-1,max(xax)])
%  xlabel('Lateral Distance (cm)')
%  ylabel('AP distance from electrode plane (cm)')





function conductivity = lung_bkgnd;
%  conductivity = 0.5;
   conductivity = 1.0;



function print_it(fname, inopt )
   if nargin==0; opt = struct(); end
   opt.horz_cut     = 50;
   opt.horz_space   = 10;
   opt.vert_cut     = 30;
   opt.vert_space   = 10;
   opt.resolution   = 300;
   opt.subfigure = '';
   opt.supersampling_factor = 1;

% Override fields  => stupid matlab puts in columns
   for fn = fieldnames( inopt )'
      fn = fn{1}; % supid again
      opt.( fn ) = inopt.( fn );
   end

   print_convert([fname, opt.subfigure, '.png'],opt);
   axis equal; ylim([-15,20]);
   print('-dpdf',[fname,'.pdf']); 

