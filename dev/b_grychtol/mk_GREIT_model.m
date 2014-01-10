function [imdl, weight]= mk_GREIT_model( fmdl, radius, weight, options )
% MK_GREIT_MODEL: make EIDORS inverse models using the GREIT approach
%   [imdl, weight]= mk_GREIT_model( mdl, radius, weight, options )
%
% Output: 
%   imdl   - GREIT inverse model
%   weight - value of the weight paramater chosed to satisfy the prescribed
%            noise figure (NF). See options.noise_figure below.
%
% Parameters:
%   mdl    - fwd model on which to do simulations, or
%          - inv model (experimental), or
%          - string specifying prepackaged models
%
%   radius - requested weighting matrix  (recommend 0.25 for 16 electrodes)
%   weight - weighting matrix (weighting of noise vs signal). Can be empty
%            options.noise_figure is specified
%   options- structure with fields:
%     imgsz         - [xsz ysz] reconstructed image size in pixels 
%                     (default: [32 32])
%     square_pixels - forces square pixels if 1 (default: 0)
%     Nsim          - number of training points (default: 1000)
%     distr         - distribution of training points:
%         0 -> original (as per GREITv1, default)
%         1 -> random, centre-heavy 
%         2 -> random, uniform
%         3 -> fixed, uniform (debug
%         Alternatively, a 3xN xyz or 4xN xyzr matrix can be provided. If
%         xyzr is given, r overrides target_size specifications.
%     target_size - size of simulated targets as proportion of mesh radius
%         (default: 0.02). Can be specified as [min_size max_size] for 
%         random variation
%     target_plane - the (mean) height z at which simulation targets are
%         placed. This controls the image plane. Default: mean electrode
%         height
%     target_offset - maximum allowed vertical displacement from the
%         target_plane (default: 0). Can be specified as
%         [down_offset up_offset].
%     noise_figure - the noise figure (NF) to achieve. Overwrites weight 
%         which will be optimised to achieve the target NF.
%     noise_figure_targets - circular target(s) to use for NF calculation
%         as an array of coordinates and radius xyzr [4xN] (default: single
%         target at the center at average electrode height with radius of
%         opt.target_size. Note that multiple targets are simultaneously
%         simulated in a single measurement, meaning they should not
%         overlap.
%     extra_noise - extra noise samples (such as electrode movement)
%     training_data - when specific target positions are provided in
%         opt.distr, training data can be provided as follows:
%         training_data.vh - a single (Mx1) reference frame
%         training_data.vi - MxNsim training frames matching the targets 
%         specified in opt.distr
%     desired_solution_fn - specify a function to calculate the desired 
%         image. It must have the signature:
%         D = my_function( xyc, radius, options); 
%         See CALC_GREIT_RM for details.
%
% NOTE
%   currently extra_noise is not supported
%   currently weighting matrix must be scalar
%               
% Examples
%   imdl =  mk_GREIT_model( 'c=1;h=2;r=.08;ce=16;bg=1;st=1;me=1;nd', 0.25, 10);
% OR
%   fmdl = mk_library_model('adult_male_16el');
%   fmdl.stimulation = stim;
%   fmdl.normalize_measurements = 1;
%   opt.noise_figure = 0.5; 
%   imdl = mk_GREIT_model(fmdl,0.25,5,opt);
%
% CITATION_REQUEST:
% AUTHOR: A Adler et al.
% TITLE: GREIT: a unified approach to 2D linear EIT reconstruction of lung
% images
% JOURNAL: Phys Meas
% YEAR: 2009
% VOL: 30
% NUM: 6
% PAGE: S35-55
% LINK: http://iopscience.iop.org/0967-3334/30/6/S03
%
% See also CALC_GREIT_RM

% (C) 2010 Andy Adler. License: GPL version 2 or version 3
% $Id$

citeme(mfilename);

if isstr(fmdl) && strcmp(fmdl,'UNIT_TEST'); do_unit_test; return; end

if nargin < 4, options = [];end
[imdl,fmdl,imgs] = parse_fmdl(fmdl);
options = parse_options(options,fmdl,imdl);

cache_obj= { fmdl, imdl, imgs, radius, weight, options};

out= eidors_obj('get-cache', cache_obj, 'mk_GREIT_model');
if ~isempty(out)
   eidors_msg('mk_GREIT_model: using cached value', 3);
   imdl= out{1}; weight=out{2}; return
end
[imdl, weight]= mk_GREIT_model_calc( fmdl, imdl, imgs, radius, weight, options);

eidors_obj('set-cache', cache_obj, 'mk_GREIT_model', {imdl,weight});
eidors_msg('mk_GREIT_model: setting cached value', 3);

function [imdl, weight]= mk_GREIT_model_calc( fmdl, imdl, imgs, radius, weight, opt)

Nsim = opt.Nsim;
[vi,vh,xy,opt]= stim_targets(imgs, Nsim, opt );

%Calculate rec_model (if absent)
if ~isfield(imdl,'rec_model');
   if size(opt.target_plane) == [1,1]
      [imdl.rec_model imdl.fwd_model] = mk_pixel_slice(fmdl,opt.target_plane,opt);
      imdl.rec_model.nodes(:,3) = []; % the third dimension complicated display
      % medical orientation
      imdl.rec_model.mdl_slice_mapper.y_pts = fliplr(imdl.rec_model.mdl_slice_mapper.y_pts);
   else
      opt.do_coarse2fine = 0; % coarse2fine is needed for calculation of
                              % solution error, but is otherwise meaningless
      r = {};
      inside = logical([]);
      c2f = [];
      for i = 1:length(opt.target_plane)
         r{i} = mk_pixel_slice(fmdl,opt.target_plane(i),opt);         
         inside = [inside; r{i}.inside];
         c2f  = blkdiag( c2f, r{i}.coarse2fine);
      end
      rmdl = merge_meshes(r{:});
      rmdl.inside = inside;
      rmdl.coarse2fine = c2f;
      imdl.rec_model = rmdl;
   end
end
% user MUST specify the inside array
inside = imdl.rec_model.inside;

% IMPORTANT: opt.meshsz must now reflect the rec_model, not the fwd_model!!
minnode = min(imdl.rec_model.nodes);
maxnode = max(imdl.rec_model.nodes);
opt.meshsz = [minnode(1) maxnode(1) minnode(2) maxnode(2)];

imdl.solve = @solve_use_matrix;
log_level = eidors_msg( 'log_level', 1);

if ~isempty(opt.noise_figure)
    target = opt.noise_figure;
    if ~isempty(weight)
        eidors_msg('mk_GREIT_model: Using weight parameter as a guess, options.noise_figure is non-empty');
    else
        weight = target;
    end
    
    xyzr = opt.noise_figure_targets;
    [jnk,vi_NF] = simulate_movement(imgs,xyzr');
    vi_NF = sum(vi_NF,2); % sum the targets
    eidors_msg('mk_GREIT_model: Finding noise weighting for given Noise Figure',1);
    eidors_msg('mk_GREIT_model: This will take a while...',1);
    f = @(X) to_optimise(vh,vi,xy, radius, X, opt, inside, imdl, target, vi_NF);
    fms_opts.TolFun = 0.01*target; %don't need higher accuracy
    if exist('OCTAVE_VERSION')
       % octave doesn't currently (2013 Apr) include an fminsearch function
       [weight, NF] = fminsearch_octave(f, weight,fms_opts);
    else
       [weight, NF] = fminsearch(f, weight,fms_opts);
    end
    eidors_msg(['mk_GREIT_model: Optimal solution gives NF=' ... 
        num2str(NF+target) ' with weight=' num2str(weight)],1);
end
eidors_msg( 'log_level', log_level);
RM= calc_GREIT_RM(vh,vi, xy, radius, weight, opt );
imdl.solve_use_matrix.RM = resize_if_reqd(RM,inside);
imdl.jacobian_bkgnd = imgs;
%imdl.solve_use_matrix.map = inside;

function out = to_optimise(vh,vi,xy,radius,weight, opt, inside, imdl, ...
    target,vi_NF)

   % calculate GREIT matrix as usual
   RM = calc_GREIT_RM(vh,vi,xy, radius, weight, opt);
   imdl.solve_use_matrix.RM = resize_if_reqd(RM,inside);
   NF = calc_noise_figure(imdl,vh, vi_NF);
   eidors_msg(['NF = ', num2str(NF), ' weight = ', num2str(weight)],1);
   out = (NF - target)^2;
%    out = (mean(NF) - target)^2 + std(NF);


function  imgs = get_prepackaged_fmdls( fmdl );
  switch fmdl
    case 'c=1;h=2;r=.08;ce=16;bg=1;st=1;me=1;nd'
      fmdl = ng_mk_cyl_models([2,1,0.18],[16,1],[0.05]); 
      fmdl.stimulation = mk_stim_patterns(16,1,[0,1],[0,1],{},1);
      fmdl = mdl_normalize(fmdl,1);
      imgs= mk_image( fmdl, 1);
    otherwise
      error('specified fmdl (%s) is not understood', fmdl);
  end

function [vi,vh,xyz,opt]= stim_targets(imgs, Nsim, opt );
    fmdl = imgs.fwd_model;
   ctr =  mean(fmdl.nodes);  
   maxx = max(abs(fmdl.nodes(:,1) - ctr(1)));
   maxy = max(abs(fmdl.nodes(:,2) - ctr(2)));
   
   if numel(opt.distr) > 1
      xyzr = opt.distr;
      xyzr(4,:) = calc_radius(0,opt,Nsim);
   else
      
      switch opt.distr
         case 0 % original
            r = linspace(0,0.9, Nsim);
            th = r*4321; % want object to jump around in radius
            xyzr = [maxx*r.*cos(th); maxy*r.*sin(th);
               opt.target_plane*ones(1,Nsim);
               0.05*mean([maxx,maxy])*ones(1,Nsim)];
            
         case 1 %centre-heavy
            F = fourier_fit(opt.contour_boundary(:,1:2));
            v = linspace(0,1,Nsim*100+1); v(end)=[];
            pts = fourier_fit(F,v);
            idx_p = floor(rand(Nsim,1)*Nsim*100);
            xyzr = pts(idx_p,:)'.*repmat(rand(Nsim,1),[1 2])';
            xyzr(3,:) = calc_offset(opt.target_plane,opt,Nsim);
            
            % TODO: What size is good here and how to figure it out?
            xyzr(4,:) = calc_radius(mean([maxx maxy]),opt,Nsim);
         case 2 %uniform
            %            F = fourier_fit(opt.contour_boundary(:,1:2));
            %            v = linspace(0,1,101); v(end)=[];
            %            pts = fourier_fit(F,v);
            pts = opt.contour_boundary(:,1:2);
            % avoid edges
            pts = 0.9*( pts - repmat(ctr(1:2),length(pts),1) ) + repmat(ctr(1:2),length(pts),1);
            % using maxx and maxy below would in general not produce a
            % uniform distribution
            lim = max(maxx, maxy);
            x = ctr(1) + (rand(Nsim*10,1)-0.5)*2*lim;
            y = ctr(2) + (rand(Nsim*10,1)-0.5)*2*lim;
            IN = inpolygon(x,y,pts(:,1),pts(:,2));
            xyzr(1,:) = x(find(IN,Nsim));
            xyzr(2,:) = y(find(IN,Nsim));
            xyzr(3,:) = calc_offset(opt.target_plane,opt,Nsim);
            % TODO: What size is good here and how to figure it out?
            xyzr(4,:) = calc_radius(mean([maxx maxy]),opt,Nsim);
         case 3 % uniform, non-random
            %            F = fourier_fit(opt.elec_loc(:,1:2));
            %            v = linspace(0,1,101); v(end)=[];
            %            pts = fourier_fit(F,v);
            xyzr = [];
            for i = 1:length(opt.contour_boundary)
               pts = opt.contour_boundary{i}(:,1:2);
               lim = max(maxx, maxy);
               frac = polyarea(pts(:,1),pts(:,2)) / (2*lim)^2;
               [x,y] = ndgrid( linspace(-lim,lim,ceil(sqrt(Nsim/frac))), ...
                  linspace(-lim,lim,ceil(sqrt(Nsim/frac))));
               
               x = x+ctr(1); y = y + ctr(2);
               IN = inpolygon(x,y,pts(:,1),pts(:,2));
               tmp(1,:) = x(find(IN));
               tmp(2,:) = y(find(IN));
               tmp(3,:) = calc_offset(opt.target_plane(i),opt,size(tmp,2));
               % TODO: What size is good here and how to figure it out?
               tmp(4,:) = calc_radius(mean([maxx maxy]),opt,size(tmp,2));
               xyzr = [xyzr tmp];
            end
      end
   end
   eidors_msg(['mk_GREIT_model: Using ' num2str(size(xyzr,2)) ' points']);
   before = size(xyzr,2);
   if isfield(opt, 'training_data')
      vh = opt.training_data.vh;
      vi = opt.training_data.vi;
   else
      [vh,vi,xyzr] = simulate_movement(imgs, xyzr);
      after = size(xyzr,2);
      if(after~=before)
         eidors_msg(['mk_GREIT_model: Now using ' num2str(after) ' points']);
      end
   end
   xyz = xyzr(1:3,:);

function z = calc_offset(z0,opt,Nsim)
    if opt.random_offset
        l_bnd = opt.target_offset(1);
        width = sum(opt.target_offset(1:2));
        z = z0 - l_bnd + rand(Nsim,1)*width;
    else
        z = z0*ones(Nsim,1);
    end

function r = calc_radius(R,opt,Nsim)
   if opt.random_size
      min_sz = opt.target_size(1);
      max_sz = opt.target_size(2);
      range = max_sz - min_sz;
      r = (min_sz + rand(Nsim,1)*range)*R;
   else
      if numel(opt.target_size) == 1
         % default
         r = opt.target_size(1)*ones(Nsim,1)*R;
      else
         % absolute sizes are provided
         r = opt.target_size;
      end
   end
           
   
   
function RM = resize_if_reqd(RM,inside);
   szRM = size(RM,1);
   if sum(inside) == szRM
      % RM is fine
   elseif size(inside,1) == szRM
      RM = RM(inside,:);
   else
      error('mismatch in size of provided RecMatrix');
   end


function [imdl,fmdl,imgs] = parse_fmdl(fmdl);
   imdl = []; 
   if isstr(fmdl)
      imgs = get_prepackaged_fmdls( fmdl );
      fmdl = imgs.fwd_model;
   elseif isfield(fmdl,'type');
     switch fmdl.type
   %  if we get a fwd_model, assume uniform conductivity backgnd of 1
       case 'fwd_model'; imgs = mk_image( fmdl, 1);
   %  if we get an image, use it. It may have a non-uniform backgnd
       case 'image';     imgs = fmdl; % fmdl was an image
                         fmdl = imgs.fwd_model; % now it's a fmdl
       case 'inv_model'; imdl = fmdl;
                         fmdl = imdl.fwd_model;
                         imgs = calc_jacobian_bkgnd(imdl);
       otherwise; error('unrecognized eidors object');
     end
   else
      error('specified parameter must be an object or a string');
   end
   % Prepare model
   if isempty(imdl)
      imdl = select_imdl( fmdl,{'Basic GN dif'});
   end

function opt = parse_options(opt,fmdl,imdl);

    if ~isfield(opt, 'imgsz'),     opt.imgsz = [32 32]; end
    if ~isfield(opt, 'square_pixels')
        opt.square_pixels = 0;
    end
    % Allow imdl.rec_model to overwrite options.imgsz
    % TODO: fix this to support multi-layer models
%     if isfield(imdl,'rec_model') && ~isempty(imdl.rec_model)
%         % this assumes rec_model is a rectangular grid, as it should
%         opt.imgsz(1) = numel(unique(imdl.rec_model.nodes(:,1)))-1;
%         opt.imgsz(2) = numel(unique(imdl.rec_model.nodes(:,2)))-1;
%     end  
    
    if ~isfield(opt, 'distr'),     opt.distr = 3;       end 
    if ~isfield(opt, 'Nsim' ),     opt.Nsim  = 1000;    end
    if ~isfield(opt, 'noise_figure'), opt.noise_figure = []; end
    if isfield(opt,'extra_noise')
      error('mk_GREIT_model: doesn''t currently support extra_noise');
    end
    
    if numel(opt.distr) > 1
       % user specified target positions
       switch size(opt.distr,1)
          case 3 %xyz
          case 4 %xyzr
             opt.target_size = opt.distr(4,:);
             opt.distr       = opt.distr(1:3,:);
             opt.random_size = false;
          otherwise
             error('opt.distr must have the size of 3xN or 4xN');
       end
       opt.Nsim = size(opt.distr,2);
    end
    
    if isfield(opt, 'training_data');
       try
          opt.training_data.vh;
          opt.training_data.vi;
       catch
          error('Both reference and target measurements must be provided');
       end
       if size(opt.training_data.vi,2) ~= size(opt.distr,2)
          error('Training data and target positions must match in lenght');
       end
    end
    
    if ~isfield(opt, 'target_size')
        opt.target_size = 0.05;
    end
    
    if ~isfield(opt, 'random_size')
       if numel(opt.target_size) == 2
          if opt.target_size(1) == opt.target_size(2);
             opt.random_size = false;
          else
             opt.random_size = true;
          end
       else
          opt.random_size = false;
       end
    end
    
    % Calculate the position of the electrodes
    Nelecs = length(fmdl.electrode);
    for i=1:Nelecs
       enodesi = fmdl.electrode(i).nodes;
       elec_loc(i,:) = mean( fmdl.nodes( enodesi,:),1 );
    end
    opt.elec_loc = elec_loc;
    
    if ~isfield(opt, 'target_plane')
          opt.target_plane = mean(elec_loc(:,3));
    else
        t = opt.target_plane;
        minnode = min(fmdl.nodes);
        maxnode = max(fmdl.nodes);
        if any(t<minnode(3)) || any(t>maxnode(3))
            warning('options.target_plane is outside the model!');
            eidors_msg('mk_GREIT_model: Resorting to default target_plane');
            opt.target_plane = mean(elec_loc(:,3));
        end
    end
    if ~isfield(opt, 'target_offset')
        opt.target_offset = 0;
    end
    if sum(size(opt.target_offset)) == 2
        if opt.target_offset < 0, opt.target_offset = 0; end
        opt.target_offset(2) = opt.target_offset(1);
    end
    if any(opt.target_offset > 0)
        opt.random_offset = true;
    else
        opt.random_offset = false;
    end

    if ~isfield(opt,'noise_figure_targets');
       R = max(max(fmdl.nodes(:,1:2)) - min(fmdl.nodes(:,1:2)));
       xyzr = mean(fmdl.nodes);
       xyzr(3) = mean(opt.target_plane);
       xyzr(4) = mean(opt.target_size)*0.5*R;
       opt.noise_figure_targets = xyzr;
    end

       
    
    
    try, opt.normalize = fmdl.normalize_measurements;
    catch, 
        opt.normalize = 0;
        eidors_msg('mk_GREIT_model: fmdl.normalize_measurements not specified, assuming 0');
    end
    opt.contour_boundary = {};
    for i = 1:length(opt.target_plane)
       % find the boundary at target level (needed in many places)
       slc = mdl_slice_mesher(fmdl,[inf inf opt.target_plane(i)]);
       bnd = find_boundary(slc.fwd_model);
       opt.contour_boundary = [opt.contour_boundary, ...
                                 order_loop(slc.fwd_model.nodes(unique(bnd),:)) ];
    end
    

function do_unit_test

do_performance_test; 
% return;
figure
% Create a 3D elliptical cylinder with 16 circular electrodes 
fmdl_1= ng_mk_ellip_models([1,1.2,0.8],[16,0.5],[0.1]); %show_fem(fmdl);
% Put two balls into the elliptical cylinder
extra={'ball','solid ball = sphere(0.5,0.5,0.5;0.1);'};
[fmdl_2,mat_idx]= ng_mk_ellip_models([1,1.2,0.8],[16,0.5],[0.1],extra); 
% Set the model to use adjacent current patterns
stim = mk_stim_patterns(16,1,[0,1],[0,1],{}); 
fmdl_1.stimulation = stim;
fmdl_2.stimulation = stim;
% Simulate homogeneous voltages (background conductivity = 0.5);
img = mk_image(fmdl_2, 0.5); vh = fwd_solve(img); %show_fem(img);
% Simulate inhomogeneous voltages (ball conductivity = 1.0);
img.elem_data(mat_idx{2})= 1.0; vi = fwd_solve(img); 
show_fem(img);
% Reconstruct the image using GREITv1
imdl= mk_common_gridmdl('GREITc1'); 
img= inv_solve(imdl,vh,vi);
figure, subplot(2,2,1);
show_slices(img)

% Create a GREIT model for the ellipse
opt.noise_figure = 0.5; opt.distr = 3;opt.square_pixels = 1; %other options are defaults
fmdl_2 = mdl_normalize(fmdl_2,0);
% use the true model (inverse crime)
imdl1 = mk_GREIT_model(mk_image(fmdl_2,0.5), 0.25, [], opt);
img1= inv_solve(imdl1,vh,vi);  
subplot(2,2,2);show_slices(img1);

% use honogenous model 
fmdl_1 = mdl_normalize(fmdl_1,0);
imdl2 = mk_GREIT_model(mk_image(fmdl_1,0.5), 0.25, [], opt);
img2= inv_solve(imdl2,vh,vi); 
subplot(2,2,3); show_slices(img2);


% specify targets for NF calc
opt.noise_figure_targets = [-.5 0 .5 .2;.5 0 .5 .2;];
imdl3 = mk_GREIT_model(mk_image(fmdl_1,0.5), 0.25, [], opt);
img3= inv_solve(imdl3,vh,vi); 
subplot(2,2,4); show_slices(img3);
% cleanup
opt = rmfield(opt,'noise_figure_targets');


%% repeat with normalized data
fmdl_2 = mdl_normalize(fmdl_2,1);
% use the true model (inverse crime)
imdl3 = mk_GREIT_model(mk_image(fmdl_2,0.5), 0.25, [], opt);
img3= inv_solve(imdl3,vh,vi); 

% use honogenous model 
fmdl_1 = mdl_normalize(fmdl_1,1);
imdl4 = mk_GREIT_model(mk_image(fmdl_1,0.5), 0.25, [], opt);
img4= inv_solve(imdl4,vh,vi); 

figure
show_slices([img1 img2 img3 img4])


%% Use a prepackaged model
fmdl = mk_library_model('adult_male_16el_lungs');
fmdl.stimulation = stim;
fmdl = mdl_normalize(fmdl,1);
img = mk_image(fmdl,1);
img.elem_data([fmdl.mat_idx{2}; fmdl.mat_idx{3}],1) = 0.3;
vh = fwd_solve(img);
img.elem_data([fmdl.mat_idx{2}; fmdl.mat_idx{3}],1) = 0.4;
vi = fwd_solve(img);


fmdl2 = mk_library_model('adult_male_16el');
fmdl2.stimulation = stim;
fmdl2 = mdl_normalize(fmdl2,1);

opt.imgsz = [50 30];
opt.square_pixels = 1;
imdl = mk_GREIT_model(fmdl2,0.25,3,opt);

img = inv_solve(imdl,vh, vi);
figure
show_slices(img);


function do_performance_test
% Reconstruct GREIT Images
imdl_v1 = mk_common_gridmdl('GREITc1');
imdl_v1.inv_solve.calc_solution_error = false;

% Reconstruct backprojection Images
imdl_bp = mk_common_gridmdl('backproj');

% Recosntruct with new GREIT
% fmdl = ng_mk_cyl_models([2,1,0.05],[16,1],[0.05]); 
fmdl = mk_library_model('cylinder_16x1el_fine');
fmdl.stimulation = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'}, 1);
opt.noise_figure = 0.88;
opt.target_size = 0.1;
opt.distr = 0;
imdl_gr = mk_GREIT_model(fmdl, 0.2, [], opt);
 

test_performance( { imdl_v1, imdl_gr},fmdl );
