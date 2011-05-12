function [imdl, weight]= mk_GREIT_model( fmdl, radius, weight, options )
% MK_GREIT_MODEL: make EIDORS inverse models using the GREIT approach
%   [imdl, weight]= mk_GREIT_model( fmdl, radius, weight, options )
%
% Output: 
%   imdl   - GREIT inverse model
%   weight - value of the weight paramater chosed to satisfy the prescribed
%            noise figure (NF). See options.noise_figure below.
%
% Parameters:
%   fmdl   - fwd model on which to do simulations, or
%          - string specifying prepackaged models
%
%   radius - requested weighting matrix  (recommend 0.25 for 16 electrodes)
%   weight - weighting matrix (weighting of noise vs signal). Can be empty
%            options.noise_figure is specified
%   options- structure with fields:
%     imgsz       - [xsz ysz] reconstructed image size in pixels 
%                   (default: [32 32])
%     Nsim        - number of training points (default: 1000)
%     distr       - distribution of training points:
%         0 -> original (as per GREITv1, default)
%         1 -> random, centre-heavy 
%         2 -> random, uniform
%         3 -> fixed, uniform (debug)
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
%     extra_noise - extra noise samples (such as electrode movement)
%
% NOTE
%   currently extra_noise is not supported
%   currently weighting matrix must be scalar
              
% Example
%   imdl =  mk_GREIT_model( 'c=1;h=2;r=.08;ce=16;bg=1;st=1;me=1;nd', 0.25, 10);

% (C) 2010 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isstr(fmdl) && strcmp(fmdl,'UNIT_TEST'); do_unit_test; return; end

if isstr(fmdl)
   imgs = get_prepackaged_fmdls( fmdl );
elseif isfield(fmdl,'type');
  switch fmdl.type
%  if we get a fwd_model, assume uniform conductivity backgnd of 1
    case 'fwd_model'; imgs = mk_image( fmdl, 1);
%  if we get an image, use it. It may have a non-uniform backgnd
    case 'image';     imgs = fmdl; % fmdl was an image
                      fmdl = imgs.fwd_model; % now it's a fmdl
    otherwise; error('unrecognized eidors object');
  end
else
   error('specified parameter must be an object or a string');
end

opt = parse_options(options,fmdl);
Nsim = opt.Nsim;
[vi,vh,xy,bound,elec_loc,opt]= stim_targets(imgs, Nsim, opt );
maxnode = max(fmdl.nodes); minnode = min(fmdl.nodes);
opt.normalize = imgs.fwd_model.normalize_measurements;
opt.meshsz = [minnode(1) maxnode(1) minnode(2) maxnode(2)];


%imdl = mk_common_gridmdl('b2c', RM);
 
% Prepare model
 xgrid = linspace(minnode(1),maxnode(1),opt.imgsz(1)+1);
 ygrid = linspace(minnode(2),maxnode(2),opt.imgsz(2)+1);
 rmdl = mk_grid_model([],xgrid,ygrid);
 x_avg = conv2(xgrid, [1,1]/2,'valid');
 y_avg = conv2(ygrid, [1,1]/2,'valid');
 [x,y] = ndgrid( x_avg, y_avg);
 inside = inpolygon(x(:),y(:),bound(:,1),bound(:,2) );
 
 ff = find(~inside);
 rmdl.elems([2*ff, 2*ff-1],:)= [];
 rmdl.coarse2fine([2*ff, 2*ff-1],:)= [];
 rmdl.coarse2fine(:,ff)= [];
 
imdl = select_imdl( fmdl,{'Basic GN dif'});
imdl.solve = @solve_use_matrix;
imdl.rec_model = rmdl;

log_level = eidors_msg( 'log_level', 1);

if ~isempty(opt.noise_figure)
    if ~isempty(weight)
        eidors_msg('mk_GREIT_model: Ignoring weight parameter, options.noise_figure is non-empty');
    end
    target = opt.noise_figure;
    xyzr = mean(fmdl.nodes);
    xyzr(3) = opt.target_plane;
    xyzr(4) = opt.target_size;
    [jnk,vi_NF] = simulate_movement(imgs,xyzr');
    eidors_msg('mk_GREIT_model: Finding noise weighting for given Noise Figure',1);
    eidors_msg('mk_GREIT_model: This will take a while...',1);
    f = @(X) to_optimise(vh,vi,xy, radius, X, opt, inside, imdl, target, vi_NF);
    fms_opts.TolFun = 0.01*target; %don't need higher accuracy
    [weight, NF] = fminsearch(f, target);
    eidors_msg(['mk_GREIT_model: Optimal solution gives NF=' ... 
        num2str(NF+target) ' with weight=' num2str(weight)],1);
end
eidors_msg( 'log_level', log_level);
RM= calc_GREIT_RM(vh,vi, xy, radius, weight, opt );
imdl.solve_use_matrix.RM = resize_if_reqd(RM,inside);
%imdl.solve_use_matrix.map = inside;

function out = to_optimise(vh,vi,xy,radius,weight, opt, inside, imdl, ...
    target,vi_NF)

   % calculate GREIT matrix as usual
   RM = calc_GREIT_RM(vh,vi,xy, radius, weight, opt);
   imdl.solve_use_matrix.RM = resize_if_reqd(RM,inside);
   NF = calc_noise_params(imdl,vh, vi_NF);
   eidors_msg(['NF = ', num2str(NF), ' weight = ', num2str(weight)],1);
   out = (NF - target)^2;
%    out = (mean(NF) - target)^2 + std(NF);
   
function  imgs = get_prepackaged_fmdls( fmdl );
  switch fmdl
    case 'c=1;h=2;r=.08;ce=16;bg=1;st=1;me=1;nd'
      fmdl = ng_mk_cyl_models([2,1,0.08],[16,1],[0.05]); 
      fmdl.stimulation = mk_stim_patterns(16,1,[0,1],[0,1],{},1);
      fmdl.normalize_measurements = 1;
      imgs= mk_image( fmdl, 1);
    otherwise
      error('specified fmdl (%s) is not understood', fmdl);
  end

function [vi,vh,xy,bound,elec_loc,opt]= stim_targets(imgs, Nsim, opt );
    fmdl = imgs.fwd_model;
   ctr =  mean(fmdl.nodes);  
   maxx = max(abs(fmdl.nodes(:,1) - ctr(1)));
   maxy = max(abs(fmdl.nodes(:,2) - ctr(2)));

   % Calculate the position of the electrodes
   Nelecs = length(imgs.fwd_model.electrode);
   for i=1:Nelecs
       enodesi =     imgs.fwd_model.electrode(i).nodes;
       elec_loc(i,:) = mean( imgs.fwd_model.nodes( enodesi,:),1 );
   end
   
   if opt.target_plane == 1i
       opt.target_plane = mean(elec_loc(:,3));
   end
   
   % calculate the boundary (for external use)
   F = fourier_fit(elec_loc(:,1:2));
   v = linspace(0,1,100+1); v(end)=[];
   bound = fourier_fit(F,v);

   
   switch opt.distr 
       case 0 % original
           r = linspace(0,0.9, Nsim);
           th = r*4321; % want object to jump around in radius
           xyzr = [maxx*r.*cos(th); maxy*r.*sin(th); 
               opt.target_plane*ones(1,Nsim);
               0.05/mean([maxx,maxy])*ones(1,Nsim)];
       
       case 1 %centre-heavy
           % Now, use elec_loc to figure out the shape. We can assume the obj is
           % extruded in z
           F = fourier_fit(elec_loc(:,1:2));
           v = linspace(0,1,Nsim*100+1); v(end)=[];
           pts = fourier_fit(F,v);
           idx_p = floor(rand(Nsim,1)*Nsim*100);
           xyzr = pts(idx_p,:)'.*repmat(rand(Nsim,1),[1 2])';
           xyzr(3,:) = calc_offset(opt.target_plane,opt,Nsim);
           
           % TODO: What size is good here and how to figure it out?
           xyzr(4,:) = calc_radius(mean([maxx maxy]),opt,Nsim);
       case 2 %uniform
           F = fourier_fit(elec_loc(:,1:2));
           v = linspace(0,1,101); v(end)=[];
           pts = fourier_fit(F,v);
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
           F = fourier_fit(elec_loc(:,1:2));
           v = linspace(0,1,101); v(end)=[];
           pts = fourier_fit(F,v);
           lim = max(maxx, maxy);
           frac = polyarea(pts(:,1),pts(:,2)) / (2*lim)^2;
           [x,y] = ndgrid( linspace(-lim,lim,ceil(sqrt(Nsim/frac))), ...
                           linspace(-lim,lim,ceil(sqrt(Nsim/frac))));
                      
           x = x+ctr(1); y = y + ctr(2);    
           IN = inpolygon(x,y,pts(:,1),pts(:,2));
           xyzr(1,:) = x(find(IN));
           xyzr(2,:) = y(find(IN));
           xyzr(3,:) = calc_offset(opt.target_plane,opt,size(xyzr,2));
           % TODO: What size is good here and how to figure it out?
           xyzr(4,:) = calc_radius(mean([maxx maxy]),opt,size(xyzr,2));
           eidors_msg(['mk_GREIT_model: Using ' num2str(size(xyzr,2)) ' points']);
   end
   before = size(xyzr,2);
   [vh,vi,xyzr] = simulate_movement(imgs, xyzr);
   after = size(xyzr,2);
   if(after~=before)
       eidors_msg(['mk_GREIT_model: Now using ' num2str(after) ' points']);
   end
   xy = xyzr(1:2,:);

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
       r = opt.target_size(1)*ones(Nsim,1)*R;
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

function opt = parse_options(opt,fmdl);
    
    maxnode = max(fmdl.nodes); minnode = min(fmdl.nodes);
    
    if ~isfield(opt, 'imgsz'),     opt.imgsz = [32 32]; end
    if ~isfield(opt, 'distr'),     opt.distr = 0; end 
    if ~isfield(opt, 'Nsim' ),     opt.Nsim  = 1000; end
    if ~isfield(opt, 'noise_figure'), opt.noise_figure = []; end
    if isfield(opt,'extra_noise')
      error('mk_GREIT_model: doesn''t currently support extra_noise');
    end
    if ~isfield(opt, 'target_size')
        opt.target_size = 0.02;
    end
    if sum(size(opt.target_size)) > 2
        if opt.target_size(1) == opt.target_size(2);
            opt.random_size = false;
        else
            opt.random_size = true;
        end
    end
    if sum(size(opt.target_size)) == 2
            opt.random_size = false;
    end
    
    if ~isfield(opt, 'target_plane')
        opt.target_plane = 1i;
    else
        t = opt.target_plane;
        if t<minnode(3) || t>maxnode(3)
            warning('options.target_plane is outside the model!');
            eidors_msg('mk_GREIT_model: Resorting to default target_plane');
            opt.target_plane = 1i;
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

function do_unit_test
   imdl =  mk_GREIT_model( 'c=1;h=2;r=.08;ce=16;bg=1;st=1;me=1;nd', 0.25, 10);
   i2 =  mk_common_model('d2c2', 16);
   i2.fwd_model.normalize_measurements= 1;
   i2 = select_imdl( i2, {'Basic GN dif','Choose NF=1.0'});
   test_performance({mk_common_gridmdl('GREITc1'), imdl, i2});
