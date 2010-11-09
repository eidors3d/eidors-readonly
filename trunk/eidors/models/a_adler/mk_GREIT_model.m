function imdl= mk_GREIT_model( fmdl, radius, weight, extra_noise )
% MK_GREIT_MODEL: make EIDORS inverse models using the GREIT approach
%   imdl= mk_GREIT_model( fmdl, radius, weight, extra_noise )
%
% Parameters
%   fmdl   - fwd model on which to do simulations, or
%          - string specifying prepackaged models
%
%   radius - requested weighting matrix  (recommend 0.25 for 16 electrodes)
%   weight - weighting matrix (weighting of noise vs signal)
%
%   extra_noise - extra noise samples (such as electrode movement)
%
% NOTE
%   currently extra_noise is not supported
%   currently weighting matrix must be scalar
              
% Example
%   imdl =  mk_GREIT_model( 'c=1;h=2;r=.08;ce=16;bg=1;st=1;me=1;nd', 0.25, 10);

% (C) 2010 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isstr(fmdl) && strcmp(fmdl,'UNIT_TEST'); do_unit_test; return; end

if nargin==4;
  error('mk_GREIT_model: doesn''t currently support extra_noise');
end

if isstr(fmdl)
   imgs = get_prepackaged_fmdls( fmdl );
elseif isfield(fmdl,'type');
  switch fmdl.type
%  if we get a fwd_model, assume uniform conductivity backgnd of 1
    case 'fwd_model'; imgs = mk_image( fmdl, 1);
%  if we get an image, use it. It may have a non-uniform backgnd
    case 'image';     imgs = fmdl; % fmdl was an image
    otherwise; error('unrecognized eidors object');
  end
else
   error('specified parameter must be an object or a string');
end

Nsim = 1000;
[vi,vh,xy]= stim_targets(imgs, Nsim );

RM= calc_GREIT_RM(vh,vi, xy, radius, weight, imgs.fwd_model.normalize_measurements );
imdl = mk_common_gridmdl('b2c', RM);
imdl.fwd_model.normalize_measurements= imgs.fwd_model.normalize_measurements;

%RM= calc_GREIT_RM(vh,vi, xyc, radius weight, normalize)

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

function [vi,vh,xy]= stim_targets(imgs, Nsim );
   fmdl = imgs.fwd_model;
   r =  linspace(0,0.9,Nsim);
   ctr = [0,0, mean(fmdl.nodes(:,3))];  % Assume x,y centre is zero
   maxx = max(abs(fmdl.nodes(:,1) - ctr(1)));
   maxy = max(abs(fmdl.nodes(:,2) - ctr(2)));

   th=  r*4321; % want object to jump around in radius
   xyzr = [maxx*r.*cos(th); maxy*r.*sin(th); ctr(3)*ones(1,Nsim); 
           0.05/mean([maxx,maxy])*ones(1, Nsim)];

   [vh,vi] = simulate_movement(imgs, xyzr);
   xy = xyzr(1:2,:);


% TODO: stim_targets will need to figure out the shape of the imgs.fwd_model
% Calculate the position of the electrodes
   Nelecs = length(imgs.fwd_model.electrode);
   for i=1:Nelecs
     enodesi =     imgs.fwd_model.electrode(:).nodes; 
     elec_loc(i,:) = mean( imgs.fwd_model.nodes( enodesi,:),1 );
   end

% Now, use elec_loc to figure out the shape. We can assume the obj is extruded in z



function do_unit_test
   imdl =  mk_GREIT_model( 'c=1;h=2;r=.08;ce=16;bg=1;st=1;me=1;nd', 0.25, 10);
   i2 =  mk_common_model('d2c2', 16);
   i2.fwd_model.normalize_measurements= 1;
   i2 = select_imdl( i2, {'Basic GN dif','Choose NF=1.0'});
   test_performance({mk_common_gridmdl('GREITc1'), imdl, i2});
