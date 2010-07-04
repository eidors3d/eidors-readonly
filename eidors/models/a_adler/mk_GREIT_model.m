function imdl= mk_GREIT_model( fmdl, radius, weight )
% MK_GREIT_MODEL: make EIDORS inverse models using the GREIT approach

% (C) 2010 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isstr(fmdl) && strcmp(fmdl,'UNIT_TEST'); do_unit_test; return; end

if isstr(fmdl)
  imgs = get_prepackaged_fmdls( fmdl );
else
  imgs = mk_image( fmdl, 1);
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


function do_unit_test
   imdl =  mk_GREIT_model( 'c=1;h=2;r=.08;ce=16;bg=1;st=1;me=1;nd', 0.25, 10);
   i2 =  mk_common_model('d2c2', 16);
   i2.fwd_model.normalize_measurements= 1;
   i2 = select_imdl( i2, {'Basic GN dif','Choose NF=1.0'});
   test_performance({mk_common_gridmdl('GREITc1'), imdl, i2});
