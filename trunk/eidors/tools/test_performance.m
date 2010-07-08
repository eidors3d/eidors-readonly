function [r, params] =  test_performance( imdls, fmdl );
% TEST_PERFORMANCE: test of difference reconstruction algorithms
%   in terms of the performance parameters defined in GREIT
%
% Input:
%   imdls = cell array of inverse models to test on
%   fmdl  = fmdl of the shape to test on (default is 16 electrode tank).
%
% Model assumes the tank is elliptical with an x,y centre at 0,0
%
% [r, params] =  test_performance( imdls );
%  r =  radii of test objects
%  param = [AR, PE, RES, SD, RNG, NOI]
%
% Example
%   i_bp = mk_common_gridmdl('backproj');
%   i_gp = mk_common_gridmdl('GREITc1');
%   test_performance( { i_gp, i_gr })

% (c) 2010 Andy Adler. Licenced under GPL v2 or v3
% $Id$

if isstr(imdls) && strcmp(imdls,'UNIT_TEST'); do_unit_test; return; end

if nargin==1
   fmdl = ng_mk_cyl_models([2,1,0.05],[16,1],[0.05]); 
   fmdl.stimulation = mk_stim_patterns(16,1,[0,1],[0,1],{},1);
   imgs= mk_image( fmdl, 1);
end

Nsim = 100;
r =  linspace(0,0.9,Nsim);
ctr = [0,0, mean(fmdl.nodes(:,3))];  % Assume x,y centre is zero
maxx = max(abs(fmdl.nodes(:,1) - ctr(1)));
maxy = max(abs(fmdl.nodes(:,2) - ctr(2)));

th=  r*4321; % want object to jump around in radius
xyzr = [maxx*r.*cos(th); maxy*r.*sin(th); ctr(3)*ones(1,Nsim); 
        0.05/mean([maxx,maxy])*ones(1, Nsim)];

[vh,vi] = simulate_movement(imgs, xyzr);

for i= 1:length(imdls);
   imgr = inv_solve(imdls{i}, vh, vi);
   imgr.calc_colours.npoints = 64;
   param_GR = eval_GREIT_fig_merit(imgr, xyzr);
   pnoise = calc_noise_params( imdls{i}, vh, vi );
   params{i} = [param_GR; pnoise];
end

clf; Nparams = 6;
for j=1:Nparams;
   p = [];
   for i= 1:length(params);
      p = [p, params{i}(j,:)'];
   end

%   subplot(5,1,j);
   hig = 0.95/Nparams;
   axes('position',[0.1,0.05 + hig*(Nparams-j),.88, hig]);
   plot(r, p);
   if j==1;     axis([0,0.9,0,2.1]);        ylabel('AR');
   elseif j==2; axis([0,0.9,-0.16,0.16]);   ylabel('PE');
   elseif j==3; axis([0,0.9,0,0.41]);       ylabel('RES');
   elseif j==4; axis([0,0.9,0,0.31]);       ylabel('SD');
   elseif j==5; axis([0,0.9,0,0.61]);       ylabel('RNG');
   elseif j==6; set(gca,'YScale','log');    ylabel('NF');
   end
   if j<Nparams; set(gca,'XTickLabel',[]);end
end

function do_unit_test;
% Reconstruct GREIT Images
imdl_gr = mk_common_gridmdl('GREITc1');

% Reconstruct backprojection Images
imdl_bp = mk_common_gridmdl('backproj');

% Reconstruct GN Images
imdl_gn = select_imdl( mk_common_model('d2c2', 16), {'Basic GN dif','Choose NF=0.5'});

test_performance( { imdl_gr, imdl_bp, imdl_gn } );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vh,vi,param]=test_performance_old(fmdl,imdl,r,N)

debug = false;

%% 1. Figure out the limits of fmdl
boundary = fmdl.boundary;
b_ver = unique(boundary);
boundary = fmdl.nodes(b_ver,:);
min_b = min(boundary)+repmat(0.5*r,1,3);
max_b = max(boundary)-repmat(0.5*r,1,3);
vol = get_elem_volume(fmdl);

%% 2. create N samples within the volume
rand('state', sum(100*clock));

n_elems = size(fmdl.elems,1);


c2f = eidors_obj('get-cache', {fmdl, N, r}, 'samples');
xyzr_pt = eidors_obj('get-cache', {fmdl, N, r}, 'points');

count = 0; %count points already created

if isempty(c2f) || isempty(xyzr_pt);
    disp('Generating sample points.');
    c2f = sparse(n_elems,N);
    while count < N
        n_pts = ceil(1.5 * (N - count) ); % more than needed to limit iterations
        vec = rand_pt ( repmat(min_b,n_pts,1), repmat(max_b,n_pts,1) );
        map = mk_c2f_circ_mapping(fmdl, [vec repmat(r,n_pts,1)]');
        frac = vol'*map/(4/3*pi*r^3);
        idx = find(frac > 0.9, N - count, 'first');
        
        xyzr_pt = [xyzr_pt, vec(idx,:)'];
        newcount = count + length(idx);
        c2f(:, (count+1):newcount ) = map(:,idx);
        count = newcount;
    end
    eidors_obj('set-cache', {fmdl, N, r}, 'samples', c2f);
    eidors_obj('set-cache', {fmdl, N, r}, 'points', xyzr_pt );
    disp('Done.');
end
 


% x=110;y=150;r=20;
% x=110;y=80;r=20;
% map = mk_c2f_circ_mapping(fmdl, [x,y,r]');
% img= mk_image(map);img.fwd_model = fmdl;
% phi= linspace(0,2*pi,20);xr=r*cos(phi)+x; yr=r*sin(phi)+y;
% show_fem(img,[0,1,1]); hold on; plot(xr,yr,'r');hold off
% vol = get_elem_volume(fmdl); vol'*map/(pi*r^2)

%% 3. reconstruct the images
% create a homegenous image
   elem_data = ones(size(fmdl.elems,1),1);
   img = mk_image(fmdl, elem_data);
   disp('Calculating Jacobian. This may take a (long) while.');
   J = calc_jacobian( img );
   disp('Done.');
   vh= fwd_solve(img);
   vh = vh.meas;
   vi= vh*ones(1,N) + J*c2f;

%% 4. calculate parameters


   
   img = inv_solve(imdl,vh,vi);
   img.calc_colours.npoints=32;
   imgs = calc_slices(img);

       param= [calc_noise_params(imdl, vh, vi); ... % noise parameters
           GREIT_sim_params(  imgs, xyzr_pt)];            % image parameters
       
       r = sqrt(sum(xyzr_pt(1:2,:).^2));
       
       names = {'Noise','Amplitude','Position Error','Resolution', ...
           'Shape Deformation', 'Ringing'};
       for j=1:size(param,1)
           figure
           scatter(r, param(j,:));
           title(names{j});
       end

       
       
       
function x = rand_pt(x_min, x_max)
x = rand(size(x_min)) .* (x_max - x_min) + x_min;
