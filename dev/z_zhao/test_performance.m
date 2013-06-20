function [r, params] =  test_performance( imdls, fmdl, obj_type);
% TEST_PERFORMANCE: test of difference reconstruction algorithms
%   in terms of the performance parameters defined in GREIT
%
% Input:
%   imdls = cell array of inverse models to test on
%   fmdl  = fmdl of the shape to test on (default is 16 electrode tank).

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

% 2013-06-19 Zhao
%   obj_type = different test objects simulation
%  1:original; 2: uniform, non-random; 3:original modified for non-circular
%  shape; 4: uniform, non-random modified for non-circular shape (default)

% $Id: test_performance.m 3808 2013-04-12 17:49:37Z bgrychtol $

if isstr(imdls) && strcmp(imdls,'UNIT_TEST'); do_unit_test; return; end

if nargin==1
   fmdl = ng_mk_cyl_models([2,1,0.05],[16,1],[0.05]); 
   fmdl.stimulation = mk_stim_patterns(16,1,[0,1],[0,1],{},1);
   
elseif nargin<3
   obj_type=4;
end
imgs= mk_image( fmdl, 1);

if ~iscell(imdls)
    imdls = {imdls};
end

switch obj_type
    case 1  %original
        Nsim = 100;
        l =  linspace(0,0.9,Nsim);
        ctr = [0,0, mean(fmdl.nodes(:,3))];  % Assume x,y centre is zero
        maxx = max(abs(fmdl.nodes(:,1) - ctr(1)));
        maxy = max(abs(fmdl.nodes(:,2) - ctr(2)));
        r=0.05;
        th=  l*4321; % want object to jump around in radius
        Zpos=ctr(3)*ones(1,Nsim);
        xyzr = [maxx*l.*cos(th); maxy*l.*sin(th); Zpos; r*mean([maxx,maxy])*ones(1, Nsim)];

    case 2   %uniform, non-random; 
        N = 16;
        bnd_nodes = unique(imdls{1}.fwd_model.boundary);
        min_bb = min(imdls{1}.fwd_model.nodes(bnd_nodes,:));
        max_bb = max(imdls{1}.fwd_model.nodes(bnd_nodes,:));
        maxx=max_bb(1);
        maxy=max_bb(2);
        h = mean([min_bb(3),max_bb(3)]);            % better way to calculate Zpos
        r = 0.025 * max(max_bb(1:2) - min_bb(1:2));
        xspace = linspace(min_bb(1),max_bb(1),N);
        yspace = linspace(min_bb(2),max_bb(2),N);
        [X Y] = meshgrid(xspace,yspace);
        imgs.calc_colours.npoints = N;
        M = calc_slices(imgs,1);
        IN = M==1;
        n_px = nnz(IN);
        OUT = isnan(M);
        Zpos=h*ones(1,n_px);
        xyzr = [X(IN)'; Y(IN)'; Zpos; r*ones(1,n_px); ];
    case 3      %original modified for non-circular shape
        Nsim = 100;
        l =  linspace(0,0.9,Nsim);
        ctr = [0,0, mean(fmdl.nodes(:,3))];  % Assume x,y centre is zero
        maxx = max(abs(fmdl.nodes(:,1) - ctr(1)));
        maxy = max(abs(fmdl.nodes(:,2) - ctr(2)));
        r=0.05;
        th=  l*4321; % want object to jump around in radius
        Zpos=ctr(3)*ones(1,Nsim);
        xyzr = [maxx*l.*cos(th); maxy*l.*sin(th); Zpos; r*mean([maxx,maxy])*ones(1, Nsim)];
        N_b=64;
        [xyzr]=del_out_map(imdls{1},imgs,N_b,xyzr);
    case 4      %uniform, non-random modified for non-circular shape (default)
        N = 16;
        bnd_nodes = unique(imdls{1}.fwd_model.boundary);
        min_bb = min(imdls{1}.fwd_model.nodes(bnd_nodes,:));
        max_bb = max(imdls{1}.fwd_model.nodes(bnd_nodes,:));
        maxx=max_bb(1);
        maxy=max_bb(2);
        h = mean([min_bb(3),max_bb(3)]);            % better way to calculate Zpos
        r = 0.025 * max(max_bb(1:2) - min_bb(1:2));
        xspace = linspace(min_bb(1),max_bb(1),N);
        yspace = linspace(min_bb(2),max_bb(2),N);
        [X Y] = meshgrid(xspace,yspace);
        imgs.calc_colours.npoints = N;
        M = calc_slices(imgs,1);
        IN = M==1;
        n_px = nnz(IN);
        OUT = isnan(M);
        Zpos=h*ones(1,n_px);
        xyzr = [0.9*X(IN)'; 0.9*Y(IN)'; Zpos; r*ones(1,n_px); ];
        N_b=64;
        [xyzr]=del_out_map(imdls{1},imgs,N_b,xyzr);
       
end
[vh,vi] = simulate_movement(imgs, xyzr);
for i= 1:length(imdls);
   imgr = inv_solve(imdls{i}, vh, vi);
   imgr.calc_colours.npoints = 64;
   param_GR = eval_GREIT_fig_merit(imgr, xyzr);
   pnoise = calc_noise_figure( imdls{i}, vh, vi );
   params{i} = [param_GR; pnoise];
end
figure
levels =[inf,inf,xyzr(3,:)];
imgr.show_slices.img_cols =10;
show_slices(imgr, levels);
figure

Nparams = 6;
for j=1:Nparams;
   p = [];
   for i= 1:length(params);
      p = [p, params{i}(j,:)'];
   end
   m = mean([maxx,maxy]);
%   subplot(5,1,j);
   hig = 0.95/Nparams;
   axes('position',[0.1,0.05 + hig*(Nparams-j),.88, hig]);
   plot(1:size(xyzr,2), p);
   if j==1;     
%        axis([0,0.9,0,2.1]);        
       ylabel('AR');
   elseif j==2; 
%        axis([0,0.9,-0.16*m,0.16*m]);   
       ylabel('PE');
   elseif j==3; 
%        axis([0,0.9,0,0.41]);       
       ylabel('RES');
   elseif j==4; 
%        axis([0,0.9,0,0.31]);       
       ylabel('SD');
   elseif j==5; 
%        axis([0,0.9,0,0.61]);       
       ylabel('RNG');
   elseif j==6; set(gca,'YScale','log');    
       ylabel('NF');
   end
   if j<Nparams; set(gca,'XTickLabel',[]);end
end

function [xyzr]=del_out_map(i_gr,imgs,N,xyzr)
% N = 16;
bnd_nodes = unique(i_gr.fwd_model.boundary);
min_bb = min(i_gr.fwd_model.nodes(bnd_nodes,:));
max_bb = max(i_gr.fwd_model.nodes(bnd_nodes,:));
h = mean([min_bb(3),max_bb(3)]);            % better way to calculate Zpos
xspace = linspace(min_bb(1),max_bb(1),N);
yspace = linspace(min_bb(2),max_bb(2),N);
[X Y] = meshgrid(xspace,yspace);
imgs.calc_colours.npoints = N;
M = calc_slices(imgs,1);
IN = M==1;
IN=bwperim(IN,8); %contur
xy = [X(IN)'; Y(IN)'];
% Q1 = xy(1,:)>0 & xy(2,:)>0;
% Q2 = xy(1,:)<0 & xy(2,:)>0;
% Q3 = xy(1,:)<0 & xy(2,:)<0;
% Q4 = xy(1,:)>0 & xy(2,:)<0;
out_ind=[];
max_x=max(xy(1,:));
max_y=max(xy(2,:));
min_x=min(xy(1,:));
min_y=min(xy(2,:));
for i=1:length(xyzr)
    % in the first Quadrant, x>0, y>0
    if i==100
        i=i;
    end
    if xyzr(1,i)>max_x || xyzr(1,i)<min_x || xyzr(2,i)>max_y || xyzr(2,i)<min_y
        out_ind=[out_ind i];
    end
    
    out_boundary = sum ((repmat(xyzr(1:2,i),1,size(xy(:,:),2))-xy(:,:)).^2) <= xyzr(4,i)^2;
    if any(out_boundary)  
        out_ind=[out_ind i];
    end
    
%     if xyzr(1,i)>0 && xyzr(2,i)>0
%         out_boundary = sum ((repmat(xyzr(1:2,i),1,size(xy(:,Q1),2))-xy(:,Q1)).^2) <= xyzr(4,i)^2;
%         if any(out_boundary)    
%             out_ind=[out_ind i];
%         end
%         % in the second quadrant, x<0, y>0
%     elseif xyzr(1,i)<0 && xyzr(2,i)>0
%         out_boundary = sum ((repmat(xyzr(1:2,i),1,size(xy(:,Q2),2))-xy(:,Q2)).^2) <= xyzr(4,i)^2;
%         if any(out_boundary)
%             out_ind=[out_ind i];
%         end
%         % in the third quadrant, x<0, y<0
%     elseif xyzr(1,i)<0 && xyzr(2,i)<0
%         out_boundary = sum ((repmat(xyzr(1:2,i),1,size(xy(:,Q3),2))-xy(:,Q3)).^2) <= xyzr(4,i)^2;
%         if any(out_boundary) 
%             out_ind=[out_ind i];
%         end
%         % in the fourth quadrant, x>0, y<0
%     else  % xyzr(1,i)>0 & xyzr(2,i)<0;
%         out_boundary = sum ((repmat(xyzr(1:2,i),1,size(xy(:,Q4),2))-xy(:,Q4)).^2) <= xyzr(4,i)^2;
%         if any(out_boundary)
%             out_ind=[out_ind i];
%         end
%         
%     end
end
xyzr(:,out_ind)=[];


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

       param= [calc_noise_figure(imdl, vh, vi); ... % noise parameters
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
