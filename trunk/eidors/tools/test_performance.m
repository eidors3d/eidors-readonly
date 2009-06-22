function [vh,vi,param]=test_performance(fmdl,imdl,r,N)

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
   img = mk_image(elem_data,'name', 'image name', 'fwd_model', fmdl);
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