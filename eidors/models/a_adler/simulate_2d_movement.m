function [vh,vi,xyr_pt]= simulate_2d_movement( n_sims, fmdl, rad_pr, movefcn )
% SIMULATE_2D_MOVEMENT simulate rotational movement in 2D
% [vh,vi,xyr_pt]= simulate_2d_movement( n_points, model, rad_pr, movefcn )
%
% the target starts at (rad_pr(1),0) and rotates around 
%  clockwise
% 
%   rad_pr = [path_radius, target_radius] = [2/3, .05] (default)
% 
%   n_points = number of points to simulate (default = 200)
%
%   model = fwd_model to simulate 
%         (default use internal, or if model= []);
%
%   movefcn = 1 (Default)  radial motion where the target starts
%     at (rad_pr(1),0) and rotates clockwise
%   movefcn = 2 movement from centre to outward on to
%     at (rad_pr(1),0) on x-axis
%
%   movefcn = FUCN NAME or FUNC HANDLE
%      the function must accept the following parameters
%      [xp,yp] = movefcn(f_frac, radius);
%
% OUTPUT:
%   vh - homogeneous measurements M x 1
%   vi - target simulations       M x n_points
%   xyr_pt - x y and radius of each point 3 x n_points

%
% $Id: simulate_2d_movement.m,v 1.16 2008-07-17 18:54:37 aadler Exp $

if nargin <1
   n_sims = 200;
end

if nargin<2 || isempty(fmdl) % create our own fmdl
   n_circles = 36;
   n_elec= 16;
   fmdl= mk_fwd_model(n_circles, n_elec);
end

if nargin<3 || isempty(rad_pr)
   radius= 2/3;
   rp= .05;
else
   radius= rad_pr(1);
   rp=     rad_pr(2);
end

if nargin<4
   movefcn = 1;
end
if isnumeric(movefcn)
   if     movefcn==1
      movefcn = @rotation_path;
   elseif movefcn==2
      movefcn = @straight_out;
   else
      error('value of movefcn not understood');
   end
else
   % assume movefcn is a function 
end

    n_elems= size(fmdl.elems,1);
    img= eidors_obj('image','simulate_movement', ...
                    'fwd_model', fmdl, ...
                    'elem_data', ones(n_elems,1) );
    vh= fwd_solve(img);

    np= 512;
    np= 256;
    [x,y]=meshgrid( linspace(-1,1,np), linspace(-1,1,np) );
    [eptr,vol]= img_mapper2(fmdl.nodes', fmdl.elems', np, np);

    % there is a faster way to do this with sparse, but it is confusing
%   eptr_n= zeros(n_elems,1);
%   for i=1:n_elems; eptr_n(i) = sum( eptr(:)==i ); end
    eptr_n= sparse( eptr(:)+1,1,1, n_elems+1, 1);
    eptr_n= full(eptr_n(2:end));
    
    target_conductivity= .2;
 
    for i=1:n_sims
       f_frac= (i-1)/n_sims;
       fprintf('simulating %d / %d \n',i,n_sims);

      [xp,yp]= feval(movefcn, f_frac, radius);
       xyr_pt(:,i)= [xp;-yp;rp]; % -y because images and axes are reversed

       ff= find( (x(:)-xp).^2 + (y(:)-yp).^2 <= rp^2 )';
       obj_n= sparse( eptr(ff)+1,1,1, n_elems+1, 1);
       obj_n= full(obj_n(2:end));

%      img.elem_data= 1 + target_conductivity * (obj_n./eptr_n);
       img.elem_data= 1 + target_conductivity * (obj_n./vol);

       vi(i)= fwd_solve( img );
%show_fem(img);drawnow; keyboard
    end

% convert to data matrix
vi= [vi(:).meas]; 
vh= [vh(:).meas];

%   movefcn = 1 (Default)  rotational motion where the target starts
%     at (rad_pr(1),0) and rotates clockwise
% calculate x,y position of point, given f_frac of path
function [xp,yp] = rotation_path(f_frac, radius);
   xp= radius * cos(f_frac*2*pi);
   yp= radius * sin(f_frac*2*pi);

function [xp,yp] = straight_out(f_frac, radius);
   xp= radius*f_frac;
   yp= 0;

% THis is the code copied from calc_slices
% Search through each element and find the points which
% are in that element
function [EPTR,VOL]= img_mapper2(NODE, ELEM, npx, npy );
  xmin = min(NODE(1,:));    xmax = max(NODE(1,:));
  xmean= mean([xmin,xmax]); xrange= xmax-xmin;

  ymin = min(NODE(2,:));    ymax = max(NODE(2,:));
  ymean= mean([ymin,ymax]); yrange= ymax-ymin;

  range= max([xrange, yrange]);
  [x y]=meshgrid( ...
      linspace( xmean - range*0.50, xmean + range*0.50, npx ), ...
      linspace( ymean + range*0.50, ymean - range*0.50, npy ) );
  v_yx= [-y(:) x(:)];
  turn= [0 -1 1;1 0 -1;-1 1 0];
  EPTR=zeros(npy,npx);
  % for each element j, we get points on the simplex a,b,c
  %   area A = abc
  %   for each candidate point d,
  %      area AA = abd + acd + bcd
  %      d is in j if AA = A
  e= size(ELEM,2);
  VOL= zeros(e,1);
  for j= 1:e
    % calculate area of three subtrianges to each candidate point.
    xy= NODE(:,ELEM(:,j))';
    % come up with a limited set of candidate points which
    % may be within the simplex
    endr=find( y(:)<=max(xy(:,2)) & y(:)>=min(xy(:,2)) ...
             & x(:)<=max(xy(:,1)) & x(:)>=min(xy(:,1)) );
    % a is determinant of matrix [i,j,k, xy]
    a= xy([2;3;1],1).*xy([3;1;2],2)- xy([3;1;2],1).*xy([2;3;1],2);
    VOL(j)= abs(sum(a));

    aa= sum(abs(ones(length(endr),1)*a'+ ...
                v_yx(endr,:)*xy'*turn)');
    endr( abs( ( VOL(j)-aa ) ./ VOL(j) ) >1e-8)=[];
    EPTR(endr)= j;
  end %for j=1:ELEM


function mdl_2d= mk_fwd_model(n_circles, n_elec)
    params= mk_circ_tank(n_circles, [], n_elec); 
    n_rings= 1;
    options= {'no_meas_current','no_rotate_meas','do_redundant'};
    [st, els]= mk_stim_patterns(n_elec, n_rings, '{ad}','{ad}', options, 10);
    params.stimulation= st;
    params.meas_select= els;
    params.solve=      'aa_fwd_solve';
    params.system_mat= 'aa_calc_system_mat';
    params.jacobian=   'aa_calc_jacobian';
    params.normalize_measurements= 1;
    params.np_fwd_solve.perm_sym= '{n}';
    mdl_2d   = eidors_obj('fwd_model', params);
