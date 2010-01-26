function tank_shape = ng_mk_extruded_model(shape, elec_pos, elec_shape)
% NG_MAKE_EXTRUDED_MODELS: create extruded models using netgen
% [fmdl,mat_idx] = ng_mk_extruded_models(cyl_shape, elec_pos, ...
%                 elec_shape, extra_ng_code);
% INPUT:
% shape = {[x,y],height,maxsz}
%   [x,y]       -> N-by-2 CLOCKWISE list of points defining the 2D shape
%   height (OPT)-> (default = 1) if height = 0 calculate a 2D model
%   maxsz       -> max size of mesh elems (default = course mesh)
%
% ELECTRODE POSITIONS:
%  elec_pos = [n_elecs_per_plane,z_planes] 
%     OR
%  elec_pos = [degrees,z] centres of each electrode (N_elecs x 2)
%
% ELECTRODE SHAPES::
%  elec_shape = [width,height, {maxsz, pem}]  % Rectangular elecs
%     OR
%  elec_shape = [radius, {0, maxsz, pem} ]  % Circular elecs
%     maxsz  (OPT)  -> max size of mesh elems (default = course mesh)
%     pem  (OPT)  -> 1: Point Electrode Model, 0: Complete Electrode Model (default)
%
% Specify either a common electrode shape or for each electrode

% (C) Bartlomiej Grychtol, 2010. Licenced under GPL v2 or v3
% $Id$

[tank_height, tank_shape, tank_maxh, is2D] = parse_shape(shape);
[elecs, centres] = parse_elecs(elec_pos, elec_shape, tank_shape, tank_height, is2D );

hold on
plot(centres(:,1),centres(:,2),'sk')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TANK SHAPE (struct):
%         vertices: [Nx2] 
%             size: 0.5* length of the diagonal of the containing rectangle
%     edge_normals: [Nx2]
%       vertex_dir: [Nx2] direction of vertex movement when scaling
%         centroid: [x y]
%   vertices_polar: [Nx2] Phi, r
% 
function [tank_height, tank_shape, tank_maxh, is2D] = parse_shape(shape)
    % parses the shape input

    %defaults
    tank_height = 1;
    is2D = false;
    tank_maxh = 0;
    tank_shape = [];

    if iscell(shape) && length(shape)>1
        points = shape{1};
        tank_height = shape{2};
        if length(shape)>2
            tank_maxh = shape{3};
        end
    else
        points = shape;
    end
    
    tank_shape.vertices = points;
    % diagonal of the containing rectangle:
    tank_shape.size = 0.5 * sqrt(sum(range(points).^2));
    
    if tank_height==0
        is2D = 1;
        % Need some width to let netgen work, but not too much so
        % that it meshes the entire region
        tank_height = tank_shape.size/5; % initial extimate
        if tank_maxh>0
            tank_height = min(tank_height,2*tank_maxh);
        end
    end


    tank_shape.edge_normals = [];
    tank_shape.vertex_dir = [];

    tmp = points;
    tmp(end+1,:) = tmp(1,:); %duplicate first vertex at the end;

    edges = diff(tmp,1);
    tmp = [];
    % Normal to vector (x y) is (-y x).
    % It points outward for clockwise definition
    tmp = circshift(edges, [0 1]); %swap coords
    %normalize
    lngth = sqrt(sum(tmp.^2, 2));
    tmp(:,1) = -tmp(:,1) ./ lngth;
    tmp(:,2) = tmp(:,2)  ./ lngth;
    tank_shape.edge_normals = tmp;

    tank_shape.vertex_dir = calc_vertex_dir(points, edges, ...
        tank_shape.edge_normals);
    
    tank_shape.centroid = calc_centroid(points);
    tmp = [];
    polar = zeros(size(points));
    for i = 1:length(points)
        tmp = points(i,:) - tank_shape.centroid;
        [polar(i,1) polar(i,2)]  = cart2pol(tmp(1),tmp(2));
    end
    tank_shape.vertices_polar = polar;
    
    % debug plot
    pts = edges./2 + points;
    plot(tank_shape.vertices(:,1),tank_shape.vertices(:,2),'-o'); hold on;
    plot(tank_shape.centroid(:,1),tank_shape.centroid(:,2),'+');
    plot(tank_shape.vertices(:,1)+0.05*tank_shape.vertex_dir(:,1),...
        tank_shape.vertices(:,2)+0.05*tank_shape.vertex_dir(:,2),'ro-')
    quiver(pts(:,1),pts(:,2),tank_shape.edge_normals(:,1),tank_shape.edge_normals(:,2));
    hold off




function out = calc_vertex_dir(points, edges, edgnrm)
%     calculate the direction of vertex movement if all edges are shifted
%     outwards by 1 unit along their normals:

%     duplicate last edge at the beginning
    edg = [edges(end,:) ; edges];
    edgnrm = [edgnrm(end,:) ; edgnrm];

    out = zeros(size(points));
    for i = 1:length(points)
        p1 = points(i,:) + edgnrm(i,:);
        p2 = points(i,:) + edgnrm(i+1,:);

        dir1(1) = edgnrm(i,2); dir1(2) = -edgnrm(i,1);
        dir2(1) = edgnrm(i+1,2); dir2(2) = -edgnrm(i+1,1);
        A = [dir1' , -dir2'];
        u = (p2 - p1)';
        x = A\u;
        out(i,:) = x(1) * dir1 + p1 - points(i,:);
    end

function out = calc_centroid(points)
% Calculates the centroid of the shape
% The algorithm identifies a middle point M within the shape and then uses it
% to divide the shape into N triangles (N=number of vertices), calculates
% the area and centroid of each traingle, and finally computes the centroid
% of the shape as a mean of the centroids of the individual traingles
% weighted by their area. 

    % it never makes sense to have less than 3 points
    n_points = size(points,1);
    if  n_points == 3
        out = mean(points); % centroid of a triangle
        return
    end

    out = 0;
    pts = [points ; points(1,:)];

    % guess a point in the middle
    m = mean(points);

    if ~inpolygon(m(1),m(2),points(:,1),points(:,2))
        f1 = figure;
        set(f1,'Name', 'Select a point within the shape');
        plot(pts(:,1),pts(:,2));
        m = ginput(1);
        close(f1)
    end

    tmp = 0;
    tot_area = 0;
    for i = 1:n_points
        a = pts(i,:);
        b = pts(i+1,:);
        cntrd = (m + a + b)/3;
        area = 0.5 * abs(det([m 1; a 1; b 1]));
        tmp = tmp + cntrd*area;
        tot_area = tot_area + area;
    end

    out = tmp./tot_area;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% ELECTRODE POSITIONS:
%  elec_pos = [n_elecs_per_plane,z_planes] 
%     OR
%  elec_pos = [degrees,z] centres of each electrode (N_elecs x 2)
%
% ELECTRODE SHAPES::
%  elec_shape = [width,height, {maxsz}]  % Rectangular elecs
%     OR
%  elec_shape = [radius, {0, maxsz} ]  % Circular elecs
%     maxsz  (OPT)  -> max size of mesh elems (default = courase mesh)
% 
% OUTPUT:
%  elecs(i).pos   = [x,y,z]
%  elecs(i).shape = 'C' or 'R'
%  elecs(i).dims  = [radius] or [width,height]
%  elecs(i).maxh  = '-maxh=#' or '';
function [elecs, centres] = parse_elecs(elec_pos, elec_shape, tank_shape, hig, is2D )

   if is2D
      elec_pos(:,2) = hig/2;
   end
   
   % temp fix
   rad = tank_shape.size;

   % It never makes sense to specify only one elec
   % So elec_pos means the number of electrodes in this case
   if size(elec_pos,1) == 1
       % Parse elec_pos = [n_elecs_per_plane,z_planes] 
      n_elecs= elec_pos(1); % per plane
      th = linspace(0,2*pi, n_elecs+1)'; th(end)=[];

      on_elecs = ones(n_elecs, 1);
      el_th = []; 
      el_z  = []; 
      for i=2:length(elec_pos)
        el_th = [el_th; th];
        el_z  = [el_z ; on_elecs*elec_pos(i)];
      end
   else
      el_th = elec_pos(:,1)*2*pi/360;
      el_z  = elec_pos(:,2);
   end
      
   n_elecs= size(el_z,1); 

   if size(elec_shape,1) == 1
      elec_shape = ones(n_elecs,1) * elec_shape;
   end

   for i= 1:n_elecs
     row = elec_shape(i,:); 
     elecs(i) = elec_spec( row, is2D, hig, rad );
   end
   
   
   %centres = [rad*sin(el_th),rad*cos(el_th),el_z];
   for i= 1:n_elecs; 
       centres(i,1:2) = calc_elec_centre(tank_shape, el_th(i));
       centres(i,3) = el_z(i);
       elecs(i).pos  = centres(i,:); 
   end

   if n_elecs == 0
      elecs= struct([]); % empty
   end

   
   
    function [pos, edg_no] = calc_elec_centre(tank_shape, th)
        % The calculation relies on the theorem that if point D lies on a
        % line between B and C, but point A is not on that line, then:
        %   |BD|    |AB| sin(<DAB)
        %   ---- = ---------------
        %   |DC|    |AC| sin(<DAC)
        % Thus, B and C are vertices of our shape, A is its centroid and D
        % is the sought center of the electrode. All quantities on RHS are
        % known.
        
        % make sure th is between -pi and pi
        if th > pi; th = th - 2*pi; end
        
        
        vert_pol = tank_shape.vertices_polar; %[th, r]
        n_vert = size(vert_pol,1);
        vert_pol = [vert_pol , (1:n_vert)'];
        % Re-order the vertices -pi to pi. Now counter-clockwise
        vert_pol = sortrows(vert_pol,1); 
        % find the edge on which the elctrode lies. (Edge 1 is between
        % verticies 1 and 2)
        idx = find(vert_pol(:,1) > th, 1, 'first');
        if isempty(idx); idx = 1; end
        edg_no = vert_pol(idx,3);
        
        v1 = edg_no;
        if edg_no == n_vert % between the last and first vertex    
            v2 = 1;
        else
            v2 = v1+1;
        end
        vert_pol = [];
        
        
        vert_pol = tank_shape.vertices_polar;
        vert = tank_shape.vertices;
        cntr = tank_shape.centroid;
        % position between vertices - see first comment
        AB = sqrt(sum( (vert(v1,:) - cntr).^2 ));
        AC = sqrt(sum( (vert(v2,:) - cntr).^2 ));
        DAB = abs(vert_pol(v1,1)-th); 
        if DAB > pi, DAB = abs( DAB - 2*pi); end; 
        DAC  = abs(vert_pol(v2,1)-th);
        if DAC > pi, DAC = abs( DAC - 2*pi); end;
        ratio = AB * sin(DAB) / (AC * sin(DAC));
        
        
        pos = vert(v1,:) + ( ratio / (1 + ratio) ) * (vert(v2,:) - vert(v1,:));
        
   
   
   
   
   
   

function elec = elec_spec( row, is2D, hig, rad )
  if     is2D
     if length(row)>=2 && row(2) == -1 % Point electrodes
        % Create rectangular electrodes with bottom, cw point where we want
        elec.shape = 'P' 
        if length(row)>=3 && row(3) > 0
           elec.dims  =  row(3);
        else
           elec.dims  =  rad; % Make big if unspecified
        end
     else
        elec.shape = 'R';
        elec.dims  = [row(1),hig];
     end
  else
     if length(row)<2 || row(2) == 0 % Circular electrodes 
        elec.shape = 'C';
        elec.dims  = row(1);
     elseif row(2) == -1 % Point electrodes
        % Create rectangular electrodes with bottom, cw point where we want
        elec.shape = 'P' 
        if length(row)>=3 && row(3) > 0
           elec.dims  =  row(3);
        else
           elec.dims  =  rad; % Make big if unspecified
        end
     elseif row(2)>0      % Rectangular electrodes
        elec.shape = 'R';
        elec.dims  = row(1:2);
     else
        error('negative electrode width');
     end
  end

  if length(row)>=3 && row(3) > 0
     elec.maxh = sprintf('-maxh=%f', row(3));
  else
     elec.maxh = '';
  end

  if length(row)<4 || row(4) == 0
     elec.model = 'cem'; % Complete Electrode Model (CEM)
  else
     elec.model = 'pem'; % Point Electrode Model (PEM)
  end
  %TODO support Shunt Electrode Model (SEM)

