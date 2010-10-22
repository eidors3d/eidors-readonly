% TODO: Implement control segments in the bit that writes the file.

function [fmdl,mat_idx] = ng_mk_extruded_model(shape, elec_pos, elec_shape, ...
    extra_ng_code)
% NG_MAKE_EXTRUDED_MODEL: create extruded models using netgen
% [fmdl,mat_idx] = ng_mk_extruded_models(trunk_shape, elec_pos, ...
%                 elec_shape, extra_ng_code);
% INPUT:
% trunk_shape = { height,[x,y],curve_type,maxsz}
%   height      -> if height = 0 calculate a 2D model
%   [x,y]       -> N-by-2 CLOCKWISE list of points defining the 2D shape
%   curve_type  -> 1 - interpret as vertices (default)
%                  2 - interpret as splines with de Boor points at even 
%                  indices (legacy)
%                  3 - interpolate points (piecewise polynomial
%                  interpolation). Syntax [3, N] also specifies the number
%                  of samples to create.
%                  4 - interpolate points with Fourier descriptor. Syntax 
%                  [4, N] also specifies the number of samples to create.
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
if isstr(shape) && strcmp(shape,'UNIT_TEST'); fmdl = do_unit_test; return; end

if nargin < 4; extra_ng_code = {'',''}; end

fnstem = 'tmp';%tempname;
geofn= [fnstem,'.geo'];
ptsfn= [fnstem,'.msz'];
meshfn= [fnstem,'.vol'];

[tank_height, tank_shape, tank_maxh, is2D] = parse_shape(shape);
[elecs, centres] = parse_elecs(elec_pos, elec_shape, tank_shape, tank_height, is2D );
write_geo_file(geofn, ptsfn, tank_height, tank_shape, ...
               tank_maxh, elecs, extra_ng_code);
           
call_netgen( geofn, meshfn, ptsfn);

[fmdl,mat_idx] = ng_mk_fwd_model( meshfn, centres, 'ng', [],0.01,...
    @ng_remove_electrodes);

%delete(geofn); delete(meshfn); delete(ptsfn); % remove temp files

hold on
plot(centres(:,1),centres(:,2),'sk')
for i = 1:size(elecs,2)
    dirn = elecs(i).normal;
    quiver(centres(i,1),centres(i,2),dirn(1),dirn(2),'k');
end
hold off
axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TANK SHAPE (struct):
%         vertices: [Nx2] 
%             size: 0.5* length of the diagonal of the containing rectangle
%     edge_normals: [Nx2]
%       vertex_dir: [Nx2] direction of vertex movement when scaling
%         centroid: [x y]
%   vertices_polar: [Nx2] Phi, r
%           convex: [N] boolean array indicating external angle >= 180 deg
%       curve_type: One of three values
%                   1 - Normal, each point is a vertex
%                   2 - Spline, all even points are de Boor points
%                   3 - Same as 1 but will be converted to smooth 
% 
function [tank_height, tank_shape, tank_maxh, is2D] = parse_shape(shape)
    % parses the shape input

    %defaults
    is2D = false;
    tank_maxh = 0;
    tank_shape = [];
    tank_shape.curve_type = 1;
    curve_info = [];

    if iscell(shape) && length(shape)>2
        tank_height = shape{1};
        if ~iscell(shape{2}) || numel(shape{2}) == 1
            points = shape{2};
        else
            c = shape{2};
            points = c{1};
            tank_shape.additional_shapes = c(2:end);
        end
        tank_shape.curve_type = shape{3};
        if max(size(tank_shape.curve_type)) > 1
            curve_info = tank_shape.curve_type;
            tank_shape.curve_type = curve_info(1);
        end
%         if length(shape) > 2
%             tank_height = shape{1};
%         end
        if length(shape) > 3
            tank_maxh = shape{4};
        end
    else
        points = shape;
    end
    

    
    
    spln_sgmnts = zeros(size(points)); %default
    if tank_shape.curve_type == 2
        [points, spln_sgmnts] = remove_linear_control_points(points);
    end
    
    % piecewise polynomial interpolation
    if tank_shape.curve_type == 3 
        if ~isempty(curve_info)
            n_samples = curve_info(2);
        else
            n_samples = 50;
        end
        points = interpolate_shape(points, n_samples);
        spln_sgmnts = zeros(size(points)); % now needs to be bigger
    end
    
    % Fourier descriptor interpolation
    if tank_shape.curve_type == 4
        if ~isempty(curve_info)
            n_samples = curve_info(2);
        else
            n_samples = 50;
        end
        points = fourier_interpolate_shape(points, n_samples);
        spln_sgmnts = zeros(size(points)); % now needs to be bigger
    end
    
    tank_shape.centroid = calc_centroid(points);
    tank_shape.spln_sgmnts = spln_sgmnts;

    tank_shape.vertices = points;
    % diagonal of the containing rectangle:
    tank_shape.size = 0.5 * sqrt(sum(range(points).^2));
    
    if tank_height==0
        is2D = 1;
        % Need some width to let netgen work, but not too much so
        % that it meshes the entire region
        tank_height = tank_shape.size/5; % initial estimate
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
    
    
    tmp = [];
    polar = zeros(size(points));
    for i = 1:length(points)
        tmp = points(i,:) - tank_shape.centroid;
        [polar(i,1) polar(i,2)]  = cart2pol(tmp(1),tmp(2));
    end
    tank_shape.vertices_polar = polar;
    
    tank_shape.convex = calc_convex(tank_shape.vertices);
    
    % debug plot
    pts = edges./2 + points;
    plot(tank_shape.vertices(:,1),tank_shape.vertices(:,2),'-o'); hold on;
    plot(tank_shape.centroid(:,1),tank_shape.centroid(:,2),'+');
    plot(tank_shape.vertices(:,1)+0.05*tank_shape.vertex_dir(:,1),...
        tank_shape.vertices(:,2)+0.05*tank_shape.vertex_dir(:,2),'ro-')
    quiver(pts(:,1),pts(:,2),tank_shape.edge_normals(:,1),tank_shape.edge_normals(:,2));
    hold off
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
% points - [2N x 2] defined vertices (odd) and control points (even)
% OUTPUT:
% points   - same as points but with linear control points removed
% spln_sgmnts - boolean array indicating which segments are splines
function [points, spln_sgmnts] = remove_linear_control_points(points)
n_points = length(points);
points(end+1,:) = points(1,:);
spln_sgmnts(1:(n_points/2)) = 1;
for i = 1:2:n_points
    a = (points(i+1,:) - points(i,:));
    a = a/norm(a);
    b = (points(i+2,:) - points(i,:));
    b = b/norm(b); 
    if a(1) == b(1) && a(2) == b(2)
        spln_sgmnts(i/2 + 0.5) = 0;
    end    
end
idx = find(spln_sgmnts==0) * 2;
points(idx,:) = [];
points(end,:) = [];

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
% points - [N x 2] defined vertices
% OUTPUT:
% out    - interpolated vertices    
function out = interpolate_shape(points, n_points)
% Quadratic spline interpolation of the points provided.


[pp m] = piece_poly_fit(points);
p = linspace(0,1,n_points+1)'; p(end) = [];
[th xy] = piece_poly_fit(pp,0,p);
tmp = [th xy];
tmp = sortrows(tmp,-1);% ensure clockwise direction
xy = tmp(:,2:3);

out = xy + repmat(m, [n_points,1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
% points - [N x 2] defined vertices
% OUTPUT:
% out    - interpolated vertices    
function out = fourier_interpolate_shape(points, n_points)
% Quadratic spline interpolation of the points provided.


pp = fourier_fit(points);
p = linspace(0,1,n_points+1)'; p(end) = [];
xy = fourier_fit(pp,p);
% [th r] = cart2pol(xy);
% tmp = [th xy];
% tmp = sortrows(tmp,-1);% ensure clockwise direction
% xy = tmp(:,2:3);

out = xy;% + repmat(m, [n_points,1]);


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
        % if the edge directions are the same (accounting for round-off
        % error), return the edge normal.
        if isempty(find(abs(dir1 - dir2) > 1e-14))
            out(i,:) = edgnrm(i,:);
        else
            A = [dir1' , -dir2'];
            u = (p2 - p1)';
            x = A\u;
            out(i,:) = x(1) * dir1 + p1 - points(i,:);
        end
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

function out = calc_convex(verts)
% Returns an array of boolean values for every vertex, true if the external
% angle at this vertex is greater or equal to 180 degrees, false otherwise.
% This marks the vertices which upset the convexity of the polygon and
% require special treatment.

n_verts = size(verts,1);
tmp = [verts(end,:); verts; verts(1,:)];
verts = tmp;

for i = 2:n_verts+1
    v1 = [verts(i-1,:) - verts(i,:), 0];
    v2 = [verts(i+1,:) - verts(i,:), 0];
    cp = cross(v1',v2);
    out(i-1) = cp(3) >= 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% ELECTRODE POSITIONS:
%  elec_pos = [n_elecs_per_plane,(0=equal angles,1=equal dist),z_planes] 
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
      elec_pos(:,3) = hig/2;
   end
   
   % temp fix
   rad = tank_shape.size;

   % It never makes sense to specify only one elec
   % So elec_pos means the number of electrodes in this case
   if size(elec_pos,1) == 1
       % Parse elec_pos = [n_elecs_per_plane,z_planes] 
      n_elecs= elec_pos(1); % per plane
      switch elec_pos(2)
          case 0
              th = linspace(0,2*pi, n_elecs+1)'; th(end)=[];
          case 1
              if tank_shape.curve_type == 4
                  pp = fourier_fit(tank_shape.vertices);
                  p = linspace(0,1,n_elecs+1)'; p(end) = [];
                  th = fourier_fit(pp,p);
              else
                  pp= piece_poly_fit(tank_shape.vertices);
                  p = linspace(0,1,n_elecs+1)'; p(end) = [];
                  th = piece_poly_fit(pp,0,p);
              end
      end

      on_elecs = ones(n_elecs, 1);
      el_th = []; 
      el_z  = []; 
      for i=3:length(elec_pos)
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
%        switch tank_shape.curve_type
%            case 1
               [centres(i,1:2), normal] = calc_elec_centre(tank_shape, el_th(i));
%            case{2, 3}
%                [centres(i,1:2), normal] = calc_elec_centre_spline(tank_shape, el_th(i));
%            otherwise
%                error('Unknown curve type');
%        end
       centres(i,3) = el_z(i);
       elecs(i).pos  = centres(i,:);
       elecs(i).normal = normal;
   end

   if n_elecs == 0
      elecs= struct([]); % empty
   end

   
   
    function [pos, normal] = calc_elec_centre(tank_shape, th)
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
        
        normal = tank_shape.edge_normals(edg_no,:);
        
        
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
        if DAC ~= 0
            ratio = AB * sin(DAB) / (AC * sin(DAC));
            pos = vert(v1,:) + ( ratio / (1 + ratio) ) * (vert(v2,:) - vert(v1,:));
        else
            pos = vert(v2,:);
        end
   
   
   function [pos, normal] = calc_elec_centre_spline(tank_shape, th)
        % The calculation proceeds by finding a common point between a line
        % from the centroid outwards and the equation of the relevant
        % quadratic spline segment defined using 3 control points
        
        % make sure th is between -pi and pi
        if th > pi; th = th - 2*pi; end 
        
        vert_pol = tank_shape.vertices_polar; %[th, r]
        
        % The number of vertices must be even, but just in case...
        if mod(size(vert_pol,1),2)
            error(['The number of points must be even. '...
                'One de Boor control point for every vertex']);
        end
        
        % if the curve is defined as splines, every second point is not
        % actually a vertex. We remove them.
        if tank_shape.curve_type == 2 || tank_shape.curve_type == 3
            vert_pol(2:2:end,:) = [];
        end
      
        n_vert = size(vert_pol,1);
   
        vert_pol = [vert_pol , (1:n_vert)']; %excludes control points
        % Re-order the vertices -pi to pi. Now counter-clockwise
        vert_pol = sortrows(vert_pol,1); 
        % find the edge on which the electrode lies. Edge 1 is between
        % vertices 1 and 2.
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
        
        % correcting for the control points
        v1 = 2 * v1 - 1;
        v2 = 2 * v2 - 1;
        
        % the spline goes from point P0 to P2 such that P1-P0 is a tangent
        % at P0 and P2-P1 is a tangent at P2 
        C = tank_shape.centroid;
        P0 = tank_shape.vertices(v1,:) - C;
        P1 = tank_shape.vertices(v1+1,:) - C; % control point
        P2 = tank_shape.vertices(v2,:) - C;
        
        
        % find the gradient of the line from centroid to electrode center:
        [x, y] = pol2cart(th, 1);
        % FIXME: This doesn't crash only because of round-off errors.
        g = y/x;
        % (because we subtracted the centroid from the vertices, the line
        % passes through the origin now)
        
        % the spline is f(t) = (1-t)^2 * P0 + 2t(1-t)P1 + t^2 * P2
        % which can also be expressed as
        f = @(t) (P2 - 2*P1 + P0)*t^2 + 2*(P1 - P0)*t + P0;
        % and it's derivative:
        df = @(t) 2*(P2 - 2*P1 + P0)*t + 2*(P1 - P0);
        % to find the value of t for which the line cross, we substitute
        % p0 = y0-ax0 for P0 and so on. 
        p0 = P0(2) - g * P0(1);
        p1 = P1(2) - g * P1(1);
        p2 = P2(2) - g * P2(1);
        
        % thus we have a quadratic equation a*t^2 + b*t + c = 0 where
        a = (p2 - 2*p1 + p0);
        b = 2* (p1 - p0);
        c = p0;
        
        if abs(a) < 1e-10
            t = -c/b;
            pos = f(t) + C;
            tmp = df(t);
            normal = [-tmp(2), tmp(1)] / sqrt(sum(tmp.^2));
            return;
        end
        
        % the determinant is
        D = b^2 - 4*a*c;
        
        % find the roots
        if D == 0
            t = -b / (2 * a);

        elseif D > 0
            t1 = (-b - sqrt(D) ) / (2 * a);
            t2 = (-b + sqrt(D) ) / (2 * a);
            if t1 >= 0 && t1 <= 1
                t = t1;
            else
                t = t2;
            end
        else
            error('Something went wrong, cannot place electrode on spline');
        end
        
        pos = f(t) + C;
        tmp = df(t);
        normal = [-tmp(2), tmp(1)]/ sqrt(sum(tmp.^2));

   
   

function elec = elec_spec( row, is2D, hig, rad )
  if     is2D
     if length(row)>=2 && row(2) == -1 % Point electrodes
        % Create rectangular electrodes with bottom, cw point where we want
        elec.shape = 'P' ;
        if length(row)>=3 && row(3) > 0
           elec.dims  =  row(3);
        else
           elec.dims  =  rad; % Make big if unspecified
        end
     else
        % create circular electrodes for now, rectangular not yet supported
%         elec.shape = 'C';
%         elec.dims = row(1);
        elec.shape = 'R';
        elec.dims  = [row(1),hig];
     end
  else
     if length(row)<2 || row(2) == 0 % Circular electrodes 
        elec.shape = 'C';
        elec.dims  = row(1);
     elseif row(2) == -1 % Point electrodes
        % Create rectangular electrodes with bottom, cw point where we want
        elec.shape = 'P'; 
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

  
  
  
  
  
  
  
  
  
  
  
  
  
  
function write_geo_file(geofn, ptsfn, tank_height, tank_shape, ...
                        tank_maxh, elecs, extra_ng_code)
    fid=fopen(geofn,'w');
    write_header(fid,tank_height,tank_shape,tank_maxh,extra_ng_code);

    n_verts = size(tank_shape.vertices,1);
    n_elecs = length(elecs);
    %  elecs(i).pos   = [x,y,z]
    %  elecs(i).shape = 'C' or 'R'
    %  elecs(i).dims  = [radius] or [width,height]
    %  elecs(i).maxh  = '-maxh=#' or '';
    %  elecs(i).edg_no = i (index of the edge on which the electrode lies)
    pts_elecs_idx = [];
    %^keyboard
    for i=1:n_elecs
        name = sprintf('elec%04d',i);
        pos = elecs(i).pos;
        dirn = elecs(i).normal;
        switch elecs(i).shape
            case 'C'
                write_circ_elec(fid,name, pos, dirn,  ...
                    elecs(i).dims, tank_shape.centroid, elecs(i).maxh);
            case 'R'
                write_rect_elec(fid,name, pos, dirn,  ...
                    elecs(i).dims, tank_shape.size/10, elecs(i).maxh);
                %        case 'P'
                %          pts_elecs_idx = [ pts_elecs_idx, i];
                %          continue; % DON'T print solid cyl

            otherwise; error('unknown electorde shape');
        end
        %       fprintf(fid,'solid cyl%04d = trunk   and %s; \n',i,name);
    end
    fprintf(fid,'solid trunk = bound');
    if isfield(tank_shape,'additional_shapes')
         for i = 1:length(tank_shape.additional_shapes)
             fprintf(fid,' and not add_obj%04d',i);
         end
    end
    fprintf(fid,';\n');
    
    if isfield(tank_shape,'additional_shapes')
        for i = 1:length(tank_shape.additional_shapes)
            fprintf(fid,'solid add_obj%04dc = add_obj%04d',i,i);
            for j = (i+1):length(tank_shape.additional_shapes)
                fprintf(fid,' and not add_obj%04d',j);
            end
            fprintf(fid,[' and plane(0,0,0;0,0,-1)\n' ...
                '      and  plane(0,0,%6.2f;0,0,1)'],tank_height);
            fprintf(fid,';\n');
        end
    end
    
    if tank_maxh ~= 0
        fprintf(fid,'tlo trunk -transparent -maxh=%f;\n',tank_maxh);
    else
        fprintf(fid,'tlo trunk -transparent;\n');
    end
    if ~isempty(extra_ng_code{1})
        fprintf(fid,'tlo %s -col=[0,1,0];\n',extra_ng_code{1});
    end

    if isfield(tank_shape,'additional_shapes')
         for i = 1:length(tank_shape.additional_shapes)
             fprintf(fid,'tlo add_obj%04dc -col=[0,1,0];\n',i);
         end
    end

    for i=1:n_elecs
        if any(i == pts_elecs_idx); continue; end
        fprintf(fid,'tlo elec%04d -col=[1,0,0];\n',i);
    end


    fclose(fid); % geofn

   
   
   function write_header(fid,tank_height,tank_shape,maxsz,extra)
   if maxsz==0; 
      maxsz = '';
   else
      maxsz = sprintf('-maxh=%f',maxsz);
   end

   if ~isempty( extra{1} )
      extra{1} = [' and not ',extra{1}];
   end

   
   fprintf(fid,'#Automatically generated by ng_mk_extruded_model\n');
   fprintf(fid,'algebraic3d\n');
   fprintf(fid,'%s\n',extra{2}); % Define extra stuff here
   
   fprintf(fid,'curve3d extrsncurve=(2; 0,0,0; 0,0,%6.2f; 1; 2,1,2);\n', ...
       tank_height+1);


   write_curve(fid,tank_shape,'outer', 1.15);
   write_curve(fid,tank_shape,'inner', 0.99);
   write_curve(fid,tank_shape,'surf', 1);
   
    fprintf(fid,['solid bound= plane(0,0,0;0,0,-1)\n' ...
                '      and  plane(0,0,%6.2f;0,0,1)\n' ...
                '      and  extrusion(extrsncurve;surf;0,1,0)'...
                '%s %s;\n'],tank_height,extra{1},maxsz);
            
   fprintf(fid,['solid inner_bound= plane(0,0,0;0,0,-1)\n' ...
                '      and  plane(0,0,%6.2f;0,0,1)\n' ...
                '      and  extrusion(extrsncurve;inner;0,1,0)'...
                '%s %s;\n'],tank_height,extra{1},maxsz);

   fprintf(fid,['solid outer_bound= plane(0,0,0;0,0,-1)\n' ...
                '      and  plane(0,0,%6.2f;0,0,1)\n' ...
                '      and  extrusion(extrsncurve;outer;0,1,0)'...
                '%s %s;\n'],tank_height,extra{1},maxsz);
           
   % EVERYTHING below this line assumes additional shapes are defined
   if ~isfield(tank_shape, 'additional_shapes'), return, end
   
   for i = 1:length(tank_shape.additional_shapes)
       name_curve = sprintf('add_curve%04d',i); 
       write_curve(fid,tank_shape.additional_shapes{i},name_curve);
       name_obj = sprintf('add_obj%04d',i); 
       fprintf(fid,['solid %s= plane(0,0,%6.2f;0,0,-1)\n' ...
           '      and  plane(0,0,%6.2f;0,0,1)\n' ...
           '      and  extrusion(extrsncurve;%s;0,1,0)'...
           '%s %s;\n'],name_obj,-i,tank_height+i,name_curve,extra{1},maxsz);
   end
                   
        
   function write_curve(fid, tank_shape, name, scale)
        if nargin <4
            scale = 1;
        end
       
        is_struct = isstruct(tank_shape);
        if ~is_struct
            vertices = tank_shape;
            STRUCT = false;
            if scale ~= 1
                warning('Scale is ignored when second input is an array');
                scale = 1;
            end
        elseif scale ~= 1
            vertices = tank_shape.vertices + ...
                (scale-1)*tank_shape.vertex_dir*tank_shape.size;
        else
            vertices = tank_shape.vertices;
        end
       n_vert = size(vertices,1);
       
       fprintf(fid,'curve2d %s=(%d; \n', name, n_vert);
       
       for i = 1:n_vert
           % because of the definitions of the local axis in extrusion, the
           % x coordinate has to be multiplied by -1. This assures the
           % object appears at the expected coordinates. To maintain
           % clockwise order (required by netget) the vertices are printed
           % in the opposite order.
           fprintf(fid,'       %6.4f, %6.4f;\n',[-1 1].*vertices(n_vert-i+1,:));
           %             fprintf(fid,'       %6.2f, %6.2f;\n',vertices(i,:));
       end
       if is_struct
           spln_sgmnts = tank_shape.spln_sgmnts;
       else
           spln_sgmnts = zeros(max(size(vertices)));
       end
       n_sgmnts = length(spln_sgmnts);
       fprintf(fid,'       %d;\n',n_sgmnts);
       cv = 1; %current vertex
       for i = 1:n_sgmnts
           if spln_sgmnts(i)
               if i == n_sgmnts
                  fprintf(fid,'       %d, %d, %d, %d );\n\n\n', 3, cv,cv+1, 1);
               else
                   fprintf(fid,'       %d, %d, %d, %d; \n', 3, cv, cv+1, cv+2);
               end
               cv = cv + 2;
           else
               if i == n_sgmnts
                   fprintf(fid,'       %d, %d, %d );\n\n\n', 2, cv, 1);
               else
                   fprintf(fid,'       %d, %d, %d; \n', 2, cv, cv+1);
               end
               cv = cv + 1;
           end
       end
       
       
function write_circ_elec(fid,name,c, dirn, rd, centroid, maxh)
% writes the specification for a netgen cylindrical rod on fid,
%  named name, centerd on c,
% in the direction given by vector dirn, radius rd 
% direction is in the xy plane

    % the direction vector
    dirn(3) = 0; dirn = dirn/norm(dirn);

    fprintf(fid,'solid %s  = ', name);
    fprintf(fid,['  outer_bound and not inner_bound and '...
        'cylinder(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f;%6.3f) '...
        'and plane(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f) '...
        'and not bound %s;\n'], ...
        c(1)-dirn(1),c(2)-dirn(2),c(3)-dirn(3),...
        c(1)+dirn(1),c(2)+dirn(2),c(3)+dirn(3), rd, ...
        centroid(1), centroid(2), 0, -dirn(1), -dirn(2), dirn(3),maxh);

function write_rect_elec(fid,name,c, dirn,wh,d,maxh)
% writes the specification for a netgen cuboid on fid, named name, centerd on c,
% in the direction given by vector dirn,
% hw = [height, width]  and depth d
% direction is in the xy plane
   w = wh(1); h= wh(2);
   dirn(3) = 0; dirn = dirn/norm(dirn);
   dirnp = [-dirn(2),dirn(1),0];
   dirnp = dirnp/norm(dirnp);

   bl = c - (d/2)* dirn + (w/2)*dirnp - [0,0,h/2];
   tr = c + (d/2)* dirn - (w/2)*dirnp + [0,0,h/2];
   fprintf(fid,'solid %s  = outer_bound and not inner_bound and', name);
   fprintf(fid,' plane (%6.3f,%6.3f,%6.3f;0, 0, -1) and\n', ...
           bl(1),bl(2),bl(3));
   fprintf(fid,' plane(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f) and\n', ...
           bl(1),bl(2),bl(3),-dirn(1),-dirn(2),0);
   fprintf(fid,' plane(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f) and\n', ...
           bl(1),bl(2),bl(3),dirnp(1),dirnp(2),0);
   fprintf(fid,' plane(%6.3f,%6.3f,%6.3f;0, 0, 1) and\n', ...
           tr(1),tr(2),tr(3));
   fprintf(fid,' plane(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f) and\n', ...
           tr(1),tr(2),tr(3),dirn(1),dirn(2),0);
   fprintf(fid,' plane(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f  )%s\n', ...
           tr(1),tr(2),tr(3),-dirnp(1),-dirnp(2),0,maxh);
   fprintf(fid,' and not bound;\n');
    
function [srf,vtx,fc,bc,simp,edg,mat_ind] = ng_remove_electrodes...
    (srf,vtx,fc,bc,simp,edg,mat_ind, N_elec)
% NG_REMOVE_ELECTRODES: cleans up matrices read from a *.vol file
% [srf,vtx,fc,bc,simp,edg,mat_ind]= ng_remove_electrodes...
%     (srf,vtx,fc,bc,simp,edg,mat_ind, N_elec)
%
% Used to clean up external objects used to force electrode meshing in
% ng_mk_extruded_model.
%


% total objects:
N_obj = max(mat_ind);

% electrode simps:
e_simp_ind = mat_ind > (N_obj - N_elec);

in = unique(simp(~e_simp_ind,:));
out = unique(simp(e_simp_ind,:));
boundary = intersect(in,out);
out = setdiff(out,boundary);

ext_srf_ind = ismember(srf,out);
ext_srf_ind = ext_srf_ind(:,1) | ext_srf_ind(:,2) | ext_srf_ind(:,3);

srf(ext_srf_ind,:) = [];
bc(ext_srf_ind,:) = [];
fc(ext_srf_ind,:) = [];
simp = simp(~e_simp_ind,:);
mat_ind = mat_ind(~e_simp_ind);

% fix bc:
n_unique = numel(unique(bc));
missing = setdiff(1:n_unique, unique(bc));
spare = setdiff(unique(bc), 1:n_unique); 
for i = 1:length(missing)
    bc( bc==spare(i) ) = missing(i);
end

% fic vtx:
v = 1:size(vtx,1);
unused_v = setdiff(v, union(unique(simp),unique(srf))); 
v(unused_v) = [];
for i = 1: length(vtx)
%     simp_ind = find(simp == i);
%     srf_ind = find( srf == i);
    new_v_ind = find(v == i);
    simp( simp == i ) = new_v_ind; 
    srf( srf  == i ) = new_v_ind;
end
vtx(unused_v,:) = [];

function [fmdl, mat_idx] = do_unit_test
fmdl = [];
mat_idx = [];
    a = [
   -0.8981   -0.7492   -0.2146    0.3162    0.7935    0.9615    0.6751    0.0565   -0.3635   -0.9745
    0.1404    0.5146    0.3504    0.5069    0.2702   -0.2339   -0.8677   -0.6997   -0.8563   -0.4668 ]';
% [fmdl, mat_idx] = ng_mk_extruded_model({2,{a,0.5*a,0.2*a},1},[16,0,1],[0.01]);
load CT2

% [fmdl, mat_idx] = ng_mk_extruded_model({150,flipud(trunk),1},[16,0,75],[0.01]);

[fmdl, mat_idx] = ng_mk_extruded_model({2,{trunk/100, lung_heart_dep/100, heart/100},1},[16,1,1],[0.1]);
% img = mk_image( fmdl, 1);
%  img.elem_data(mat_idx{2}) = 1.1; 

% trunk = [    -4    -2     2     4     4     2    -2    -4
%               2     4     4     2    -2    -4    -4    -2]';
% heart_lung = [    -2    -1    -0.8  0.8  1     2     2    -2
%                    1     2     1.8  1.8  2     1    -2    -2]';    
% heart = [    -1    -1     1     1
%               0     2     2     0]';

% [fmdl, mat_idx] = ng_mk_extruded_model({2,{trunk, heart_lung, heart},1},[16,1,1],[0.1]);

 figure, show_fem( fmdl );
 
%%
xx=[
  -88.5777  -11.4887    4.6893   49.8963  122.7033  150.3033  195.5103 249.7573 ...
  258.8013  279.7393  304.9623  309.2443  322.0923  337.7963  340.6503 348.2633 ...
  357.3043  358.7333  361.5873  364.9183  365.3943  366.3453  366.3453 365.3943 ...
  362.5393  351.5943  343.5053  326.8513  299.2503  288.3073  264.9923 224.0703 ...
  206.4633  162.6833  106.5313   92.2543   57.5153    7.0733   -8.6297 -42.4167 ...
  -90.9547 -105.7057 -134.2577 -178.0367 -193.2647 -222.7687 -265.5957 -278.9197 ...
 -313.1817 -355.5337 -363.6237 -379.3267 -397.8857 -400.7407 -401.6927 -398.8377 ...
 -395.0307 -384.0867 -368.3837 -363.6247 -351.7277 -334.1217 -328.4117 -314.1357 ...
 -291.2947 -282.7297 -267.0257 -236.5707 -221.8187 -196.5977 -159.4807 -147.5837];

yy=[
 -385.8513 -386.8033 -386.3273 -384.8993 -368.7193 -353.9673 -323.0363 -283.5403 ...
 -274.9743 -254.0363 -225.4843 -217.8703 -187.4153 -140.7813 -124.6013  -86.0573 ...
  -38.4703  -29.4273   -9.9173   21.0137   32.4347   53.3727   83.8257   93.3437 ...
  114.7587  149.0237  161.8717  187.5677  222.3037  231.3447  247.5237  267.5087 ...
  271.3177  277.0297  281.3127  279.4097  274.6507  273.2227  276.5547  284.6447 ...
  295.1127  297.4927  301.7757  304.1557  302.2537  297.4947  287.5017  282.2667 ...
  259.9017  225.6387  213.7427  185.6677  141.4127  125.2337   88.5917   34.8187 ...
   17.6897  -22.2803  -73.6723  -85.0923 -117.9263 -163.6083 -176.4573 -205.9613 ...
 -245.9343 -256.4023 -275.4373 -304.9403 -315.4083 -332.0623 -352.0473 -355.3783];

a = [xx; yy]';
a = flipud(a);
th=linspace(0,2*pi,33)'; th(end)=[];
%a=[sin(th)*0.3,cos(th)];


%      fmdl = ng_mk_extruded_model({300,a,[3,50]},[16,150,1],[0.01]);
%      figure
%      show_fem(fmdl);
