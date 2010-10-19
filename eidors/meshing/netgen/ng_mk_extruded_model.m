% TODO: Implement control segments in the bit that writes the file.

function [fmdl,mat_idx] = ng_mk_extruded_model(shape, elec_pos, elec_shape, ...
    extra_ng_code)
% NG_MAKE_EXTRUDED_MODEL: create extruded models using netgen
% [fmdl,mat_idx] = ng_mk_extruded_models(trunk_shape, elec_pos, ...
%                 elec_shape, extra_ng_code);
% INPUT:
% trunk_shape = { height,[x,y],curve_type,maxsz}
%   height (OPT)-> if height = 0 calculate a 2D model
%   [x,y]       -> N-by-2 CLOCKWISE list of points defining the 2D shape
%   curve_type  -> (default = 1) 1 - interpret as vertices 
%                                2 - interpret as splines with de Boor
%                                points at even indices
%                                3 - create a smooth shape from points
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
if nargin < 4; extra_ng_code = {'',''}; end

fnstem = 'new';%tempname;
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

    if iscell(shape) && length(shape)>2
        tank_height = shape{1};
        points = shape{2};
        tank_shape.curve_type = shape{3};
%         if length(shape) > 2
%             tank_height = shape{1};
%         end
        if length(shape) > 3
            tank_maxh = shape{4};
        end
    else
        points = shape;
    end
    
    tank_shape.centroid = calc_centroid(points);
    
    
    spln_sgmnts = zeros(size(points)); %default
    if tank_shape.curve_type == 2
        [points, spln_sgmnts] = remove_linear_control_points(points);
    end 
    
    if tank_shape.curve_type == 3
        points = calc_spline_from_shape(points,tank_shape.centroid);
    end
    
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
% out    - [2N x 2] vertices (odd) and control points (even)    
function out = calc_spline_from_shape(points, centroid)
% Quadratic spline interpolation of the points provided.

% The problem is to obtain control points for the spline segments, such
% that Netgen can draw a smooth curve.
% Every spline segment S_i(t) can be expressed us:
%   S_i(t) = (1-t^2)*P_i + 2t(1-t)*C_i + t^2*P_(i+1)
% OR
%   S_i(t) = (P_(i+1) - 2C_1 + P_i)t^2 + 2(C_i - P_i)t + P_i
% where Pi is the i-th defined vertex, and Ci is the control point of
% the i-th spline segment.
% The first derivative can be expressed as:
%   S'_i(t) = 2t*P_(i+1) + 2(t-1)*P_i + 2(1-2t)*C_i
%           = 2(P_(i+1) - 2C_1 + P_i)t + 2(C_i - P_i)

% New approach:
% 0. Subtract the centroid and convert to polar coords
% 1. Calculate a cubic spline interpolation
% 2. Find inflection points of the spline segments
% 3. If inflection point is found within the segment, put a new node there
% 4. Use the gradients at each node to find the control points for
% quadratic splines

% 0. Subtract the centroid and convert to polar coords
% 1. Calculate a cubic spline interpolation
% 2. Find inflection points of the spline segments
% 3. If inflection point is found within the segment, put a new node there
% 4. Use the gradients at each node to find the control points for
% quadratic splines

    if nargin < 2
        centroid = [0 , 0];
    end

    % 0. Subtract the centroid and convert to polar coords
    n_points = size(points,1);
    points = points - repmat(centroid, [n_points,1]);
    [ppoints(:,1), ppoints(:,2)] = cart2pol(points(:,1), points(:,2));

    % 1. Calculate a cubic spline interpolation
    % First, close the loop:
    ppoints = sortrows(ppoints,1);
    r = [ppoints(:,2); ppoints(1,2)];
    rho = [ppoints(:,1) ; 2*pi+ppoints(1,1)];
    df = (r(2) - r(end-1)) / ( (2*pi - rho(2)) - rho(end-1));
    pp=spline(rho,[df; r; df]);

    % 2. Find inflection points of the spline segments
    % For a polynomial ax^3 + bx^2 + cx + d, the second derivative is
    % 6ax + 2b.
    % thus, to find if an inflection point occurs within any of the spline
    % segments, we will check if x = -b/3a belongs to that segment

    count = 1;
    newrho = [];
    newr = [];
    newgradients = [];
    for i = 1:n_points

        df = @(z) 3*pp.coefs(i,1)*z^2 + 2*pp.coefs(i,2)*z + pp.coefs(i,3);
        ddf = @(z) 6*pp.coefs(i,1)*z + 2*pp.coefs(i,2);

        gradients(i) = df(0);

        if pp.coefs(i,1) == 0,
            continue; %cannot have an inflection point
        end
        x = -pp.coefs(i,2)/(3* pp.coefs(i,1));
        if x > 0 && x < (pp.breaks(i+1) - pp.breaks(i))

            %check if the second derivative changes sign
            if sign(ddf(x-0.001)) ~= sign(ddf(x+0.001))
    %             disp(['Inflection point in segment ' num2str(i)]);
                newrho(count) = pp.breaks(i)+x;
                newr(count) = ppval(pp,newrho(count));
                newgradients(count) = df(x);
                count = count+1;
    %         else
    %             disp(['Non-Inflection point in segment ' num2str(i)]);
            end

        end
    end
    % figure
    % plot(ppoints(:,1), ppoints(:,2), 'p');
    % hold on

    % 3. If inflection point is found within the segment, put a new node there
    if ~isempty(newrho)
        ppoints = [ppoints; newrho' newr'];
        ppoints = [ppoints [gradients newgradients]'];
    else 
        ppoints = [ppoints gradients'];
    end
    ppoints = sortrows(ppoints,1);

    % 4. Using the gradients at each node, compute the control points

    % figure
    % polar(ppoints(:,1), ppoints(:,2), 'o');
    % hold on

    ppoints(end+1,:) = ppoints(1,:);
    for i = 1: size(ppoints,1)-1
        % gradinent in polar coordinates:
        a1 = ppoints(i,3);
        % offset:
        b1 = ppoints(i,2) - a1*ppoints(i,1);
        rho1 = ppoints(i,1);
        r1   = ppoints(i,2);
        % the function r = a1 * Phi + b1 in polar coordinates is a curve in
        % cartesian coordinates. Find a tangent at a vertex (b1):
        % dy = (a*rho + b)cos(rho) + a*sin(rho)
        dy1 = ( a1*rho1 + b1) * cos(rho1) + a1*sin(rho1);
        % dx = (a*rho +b)(-sin(rho)) + a*cos(rho)
        dx1 = ( a1*rho1 + b1) * (- sin(rho1)) + a1*cos(rho1);
        % the vertex in cartesian coords
        [x1, y1] = pol2cart(rho1,r1);


        % the same for the next vertex:
        a2 = ppoints(i+1,3);
        b2 = ppoints(i+1,2) - a2*ppoints(i+1,1);
        rho2 = ppoints(i+1,1);
        r2   = ppoints(i+1,2);
        dy2 = ( a2*rho2 + b2) * cos(rho2) + a2*sin(rho2);
        dx2 = ( a2*rho2 + b2) * (- sin(rho2)) + a2*cos(rho2);
        [x2, y2] = pol2cart(rho2,r2);

        A = [dx1, -dx2; dy1, -dy2];
        u = [x2 - x1 ; y2 - y1];
        t = A\u;

    %     x = ppoints(i,1):0.01:ppoints(i+1,1);
    %     yl = a1*x+b1;
    %     yr = a2*x+b2;
    %     polar(x, yl, '-b')
    %     polar(x, yr, '-k');
    %     
    %     v = 0: 0.01 : 1;
    %     ylc = dy1*v + y1;
    %     xlc = dx1*v + x1;
    %     plot(xlc, ylc);

    %     control(i,1) = x1 + t(1)*dx1;
    %     control(i,2) = y1 + t(2)*dy1;
        % store as polar coords (for sorting later on)
        [control(i,1) , control(i,2)] = cart2pol(x1 + t(1)*dx1, y1 + t(1)*dy1);
    end

    %  polar(control(:,1), control(:,2), 'ro');

    ppoints(end,:) = [];

    % a1 = ppoints(end,3);
    % b1 = ppoints(end,2) - a1*(ppoints(end,1));
    % a2 = ppoints(1,3);
    % b2 = ppoints(1,2) - a2*(2*pi + ppoints(1,1));
    % 
    % i = size(ppoints,1);
    % control(i,1) = (b2 - b1 ) / (a1 - a2);
    % control(i,2) = a1 * control(i,1) + b1;


    % integrate the control points between the defined vertices:
    n = 2*size(ppoints,1);
    new = zeros(n,2);
    new(1:2:end,:) = ppoints(:,1:2);
    new(2:2:end,:) = control;
    new = flipud(new);
    %move first point to the end
    new = circshift(new,[-1 0]);

    % convert to cartesian coords:
    [out(:,1), out(:,2)] = pol2cart(new(:,1),new(:,2));
    % shift back to centroid
    out = out + repmat(centroid, [n,1]);




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
       switch tank_shape.curve_type
           case 1
               [centres(i,1:2), normal] = calc_elec_centre(tank_shape, el_th(i));
           case{2, 3}
               [centres(i,1:2), normal] = calc_elec_centre_spline(tank_shape, el_th(i));
           otherwise
               error('Unknown curve type');
       end
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
                %        case 'R'
                %          write_rect_elec(fid,name, pos, pos,  ...
                %                elecs(i).dims, tank_radius, elecs(i).maxh);
                %        case 'P'
                %          pts_elecs_idx = [ pts_elecs_idx, i];
                %          continue; % DON'T print solid cyl

            otherwise; error('huh? shouldn`t get here');
        end
        %       fprintf(fid,'solid cyl%04d = trunk   and %s; \n',i,name);
    end
    %
    if tank_maxh ~= 0
        fprintf(fid,'tlo bound -transparent -maxh=%f;\n',tank_maxh);
    else
        fprintf(fid,'tlo bound -transparent;\n');
    end


    for i=1:n_elecs
        if any(i == pts_elecs_idx); continue; end
        fprintf(fid,'tlo elec%04d -col=[1,0,0];\n',i);
    end

    if ~isempty(extra_ng_code{1})
        fprintf(fid,'tlo %s -col=[0,1,0];\n',extra_ng_code{1});
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
           
                   
        
   function write_curve(fid, tank_shape, name, scale)
        if nargin <4
            scale = 1;
        end
       
        if scale ~= 1
            vertices = tank_shape.vertices + (scale-1)*tank_shape.vertex_dir;
        else
            vertices = tank_shape.vertices;
        end
       n_vert = size(tank_shape.vertices,1);
       fprintf(fid,'curve2d %s=(%d; \n', name, n_vert);
       
       for i = 1:n_vert
           % because of the definitions of the local axis in extrusion, the
           % x coordinate has to be multiplied by -1. This assures the
           % object appears at the expected coordinates. To maintain
           % clockwise order (required by netget) the vertices are printed
           % in the opposite order.
           fprintf(fid,'       %6.2f, %6.2f;\n',[-1 1].*vertices(n_vert-i+1,:));
           %             fprintf(fid,'       %6.2f, %6.2f;\n',vertices(i,:));
       end
       spln_sgmnts = tank_shape.spln_sgmnts;
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
     
function [srf,vtx,fc,bc,simp,edg,mat_ind] = ng_remove_electrodes...
    (srf,vtx,fc,bc,simp,edg,mat_ind, N_elec)
% NG_REMOVE_ELECTRODES: cleans up matrices read from a *.vol file
% [srf,vtx,fc,bc,simp,edg,mat_ind]= ng_remove_electrodes...
%     (srf,vtx,fc,bc,simp,edg,mat_ind, N_elec)
%
% Used to clean up external objects used to force electrode meshing in
% ng_mk_extruded_model.
%
% (C) Bartlomiej Grychtol, 2010. Licenced under GPL v2 or v3
% $Id$

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



     