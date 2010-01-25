function tank_shape = ng_mk_extruded_model(shape)
% NG_MAKE_EXTRUDED_MODELS: create extruded models using netgen
% [fmdl,mat_idx] = ng_mk_extruded_models(cyl_shape, elec_pos, ...
%                 elec_shape, extra_ng_code);
% INPUT:
% shape = {[x,y],height,maxsz}
%   [x,y]       -> N-by-2 CLOCKWISE list of points defining the 2D shape
%   height (OPT)-> (default = 1) if height = 0 calculate a 2D model
%   maxsz       -> max size of mesh elems (default = course mesh)

% (C) Bartlomiej Grychtol, 2010. Licenced under GPL v2 or v3
% $Id$

[tank_height, tank_shape, tank_maxh, is2D] = parse_shape(shape);






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
    tank_shape.size = sqrt(sum(range(points).^2));
    
    if tank_height==0
        is2D = 1;
        % Need some width to let netgen work, but not too much so
        % that it meshes the entire region
        tank_height = tank_shape.size/10; % initial extimate
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


