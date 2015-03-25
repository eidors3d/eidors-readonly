function hh = show_fem(mdl, options)
% SHOW_FEM: show the EIDORS3D finite element model
% hh = show_fem(mdl, options)
% mdl is an EIDORS3D 'model' or 'image' structure
% hh = handle to the plotted model
%
% options may be specified by a list (for compatibility purposes)
%
% options specifies a set of options
%   options(1) => show colourbar
%   options(2) => show numbering on electrodes
%   options(3) => number elements (==1) or nodes (==2);
%
% for detailed control of colours, use the following fields (default values
% are indicated in parenthesis)
%
%     options.viewpoint.az (-37.5)
%     options.viewpoint.el (30.0)
% 
%     options.edge.color ([0 0 0])
%     options.edge.width (1)
%     options.edge.significant.show (true)
%     options.edge.significant.color ([0 0 0])
%     options.edge.significant.width (2)
%     options.edge.significant.angle (45)
%     options.edge.significant.viewpoint_dependent.show (true)
%     options.edge.significant.viewpoint_dependent.color ([0 0 0])
%     options.edge.significant.viewpoint_dependent.width (2)
%     options.edge.significant.viewpoint_dependent.callback (true)
%  
%     options.colorbar.show (false)
% 
%     options.electrode.number.show (false)
%     options.electrode.number.font_size (10)
%     options.electrode.number.font_weight ('bold')
%     options.electrode.number.color ([.6 0 0])
%     options.electrode.number.background_color ('none')
%     options.electrode.number.scale_position (1.15)
% 
%     options.element.number.show (false)
%     options.element.number.font_size (7)
%     options.element.number.font_weight ('normal')
%     options.element.number.color ([0 0 0])
%     options.element.number.background_color ('none')
%
%     options.node.number.show (false)
%     options.node.number.font_size (7)
%     options.node.number.font_weight ('normal')
%     options.node.number.color ([0 0 0.5])
%     options.node.number.background_color ([1 1 1])
%

% (C) 2015 Herve Gagnon. License: GPL version 2 or version 3

% TESTS
switch nargin
    case 0
        error('Insufficient parameters for show_fem');
    case 1
        if ischar(mdl) && strcmp(mdl,'UNIT_TEST'); 
            do_unit_test; 
            return; 
        end
        if (isstruct(mdl) && isfield(mdl, 'show_fem'))
            options = mdl.show_fem;
        else
            options = [];
        end
    case 2
    otherwise
        error('Too many parameters for show_fem');
end

[img, mdl, opts] = proc_params(mdl, options);

if ~ishold
    cla;
end

mdl = find_sub_elements(mdl);

% Convert nodes to 3D if necessary.
mdl.nodes = mdl.nodes*eye(size(mdl.nodes, 2), 3);

hh = draw_all(img, mdl, opts);

% Setup callback function.
if (opts.edge.significant.viewpoint_dependent.callback)
    h3d = rotate3d;
    set(gca, 'UserData', struct('name', 'show_fem_data', 'img', img, 'mdl', mdl, 'opts', opts));
    set(h3d, 'ActionPostCallback' , @refresh_current_axis)
end

if (nargout == 0)
    clear hh; 
end

function [img, mdl, opts] = proc_params(mdl, src_opts)

    % Assign default viewpoint option values.
    if size(mdl.nodes,2) == 2
       opts.viewpoint.az = 0;
       opts.viewpoint.el = 90;
    else
       opts.viewpoint.az = -37.5;
       opts.viewpoint.el =  30.0;
    end

    % Assign default edge option values.
    opts.edge.color   = [0 0 0];
    opts.edge.width   = 1;
    opts.edge.significant.show  = true;
    opts.edge.significant.color = [0 0 0];
    opts.edge.significant.width = 2;
    opts.edge.significant.angle = 45;
    opts.edge.significant.viewpoint_dependent.show     = true;
    opts.edge.significant.viewpoint_dependent.color    = [0 0 0];
    opts.edge.significant.viewpoint_dependent.width    = 2;
    opts.edge.significant.viewpoint_dependent.callback = true;
 
    % Assign default colorbar option values.
    opts.colorbar.show = false;

    % Assign default electrode option values.
    opts.electrode.number.show             = false;
    opts.electrode.number.font_size        = 10;
    opts.electrode.number.font_weight      = 'bold';
    opts.electrode.number.color            = [.6 0 0];
    opts.electrode.number.background_color = 'none';
    opts.electrode.number.scale_position   = 1.15;

    % Assign default element option values.
    opts.element.number.show             = false;
    opts.element.number.font_size        = 7;
    opts.element.number.font_weight      = 'normal';
    opts.element.number.color            = [0 0 0];
    opts.element.number.background_color = 'none';

    % Assign default node option values.
    opts.node.number.show             = false;
    opts.node.number.font_size        = 7;
    opts.node.number.font_weight      = 'normal';
    opts.node.number.color            = [0 0 0.5];
    opts.node.number.background_color = [1 1 1];
    
    if (nargin == 2)
        if (isstruct(src_opts))
            opts = copy_field(opts, src_opts, 'viewpoint', 'az');
            opts = copy_field(opts, src_opts, 'viewpoint', 'el');
            opts = copy_field(opts, src_opts, 'edge', 'color');
            opts = copy_field(opts, src_opts, 'edge', 'width');
            
            if (isfield(src_opts, 'edge'))
                opts.edge = copy_field(opts.edge, src_opts.edge, 'significant', 'show');
                opts.edge = copy_field(opts.edge, src_opts.edge, 'significant', 'color');
                opts.edge = copy_field(opts.edge, src_opts.edge, 'significant', 'width');
                opts.edge = copy_field(opts.edge, src_opts.edge, 'significant', 'angle');
                if (isfield(src_opts.edge, 'significant'))
                    opts.edge.significant = copy_field(opts.edge.significant, src_opts.edge.significant, 'viewpoint_dependent', 'show');
                    opts.edge.significant = copy_field(opts.edge.significant, src_opts.edge.significant, 'viewpoint_dependent', 'color');
                    opts.edge.significant = copy_field(opts.edge.significant, src_opts.edge.significant, 'viewpoint_dependent', 'width');
                    opts.edge.significant = copy_field(opts.edge.significant, src_opts.edge.significant, 'viewpoint_dependent', 'callback');
                end
            end

            opts = copy_field(opts, src_opts, 'colorbar', 'show');

            if (isfield(src_opts, 'electrode'))
                opts.electrode = copy_field(opts.electrode, src_opts.electrode, 'number', 'show');
                opts.electrode = copy_field(opts.electrode, src_opts.electrode, 'number', 'font_size');
                opts.electrode = copy_field(opts.electrode, src_opts.electrode, 'number', 'font_weight');
                opts.electrode = copy_field(opts.electrode, src_opts.electrode, 'number', 'color');
                opts.electrode = copy_field(opts.electrode, src_opts.electrode, 'number', 'background_color');
                opts.electrode = copy_field(opts.electrode, src_opts.electrode, 'number', 'scale_position');                
            end
            
            if (isfield(src_opts, 'element'))
                opts.element = copy_field(opts.element, src_opts.element, 'number', 'show');
                opts.element = copy_field(opts.element, src_opts.element, 'number', 'font_size');
                opts.element = copy_field(opts.element, src_opts.element, 'number', 'font_weight');
                opts.element = copy_field(opts.element, src_opts.element, 'number', 'color');
                opts.element = copy_field(opts.element, src_opts.element, 'number', 'background_color');               
            end

            if (isfield(src_opts, 'node'))
                opts.node = copy_field(opts.node, src_opts.node, 'number', 'show');
                opts.node = copy_field(opts.node, src_opts.node, 'number', 'font_size');
                opts.node = copy_field(opts.node, src_opts.node, 'number', 'font_weight');
                opts.node = copy_field(opts.node, src_opts.node, 'number', 'color');
                opts.node = copy_field(opts.node, src_opts.node, 'number', 'background_color');               
            end
            
            
        else % Support options from old version of show_fem for compatibility.
            % fill in default options
            optionstr = zeros(1, 100);
            optionstr(1:length(src_opts)) = src_opts;

            opts.colorbar.show         = optionstr(1);
            opts.electrode.number.show = optionstr(2);
            switch optionstr(3)
                case 1
                    opts.element.number.show = 1;
                case 2
                    opts.node.number.show = 1;
                case 3
                    opts.element.number.show = 1;
                    opts.node.number.show = 1;
            end
        end
    end

    % Display first image if several images are available.
    if (numel(mdl) > 1)
        eidors_msg('warning: show_fem only shows first image',1);
        mdl = mdl(1);
    end
 
    % if we have an img input, define mdl
    if strcmp(mdl(1).type , 'image' )
        img = mdl;
        mdl = img.fwd_model;
    else
        img = [];
    end

function refresh_current_axis(obj, evd)
UserData = get(gca, 'UserData');
if (all(isfield(UserData, {'name', 'img', 'mdl', 'opts'})) && ...
        strcmp(UserData.name, 'show_fem_data'))
    [az, el] = view(gca);
    if (az ~= UserData.opts.viewpoint.az || el ~= UserData.opts.viewpoint.el)
        UserData.opts.viewpoint.az = az;
        UserData.opts.viewpoint.el = el;
        cla;
        draw_all(UserData.img, UserData.mdl, UserData.opts);
        set(gca, 'UserData', UserData);
    end
end

function hh = draw_all(img, mdl, opts)

    hh = draw_fem(img, mdl, opts);

    % Number elements if necessary.
    if (opts.element.number.show)
        draw_numbers(interp_mesh(mdl), opts.element.number);
    end

    % Number nodes if necessary.
    if (opts.node.number.show)
        draw_numbers(mdl.nodes, opts.node.number);
    end

    % Number electrodes if necessary.
    if (opts.electrode.number.show && isfield(mdl, 'electrode'))
        n_elec = numel(mdl.electrode);
        mesh_center = ones(n_elec, 1)*mean(mdl.nodes, 1);

        elec_pos = zeros(n_elec, size(mesh_center, 2));

        for i = 1:n_elec
            if (isfield(mdl.electrode(i), 'nodes'))
                elec_pos(i, :) = mean(mdl.nodes(mdl.electrode(i).nodes, :), 1);
            elseif (isfield(mdl.electrode(i), 'pos'))
                elec_pos(i, :) = mean(mdl.electrode(i).pos, 1)*eye(size(...
                                                        elec_pos(i, :), 2), 3);
            end
        end

        draw_numbers(opts.electrode.number.scale_position*(elec_pos - mesh_center) + mesh_center, opts.electrode.number);
    end

    if (size(mdl.elems, 2) == 4)  
        view(opts.viewpoint.az, opts.viewpoint.el);
    end

    axis image;
    
function hh = draw_fem(img, mdl, opts)
    
    % Assign colors and alpha to elements.
    if (size(mdl.elems, 2) == 4)
        if ~isempty(img)
            elems = mdl.sub_elements;
            electrode_field = 'sub_elements';
        else
            elems = mdl.boundary;
            electrode_field = 'boundary';
        end
        triangle_alpha = 0.5*ones(size(elems, 1), 1);
        triangle_edge_color = 'none';
        triangle_edge_width = 0.5;
        
        % Color mesh edges.
        edge_edges = mdl.boundary_edges;
        edge_width = ones(size(edge_edges, 1), 1)*opts.edge.width;
        edge_color = ones(size(edge_edges, 1), 1)*opts.edge.color;
        
        if opts.edge.significant.show
            if opts.edge.significant.viewpoint_dependent.show 
                % Highlight profile of boundary according to viewing angle.
                T = viewmtx(opts.viewpoint.az, opts.viewpoint.el);
                transformed_boundary_normal_vector = mdl.boundary_normal_vector*T(1:3,1:3)';
                boundary_edges_idx2 = (transformed_boundary_normal_vector(mdl.edge2boundary(:, 1), 3).*transformed_boundary_normal_vector(mdl.edge2boundary(:, 2), 3) <= 0);

                edge_width(boundary_edges_idx2)    = opts.edge.significant.viewpoint_dependent.width;
                edge_color(boundary_edges_idx2, :) = ones(sum(boundary_edges_idx2), 1)*opts.edge.significant.viewpoint_dependent.color;
            end

            % Highlight significant edges where boundary shape bends more than
            % 45 degrees.
            boundary_edges_idx = dot(mdl.boundary_normal_vector(...
                                     mdl.edge2boundary(:, 1), :), ...
                                     mdl.boundary_normal_vector(...
                                     mdl.edge2boundary(:, 2), :), 2) ...
                                     <= cosd(opts.edge.significant.angle);
            edge_width(boundary_edges_idx)    = opts.edge.significant.width;
            edge_color(boundary_edges_idx, :) = ones(sum(boundary_edges_idx), 1)*opts.edge.significant.color;
        end
    else
        elems = mdl.elems;
        electrode_field = 'boundary';
        triangle_alpha = ones(size(elems, 1), 1);
        triangle_edge_color = opts.edge.color;
        triangle_edge_width = opts.edge.width;
        
        % Highlight mesh boundary.
        if opts.edge.significant.show
            edge_edges = mdl.boundary;
            edge_width = opts.edge.significant.width*ones(size(edge_edges, 1), 1);
            edge_color = ones(size(edge_edges, 1), 1)*opts.edge.significant.color;
        else
            edge_edges = [];
        end
    end
    
    if ~isempty(img)
        if (size(mdl.elems, 2) == 4)
            triangle_color = calc_colours(mdl.element2sub_element*...
                                get_img_data(img), img, opts.colorbar.show);
                            
            factor = 3/8;
            
            alpha_map = zeros(size(colormap, 1), 1);
            alpha = linspace(0, 1, round(size(colormap, 1)*factor))';
            alpha_map(1:size(alpha, 1)) = alpha(end:-1:1);
            alpha_map(end:-1:(end-size(alpha, 1)+1)) = alpha(end:-1:1);
            
            factor = 3/8;
            alpha_map2 = ones(size(colormap, 1), 1);
            alpha2 = linspace(0, 1, round(size(colormap, 1)*factor))';
            alpha_map2((1:size(alpha2, 1)) + size(colormap, 1)/2) = alpha2;
            alpha_map2((end:-1:(end-size(alpha2, 1)+1)) - size(colormap, 1)/2) = alpha2;
            
            factor = 3/8;
            alpha_map3 = ones(size(colormap, 1), 1);
            alpha3 = linspace(0, 1, round(size(colormap, 1)*factor))';
            idx1 = (size(colormap, 1)/2 - size(alpha3, 1))/2;
            alpha_map3((-idx1:idx1) + size(colormap, 1)/2) = 0;
            alpha_map3((1:size(alpha3, 1)) + size(colormap, 1)/2 + (size(colormap, 1)/2 - size(alpha3, 1))/2) = alpha3;
            alpha_map3((end:-1:(end-size(alpha3, 1)+1)) - size(colormap, 1)/2 - (size(colormap, 1)/2 - size(alpha3, 1))/2) = alpha3;
    
            triangle_alpha = alpha_map3(triangle_color);
            
            % Restore transparency for boundary elements;
            alpha_idx =  triangle_alpha(mdl.boundary_to_sub_element_idx, :) < 0.5;
            triangle_alpha(mdl.boundary_to_sub_element_idx(alpha_idx), :) = 0.5;
        else
            triangle_color = calc_colours(img, [], opts.colorbar.show);
        end
        triangle_color = squeeze(triangle_color(:, 1, :));
    else
        triangle_color = ones(size(elems, 1), 1)*[1 1 1];
    end
    
    % Assign colors for electrodes.
    marker_position = [];
    marker_color    = [];
    
    if (isfield(mdl, 'electrode'))
        for i = 1:numel(mdl.electrode)
            if (isfield(mdl.electrode(i), electrode_field) && ...
                    ~isempty(mdl.electrode(i).(electrode_field)))
                
                if (size(mdl.elems, 2) == 4)
                    triangle_color(...
                        mdl.electrode(i).(electrode_field), :) = ...
                        ones(numel(mdl.electrode(i).(electrode_field)), ...
                        1)*electr_colour(i, size(triangle_color, 2));

                    triangle_alpha(mdl.electrode(i).(electrode_field), ...
                                                                    :) = 1;
                    % Highlight electrode edges.
                    boundary_edges = [mdl.electrode(i).boundary_inner_edges;
                                      mdl.electrode(i).boundary_outer_edges];
                    edge_color(boundary_edges, :) = 0.5*ones(size(...
                        boundary_edges, 1), 1)*electr_colour(i, 3);
                    edge_width(mdl.electrode(i).boundary_outer_edges) = 1;
                else

                    edge_width(mdl.electrode(i).(electrode_field), :) = 3;
                    edge_color(mdl.electrode(i).(electrode_field), :) = ...
                        ones(numel(mdl.electrode(i).(electrode_field)), ...
                        1)*electr_colour(i, 3);
                end
            else
                if (isfield(mdl.electrode(i), 'nodes'))
                    marker_position(end + 1, :) = mean(mdl.nodes(...
                                            mdl.electrode(i).nodes, :), 1);
                elseif (isfield(mdl.electrode(i), 'pos'))
                    marker_position(end + 1, :) = mean(...
                                               mdl.electrode(i).pos, 1)*...
                                  eye(size(marker_position(end, :), 2), 3);
                end
                marker_color(end + 1, :) = electr_colour(i, 3);
            end
        end
    end

    % Draw triangles
    hh = draw_triangles(elems, mdl.nodes, ...
                      triangle_color, triangle_alpha, ...
                      triangle_edge_color, triangle_edge_width);

    % Draw edges if necessary
    if (~isempty(edge_edges))
        draw_edges(edge_edges, mdl.nodes, edge_width, edge_color);
    end
                  
    % Draw markers if necessary.
    if (~isempty(marker_position))
        draw_markers(marker_position, marker_color, 9);
    end
    
function draw_numbers(positions, opts)
    text(positions(:,1), positions(:,2), positions(:,3), ...
        arrayfun(@(x){int2str(x)}, 1:size(positions, 1)), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   'middle', ...
        'FontSize',            opts.font_size, ...
        'FontWeight',          opts.font_weight, ...
        'Color',               opts.color, ...
        'BackgroundColor',     opts.background_color);

function hh = draw_markers(points, colour, marker_size)
    [unique_colour, unused, colour_idx] = unique(colour, 'rows');
    
    for i = 1:size(unique_colour, 1)
        points_idx = colour_idx == i;
        hh = line(points(points_idx, 1), points(points_idx, 2), ...
                  points(points_idx, 3), ...
                  'LineStyle', 'none', ...
                  'Marker', 'o', 'MarkerSize', marker_size, ...
                  'MarkerFaceColor', unique_colour(i, :), ...
                  'MarkerEdgeColor', unique_colour(i, :));
    end
        
function hh = draw_edges(edges, vertices, width_data, color_data)
    [unique_width_color_data, unused, width_color_idx] = ... 
                                   unique([width_data color_data], 'rows');                            
    for i = 1:size(unique_width_color_data, 1)
        if (unique_width_color_data(i, 1) > 0)
            edge_idx = (width_color_idx == i);
            points_list = [];
            points_list(1:3:3*size(edges(edge_idx, :), 1), :) = ...
                                           vertices(edges(edge_idx, 1), :);
            points_list(2:3:3*size(edges(edge_idx, :), 1), :) = ...
                                           vertices(edges(edge_idx, 2), :);
            points_list(3:3:3*size(edges(edge_idx, :), 1), :) = ...
                                       nan(size(edges(edge_idx, :), 1), 3);
            hh = line(points_list(:, 1), points_list(:, 2), ...
                      points_list(:, 3), 'LineWidth', ...
                      unique_width_color_data(i, 1), 'LineStyle', '-', ...
                      'Color', unique_width_color_data(i, 2:end));
        end
    end

function hh = draw_triangles(faces, vertices, color_data, alpha_data, ...
                             edge_color, edge_width)
                         
    if (size(color_data, 1) == 1)
        color_data = ones(size(faces, 1), 1)*color_data;
        face_color = 'flat';
    elseif (size(color_data, 1) == size(faces, 1))
        face_color = 'flat';
    elseif (size(color_data, 1) == size(vertices, 1))
        face_color = 'interp';        
    else
        eidors_msg('warning: color data and mesh do not match. Showing grey', 1);
        color_data = 0.5*ones(size(faces, 1), 3);
        face_color = 'flat';      
    end
    
    if (size(alpha_data, 1) == 1)
        alpha_data = ones(size(faces, 1), 1)*alpha_data;
        face_color = 'flat';
    elseif (size(alpha_data, 1) == size(faces, 1))
        face_alpha = 'flat';
    elseif (size(alpha_data, 1) == size(vertices, 1))
        face_alpha = 'interp';          
    else
        eidors_msg('warning: alpha data and mesh do not match. Showing opaque', 1);
        alpha_data = 1;
        face_alpha = 'flat';        
    end

    hh = patch('Faces', faces, 'Vertices', vertices, ...
               'AlphaDataMapping', 'none', ...
               'CDataMapping', 'direct', ...
               'FaceVertexCData', color_data, ...
               'FaceVertexAlphaData', alpha_data, ...
               'FaceColor', face_color, 'FaceAlpha', face_alpha, ...
               'EdgeColor', edge_color, 'FaceLighting', 'flat', ...
               'LineWidth', edge_width);

function colour = electr_colour(e, colormap_width)
    switch (e)
        case 1
            desired_colour = [0 .7 0]; % light green electrode #1
        case 2
            desired_colour = [0 .5 0]; % mid-green electrode #2
        otherwise
            desired_colour = [0 .3 0]; % dark green
    end

    switch(colormap_width)
        case 1
            map = colormap;
            colormap_idx = find(all(map == ones(size(map, 1), 1)* ...
                                                       desired_colour, 2));
            if ~isempty(colormap_idx)
                colour = colormap_idx;
            else
                map = [map; desired_colour];
                colormap(map);
                colour = size(map, 1);
            end
        case 3 
            colour = desired_colour;
        otherwise
            error('Invalid colormap width.');
    end
    
function mdl = find_sub_elements(mdl)
    if (~isfield(mdl, 'sub_elements'))
        % Number of nodes per elements.
        n_nodes_per_elem = size(mdl.elems, 2);
        
        % Find sub-elements.
        combos= combnk(1:n_nodes_per_elem, n_nodes_per_elem - 1);
        mdl.sub_elements = sort( ...
                           reshape(mdl.elems(:, combos')', ...
                                   n_nodes_per_elem - 1, []), ...
                                 1)';
                       
        % Vector that associates each sub-element with
        % corresponding element.
        sub_element_idx = reshape(ones(n_nodes_per_elem, 1) ...
                                           *(1:size(mdl.elems, 1)),[] , 1);
                                       
        sub_element_other_node_idx = reshape((n_nodes_per_elem:-1:1)'...
                                     *ones(1, size(mdl.elems, 1)), [] , 1);
                                 
        sub_element_other_node = mdl.elems(sub2ind(size(mdl.elems), sub_element_idx, sub_element_other_node_idx));
                       
        % Remove duplicate sub-elements. 
        [mdl.sub_elements, ia, ic] = unique(...
                                       mdl.sub_elements, 'rows', 'sorted');
                                   
         sub_element_other_node = sub_element_other_node(ia, :);                    
                                
        % Compute a matrix that converts a property defined on the elements
        % to a property defined on the sub-elements.
        mdl.element2sub_element = sparse(ic, sub_element_idx, ...
                             ones(n_nodes_per_elem*size(mdl.elems, 1), 1));
                         
        % Normalize each row to one to account for boundary sub-elements 
        % and internal sub-elements.
        mdl.element2sub_element = spdiags(1./sum(...
                                     mdl.element2sub_element, 2), ...
                                    0, size(mdl.sub_elements, 1), ...
                      size(mdl.sub_elements, 1))*mdl.element2sub_element;
                  
        % Extract boundary from sub-elements.
        mdl.boundary = mdl.sub_elements;
        boundary_other_node = sub_element_other_node;
        mdl.boundary_to_sub_element_idx = (1:size(mdl.sub_elements, 1))';

        sorted_ic = sort(ic);

        mdl.boundary(sorted_ic(diff(sorted_ic) == 0), :) = [];
        boundary_other_node(sorted_ic(diff(sorted_ic) == 0), :) = [];
        mdl.boundary_to_sub_element_idx(sorted_ic(diff(sorted_ic) == 0)) = [];
      
        if (size(mdl.boundary, 2) == 3)
            % Compute normal vectors.
            Node1 = mdl.nodes(mdl.boundary(:, 1), :);
            Node2 = mdl.nodes(mdl.boundary(:, 2), :);
            Node3 = mdl.nodes(mdl.boundary(:, 3), :);
            Node4 = mdl.nodes(boundary_other_node, :);
            mdl.boundary_normal_vector = cross(Node1 - Node2, ...
                                               Node3 - Node2, 2);
 
            % Normalize normal vectors.
            norm = 1./sqrt(sum(mdl.boundary_normal_vector ...
                                         .*mdl.boundary_normal_vector, 2));
            mdl.boundary_normal_vector = mdl.boundary_normal_vector ...
                                                        .*[norm norm norm];
              
            % Check if normal is outward. If not, invert direction.
            normal_vector_idx = dot(mdl.boundary_normal_vector, ...
                                                     Node4 - Node2, 2) > 0;
            mdl.boundary_normal_vector(normal_vector_idx, :) = ...
                         -mdl.boundary_normal_vector(normal_vector_idx, :);
                                                    
            % Find boundary edges.
            mdl.boundary_edges = sort(reshape(mdl.boundary(:, ...
                                            combnk(1:3, 2)')', 2, []), 1)';

            % Vector that associates each edge with its 
            % corresponding two boundary elements.
            sub_boundary_edges_idx = reshape(ones(3, 1) ...
                                       *(1:size(mdl.boundary, 1)), [] , 1);

            [mdl.boundary_edges, sorted_idx] = sortrows(mdl.boundary_edges);

            mdl.boundary_edges = mdl.boundary_edges(1:2:end, :);

            mdl.edge2boundary = reshape(sub_boundary_edges_idx(...
                                                      sorted_idx), 2, [])';               
        end

        % Extract electrode boundary if necessary.
        if (isfield(mdl, 'electrode'))
            for i = 1:numel(mdl.electrode)
                mdl.electrode(i).boundary = find(all(ismember(...
                                mdl.boundary, mdl.electrode(i).nodes), 2));
                mdl.electrode(i).sub_elements = ...
                                    mdl.boundary_to_sub_element_idx( ...
                                                mdl.electrode(i).boundary);
                
                % Find electrode edges.
                if (isfield(mdl, 'boundary_edges'))
                    edge_candidates = sum(ismember(mdl.edge2boundary, ...
                                            mdl.electrode(i).boundary), 2);
                    mdl.electrode(i).boundary_inner_edges = find(...
                                                     edge_candidates == 2);
                    mdl.electrode(i).boundary_outer_edges = find(...
                                                     edge_candidates == 1);                                           
                end
            end
        end
    end
   
function dest_struct = copy_field(dest_struct, src_struct, parent_field_name, child_field_name)
    if (isfield(src_struct, parent_field_name) && ...
        isfield(dest_struct, parent_field_name) && ...
        isfield(src_struct.(parent_field_name), child_field_name))
            dest_struct.(parent_field_name).(child_field_name) = ...
                src_struct.(parent_field_name).(child_field_name);
    end

% TESTS:
function do_unit_test
   ver= eidors_obj('interpreter_version');
   clf

   img=calc_jacobian_bkgnd(mk_common_model('a2c0',8)); 
   img.elem_data=rand(size(img.fwd_model.elems,1),1);
   subplot(3,4,1); show_fem(img.fwd_model,[0,0,1]) 
   title('regular mesh numbered');

if ~ver.isoctave 
   imgn = rmfield(img,'elem_data');
   imgn.node_data=rand(size(img.fwd_model.nodes,1),1);
   subplot(3,4,9); show_fem(imgn) 
   title('interpolated node colours');
end

   img2(1) = img; img2(2) = img;
   subplot(3,4,2); show_fem(img,[1]);
   title('colours with legend');
   subplot(3,4,3); show_fem(img2,[0,1]);
   title('colours with legend');
   img.calc_colours.mapped_colour = 0; % USE RGB colours
   subplot(3,4,4); show_fem(img,[0,1,1]);
   title('RGB colours');
   subplot(3,4,4); show_fem(img);
   title('RGB colours');

   img.elem_data = [1:10];
   subplot(3,4,12);show_fem(img); %Should show grey
   title('error -> show grey');

if ~ver.isoctave
   imgn.calc_colours.mapped_colour = 0; % USE RGB colours
   subplot(3,4,10);show_fem(imgn,[0,1]) 
   title('interpolated node colours');


   subplot(3,4,11);hh=show_fem(imgn); set(hh,'EdgeColor',[0,0,1]);
   title('with edge colours');

end

   img3=calc_jacobian_bkgnd(mk_common_model('n3r2',[16,2]));
   img3.elem_data= randn(828,1);                       
   subplot(3,4,5); show_fem(img3.fwd_model) 
   subplot(3,4,6); show_fem(img3,[1])
   subplot(3,4,7); show_fem(img3,[1,1])
   subplot(3,4,8); show_fem(img3,[1,1,1])
