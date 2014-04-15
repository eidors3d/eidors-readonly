function [fmdl, mat_idx] = ng_mk_geometric_models(body_geometry, electrode_geometry)
% NG_MK_GEOMETRIC_MODELS: create geometric mesh models using Netgen. Body
% and electrodes are defined as combinations of solid primitives. The 3-D
% surface intersection of the electrode and body volumes define the
% electrode surfaces.
%
%[fmdl, mat_idx] = ng_mk_geometric_models(body_geometry, electrode_geometry)
%
% INPUT:
%  body_geometry      - Structure whose fields describe body geometries as
%                       combinations of unions, intersections and 
%                       complements of solid primitives. Cell array
%                       should be used when describing several body 
%                       geometries.
%
%  electrode_geometry - Structure whose fields describe electrode
%                       geometries as combinations of unions, intersections
%                       and complements of solid primitives. Cell array
%                       should be used when describing several electrode 
%                       geometries.
%
% The following field names are available for geometry descriptions. Each 
% field is followed by the available subfields whose default values if left 
% unspecified are indicated between parentheses. 
%
% complement_flag:    If true, the desired geometry description is the
%                     complement of the geometry description. (false)
%
% body_of_revolution: A body of revolution is described with the following
%                     subfields: axis_point_a ([0; 0; 0]), axis_point_b 
%                     ([0; 0; 1]), points (1 1; 1 2; 2 2; 2 1]),
%                     segments ([1 2; 2 3; 3 4; 4 1]), complement_flag
%                     (false).    
%
% cone:               A cone is described with the following subfields:
%                     bottom_center ([0; 0; 0]), bottom_radius (1),
%                     top_center ([0; 0; 1]), top_radius (0.5),
%                     complement_flag (false).                    
%
% cylinder:           A cylinder is described with the following subfields:
%                     bottom_center ([0; 0; 0]), top_center ([0; 0; 1]),
%                     radius (1), complement_flag (false).
%
% ellipsoid:          An ellipsoid is described with the following
%                     subfields: center ([0; 0; 0]), axis_a ([1; 0; 0]),
%                     axis_b ([0; 1; 0]), axis_c ([0; 0; 1]),
%                     complement_flag (false).
%
% elliptic_cylinder:  An elliptic cylinder is described with the following
%                     subfields: bottom_center ([0; 0; 0]),
%                     top_center ([0; 0; 1]), axis_a ([1; 0; 0]),
%                     axis_b ([0; 1; 0]), complement_flag (false). 
%
% enter_body_flag:    This flag can be used only for electrode geometry
%                     descriptions to indicate that the associated
%                     electrode solid enters the body solids. It can only
%                     be defined at the top level of each geometry
%                     description. If this flag is true, it means that the
%                     volume of the electrode intersecting with any body is
%                     part of the electrode, otherwise it is part of a
%                     body. (false)
%
% half_space          A half-space is described by the following subfields:
%                     point ([0; 0; 0]), outward_normal_vector ([0; 0; 1]),
%                     complement_flag (false).
%
% intersection:       This fields indicates to perform the intersection of
%                     all subfields. Subfields: complement_flag (false).
%
% keep_material_flag: This flag can be used only for electrode geometry 
%                     descriptions to indicate that the associated
%                     electrode material should be kept in the final mesh.
%                     It can only be defined at the top level of each
%                     geometry description. If true, it means that the
%                     volume of the electrode is meshed. Volume elements
%                     that are part of the mesh are indicated in mat_idx
%                     output argument. (false)
%
% max_edge_length:    This parameter is used to adjust the maximum size of
%                     the element composing the mesh. It can only be used
%                     at the top level of each geometry description. (inf)
%
% name:               This parameter is used to name the geometry
%                     description.
%
% ortho_brick:        An ortho-brick is described by the following
%                     subfields: opposite_corner_a ([0; 0; 0]),
%                     opposite_corner_b ([1; 1; 1]), complement_flag
%                     (false).
%
% parallelepiped:     A parallelepiped is described by the following
%                     subfields: vertex ([0; 0; 0]), vector_a ([1; 0; 0]),
%                     vector_b ([0; 1; 0]), vector_c ([0; 0; 1]),
%                     complement_flag (false).
%
% point:              This parameter describes a point. It can only be used
%                     at the top level of an electrode geometry
%                     description. It must be the only field in the
%                     structure. ([])
%
% sphere:             A sphere is described with the following subfields:
%                     center ([0; 0; 0]), radius (1), complement_flag
%                     (false).
%
% union:              This will perform the union of all its subfields.
%                     There is an implicit union operator at the top level
%                     of the geometry structure. Subfields: complement_flag
%                     (false)
%
% OUTPUT:
%  fmdl               - EIDORS forward model object.
%
%  mat_idx            - Vector indicating for each mesh element the indices
%                       of materials corresponding as separately defined
%                       by input argument body_geometry.
%
% USAGE EXAMPLES:
% % 3D cylinder with radius 1. One plane of 16 electrodes with radius 0.1
%   body_geometry.cylinder = struct;
%   n_elect = 16;
%   theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
%   for i = 1:n_elect
%     electrode_geometry{i}.sphere.center = [cos(theta(i)) sin(theta(i)) 0.5];
%     electrode_geometry{i}.sphere.radius = 0.1;
%   end
%   fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);

% (C) Hervé Gagnon, 2013. Licenced under GPL v2 or v3

% Check if function is called in UNIT_TEST mode.
if (ischar(body_geometry) && strcmp(body_geometry, 'UNIT_TEST'))
    do_unit_test; 
    return; 
end

% Validate function parameters.
if (isempty(body_geometry) || ~isstruct(body_geometry) && ~iscell(body_geometry))
   error('Parameter body_geometry must be a structure or a cell.');
end

% Check if parameter electrode_geometry is specified.
if (nargin < 2 || isempty(electrode_geometry))
    electrode_geometry = {};
end

if (~isstruct(electrode_geometry) && ~iscell(electrode_geometry))
   error('Parameter electrode_geometry must be a structure or a cell.');
end

% If parameters are not of cell type, convert them to cell.
if (~iscell(body_geometry))
    body_geometry = {body_geometry};
end

if (~iscell(electrode_geometry))
    electrode_geometry = {electrode_geometry};
end

% Check if result is already in cache. Otherwise, compute and store in cache. 
cache_obj = {body_geometry, electrode_geometry};

fmdl_mat_idx = eidors_obj('get-cache', cache_obj, 'ng_mk_geometric_models');

if (isempty(fmdl_mat_idx))
   fmdl_mat_idx = mk_geometric_models(body_geometry, electrode_geometry);
   eidors_cache('boost_priority', -2); % netgen objs are low priority
   eidors_obj('set-cache', cache_obj, 'ng_mk_geometric_models', fmdl_mat_idx);
   eidors_cache('boost_priority', +2); % restore former priority
end

% Reformat output arguments. 
fmdl    = fmdl_mat_idx{1};
mat_idx = fmdl_mat_idx{2};

function [fmdl_mat_idx] = mk_geometric_models(body_geometry, electrode_geometry)     
    % Find number of subdomains
    n_body      = numel(body_geometry);
    n_electrode = numel(electrode_geometry);

    % Allocate cell memory.
    body_solid_code       = cell(size(body_geometry));
    body_extra_code       = cell(size(body_geometry));
    body_extra_param      = cell(size(body_geometry));
    electrode_solid_code  = cell(size(electrode_geometry));
    electrode_extra_code  = cell(size(electrode_geometry));
    electrode_extra_param = cell(size(electrode_geometry));
    
    % Parse geometry for each body subdomain.
    for i = 1:n_body
        if (isfield(body_geometry{i}, 'point'))
            error('Field name "point" is not allowed for body geometry.');
        end
        if (isfield(body_geometry{i}, 'enter_body_flag'))
            error('Field name "enter_body_flag" is not allowed for body geometry.');
        end
        if (isfield(body_geometry{i}, 'keep_material_flag'))
            error('Field name "keep_material_flag" is not allowed for body geometry.');
        end
        [body_solid_code{i} body_extra_code{i} body_extra_param{i}] = parse_geometry(body_geometry{i});  
    end
    
    % Parse geometry for each body subdomain.
    for i = 1:n_electrode
        [electrode_extra_code{i} electrode_extra_param{i}] = parse_geometry_point(electrode_geometry{i});
        if (isempty(electrode_extra_param{i}.point))
            [electrode_solid_code{i} electrode_extra_code{i} electrode_extra_param{i}] = parse_geometry(electrode_geometry{i}); 
        end
    end

    % Define temporary unique filenames. 
    fn_prefix = tempname;
    geo_fn    = [fn_prefix, '.geo'];
    vol_fn    = [fn_prefix, '.vol'];
   
     % Add body names if unspecified.   
    for i = 1:numel(body_solid_code)
        if (isempty(body_extra_param{i}.name))
            body_extra_param{i}.name = sprintf('body%04d', i);
        end
    end
    
    % Add electrode names if unspecified.
    for i = 1:numel(electrode_solid_code)
        if (isempty(electrode_extra_param{i}.name))
            electrode_extra_param{i}.name = sprintf('electrode%04d', i);
        end
    end

    % Write geo file for netgen.
    write_geo_file(geo_fn, body_solid_code, electrode_solid_code, body_extra_code, electrode_extra_code, body_extra_param, electrode_extra_param);
   
    % Call netgen.
    call_netgen(geo_fn, vol_fn);
 
    % Read vol file generated by netgen.
    fmdl_mat_idx{1} = read_vol_file(vol_fn, electrode_extra_param);
    
    % Delete temporary files.
    %delete(geo_fn);
    %delete(vol_fn);

    % Complete fmdl object.
    fmdl_mat_idx{1} = complete_fmdl(fmdl_mat_idx{1}, electrode_extra_param);

    % Assign mat_idx value from fmdl.
    fmdl_mat_idx{2} = fmdl_mat_idx{1}.mat_idx;

function radius = assign_radius(struct, n_structs, struct_name, field_name, default_value)

    radius = default_value;
    
    for i = 1:n_structs
        value = struct(i).(field_name);

        if (~isempty(value))
            if (isscalar(value) && isnumeric(value) && isreal(value) && value > 0)
                radius(i) = value;
            else 
                error('%s(%d).%s value is not valid.', struct_name, i, field_name);
            end
        end
    end
    
function point = assign_point(struct, n_structs, struct_name, field_name, default_value)
        
    point = default_value;
    
    for i = 1:n_structs
        value =  struct(i).(field_name);

        if (~isempty(value))
            if (numel(value) == 3 && isnumeric(value) && isreal(value))
                point(:, i) = value;
            else
                error('%s(%d).%s value is not valid.', struct_name, i, field_name);
            end
        end
    end
 
function point_list = assign_list_of_3D_points(struct, n_structs, struct_name, field_name, default_value)
        
    point_list = cell(n_structs, 1);
    
    for i = 1:n_structs
        
        point_list{i} = default_value;
        
        value = struct(i).(field_name);

        if (~isempty(value))
            if (size(value, 2) == 3 && isnumeric(value) && isreal(value))
                point_list{i} = value;
            else
                error('%s(%d).%s value is not valid.', struct_name, i, field_name);
            end
        end
    end
    
function segment_list = assign_list_of_2D_points(struct, n_structs, struct_name, field_name, default_value)
        
    segment_list = cell(n_structs, 1);
    
    for i = 1:n_structs
        
        segment_list{i} = default_value;
        
        value = struct(i).(field_name);

        if (~isempty(value))
            if (size(value, 2) == 2 && isnumeric(value) && isreal(value))
                segment_list{i} = value;
            else
                error('%s(%d).%s value is not valid.', struct_name, i, field_name);
            end
        end
    end 
    
 function segment_list = assign_list_of_2D_or_3D_points(struct, n_structs, struct_name, field_name, default_value)
        
    segment_list = cell(n_structs, 1);
    
    for i = 1:n_structs
        
        segment_list{i} = default_value;
        
        value = struct(i).(field_name);

        if (~isempty(value))
            if ((size(value, 2) == 2 || size(value, 2) == 3) && isnumeric(value) && isreal(value))
                segment_list{i} = value;
            else
                error('%s(%d).%s value is not valid.', struct_name, i, field_name);
            end
        end
    end    
    
function flag = assign_flag(struct, n_structs, struct_name, field_name, default_value)

    flag = default_value;

    for i = 1:n_structs
        value = struct(i).(field_name);
 
        if (~isempty(value))
            if (isscalar(value) && (islogical(value) || (isnumeric(value) && ...
                                     isreal(value) && (value == 0 || value == 1))))
                flag(i) = value;
            else
                error('%s(%d).%s value is not valid.', struct_name, i, field_name);
            end
        end
    end

function [extra_code extra_param] = parse_geometry_point(geometry)

    % Initialize extra param values to default values.
    extra_code                     = '';
    extra_param.point              = [];
    extra_param.max_edge_length    = inf;
    extra_param.enter_body_flag    = false;
    extra_param.keep_material_flag = false;
    extra_param.name               = '';

    % Check if geometry is a point, return otherwise. 
    if (isfield(geometry, 'point'))
        % Check if a single point is defined.
        if (numel(geometry) ~= 1)
            error('Field name "point" must define only a single point.');
        end
        
        % Get structure field names.
        field_names = fieldnames(geometry);
        n_fields = numel(field_names);
        
        % Check if it is the only field names.
        if (n_fields ~= 1)
            if (isfield(geometry, 'name'))
                extra_param.name = geometry.name;
            else
                error('Field name "point" must be used as a single field.');
            end
        end
        
        % Check point is 3D.
        if (numel(geometry.point) ~= 3)
            error('geometry.point value is not valid.');
        end
        
        % Set returning values.
        extra_param.point = geometry.point(:);
        
        extra_code = sprintf('point(%g, %g, %g);\n', extra_param.point);
    end
    
function [geo_code extra_code extra_param] = parse_geometry(geometry, field_operator_string, element_operator_string)

    % Extra code is initialized empty;
    extra_code = '';
    
    % Initialize extra param values to default values.
    extra_param.max_edge_length     = inf;
    extra_param.enter_body_flag     = false;
    extra_param.keep_material_flag  = false;
    extra_param.name                = '';
    
    % If called from top_level, operators default to or.
    if (nargin == 1)
        top_level_flag = 1;
        field_operator_string = ' or ';
        element_operator_string = ' or ';
    else
        top_level_flag = 0;
    end

    eidors_msg('@@@ called with (field: "%s", element: "%s").', field_operator_string, element_operator_string, 3);
    
    % Validate that geometry is a non-empty structure.
    if (~isstruct(geometry) || isempty(geometry))
        error('Parameter geometry must be a valid structure.');        
    else
        % Get number of geometries.
        n_geometries = numel(geometry);
        
        % Get structure field names.
        field_names = fieldnames(geometry);
        n_fields = numel(field_names);

        % Recursively parse all geometry fields.
        geo_code = '(';
        for i = 1:n_geometries
            % complement_flag field has to be processed first.
            if (isfield(geometry(i), 'complement_flag') && ~isempty(geometry(i).complement_flag) && geometry(i).complement_flag)
                 geo_code = [geo_code '(not('];
            else
                 geo_code = [geo_code '('];
            end
            % Process name field.
            if (isfield(geometry(i), 'name'))
                if (~top_level_flag)
                    error('Field "name" can only be specified at the top level of the geometry description');
                end
                extra_param.name = geometry(i).name;

                if (isempty(extra_param.name) || ~ischar(extra_param.name))
                    error('name value is not valid.');
                end
            end
           % Process max_edge_length field.
            if (isfield(geometry(i), 'max_edge_length'))
                if (~top_level_flag)
                    error('Field "max_edge_length" can only be specified at the top level of the geometry description');
                end
                extra_param.max_edge_length = geometry(i).max_edge_length;

                if (isempty(extra_param.max_edge_length) || ~isscalar(extra_param.max_edge_length) || ~isnumeric(extra_param.max_edge_length) || ~isreal(extra_param.max_edge_length) || extra_param.max_edge_length <= 0)
                    error('max_edge_length value is not valid.');
                end
            end
            % Process enter_body_flag field.
            if (isfield(geometry(i), 'enter_body_flag'))
                if (~top_level_flag)
                    error('Field "enter_body_flag" can only be specified at the top level of the geometry description');
                end
                extra_param.enter_body_flag = geometry(i).enter_body_flag;
                
                if (isempty(extra_param.enter_body_flag) || ~isscalar(extra_param.enter_body_flag) || (~islogical(extra_param.enter_body_flag) && ...
                        (~isnumeric(extra_param.enter_body_flag) || ~isreal(extra_param.enter_body_flag) || (extra_param.enter_body_flag ~= 0 && extra_param.enter_body_flag ~= 1))))
                    error('Field "enter_body_flag value" is not valid.');
                end
            end
            % Process keep_material_flag field.
            if (isfield(geometry(i), 'keep_material_flag'))
                if (~top_level_flag)
                    error('Field "keep_material_flag" can only be specified at the top level of the geometry description');
                end
                extra_param.keep_material_flag = geometry(i).keep_material_flag;
                
                if (isempty(extra_param.keep_material_flag) || ~isscalar(extra_param.keep_material_flag) || (~islogical(extra_param.keep_material_flag) && ...
                        (~isnumeric(extra_param.keep_material_flag) || ~isreal(extra_param.keep_material_flag) || (extra_param.keep_material_flag ~= 0 && extra_param.keep_material_flag ~= 1))))
                    error('Field "keep_material_flag" value is not valid.');
                end
            end  
            first_internal_term = 1;
            for j = 1:n_fields
                if (~isempty(geometry(i).(field_names{j})) && ~strcmp(field_names{j}, 'complement_flag') && ~strcmp(field_names{j}, 'name') && ~strcmp(field_names{j}, 'max_edge_length') && ~strcmp(field_names{j}, 'keep_material_flag') && ~strcmp(field_names{j}, 'enter_body_flag'))
                    if (first_internal_term)
                        first_internal_term = 0;
                    else
                        geo_code = [geo_code field_operator_string];
                    end
                    switch (field_names{j})
                        case 'body_of_extrusion'
                            [geo_code_temp extra_code_temp] = ...
                                        parse_geometry_body_of_extrusion(geometry(i).(field_names{j}), field_operator_string);
                            geo_code   = [geo_code geo_code_temp];        
                            extra_code = [extra_code extra_code_temp];
                        case 'body_of_revolution'
                            [geo_code_temp extra_code_temp] = ...
                                        parse_geometry_body_of_revolution(geometry(i).(field_names{j}), field_operator_string);
                            geo_code   = [geo_code geo_code_temp];        
                            extra_code = [extra_code extra_code_temp];
                        case 'cone'
                            geo_code = [geo_code ...
                                        parse_geometry_cone(geometry(i).(field_names{j}), field_operator_string)];
                        case 'cylinder'
                            geo_code = [geo_code ...
                                        parse_geometry_cylinder(geometry(i).(field_names{j}), field_operator_string)];
                        case 'ellipsoid'
                            geo_code = [geo_code ...
                                        parse_geometry_ellipsoid(geometry(i).(field_names{j}), field_operator_string)];
                        case 'elliptic_cylinder'
                            geo_code = [geo_code ...
                                        parse_geometry_elliptic_cylinder(geometry(i).(field_names{j}), field_operator_string)];
                        case 'half_space'
                            geo_code = [geo_code ...
                                        parse_geometry_half_space(geometry(i).(field_names{j}), field_operator_string)];    
                        case 'ortho_brick'
                            geo_code = [geo_code ...
                                        parse_geometry_ortho_brick(geometry(i).(field_names{j}), field_operator_string)];
                        case 'parallelepiped'
                            geo_code = [geo_code ...
                                        parse_geometry_parallelepiped(geometry(i).(field_names{j}), field_operator_string)];       
                        case 'sphere'
                            geo_code = [geo_code ...
                                        parse_geometry_sphere(geometry(i).(field_names{j}), field_operator_string)];
                        case 'intersection'
                            [geo_code_temp extra_code_temp] = ...
                                        parse_geometry(geometry(i).(field_names{j}), ' and ', field_operator_string);
                            geo_code   = [geo_code geo_code_temp];        
                            extra_code = [extra_code extra_code_temp];
                        case 'union'
                            [geo_code_temp extra_code_temp] = ...
                                        parse_geometry(geometry(i).(field_names{j}), ' or ', field_operator_string);
                            geo_code   = [geo_code geo_code_temp];        
                            extra_code = [extra_code extra_code_temp];
                        otherwise
                            error(['Field name "%s" is not valid for a geometry.\nAvailable field names for a geometry are: '...
                                   'complement_flag, intersection, union, body_of_extrusion, body_of_revolution, cone, cylinder, ellipsoid, elliptic_cylinder, half_space, ortho_brick, parallelepiped, point, sphere, keep_material_flag, enter_body_flag, name, and max_edge_length.'], field_names{j});
                    end
                end
            end
            if (isfield(geometry(i), 'complement_flag') && ~isempty(geometry(i).complement_flag) && geometry(i).complement_flag)
                geo_code = [geo_code '))'];
            else
                geo_code = [geo_code ')'];  
            end
           
            if (i < n_geometries)
                geo_code = [geo_code element_operator_string];         
            end           
        end
        geo_code = [geo_code ')'];  
    end
    
    eidors_msg('@@@ returned with (field: "%s", element: "%s").', field_operator_string, element_operator_string, 3);
 
function [geo_code extra_code] = parse_geometry_body_of_extrusion(body_of_extrusion, operator_string)

    eidors_msg('@@@ called with "%s" starting.', operator_string, 3);

    if (~isstruct(body_of_extrusion) || isempty(body_of_extrusion))
        error('Parameter body_of_extrusion must be a valid structure.');        
    else
        % Get number of body_of_extrusion.
        n_body_of_extrusions = numel(body_of_extrusion);
        
        % Get structure field names.
        field_names = fieldnames(body_of_extrusion);
        n_fields = numel(field_names);
        
        % Assign default values.
        vector_d            = [0; 1; 0]*ones(1, n_body_of_extrusions);
        profile_points{1}   = [1 1; 1 2; 2 2; 2 1];
        profile_segments{1} = [1 2; 2 3; 3 4; 4 1];
        path_points{1}      = [0 0 0; 0 0 1; 0 0 2; 0 0 3];
        path_segments{1}    = [1 2; 2 3; 3 4];
        complement_flag     = false(1, n_body_of_extrusions);
        
        % Parse all structure fields.
        for i = 1:n_fields
            switch (field_names{i})
                case 'vector_d'
                    vector_d = assign_point(body_of_extrusion, n_body_of_extrusions, 'body_of_extrusion', field_names{i}, vector_d);
                case 'profile_points'
                    profile_points = assign_list_of_2D_points(body_of_extrusion, n_body_of_extrusions, 'body_of_extrusion', field_names{i}, profile_points{1});
                case 'profile_segments'
                    profile_segments = assign_list_of_2D_or_3D_points(body_of_extrusion, n_body_of_extrusions, 'body_of_extrusion', field_names{i}, profile_segments{1});
                case 'path_points'
                    path_points = assign_list_of_3D_points(body_of_extrusion, n_body_of_extrusions, 'body_of_extrusion', field_names{i}, path_points{1});
                case 'path_segments'
                    path_segments = assign_list_of_2D_or_3D_points(body_of_extrusion, n_body_of_extrusions, 'body_of_extrusion', field_names{i}, path_segments{1});
                case 'complement_flag'
                    complement_flag = assign_flag(body_of_extrusion, n_body_of_extrusions, 'body_of_extrusion', field_names{i}, complement_flag);
                otherwise
                    error(['Field name "%s" is not allowed for a body_of_extrusion!\nAllowed field names for a body_of_extrusion are: ' ...
                           'path_points, path_segments, profile_points, profile_segments, vector_d, and complement_flag.'], field_names{i});
            end
        end
        
        % Start geo code with an opening parenthesis.
        geo_code = '(';
        extra_code = '';
        
        % Add geo code for each body_of_extrusion.
        for i = 1:n_body_of_extrusions
            
            for j = 1:size(path_segments{i}, 1)
                if (dot(vector_d(:,i), path_points{i}(path_segments{i}(j, 1), :) - path_points{i}(path_segments{i}(j, end), :)) ~= 0)
                    error('vector_d and path must be perpendicular for a body of extrusion.');
                end
            end
            
            n_points = size(profile_points{i}, 1);
            n_segments = size(profile_segments{i}, 1);
            
            if (size(profile_segments{i}, 2) == 2)
                extra_code = [extra_code sprintf('curve2d Extrusion2DProfileCurve%d = (%g ; ', i, n_points) ...
                                         sprintf('%g, %g ; ', profile_points{i}') ...
                                         sprintf(' %g ', n_segments) ...
                                         sprintf('; 2, %g, %g ', profile_segments{i}') ...
                                         sprintf(');\n\n')];
            else
                extra_code = [extra_code sprintf('curve2d Extrusion2DProfileCurve%d = (%g ; ', i, n_points) ...
                                         sprintf('%g, %g ; ', profile_points{i}') ...
                                         sprintf(' %g ', n_segments) ...
                                         sprintf('; 3, %g, %g, %g ', profile_segments{i}') ...
                                         sprintf(');\n\n')];
            end
  
            n_points = size(path_points{i}, 1);
            n_segments = size(path_segments{i}, 1);
            
            if (size(path_segments{i}, 2) == 2)
                extra_code = [extra_code sprintf('curve3d Extrusion3DPathCurve%d = (%g ; ', i, n_points) ...
                                         sprintf('%g, %g, %g ; ', path_points{i}') ...
                                         sprintf(' %g ', n_segments) ...
                                         sprintf('; 2, %g, %g ', path_segments{i}') ...
                                         sprintf(');\n\n')];
            else
                 extra_code = [extra_code sprintf('curve3d Extrusion3DPathCurve%d = (%g ; ', i, n_points) ...
                                         sprintf('%g, %g, %g ; ', path_points{i}') ...
                                         sprintf(' %g ', n_segments) ...
                                         sprintf('; 3, %g, %g, %g ', path_segments{i}') ...
                                         sprintf(');\n\n')];               
            end
                                       
            if (complement_flag(i))
                geo_code = [geo_code 'not '];
            end
            
            % Check if path is closed.
            if (path_segments{i}(end) == path_segments{i}(1))
                geo_code = [geo_code sprintf('extrusion(Extrusion3DPathCurve%d ; Extrusion2DProfileCurve%d ; %g, %g, %g)', i, i, vector_d(:,i))];
            else
                %error('Unclosed path are not yet supported for body of extrusion!');
                %for j = 1: [1 1; 1 2; 2 2; 2 1];
                %polyhedron1 = 'polyhedron(-1, 1, 0; -1, 2, 0; -2, 2, 0; -2, 1, 0; -1.5, 1.5, 0; -1.5, 1.5, 1 ;; 2, 1, 5; 3, 2, 5; 4, 3, 5; 1, 4, 5; 1, 2, 6; 2, 3, 6; 3, 4, 6; 4, 1, 6)';
                %polyhedron2 = 'polyhedron(-1, 1, 3; -1, 2, 3; -2, 2, 3; -2, 1, 3; -1.5, 1.5, 3; -1.5, 1.5, 2 ;; 1, 2, 5; 2, 3, 5; 3, 4, 5; 4, 1, 5; 2, 1, 6; 3, 2, 6; 4, 3, 6; 1, 4, 6)'; 
                %polyhedron1 = 'plane(0, 0, 3; 0, 0, 1)';
                %polyhedron2 = 'plane(0, 0, 0; 0, 0, -1)';
                first_point  = path_points{i}(path_segments{i}(1, 1), :);
                first_vector = first_point - path_points{i}(path_segments{i}(1, end), :);
                last_point   = path_points{i}(path_segments{i}(end, end), :);
                last_vector  = last_point - path_points{i}(path_segments{i}(end, 1), :);
                geo_code = [geo_code sprintf('(extrusion(Extrusion3DPathCurve%d ; Extrusion2DProfileCurve%d ; %g, %g, %g) and plane(%g, %g, %g; %g, %g, %g) and plane(%g, %g, %g; %g, %g, %g))', ...
                            i, i, vector_d(:,i), first_point, first_vector, last_point, last_vector)];
                %geo_code = [geo_code sprintf('(%s or %s)', polyhedron1, polyhedron2)];
            end

            if (i < n_body_of_extrusions)
                geo_code = [geo_code operator_string];
            else
                geo_code = [geo_code ')'];             
            end
        end
    end
    
    eidors_msg('@@@ called with "%s" returning.', operator_string, 3);

function [geo_code extra_code] = parse_geometry_body_of_revolution(body_of_revolution, operator_string)

    eidors_msg('@@@ called with "%s" starting.', operator_string, 3);

    if (~isstruct(body_of_revolution) || isempty(body_of_revolution))
        error('Parameter body_of_revolution must be a valid structure.');        
    else
        % Get number of body_of_revolution.
        n_body_of_revolutions = numel(body_of_revolution);
        
        % Get structure field names.
        field_names = fieldnames(body_of_revolution);
        n_fields = numel(field_names);
        
        % Assign default values.
        axis_point_a   = [0;0;0]*ones(1, n_body_of_revolutions);
        axis_point_b   = [0;0;1]*ones(1, n_body_of_revolutions);
        points{1}      = [1 1; 1 2; 2 2; 2 1];
        segments{1}    = [1 2; 2 3; 3 4; 4 1];
        complement_flag = false(1, n_body_of_revolutions);
        
        % Parse all structure fields.
        for i = 1:n_fields
            switch (field_names{i})
                case 'axis_point_a'
                    axis_point_a = assign_point(body_of_revolution, n_body_of_revolutions, 'body_of_revolution', field_names{i}, axis_point_a);
                case 'axis_point_b'
                    axis_point_b = assign_point(body_of_revolution, n_body_of_revolutions, 'body_of_revolution', field_names{i}, axis_point_b);
                case 'points'
                    points = assign_list_of_2D_points(body_of_revolution, n_body_of_revolutions, 'body_of_revolution', field_names{i}, points{1});
                case 'segments'
                    segments = assign_list_of_2D_or_3D_points(body_of_revolution, n_body_of_revolutions, 'body_of_revolution', field_names{i}, segments{1});
                case 'complement_flag'
                    complement_flag = assign_flag(body_of_revolution, n_body_of_revolutions, 'body_of_revolution', field_names{i}, complement_flag);
                otherwise
                    error(['Field name ''%s'' is not valid for a body_of_revolution.\Available field names for a body_of_revolution are: ' ...
                           'axis_point_a, axis_point_b, points, segments, and complement_flag.'], field_names{i});
            end
        end
        
        % Start geo code with an opening parenthesis.
        geo_code = '(';
        extra_code = '';
        
        % Add geo code for each body_of_revolution.
        for i = 1:n_body_of_revolutions
            
            n_points = size(points{i}, 1);
            n_segments = size(segments{i}, 1);
            
            if (size(segments{i}, 2) == 2)
                extra_code = [extra_code sprintf('curve2d Revolution2DCurve%d = (%g ; ', i, n_points) ...
                                               sprintf('%g, %g ; ', points{i}') ...
                                               sprintf(' %g ', n_segments) ...
                                               sprintf('; 2, %g, %g ', segments{i}') ...
                                               sprintf(');\n\n')];
            else
                extra_code = [extra_code sprintf('curve2d Revolution2DCurve%d = (%g ; ', i, n_points) ...
                                               sprintf('%g, %g ; ', points{i}') ...
                                               sprintf(' %g ', n_segments) ...
                                               sprintf('; 3, %g, %g, %g ', segments{i}') ...
                                               sprintf(');\n\n')];
            end
            
            if (complement_flag(i))
                geo_code = [geo_code 'not '];
            end
            
            geo_code = [geo_code sprintf('revolution(%g, %g, %g ; %g, %g, %g ; Revolution2DCurve%d)', ...
                        axis_point_a(:, i), axis_point_b(:, i), i)];

            if (i < n_body_of_revolutions)
                geo_code = [geo_code operator_string];
            else
                geo_code = [geo_code ')'];             
            end
        end
    end
    
    eidors_msg('@@@ called with "%s" returning.', operator_string, 3);
    
function geo_code = parse_geometry_cone(cone, operator_string)

    eidors_msg('@@@ called with "%s" starting.', operator_string, 3);

    if (~isstruct(cone) || isempty(cone))
        error('Parameter cone must be a valid structure.');        
    else
        % Get number of cones.
        n_cones = numel(cone);
        
        % Get structure field names.
        field_names = fieldnames(cone);
        n_fields = numel(field_names);
        
        % Assign default values.
        top_radius      = 0.5*ones(1, n_cones);
        bottom_radius   = ones(1, n_cones);
        top_center      = [0;0;1]*ones(1, n_cones);
        bottom_center   = [0;0;0]*ones(1, n_cones);
        complement_flag = false(1, n_cones);
        
        % Parse all structure fields.
        for i = 1:n_fields
            switch (field_names{i})
                case 'top_radius'
                    top_radius = assign_radius(cone, n_cones, 'cone', field_names{i}, top_radius);
                case 'bottom_radius'
                    bottom_radius = assign_radius(cone, n_cones, 'cone', field_names{i}, bottom_radius);
                case {'top_center', 'top_centre'}
                    top_center = assign_point(cone, n_cones, 'cone', field_names{i}, top_center);
                case {'bottom_center', 'bottom_centre'}
                    bottom_center = assign_point(cone, n_cones, 'cone', field_names{i}, bottom_center);
                case 'complement_flag'
                    complement_flag = assign_flag(cone, n_cones, 'cone', field_names{i}, complement_flag);
                otherwise
                    error(['Field name ''%s'' is not valid for a cone.\Available field names for a cone are: ' ...
                           'bottom_center, bottom_radius, top_center, top_radius, and complement_flag.'], field_names{i});
            end
        end
        
        % Start geo code with an opening parenthesis.
        geo_code = '(';

        % Add geo code for each cone.
        for i = 1:n_cones
            if (complement_flag(i))
                geo_code = [geo_code 'not '];
            end
            
            n_vector = top_center(:,i) - bottom_center(:,i); 
            
            geo_code = [geo_code sprintf('(cone(%g, %g, %g ; %g ; %g, %g, %g ; %g) and plane(%g, %g, %g ; %g, %g, %g) and plane(%g, %g, %g ; %g, %g, %g))', ...
                        bottom_center(:,i), bottom_radius(i), top_center(:,i), top_radius(i), bottom_center(:,i), -n_vector, top_center(:,i), n_vector)];

            if (i < n_cones)
                geo_code = [geo_code operator_string];
            else
                geo_code = [geo_code ')'];             
            end
        end
    end
    
    eidors_msg('@@@ called with "%s" returning.', operator_string, 3);
    
function geo_code = parse_geometry_cylinder(cylinder, operator_string)

    eidors_msg('@@@ called with "%s" starting.', operator_string, 3);

    if (~isstruct(cylinder) || isempty(cylinder))
        error('Parameter cylinder must be a valid structure.');        
    else
        % Get number of cylinders.
        n_cylinders = numel(cylinder);
        
        % Get structure field names.
        field_names = fieldnames(cylinder);
        n_fields = numel(field_names);
        
        % Assign default values.
        radius          = ones(1, n_cylinders);
        top_center      = [0;0;1]*ones(1, n_cylinders);
        bottom_center   = [0;0;0]*ones(1, n_cylinders);
        complement_flag = false(1, n_cylinders);
        
        % Parse all structure fields.
        for i = 1:n_fields
            switch (field_names{i})
                case 'radius'
                    radius = assign_radius(cylinder, n_cylinders, 'cylinder', field_names{i}, radius);     
                case {'top_center', 'top_centre'}
                    top_center = assign_point(cylinder, n_cylinders, 'cylinder', field_names{i}, top_center);
                case {'bottom_center', 'bottom_centre'}
                    bottom_center = assign_point(cylinder, n_cylinders, 'cylinder', field_names{i}, bottom_center);
                case 'complement_flag'
                    complement_flag = assign_flag(cylinder, n_cylinders, 'cylinder', field_names{i}, complement_flag);
                otherwise
                    error(['Field name ''%s'' is not valid for a cylinder!\nAvailable field names for a cylinder are: '...
                           'bottom_center, top_center, radius, and complement_flag.'], field_names{i});
            end
        end
        
        % Start geo code with an opening parenthesis.
        geo_code = '(';

        % Add geo code for each cylinder.
        for i = 1:n_cylinders
            if (complement_flag(i))
                geo_code = [geo_code 'not '];
            end
            
            n_vector = top_center(:,i) - bottom_center(:,i); 
            
            geo_code = [geo_code sprintf('(cylinder(%g, %g, %g ; %g, %g, %g ; %g) and plane(%g, %g, %g ; %g, %g, %g) and plane(%g, %g, %g ; %g, %g, %g))', ...
                        bottom_center(:,i), top_center(:,i), radius(i), bottom_center(:,i), -n_vector, top_center(:,i), n_vector)];
                    
            if (i < n_cylinders)
                geo_code = [geo_code operator_string];
            else
                geo_code = [geo_code ')'];             
            end
        end
    end
    
    eidors_msg('@@@ called with "%s" returning.', operator_string, 3);
    
function geo_code = parse_geometry_ellipsoid(ellipsoid, operator_string)

    eidors_msg('@@@ called with "%s" starting.', operator_string, 3);

    if (~isstruct(ellipsoid) || isempty(ellipsoid))
        error('Parameter ellipsoid must be a valid structure.');        
    else
        % Get number of ellipsoids.
        n_ellipsoids = numel(ellipsoid);
        
        % Get structure field names.
        field_names = fieldnames(ellipsoid);
        n_fields = numel(field_names);
  
        % Assign default values.
        axis_a          = [1;0;0]*ones(1, n_ellipsoids);
        axis_b          = [0;1;0]*ones(1, n_ellipsoids);
        axis_c          = [0;0;1]*ones(1, n_ellipsoids);
        center          = [0;0;0]*ones(1, n_ellipsoids);
        complement_flag = false(1, n_ellipsoids);
        
        % Parse all structure fields.
        for i = 1:n_fields
            switch (field_names{i})   
                case {'center', 'centre'}
                    center = assign_point(ellipsoid, n_ellipsoids, 'ellipsoid', field_names{i}, center);
                case 'axis_a'
                    axis_a = assign_point(ellipsoid, n_ellipsoids, 'ellipsoid', field_names{i}, axis_a);
                case 'axis_b'
                    axis_b = assign_point(ellipsoid, n_ellipsoids, 'ellipsoid', field_names{i}, axis_b);
                case 'axis_c'
                    axis_c = assign_point(ellipsoid, n_ellipsoids, 'ellipsoid', field_names{i}, axis_c);
                case 'complement_flag'
                    complement_flag = assign_flag(ellipsoid, n_ellipsoids, 'ellipsoid', field_names{i}, complement_flag);
                otherwise
                     error(['Field name ''%s'' is not valid for an ellipsoid!\nAvailable field names for an ellipsoid are: '...
                           'center, axis_a, axis_b, axis_c, and complement_flag.'], field_names{i});
            end
        end
        
        % Start geo code with an opening parenthesis.
        geo_code = '(';

        % Add geo code for each ellipsoid.
        for i = 1:n_ellipsoids
            if (dot(axis_a(:,i), axis_b(:,i)) ~= 0)
                error('axis_a and axis_b have to be perpendicular for an ellipsoid.');
            elseif (dot(axis_a(:,i), axis_c(:,i)) ~= 0)
                error('axis_a and axis_c have to be perpendicular for an ellipsoid.');
            elseif (dot(axis_b(:,i), axis_c(:,i)) ~= 0)
                error('axis_b and axis_c have to be perpendicular for an ellipsoid.');
            end
            
            if (complement_flag(i))
                geo_code = [geo_code 'not '];
            end
                
            geo_code = [geo_code sprintf('ellipsoid(%g, %g, %g ; %g, %g, %g ; %g, %g, %g ; %g, %g, %g)', ...
                        center(:,i), axis_a(:, i), axis_b(:,i) , axis_c(:,i))];
                    
            if (i < n_ellipsoids)
                geo_code = [geo_code operator_string];
            else
                geo_code = [geo_code ')'];             
            end
        end
    end
    
    eidors_msg('@@@ called with "%s" returning.', operator_string, 3);
      
function geo_code = parse_geometry_elliptic_cylinder(elliptic_cylinder, operator_string)

    eidors_msg('@@@ called with "%s" starting.', operator_string, 3);

    if (~isstruct(elliptic_cylinder) || isempty(elliptic_cylinder))
        error('Parameter elliptic_cylinder must be a valid structure.');        
    else
        % Get number of elliptic_cylinders.
        n_elliptic_cylinders = numel(elliptic_cylinder);
        
        % Get structure field names.
        field_names = fieldnames(elliptic_cylinder);
        n_fields = numel(field_names);
        
        % Assign default values.
        top_center      = [0;0;1]*ones(1, n_elliptic_cylinders);
        bottom_center   = [0;0;0]*ones(1, n_elliptic_cylinders);
        axis_a          = [1;0;0]*ones(1, n_elliptic_cylinders);
        axis_b          = [0;1;0]*ones(1, n_elliptic_cylinders);
        complement_flag = false(1, n_elliptic_cylinders);
        
        % Parse all structure fields.
        for i = 1:n_fields
            switch (field_names{i})   
                case {'top_center', 'top_centre'}
                    top_center = assign_point(elliptic_cylinder, n_elliptic_cylinders, 'elliptic_cylinder', field_names{i}, top_center);
                case {'bottom_center', 'bottom_centre'}
                    bottom_center = assign_point(elliptic_cylinder, n_elliptic_cylinders, 'elliptic_cylinder', field_names{i}, bottom_center);
                case 'axis_a'
                    axis_a = assign_point(elliptic_cylinder, n_elliptic_cylinders, 'elliptic_cylinder', field_names{i}, axis_a);
                case 'axis_b'
                    axis_b = assign_point(elliptic_cylinder, n_elliptic_cylinders, 'elliptic_cylinder', field_names{i}, axis_b);
                case 'complement_flag'
                    complement_flag = assign_flag(elliptic_cylinder, n_elliptic_cylinders, 'elliptic_cylinder', field_names{i}, complement_flag);
                otherwise
                    error(['Field name ''%s'' is not valid for an elliptic cylinder!\nAvailable field names for an elliptic cylinder are: '...
                           'bottom_center, top_center, axis_a, axis_b, and complement_flag.'], field_names{i});
            end
        end
        
        % Start geo code with an opening parenthesis.
        geo_code = '(';

        % Add geo code for each cylinder.
        for i = 1:n_elliptic_cylinders
            if (complement_flag(i))
                geo_code = [geo_code 'not '];
            end
            
            central_axis = top_center(:,i) - bottom_center(:,i);
            
            if (dot(axis_a(:,i), axis_b(:,i)) ~= 0)
                error('axis_a and axis_b have to be perpendicular for an elliptic cylinder.');
            elseif (dot(axis_a(:,i), central_axis(:,i)) ~= 0)
                error('axis_a and the central axis have to be perpendicular for an elliptic cylinder.');
            elseif (dot(axis_b(:,i), central_axis(:,i)) ~= 0)
                error('axis_b and the central axis have to be perpendicular for an elliptic cylinder.');
            end
            
            geo_code = [geo_code sprintf('(ellipticcylinder(%g, %g, %g ; %g, %g, %g ; %g, %g, %g) and plane(%g, %g, %g ; %g, %g, %g) and plane(%g, %g, %g ; %g, %g, %g))', ...
                        bottom_center(:,i), axis_a(:,i), axis_b(:,i), bottom_center(:,i), -central_axis, top_center(:,i), central_axis)];
                    
            if (i < n_elliptic_cylinders)
                geo_code = [geo_code operator_string];
            else
                geo_code = [geo_code ')'];             
            end
        end
    end
    
    eidors_msg('@@@ called with "%s" returning.', operator_string, 3);
    
function geo_code = parse_geometry_half_space(half_space, operator_string)

    eidors_msg('@@@ called with "%s" starting.', operator_string, 3);

    if (~isstruct(half_space) || isempty(half_space))
        error('Parameter half_space must be a valid structure.');        
    else
        % Get number of half_spaces.
        n_half_spaces = numel(half_space);
        
        % Get structure field names.
        field_names = fieldnames(half_space);
        n_fields = numel(field_names);
        
        % Assign default values.
        point           = [0;0;0]*ones(1, n_half_spaces);
        outward_normal_vector  = [0;0;1]*ones(1, n_half_spaces);
        complement_flag = false(1, n_half_spaces);
        
        % Parse all structure fields.
        for i = 1:n_fields
            switch (field_names{i})  
                case 'point'
                    point = assign_point(half_space, n_half_spaces, 'half_space', field_names{i}, point);
                case 'outward_normal_vector'
                    outward_normal_vector = assign_point(half_space, n_half_spaces, 'half_space', field_names{i}, outward_normal_vector);
                case 'complement_flag'
                    complement_flag = assign_flag(half_space, n_half_spaces, 'half_space', field_names{i}, complement_flag);
                otherwise
                    error(['Field name ''%s'' is not valid for a half_space!\Available field names for a half_space are: ' ...
                           'point, outward_normal_vector, and complement_flag.'], field_names{i});
            end
        end
        
        % Start geo code with an opening parenthesis.
        geo_code = '(';

        % Add geo code for each half_space.
        for i = 1:n_half_spaces
            if (complement_flag(i))
                geo_code = [geo_code 'not '];
            end
            
            geo_code = [geo_code sprintf('plane(%g, %g, %g ; %g, %g, %g)', ...
                        point(:,i), outward_normal_vector(:,i))];
                    
            if (i < n_half_spaces)
                geo_code = [geo_code operator_string];
            else
                geo_code = [geo_code ')'];             
            end
        end
    end
    
    eidors_msg('@@@ called with "%s" returning.', operator_string, 3);

function geo_code = parse_geometry_ortho_brick(ortho_brick, operator_string)

    eidors_msg('@@@ called with "%s" starting.', operator_string, 3);

    if (~isstruct(ortho_brick) || isempty(ortho_brick))
        error('Parameter ortho_brick must be a valid structure.');        
    else
        % Get number of ortho_bricks.
        n_ortho_bricks = numel(ortho_brick);
        
        % Get structure field names.
        field_names = fieldnames(ortho_brick);
        n_fields = numel(field_names);
        
        % Assign default values.
        opposite_corner_a = [0;0;0]*ones(1, n_ortho_bricks);
        opposite_corner_b = [1;1;1]*ones(1, n_ortho_bricks);
        complement_flag   = zeros(1, n_ortho_bricks);
        
        % Parse all structure fields.
        for i = 1:n_fields
            switch (field_names{i})  
                case {'opposite_corner_a'}
                    opposite_corner_a = assign_point(ortho_brick, n_ortho_bricks, 'ortho_brick', field_names{i}, opposite_corner_a);
                case {'opposite_corner_b'}
                    opposite_corner_b = assign_point(ortho_brick, n_ortho_bricks, 'ortho_brick', field_names{i}, opposite_corner_b);
                case 'complement_flag'
                    complement_flag = assign_flag(ortho_brick, n_ortho_bricks, 'ortho_brick', field_names{i}, complement_flag);
                otherwise
                    error(['Field name "%s" is not allowed for an ortho_brick!\nAllowed field names for an ortho_brick are: ' ...
                           'opposite_corner_a, opposite_corner_b, and complement_flag.'], field_names{i});
            end
        end
        
        % Start geo code with an opening parenthesis.
        geo_code = '(';

        % Add geo code for each ortho_brick.
        for i = 1:n_ortho_bricks
            if (complement_flag(i))
                geo_code = [geo_code 'not '];
            end
            
            geo_code = [geo_code sprintf('orthobrick(%g, %g, %g ; %g, %g, %g)', ...
                        min([opposite_corner_a(:, i) opposite_corner_b(:, i)], [], 2), ...
                        max([opposite_corner_a(:, i) opposite_corner_b(:, i)], [], 2))];
                    
            if (i < n_ortho_bricks)
                geo_code = [geo_code operator_string];
            else
                geo_code = [geo_code ')'];             
            end
        end
    end
    
    eidors_msg('@@@ called with "%s" returning.', operator_string, 3);
    
function geo_code = parse_geometry_parallelepiped(parallelepiped, operator_string)

    eidors_msg('@@@ called with "%s" starting.', operator_string, 3);

    if (~isstruct(parallelepiped) || isempty(parallelepiped))
        error('Parameter parallelepiped must be a valid structure.');        
    else
        % Get number of parallelepiped.
        n_parallelepipeds = numel(parallelepiped);
        
        % Get structure field names.
        field_names = fieldnames(parallelepiped);
        n_fields = numel(field_names);
        
        % Assign default values.       
        vertex          = [0;0;0]*ones(1, n_parallelepipeds);
        vector_a        = [1;0;0]*ones(1, n_parallelepipeds);
        vector_b        = [0;1;0]*ones(1, n_parallelepipeds);
        vector_c        = [0;0;1]*ones(1, n_parallelepipeds);
        complement_flag = zeros(1, n_parallelepipeds);
        
        % Parse all structure fields.
        for i = 1:n_fields
            switch (field_names{i})  
                case 'vertex'
                    vertex = assign_point(parallelepiped, n_parallelepipeds, 'parallelepiped', field_names{i}, vertex);
                case 'vector_a'
                    vector_a = assign_point(parallelepiped, n_parallelepipeds, 'parallelepiped', field_names{i}, vector_a);
                case 'vector_b'
                    vector_b = assign_point(parallelepiped, n_parallelepipeds, 'parallelepiped', field_names{i}, vector_b);
                case 'vector_c'
                    vector_c = assign_point(parallelepiped, n_parallelepipeds, 'parallelepiped', field_names{i}, vector_c);
                case 'complement_flag'
                    complement_flag = assign_flag(parallelepiped, n_parallelepipeds, 'parallelepiped', field_names{i}, complement_flag);
                otherwise
                    error(['Field name "%s" is not allowed for a parallelepiped!\nAllowed field names for a parallelepiped are: ' ...
                           'vertex, vector_a, vector_b, vector_c, and complement_flag.'], field_names{i});
            end
        end
        
        % Start geo code with an opening parenthesis.
        geo_code = '(';

        % Add geo code for each parallelepiped.
        for i = 1:n_parallelepipeds
            if (complement_flag(i))
                geo_code = [geo_code 'not '];
            end
             
            % Make sure all vectors are not coplanar.
            if (abs(dot(vector_a(:,i), cross(vector_b(:,i), vector_c(:,i)))) < eps)
                error('parallelepiped(%d) description includes coplanar vectors.', i);
            end
            
            % Compute opposite vertex.
            opposite_vertex = vertex(:,i) + vector_a(:,i) + vector_b(:,i) + vector_c(:,i);
            
            % Compute normal vectors.
            n_vector_ab = cross(vector_a(:,i), vector_b(:,i));
            n_vector_ac = cross(vector_a(:,i), vector_c(:,i));
            n_vector_bc = cross(vector_b(:,i), vector_c(:,i));
            
            % Check normal vectors directions.
            if (dot(n_vector_ab, vector_c(:,i)) < 0)
                n_vector_ab = -n_vector_ab;
            end
            if (dot(n_vector_ac, vector_b(:,i)) < 0)
                n_vector_ac = -n_vector_ac;
            end
            if (dot(n_vector_bc, vector_a(:,i)) < 0)
                n_vector_bc = -n_vector_bc;
            end
            
            geo_code = [geo_code sprintf(['(plane(%g, %g, %g ; %g, %g, %g) and' ...
                                          ' plane(%g, %g, %g ; %g, %g, %g) and' ...
                                          ' plane(%g, %g, %g ; %g, %g, %g) and' ...
                                          ' plane(%g, %g, %g ; %g, %g, %g) and' ...
                                          ' plane(%g, %g, %g ; %g, %g, %g) and' ...
                                          ' plane(%g, %g, %g ; %g, %g, %g))'], ...
                        vertex(:,i), -n_vector_ab, ...
                        vertex(:,i), -n_vector_ac, ...
                        vertex(:,i), -n_vector_bc, ...
                        opposite_vertex, n_vector_ab, ...
                        opposite_vertex, n_vector_ac, ...
                        opposite_vertex, n_vector_bc)];
                    
            if (i < n_parallelepipeds)
                geo_code = [geo_code operator_string];
            else
                geo_code = [geo_code ')'];             
            end
        end
    end
    
    eidors_msg('@@@ called with "%s" returning.', operator_string, 3);
    
function geo_code = parse_geometry_sphere(sphere, operator_string)

    eidors_msg('@@@ called with "%s" starting.', operator_string, 3);

    if (~isstruct(sphere) || isempty(sphere))
        error('Parameter sphere must be a valid structure.');        
    else
        % Get number of spheres.
        n_spheres = numel(sphere);
        
        % Get structure field names.
        field_names = fieldnames(sphere);
        n_fields = numel(field_names);
  
        % Assign default values.
        radius          = ones(1, n_spheres);
        center          = [0;0;0]*ones(1, n_spheres);
        complement_flag = false(1, n_spheres);
        
        % Parse all structure fields.
        for i = 1:n_fields
            switch (field_names{i})
                case {'center', 'centre'}
                    center = assign_point(sphere, n_spheres, 'sphere', field_names{i}, center);
                case 'radius'
                    radius = assign_radius(sphere, n_spheres, 'sphere', field_names{i}, radius);     
                case 'complement_flag'
                    complement_flag = assign_flag(sphere, n_spheres, 'sphere', field_names{i}, complement_flag);
                otherwise
                    error(['Field name ''%s'' is not valid for a sphere.\nAvailable field names for a sphere are: ' ...
                           'center, radius, and complement_flag.'], field_names{i});
            end
        end
        
        % Start geo code with an opening parenthesis.
        geo_code = '(';

        % Add geo code for each sphere.
        for i = 1:n_spheres
            if (complement_flag(i))
                geo_code = [geo_code 'not '];
            end
                
            geo_code = [geo_code sprintf('sphere(%g, %g, %g ; %g)', ...
                        center(:,i), radius(i))];
                    
            if (i < n_spheres)
                geo_code = [geo_code operator_string];
            else
                geo_code = [geo_code ')'];             
            end
        end
    end
    
    eidors_msg('@@@ called with "%s" returning.', operator_string, 3);
                             
function write_geo_file(geo_fn, body_solid_code, electrode_solid_code, body_extra_code, electrode_extra_code, body_extra_param, electrode_extra_param)
    
    % Open geo file for writing.
    fid = fopen(geo_fn, 'w');
    
    if (fid == -1)
        error('Unable to open file %s for writing.', geo_fn);
    end
    
    % Write header for geo file.
    fprintf(fid, '#Automatically generated by ng_mk_geometric_models\n\n');
    fprintf(fid, 'algebraic3d\n\n');
    
    % Assemble a string to represent the union of all bodies.
    total_body_solid = '(';
   
    for i = 1:numel(body_solid_code)
        total_body_solid = [total_body_solid body_extra_param{i}.name];

        if (i < numel(body_solid_code))
            total_body_solid = [total_body_solid ' or '];
        else
            total_body_solid = [total_body_solid ')'];             
        end
    end
    
    % Assemble a string to represent the union of all electrodes entering the body.
    total_electrode_solid = '(';
    n_total_electrode_solid = 0;
   
    for i = 1:numel(electrode_solid_code)
        if (electrode_extra_param{i}.enter_body_flag)
            if (n_total_electrode_solid > 0)
                total_electrode_solid = [total_electrode_solid ' or '];
            end
            total_electrode_solid = [total_electrode_solid electrode_extra_param{i}.name];
            n_total_electrode_solid = n_total_electrode_solid + 1;
        end
    end
    total_electrode_solid = [total_electrode_solid ')'];   
    
    % Write body_extra_code and electrode_extra_code in geo file
    for i = 1:numel(body_extra_code)
        if (~isempty(body_extra_code{i}))
            fprintf(fid, body_extra_code{i});
        end
    end
    for i = 1:numel(electrode_extra_code)
        if (~isempty(electrode_extra_code{i}))
            fprintf(fid, electrode_extra_code{i});
        end
    end
    fprintf(fid, '\n');
 
    % Write electrode solids that enter the body in geo file.
    for i = 1:numel(electrode_solid_code)
        if (~isempty(electrode_solid_code{i}) && electrode_extra_param{i}.enter_body_flag)
            fprintf(fid, 'solid %s = %s;\n\n', electrode_extra_param{i}.name, electrode_solid_code{i});
        end
    end
    
    % Write body solids in geo file.
    for i = 1:numel(body_solid_code)
        if (n_total_electrode_solid == 0)
            fprintf(fid, 'solid %s = %s;\n\n', body_extra_param{i}.name, body_solid_code{i});
        else
            fprintf(fid, 'solid %s = not %s and %s;\n\n', body_extra_param{i}.name, total_electrode_solid, body_solid_code{i});            
        end
    end
 
    % Write electrode solids that do not enter the body in geo file.
    for i = 1:numel(electrode_solid_code)
        if (~isempty(electrode_solid_code{i}) && ~electrode_extra_param{i}.enter_body_flag)
            fprintf(fid, 'solid %s = not %s and %s;\n\n', electrode_extra_param{i}.name, total_body_solid, electrode_solid_code{i});
        end
    end
    
    % Write electrode tlos in geo file.
    for i = 1:numel(electrode_solid_code)
        if (~isempty(electrode_solid_code{i}))
            if (isinf(electrode_extra_param{i}.max_edge_length))
                fprintf(fid, 'tlo %s -col=[1,0,0] -material=%s;\n', electrode_extra_param{i}.name, electrode_extra_param{i}.name);
            else
                fprintf(fid, 'tlo %s -col=[1,0,0] -material=%s -maxh=%g;\n', electrode_extra_param{i}.name, electrode_extra_param{i}.name, electrode_extra_param{i}.max_edge_length);          
            end
        end
    end
    
    % Write body tlos in geo file.
    for i = 1:numel(body_solid_code)
        if (isinf(body_extra_param{i}.max_edge_length))
            fprintf(fid, 'tlo %s -col=[0,1,0] -material=%s;\n', body_extra_param{i}.name, body_extra_param{i}.name);
        else
            fprintf(fid, 'tlo %s -col=[0,1,0] -material=%s -maxh=%g;\n', body_extra_param{i}.name, body_extra_param{i}.name, body_extra_param{i}.max_edge_length);            
        end
    end
    
    % Close file.
    fclose(fid);

function mat = read_mat_from_file(fid, nrows, ncols)
    mat = fscanf(fid, '%g', [ncols, nrows])';

    % Skip to next line.
    if (~isempty(fgetl(fid)))
        error('Last line was only partialy read.');
    end
    
function fmdl = read_vol_file(vol_fn, electrode_extra_param)

    % Open file for reading.
    fid = fopen(vol_fn, 'r');

    if (fid == -1)
        error('Unable to open file %s for reading.', vol_fn);
    end
    
    % Read a first line in vol file.
    line = fgetl(fid);
   
    % While no EOF or "endmesh" keyword is found.
    while (ischar(line) && ~strcmp(line, 'endmesh'))
        
        % Parse every line if not comment or empty line
        if (~isempty(line) && line(1) ~= '#') % Supposing '#' is always the first character of a comment line.
            switch(line)
                case 'mesh3d'   % Nothing to do.
                case 'dimension'
                    dimension = read_mat_from_file(fid, 1, 1);
                    if (dimension ~= 3)
                        error('unknown dimension %g in vol file.', dimension);
                    end
                case 'geomtype'
                    geomtype = read_mat_from_file(fid, 1, 1);
                    if (geomtype ~= 0)
                        error('unknown %g geomtype in vol file.', geomtype);
                    end
                case 'surfaceelements'
                    %# surfnr    bcnr   domin  domout      np      p1      p2      p3
                    n_surface_elements = read_mat_from_file(fid, 1, 1);
                    if (n_surface_elements)
                        surface_elements   = read_mat_from_file(fid, n_surface_elements, 8);
                    else
                        error('vol file contains no surface elements. There is probably something wrong with the provided geometry description.');    
                    end
                case 'volumeelements'
                    %#  matnr      np      p1      p2      p3      p4
                    n_volume_elements = read_mat_from_file(fid, 1, 1);
                    if (n_volume_elements)
                        volume_elements   = read_mat_from_file(fid, n_volume_elements, 6);
                    else
                        error('vol file contains no volume elements. There is probably something wrong with the provided geometry description.');    
                    end
                case 'edgesegmentsgi2'
                    %# surfid  0   p1   p2   trignum1    trignum2   domin/surfnr1    domout/surfnr2   ednr1   dist1   ednr2   dist2 
                    n_edge_segments_sgi2 = read_mat_from_file(fid, 1, 1);
                    edge_segments_sgi2   = read_mat_from_file(fid, n_edge_segments_sgi2, 12);
                case 'points'
                    %#          X             Y             Z
                    n_points = read_mat_from_file(fid, 1, 1);
                    if (n_points)
                        points   = read_mat_from_file(fid, n_points, 3);
                    else
                        error('vol file contains no points. There is probably something wrong with the provided geometry description.');                       
                    end
                case 'materials'
                    n_materials = read_mat_from_file(fid, 1, 1);
                    if (n_materials)
                        materials   = cell(n_materials, 2);
                        % Read and parse each material line.
                        for i = 1:n_materials
                            material_line = fgetl(fid);
                            sscanf_result = sscanf(material_line, '%g%c%s')';
                            materials{i, 1} = sscanf_result(1);
                            materials{i, 2} = char(sscanf_result(3:end));
                        end
                    else
                        error('vol file contains no materials. There is probably something wrong with the provided geometry description.');                             
                    end
                case 'face_colours'
                    %#   Surfnr     Red     Green     Blue
                    n_face_colours = read_mat_from_file(fid, 1, 1);
                    face_colours   = read_mat_from_file(fid, n_face_colours, 4);
                otherwise
                    error('unknown "%s" line in vol file.', line);
            end
        end
        
        % Read next line in vol file.
        line = fgetl(fid);
    end
    
    % Close file.
    fclose(fid);
    
    if (~exist('points', 'var'))
        error('Point description is missing from vol file.');
    end
 
    if (~exist('volume_elements', 'var'))
        error('Volume element description is missing from vol file.');
    end
    
    if (~exist('surface_elements', 'var'))
        error('Surface element description is missing from vol file.');
    end
    
    if (~exist('materials', 'var'))
        error('Material description is missing from vol file.');
    end
    
    % Find electrode and body material indices.
    electrode_material = [];
    for i = 1:n_materials
        material_name   = materials{i, 2};
        material_number = materials{i, 1};
        
%         if (strncmp(material_name, 'electrode', 9))
%             % Extract electrode number from material_name
%             electrode_number = str2double(material_name(10:end));
%             electrode_material(electrode_number) = material_number;
%         end
        for j = 1:numel(electrode_extra_param)
            if (strcmp(material_name, electrode_extra_param{j}.name))
                electrode_material(j) = material_number;
            end
        end
    end
   
    % Remove electrode material if necessary
    original_n_nodes     = size(points, 1);
    original_n_elements  = size(volume_elements, 1);
    original_n_surfaces  = size(surface_elements, 1);
    original_n_materials = size(materials, 1);

    for i = 1:numel(electrode_material)
        if (~electrode_extra_param{i}.keep_material_flag)
            % Remove unwanted volume elements
            volume_elements(volume_elements(:, 1) == electrode_material(i), :) = [];

            % Remove unwanted surface elements
            surface_elements(surface_elements(:, 3) == electrode_material(i) & ...
                             surface_elements(:, 4) == 0 | ...
                             surface_elements(:, 4) == electrode_material(i) & ...
                             surface_elements(:, 3) == 0, :) = [];
        end
    end

    % Find nodes that are now unused.
    unused_nodes = true(1, size(points, 1));
    unused_nodes(volume_elements(:, 3:6))  = false;
    unused_nodes(surface_elements(:, 6:8)) = false;     

    % Remove unused points.
    points(unused_nodes, :) = [];

    % Compute new node indices after node removal. 
    new_node_index = (1:original_n_nodes) - cumsum(unused_nodes);   

    % Update node indices for surface and volume elements.
    surface_elements(:, 6:8) = new_node_index(surface_elements(:, 6:8));
    volume_elements(:, 3:6)  = new_node_index(volume_elements(:, 3:6));

    % Find materials that are now unused.
    unused_materials = true(1, size(materials, 1));
    unused_materials(volume_elements(:, 1)) = false;

    % Remove unused materials.
    materials(unused_materials, :) = [];

    % Compute new material indices after material removal. 
    new_material_index = (1:original_n_materials) - cumsum(unused_materials);   

    % Update material indices for volume elements.
    volume_elements(:, 1)  = new_material_index(volume_elements(:, 1));

    eidors_msg('@@@ Removed %d nodes, %d elements, %d surfaces and %d materials', ...
        original_n_nodes     - size(points, 1), ...
        original_n_elements  - size(volume_elements, 1), ...
        original_n_surfaces  - size(surface_elements, 1), ...
        original_n_materials - size(materials, 1), 3);
    
    % Assign mesh data to fmdl structures.
    fmdl.nodes            = points;
    fmdl.elems            = volume_elements(:, 3:6);
    fmdl.boundary         = surface_elements(:, 6:8);
    fmdl.boundary_numbers = surface_elements(:, 2);
    fmdl.mat_idx          = volume_elements(:, 1);
    fmdl.mat_name         = materials(:, 2);
    
    % Find electrode surfaces and nodes.
    for i = 1:numel(electrode_material)
        % Find surfaces that are part of the electrodes.
        if (electrode_extra_param{i}.keep_material_flag)
            electrode_boundary = ...
               sort(find(surface_elements(:, 3) == 0 & ...
                         surface_elements(:, 4) == electrode_material(i) | ...
                         surface_elements(:, 4) == 0 & ...
                         surface_elements(:, 3) == electrode_material(i)))';
        else
            electrode_boundary = ...
                sort(find(surface_elements(:, 3) == electrode_material(i) | ...
                          surface_elements(:, 4) == electrode_material(i)))';
        end
     
        if (isempty(electrode_boundary))
            eidors_msg('WARNING: Electrode #%04d has been removed since it does not contact any body.', i, 2);
        else
            fmdl.electrode(i).boundary = electrode_boundary;
            
            % Find nodes that are part of the electrodes.
            fmdl.electrode(i).nodes = ...
                unique(fmdl.boundary(fmdl.electrode(i).boundary(:), :))';                          

            % Assign default contact impedance.
            fmdl.electrode(i).z_contact = 0.01;
            
            if (~isempty(electrode_extra_param{i}.name))
                fmdl.electrode(i).name = electrode_extra_param{i}.name;
            end
        end
    end

function fmdl = complete_fmdl(fmdl, electrode_extra_param)
 
    % Find center point of domain.
    domain_center  = (max(fmdl.nodes)-min(fmdl.nodes))/2 + min(fmdl.nodes);
    domain_centers = ones(size(fmdl.nodes, 1), 1)*domain_center;
    
    % Find node closest to center for ground node.
    [unused, min_idx] = min(sum((fmdl.nodes - domain_centers).^2, 2));
    fmdl.gnd_node     = min_idx(1);

    fmdl.np_fwd_solve.perm_sym = '{n}';

    fmdl.name = 'ng_mk_geometric_models';

    fmdl.solve=      'eidors_default';
    fmdl.jacobian=   'eidors_default';
    fmdl.system_mat= 'eidors_default';

    fmdl.normalize_measurements = 0;
    
    for i = 1:numel(electrode_extra_param)
        if (isfield(electrode_extra_param{i}, 'point'))
            % Find center point of domain.
            electrode_points = ones(size(fmdl.nodes, 1), 1)*electrode_extra_param{i}.point';

            % Find node closest to the electrode point.
            [unused, min_idx]       = min(sum((fmdl.nodes - electrode_points).^2, 2));
            fmdl.electrode(i).nodes = min_idx(1);
            fmdl.electrode(i).boundary = [];

            % Assign default contact impedance.
            fmdl.electrode(i).z_contact = 0.01;
        end
    end

    fmdl = eidors_obj('fwd_model', fmdl);

function do_unit_test
    for tn = 1:do_test_number(0)
        eidors_msg('ng_mk_geometric_models: unit_test %02d', tn, 1);
        fmdl = do_test_number(tn);
        show_fem(fmdl);
        drawnow;
    end

function fmdl = do_test_number(tn)
    switch tn
        % Simple 3D cylinder. Radius = 1 with no electrodes
        case 1;
            body_geometry.cylinder = struct;
            fmdl = ng_mk_geometric_models(body_geometry);
        % Simple 3D cylinder. Radius = 1 with 16 spherical electrodes.
        case 2;
            body_geometry.cylinder = struct;
            n_elect = 16;
            theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
            for i = 1:n_elect
                electrode_geometry{i}.sphere.center = [cos(theta(i)) sin(theta(i)) 0.5];
                electrode_geometry{i}.sphere.radius = 0.1;
            end
            fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
        % Simple 3D cylinder. Radius = 1 with 16 cylindrical electrodes.
        case 3;
            body_geometry.cylinder = struct;
            n_elect = 16;
            theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
            for i = 1:n_elect
                electrode_geometry{i}.cylinder.top_center    = [1.03*cos(theta(i)) 1.03*sin(theta(i)) 0.5];
                electrode_geometry{i}.cylinder.bottom_center = [0.97*cos(theta(i)) 0.97*sin(theta(i)) 0.5];
                electrode_geometry{i}.cylinder.radius = 0.1;
            end
            fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
        case 4;
            body_geometry.cylinder = struct;
            n_elect = 16;
            theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
            for i = 1:n_elect
                electrode_geometry{i}.cylinder.top_center    = [1.03*cos(theta(i)) 1.03*sin(theta(i)) 0.5];
                electrode_geometry{i}.cylinder.bottom_center = [0.97*cos(theta(i)) 0.97*sin(theta(i)) 0.5];
                electrode_geometry{i}.cylinder.radius = 0.1;
                electrode_geometry{i}.keep_material_flag = 1;
            end
            fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
        case 5;
            body_geometry.cylinder = struct;
            n_elect = 16;
            theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
            for i = 1:n_elect
                electrode_geometry{i}.cylinder.top_center    = [1.03*cos(theta(i)) 1.03*sin(theta(i)) 0.5];
                electrode_geometry{i}.cylinder.bottom_center = [0.97*cos(theta(i)) 0.97*sin(theta(i)) 0.5];
                electrode_geometry{i}.cylinder.radius = 0.1;
                electrode_geometry{i}.enter_body_flag = 1;
            end
            fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
        case 6;
            body_geometry.cylinder = struct;
            n_elect = 16;
            theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
            for i = 1:n_elect
                electrode_geometry{i}.cylinder.top_center    = [1.03*cos(theta(i)) 1.03*sin(theta(i)) 0.5];
                electrode_geometry{i}.cylinder.bottom_center = [0.97*cos(theta(i)) 0.97*sin(theta(i)) 0.5];
                electrode_geometry{i}.cylinder.radius = 0.1;
                electrode_geometry{i}.keep_material_flag = 1;
                electrode_geometry{i}.enter_body_flag = 1;                
            end
            fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
        case 7;
            body_geometry.cylinder = struct;
            body_geometry.sphere.center = [0 0 1];
            n_elect = 16;
            theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
            for i = 1:n_elect
                electrode_geometry{i}.cylinder.top_center    = [1.03*cos(theta(i)) 1.03*sin(theta(i)) 0.5];
                electrode_geometry{i}.cylinder.bottom_center = [0.97*cos(theta(i)) 0.97*sin(theta(i)) 0.5];
                electrode_geometry{i}.cylinder.radius = 0.1;
            end
            fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
        case 8;
            body_geometry.cylinder  = struct;
            body_geometry.sphere(1) = struct;  
            body_geometry.sphere(2).center = [0 0 1];         
            n_elect = 16;
            theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
            for i = 1:n_elect
                electrode_geometry{i}.cylinder.top_center    = [1.03*cos(theta(i)) 1.03*sin(theta(i)) 0.5];
                electrode_geometry{i}.cylinder.bottom_center = [0.97*cos(theta(i)) 0.97*sin(theta(i)) 0.5];
                electrode_geometry{i}.cylinder.radius = 0.1;
            end
            fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);   
        case 9;
            body_geometry.intersection.cylinder(1) = struct;
            body_geometry.intersection.cylinder(2).radius     = 0.5;
            body_geometry.intersection.cylinder(2).complement_flag = 1;   
            n_elect = 16;
            theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
            for i = 1:n_elect
                electrode_geometry{i}.cylinder.top_center    = [1.03*cos(theta(i)) 1.03*sin(theta(i)) 0.5];
                electrode_geometry{i}.cylinder.bottom_center = [0.97*cos(theta(i)) 0.97*sin(theta(i)) 0.5];
                electrode_geometry{i}.cylinder.radius = 0.1;
            end
            fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
        case 10;
            body_geometry.intersection(1).sphere(1).radius     = 0.5;
            body_geometry.intersection(1).sphere(1).center     = [0 0 2];
            body_geometry.intersection(1).sphere(1).complement_flag = 1;
            body_geometry.intersection(1).sphere(2).center     = [0 0 2];
            body_geometry.intersection(2).cylinder(1).top_center = [0 0 2];
            body_geometry.intersection(2).cylinder(2).radius     = 0.5;
            body_geometry.intersection(2).cylinder(2).top_center = [0 0 2];
            body_geometry.intersection(2).cylinder(2).complement_flag = 1;   
            n_elect = 16;
            theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
            for i = 1:n_elect
                electrode_geometry{i}.cylinder.top_center    = [1.03*cos(theta(i)) 1.03*sin(theta(i)) 0.5];
                electrode_geometry{i}.cylinder.bottom_center = [0.97*cos(theta(i)) 0.97*sin(theta(i)) 0.5];
                electrode_geometry{i}.cylinder.radius = 0.1;
            end
            fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
        case 11;
            body_geometry.intersection.union(1).sphere.radius = 0.5;
            body_geometry.intersection.union(1).sphere.center = [0 0 2];
            body_geometry.intersection.union(1).cylinder.radius = 0.5;
            body_geometry.intersection.union(1).cylinder.top_center = [0 0 2];
            body_geometry.intersection.union(1).complement_flag = 1;
            body_geometry.intersection.union(2).sphere.center = [0 0 2];
            body_geometry.intersection.union(2).cylinder.top_center = [0 0 2]; 
            n_elect = 16;
            theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
            for i = 1:n_elect
                electrode_geometry{i}.cylinder.top_center    = [1.03*cos(theta(i)) 1.03*sin(theta(i)) 0.5];
                electrode_geometry{i}.cylinder.bottom_center = [0.97*cos(theta(i)) 0.97*sin(theta(i)) 0.5];
                electrode_geometry{i}.cylinder.radius = 0.1;
            end
            fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
        case 12;
            body_geometry.cone = struct; 
            n_elect = 16;
            theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
            for i = 1:n_elect
                electrode_geometry{i}.cylinder.top_center    = [0.85*cos(theta(i)) 0.85*sin(theta(i)) 0.5];
                electrode_geometry{i}.cylinder.bottom_center = [0.65*cos(theta(i)) 0.65*sin(theta(i)) 0.5];
                electrode_geometry{i}.cylinder.radius = 0.1;
            end
            fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
        case 13;
            body_geometry.cone = struct; 
            n_elect = 16;
            theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
            for i = 1:n_elect
                electrode_geometry{i}.sphere.center = [0.75*cos(theta(i)) 0.75*sin(theta(i)) 0.5];
                electrode_geometry{i}.sphere.radius = 0.1;
            end
            fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
        case 14;
            body_geometry.cone(1).top_center = [0 0 1.5];
            body_geometry.cone(1).bottom_center = [0 0 0.5];
            body_geometry.cone(2).top_center = [0 0 -1.5];
            body_geometry.cone(2).bottom_center = [0 0 -0.5];
            body_geometry.cylinder.top_center    = [0, 0, 0.5];
            body_geometry.cylinder.bottom_center = [0, 0, -0.5];
            n_elect = 16;
            theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
            for i = 1:n_elect
                electrode_geometry{i}.sphere.center = [0.75*cos(theta(i)) 0.75*sin(theta(i)) 1.0];
                electrode_geometry{i}.sphere.radius = 0.1;
                electrode_geometry{i + n_elect}.sphere.center = [cos(theta(i)) sin(theta(i)) 0];
                electrode_geometry{i + n_elect}.sphere.radius = 0.15;
                electrode_geometry{i + 2*n_elect}.sphere.center = [0.75*cos(theta(i)) 0.75*sin(theta(i)) -1.0];
                electrode_geometry{i + 2*n_elect}.sphere.radius = 0.1;
            end
            fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
        case 15
            body_geometry.ortho_brick = struct;
            fmdl = ng_mk_geometric_models(body_geometry);
        case 16
            body_geometry.intersection.ortho_brick.opposite_corner_a = [0 0 0];
            body_geometry.intersection.ortho_brick.opposite_corner_b = [5 5 4];
            for i = 1:4; 
                for j = 1:4; 
                    body_geometry.intersection.cylinder(i,j).radius = 0.15;
                    body_geometry.intersection.cylinder(i,j).top_center = [i, j, 4];
                    body_geometry.intersection.cylinder(i,j).bottom_center = [i, j, 2];
                    body_geometry.intersection.cylinder(i,j).complement_flag = 1;
                end; 
            end;
            fmdl = ng_mk_geometric_models(body_geometry);    
        case 17
            body_geometry.intersection.ortho_brick.opposite_corner_a = [0 0 0];
            body_geometry.intersection.ortho_brick.opposite_corner_b = [5 5 4];
            for i = 1:4; 
                for j = 1:4; 
                    body_geometry.intersection.cylinder(i, j).radius = 0.15;
                    body_geometry.intersection.cylinder(i, j).top_center    = [i, j, 4];
                    body_geometry.intersection.cylinder(i, j).bottom_center = [i, j, 2];
                    body_geometry.intersection.cylinder(i, j).complement_flag = 1;
                    electrode_geometry{i, j, 1}.cylinder.radius        = 0.2;
                    electrode_geometry{i, j, 1}.cylinder.top_center    = [i, j, 3.1];
                    electrode_geometry{i, j, 1}.cylinder.bottom_center = [i, j, 2.9];
                    electrode_geometry{i, j, 2}.cylinder.radius        = 0.2;
                    electrode_geometry{i, j, 2}.cylinder.top_center    = [i, j, 2.2];
                    electrode_geometry{i, j, 2}.cylinder.bottom_center = [i, j, 2.0];
                end; 
            end;
            fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
        case 18
            body_geometry.parallelepiped  = struct;
            body_geometry.max_edge_length = 0.15;
            fmdl = ng_mk_geometric_models(body_geometry);
        case 19
            body_geometry.parallelepiped.vertex   = [ 0;  0;  0];
            body_geometry.parallelepiped.vector_a = [ 1;  1;  0];
            body_geometry.parallelepiped.vector_b = [ 0;  1;  1];
            body_geometry.parallelepiped.vector_c = [ 1;  0;  1];
            body_geometry.max_edge_length = 0.15;
            fmdl = ng_mk_geometric_models(body_geometry);
        case 20
            body_geometry.intersection.ortho_brick.opposite_corner_a = [-15, -15, 0];
            body_geometry.intersection.ortho_brick.opposite_corner_b = [15, 15, 5];
            body_geometry.intersection.half_space.point = [0, 0, 5];
            body_geometry.intersection.half_space.outward_normal_vector = [-1, -1, 5];
            
            fmdl = ng_mk_geometric_models(body_geometry);
        case 21
            body_geometry.ellipsoid.axis_a = [1 0 0];
            body_geometry.ellipsoid.axis_b = [0 2 0];
            body_geometry.ellipsoid.axis_c = [0 0 3];
            fmdl = ng_mk_geometric_models(body_geometry);   
        case 22
            body_geometry.ellipsoid.axis_a = [1 0 0];
            body_geometry.ellipsoid.axis_b = [0 1 1];
            body_geometry.ellipsoid.axis_c = [0 -2 2];
            fmdl = ng_mk_geometric_models(body_geometry);   
        case 23
            body_geometry.elliptic_cylinder.top_center = [0, 0, 10];
            body_geometry.elliptic_cylinder.bottom_center = [0, 0, 0];           
            body_geometry.elliptic_cylinder.axis_a = [1 0 0];
            body_geometry.elliptic_cylinder.axis_b = [0 2 0];  
            fmdl = ng_mk_geometric_models(body_geometry);
        case 24
            body_geometry.elliptic_cylinder.top_center = [0, 5, 5];
            body_geometry.elliptic_cylinder.bottom_center = [0, 0, 0];           
            body_geometry.elliptic_cylinder.axis_a = [1 0 0];
            body_geometry.elliptic_cylinder.axis_b = [0 -2 2];  
            fmdl = ng_mk_geometric_models(body_geometry);
        case 25
            body_geometry.body_of_revolution = struct;
            body_geometry.max_edge_length = 0.15;
            fmdl = ng_mk_geometric_models(body_geometry);
        case 26
            body_geometry.body_of_revolution.points   = [1 1; 1 2; 2 1.5; 2 1];
            body_geometry.body_of_revolution.segments = [1 2; 2 3; 3 4; 4 1];
            body_geometry.max_edge_length = 0.15;
            fmdl = ng_mk_geometric_models(body_geometry);
        case 27
            n_points = 24;
            theta = linspace(0, 2*pi, n_points+1)'; theta(end) = [];
            body_geometry.body_of_revolution.points   = 2 + [sin(theta) cos(theta)];
            body_geometry.body_of_revolution.segments = [(1:n_points)' [(2:n_points) 1]'];
            body_geometry.max_edge_length = 0.15;
            fmdl = ng_mk_geometric_models(body_geometry);
        case 28
            n_points = 24;
            theta = linspace(0, 2*pi, n_points+1)'; theta(end) = [];
            body_geometry.body_of_revolution.points   = 2 + [sin(theta) cos(theta)];
            body_geometry.body_of_revolution.segments = [(1:2:n_points)' (2:2:n_points)' [(3:2:n_points) 1]'];
            body_geometry.max_edge_length = 0.15;
            fmdl = ng_mk_geometric_models(body_geometry);
        case 29
            body_geometry{1}.cylinder(1).radius        = 0.5;
            body_geometry{1}.cylinder(1).top_center    = [0 0 0.75];
            body_geometry{1}.cylinder(1).bottom_center = [0 0 0.25];
            body_geometry{1}.name                      = 'Object';           
            body_geometry{2}.cylinder(2).radius        = 1;
            body_geometry{2}.name                      = 'Tank';
            n_elect = 16;
            theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
            for i = 1:n_elect
                electrode_geometry{i}.cylinder.top_center    = [1.03*cos(theta(i)) 1.03*sin(theta(i)) 0.5];
                electrode_geometry{i}.cylinder.bottom_center = [0.97*cos(theta(i)) 0.97*sin(theta(i)) 0.5];
                electrode_geometry{i}.cylinder.radius = 0.1;
            end
            fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
        case 30
            body_geometry{1}.sphere.radius     = 0.25;
            body_geometry{1}.sphere.center     = [0 0 0.5];
            body_geometry{1}.name              = 'Sphere';
            body_geometry{2}.cylinder.radius   = 1;
            body_geometry{2}.name              = 'Tank';           
            n_elect = 16;
            theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
            for i = 1:n_elect
                electrode_geometry{i}.sphere.center = [cos(theta(i)) sin(theta(i)) 0.5];
                electrode_geometry{i}.sphere.radius = 0.1;
            end
            fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
       case 31
            n_sphere = 8;
            theta = linspace(0, 2*pi, n_sphere+1); theta(end) = [];   
            for i = 1:n_sphere
                body_geometry{i}.sphere.radius   = 0.2;
                body_geometry{i}.sphere.center   = [0.65*cos(theta(i)) 0.65*sin(theta(i)) 0.5];  
                body_geometry{i}.max_edge_length = 0.025*(1 + rem(i,2));
                body_geometry{i}.name            = sprintf('Sphere%d', i);  
            end        
            body_geometry{n_sphere+1}.cylinder.radius = 1;
            body_geometry{n_sphere+1}.name            = 'Tank';  
            n_elect = 16;
            theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
            for i = 1:n_elect
                electrode_geometry{i}.sphere.center = [cos(theta(i)) sin(theta(i)) 0.5];
                electrode_geometry{i}.sphere.radius = 0.1;
                electrode_geometry{i}.max_edge_length = 0.025*(1 + rem(i,2));
            end
            fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
       case 32
            body_geometry.cylinder = struct;
            n_elect = 16;
            theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
            for i = 1:n_elect
                electrode_geometry{i}.point = [cos(theta(i)) sin(theta(i)) 0.5];
            end
            fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
       case 33     
            body_geometry.cylinder = struct;
            n_elect = 16;
            theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
            for i = 1:n_elect
                if (rem(i,2))
                    electrode_geometry{i}.point = [cos(theta(i)) sin(theta(i)) 0.5];
                    electrode_geometry{i}.name  = sprintf('Point_Electrode%d', ceil(i/2));
                else
                    electrode_geometry{i}.sphere.center = [cos(theta(i)) sin(theta(i)) 0.5];
                    electrode_geometry{i}.sphere.radius = 0.1;
                    electrode_geometry{i}.name          = sprintf('Circular_Electrode%d', floor(i/2));
                end
            end
            fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
       case 34
            body_geometry.body_of_extrusion = struct;
            body_geometry.max_edge_length = 0.15;
            fmdl = ng_mk_geometric_models(body_geometry);
        case 35
            body_geometry.body_of_extrusion.path_points   = [0 0 0; 0.25 0 1; 0.25 0 2; 0.25 0 3; 0 0 4];
            body_geometry.body_of_extrusion.path_segments = [1 2; 2 3; 3 4; 4 5];
            body_geometry.max_edge_length = 0.15;
            fmdl = ng_mk_geometric_models(body_geometry);
       case 36
            n_points = 16;
            theta = linspace(0, 2*pi, n_points+1)'; theta(end) = [];
            body_geometry.body_of_extrusion.profile_points   = 0.2*(2 + [0.75*sin(theta) cos(theta)]);
            body_geometry.body_of_extrusion.profile_segments = [(1:n_points)' [(2:n_points) 1]'];
            n_points = 32;
            theta = linspace(0, 2*pi, n_points+1)'; theta(end) = [];          
            body_geometry.body_of_extrusion.path_points   = 1*(2 + [sin(theta) 1.5*cos(theta) zeros(n_points, 1)]);
            body_geometry.body_of_extrusion.path_segments = [(1:n_points)' [(2:n_points) 1]'];
            body_geometry.body_of_extrusion.vector_d      = [0; 0; 1];
            body_geometry.max_edge_length = 0.15;
            fmdl = ng_mk_geometric_models(body_geometry); 
        case 0; fmdl = 36; % Return maximum number of tests.
        otherwise;
            error('Invalid test number.')
    end
