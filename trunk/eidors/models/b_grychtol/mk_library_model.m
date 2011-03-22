function out = mk_library_model(shape,elec_pos,elec_shape,maxsz)


if ischar(shape) && strcmp(shape,'LIBRARY_PATH')
        switch nargin
            case 1
                out = get_path;
            case 2
                set_path(elec_pos);
        end
elseif ischar(shape) && strcmp(shape, 'UNIT_TEST')
    out = do_unit_test; return;
else
    if ~iscell(shape) 
        shape = {shape, 'boundary'}; 
    elseif numel(shape) == 1
        shape{2} = 'boundary';
    end
    fname = make_filename(shape,elec_pos,elec_shape,maxsz);
    fname = [get_path '/' fname '.mat'];
    if exist(fname,'file')
        eidors_msg('MK_LIBRARY_MODEL: Using stored model');
        load(fname);
        out = fmdl;
    else
       s_shape = split_var_strings(shape(2:end));
       shapes = shape_library('get',shape{1},s_shape(1,:));
       if ~iscell(shapes), shapes = {shapes}; end
       %apply any indeces specified
       for i = 1:numel(shapes)
           eval(['shapes{i} = shapes{i}' s_shape{2,i} ';']);
       end
       [fmdl, mat_idx] = ng_mk_extruded_model({1,shapes,[4,50],maxsz},...
           elec_pos,elec_shape);
       fmdl.mat_idx = mat_idx;
       save(fname,'fmdl');
       out = fmdl;
    end
end

function str = make_filename(shape, elec_pos, elec_shape, maxsz)
%at this point, shape is a cell array of strings e.g. {'pig_23kg','lungs')
str = shape{1};
shape(1) = []; %remove the first element
shape = sort(shape); %sort the others
for i = 1:numel(shape)
    str = [str '_' shape{i}];
end
str = [str '_EP'];
for i = 1:numel(elec_pos)
    str = [str '_' num2str(elec_pos(i))];
end
str = [str '_ES'];
for i = 1:numel(elec_shape)
    str = [str '_' num2str(elec_shape(i))];
end
if ~isempty(maxsz)
    str = [str '_maxsz_' num2str(maxsz)];
end

%remove illegal characters
str = genvarname(str);

function clean = split_var_strings(strc)
for i = 1:numel(strc)
    [clean{1,i} clean{2,i}] = strtok(strc{i},'([{');
end


function out = get_path
global eidors_paths
out = eidors_paths.model_cache;

function set_path(val)
global eidors_paths
eidors_paths.model_cache = val;
%if folder doesn't exist, create it
if ~exist(val,'dir')
    mkdir(val);
end

function out = do_unit_test

out = mk_library_model({'pig_23kg','boundary','lungs(1:3:end,:)'},[16 0.05 0.5],[0.1],0.1);
