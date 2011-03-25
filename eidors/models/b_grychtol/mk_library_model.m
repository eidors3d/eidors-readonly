function out = mk_library_model(shape,elec_pos,elec_shape,maxsz)
%MK_LIBRARY_MODEL [Experimental] FEM models based on library shapes 
%
% MK_LIBRARY_MODEL(shape,elec_pos,elec_shape,maxsz) where:
%   shape -  a cell array of strings and
%     shape{1} is the shape_library model used (run shape_libary('list') to
%     get a list
%     shape{2} is the the boundary. If absent, 'boundary' is assumed.
%     shape{3..} are strings specifying additional inclusions (such as
%     lungs)
%   elec_pos - a vector specifying electrode positions. See
%     NG_MK_EXTRUDED_MODEL for details. To use the electrode positions
%     stored in the 'electrode' field in the shape_libary, specify elec_pos
%     as 'original'
%   elec_shape - a vector specifying electrode shapes. See
%     NG_MK_EXTRUDED_MODEL for details.
%   maxsz - maximum FEM size (default: course mesh)
%
% QUICK ACCESS TO COMMON MODELS:
%   MK_LIBRARY_MODEL(str) where str is a single string specifying a model.
%   Use MK_LIBRARY_MODEL('list') to obtain a list of available models.

% (C) 2011 Bartlomiej Grychtol. License: GPL version 2 or version 3
% $Id$

if ischar(shape)
    switch shape
        case 'LIBRARY_PATH'
            switch nargin
                case 1
                    out = get_path;
                case 2
                    set_path(elec_pos);
            end
        case 'list'
            out = list_predef_model_strings;
        case 'UNIT_TEST'
            out = do_unit_test; return;
        otherwise
            out = predef_model(shape);
    end
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
       if ischar(elec_pos) && strcmp(elec_pos,'original')
          el = shape_library('get',shape{1},'electrodes');
          electh= atan2(el(:,2),el(:,1))*180/pi;
          elec_pos = [electh,0.5*ones(size(electh))];
       end
       [fmdl, mat_idx] = ng_mk_extruded_model({1,shapes,[4,50],maxsz},...
           elec_pos,elec_shape);
       fmdl.mat_idx = mat_idx;
       save(fname,'fmdl');
       out = fmdl;
    end
end

%%%%%
% Lists predefined models (append when adding)
function out = list_predef_model_strings
out = {'pig_23kg_16el';
    'pig_23kg_32el';
    'pig_23kg_16el_lungs';
    'pig_23kg_32el_lungs'};

%%%%%
% Use predefined model
function out = predef_model(str)
switch str
    case 'pig_23kg_16el'
        out = mk_library_model({'pig_23kg','boundary'},...
            [16 1 0.5],[0.05],0.08);
    case 'pig_23kg_32el'
        out = mk_library_model({'pig_23kg','boundary'},...
            [32 1 0.5],[0.05],0.08);
    case 'pig_23kg_16el_lungs'
        out = mk_library_model({'pig_23kg','boundary','lungs(1:2:end,:)'},...
            [16 1 0.5],[0.05],0.08);
    case 'pig_23kg_32el_lungs'
        out = mk_library_model({'pig_23kg','boundary','lungs(1:2:end,:)'},...
            [32 1 0.5],[0.05],0.08);
    otherwise
        error('No such model');
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
if ischar(elec_shape)
    str = [str '_' elec_shape];
else
    for i = 1:numel(elec_shape)
        str = [str '_' num2str(elec_shape(i))];
    end
end
if ~isempty(maxsz)
    str = [str '_maxsz_' num2str(maxsz)];
end

%remove colons
str = strrep(str,':','-');

function clean = split_var_strings(strc)
for i = 1:numel(strc)
    [clean{1,i} clean{2,i}] = strtok(strc{i},'([{');
end


function out = get_path
global eidors_objects
out = eidors_objects.model_cache;

function set_path(val)
global eidors_objects
eidors_objects.model_cache = val;
%if folder doesn't exist, create it
if ~exist(val,'dir')
    mkdir(val);
end

function out = do_unit_test
models = mk_library_model('list');
for i = 1:numel(models)
    mdl = mk_library_model(models{i});
    img = mk_image(mdl,1);
    if numel(mdl.mat_idx) >1
        img.elem_data(mdl.mat_idx{2:end}) = 0.25;
    end
    figure
    show_fem(img,[0,1,0]);
    title(models{i},'Interpreter','none');
end


out = mk_library_model({'pig_23kg','boundary','lungs(1:2:end,:)'},[32 0 0.5],[0.05],0.1);

