function out = mk_library_model(shape,elec_pos,elec_shape,maxsz,nfft)
%MK_LIBRARY_MODEL - FEM models based on library shapes 
%
% MK_LIBRARY_MODEL(shape,elec_pos,elec_shape,maxsz,nfft) where:
%   shape -  a cell array of strings and
%     shape{1} is the shape_library model used
%          (run shape_libary('list') to get a list
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
%   nfft  - number of points to create along the boundary (default: 50)
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
       if nargin < 5
           nfft = 50;
       end
       [fmdl, mat_idx] = ng_mk_extruded_model({1,shapes,[4,nfft],maxsz},...
           elec_pos,elec_shape);
       fmdl.mat_idx = mat_idx;
       save(fname,'fmdl');
       out = fmdl;
    end
end

%%%%%
% Lists predefined models (append when adding)
function out = list_predef_model_strings
out = {
    'adult_male_16el';
    'adult_male_32el';
    'adult_male_16el_lungs';
    'adult_male_32el_lungs';
    'cylinder_16x1el_coarse';
    'cylinder_16x1el_fine';
    'cylinder_16x1el_vfine';   
    'cylinder_16x2el_coarse';
    'cylinder_16x2el_fine';
    'cylinder_16x2el_vfine';
    'neonate_16el';
    'neonate_32el';
    'neonate_16el_lungs';
    'neonate_32el_lungs';
    'pig_23kg_16el';
    'pig_23kg_32el';
    'pig_23kg_16el_lungs';
    'pig_23kg_32el_lungs';
    };

%%%%%
% Use predefined model
function out = predef_model(str)
switch str
    case 'adult_male_16el'
        out = mk_library_model({'adult_male','boundary'},...
            [16 1 0.5],[0.05],0.08);
    case 'adult_male_32el'
        out = mk_library_model({'adult_male','boundary'},...
            [32 1 0.5],[0.05],0.08);
    case 'adult_male_16el_lungs'
        out = mk_library_model({'adult_male','boundary','left_lung','right_lung'},...
            [16 1 0.5],[0.05],0.08);
    case 'adult_male_32el_lungs'
        out = mk_library_model({'adult_male','boundary','left_lung','right_lung'},...
            [32 1 0.5],[0.05],0.08);
        
    case 'cylinder_16x1el_coarse'
        out = ng_mk_cyl_models([10,15],[16,5],[0.5,0,0.18]);
    case 'cylinder_16x1el_fine'  
        out = ng_mk_cyl_models([10,15,1.1],[16,5],[0.5,0,0.15]);
    case 'cylinder_16x1el_vfine' 
        out = ng_mk_cyl_models([10,15,0.8],[16,5],[0.5,0,0.08]);
    case 'cylinder_16x2el_coarse' 
        out = ng_mk_cyl_models([30,15],[16,10,20],[0.5,0,0.18]);
    case 'cylinder_16x2el_fine'  
        out = ng_mk_cyl_models([30,15,1.5],[16,10,20],[0.5,0,0.15]);
    case 'cylinder_16x2el_vfine' 
        out = ng_mk_cyl_models([30,15,0.8],[16,10,20],[0.5,0,0.08]);
        
        
    case 'neonate_16el'
        out = mk_library_model({'neonate','boundary'},[16 1 0.5],0.1,0.08,49);
    case 'neonate_32el'
        out = mk_library_model({'neonate','boundary'},[32 1 0.5],0.07,0.08,49);
    case 'neonate_16el_lungs'
        out = mk_library_model({'neonate','boundary','left_lung','right_lung'},[16 1 0.5],0.1,0.08,49);
    case 'neonate_32el_lungs'
        out = mk_library_model({'neonate','boundary','left_lung','right_lung'},[32 1 0.5],0.07,0.08,49);
        
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
            [32 1 0.5],[0.05],0.08);    case 'pig_23kg_16el'
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
%give the model a name
out.name = str;


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
    ver= eidors_obj('interpreter_version');
    if ver.ver<4 || ver.ver>=7
       mkdir(val);
    else
 % matlab 6.x has a stupid mkdir function
       system(['mkdir ',val]); 
    end
end

function out = do_unit_test
models = mk_library_model('list');
for i = 1:numel(models)
    mdl = mk_library_model(models{i});
    img = mk_image(mdl,1);
    try   n = numel(mdl.mat_idx); catch n =1; end
    if n >1
        for j = 2:n
            img.elem_data(mdl.mat_idx{j}) = 0.25;
        end
    end
    figure
    show_fem(img,[0,1,0]);
    title(models{i},'Interpreter','none');
end


out = mk_library_model({'pig_23kg','boundary','lungs(1:2:end,:)'},[32 1 0.5],[0.05],0.08);

