function out = mk_library_model(shape,elec_pos,elec_shape,maxsz,nfft,scale)
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
%     If nfft==0, no interpolation takes place.
%   scale - avoids some Netgen issues by scaling the contours before 
%     calling netgen and scaling the resulting model back afterwards
%     (default: 1). Note that electrode and maxh specifications are not
%     scaled.
%
% QUICK ACCESS TO COMMON MODELS:
%   MK_LIBRARY_MODEL(str) where str is a single string specifying a model.
%   Use MK_LIBRARY_MODEL('list') to obtain a list of available models.
%
% PATH TO LIBRARY MODELS
%   'LIBRARY_PATH' - get or set library path

% (C) 2011 Bartlomiej Grychtol. License: GPL version 2 or version 3
% $Id$

% Fill in defaults:
if nargin < 6; scale = 1;          end
if nargin < 5; nfft = 50;          end

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
            out = mdl_normalize(out, 0); % not normalized by default
    end
else
   if ~iscell(shape)
      shape = {shape, 'boundary'};
   elseif numel(shape) == 1
      shape{2} = 'boundary';
   end
   fname = make_filename(shape,elec_pos,elec_shape,maxsz, nfft, scale);
   out = load_stored_model(fname);
   if ~isempty(out)
      return
   end
   s_shape = split_var_strings(shape(2:end));
   shapes = shape_library('get',shape{1},s_shape(1,:));
   if ~iscell(shapes), shapes = {shapes}; end
   %apply any indeces specified
   for i = 1:numel(shapes)
      eval(sprintf('shapes{i} = %f*shapes{i}%s;',scale,s_shape{2,i}));
   end
   if ischar(elec_pos) && strcmp(elec_pos,'original')
      el = shape_library('get',shape{1},'electrodes');
      electh= atan2(el(:,2),el(:,1))*180/pi;
      elec_pos = [electh,0.5*ones(size(electh))];
   end
   
   if nfft > 0
      [fmdl, mat_idx] = ng_mk_extruded_model({scale,shapes,[4,nfft],maxsz},...
         elec_pos,elec_shape);
   else
      [fmdl, mat_idx] = ng_mk_extruded_model({scale,shapes,0,maxsz},...
         elec_pos,elec_shape);
   end
   fmdl.nodes = fmdl.nodes/scale;
   fmdl.mat_idx = mat_idx;
   store_model(fmdl,fname)
   out = fmdl;
   out = mdl_normalize(out, 0); % not normalized by default
end




function out = load_stored_model(fname)
out = [];
fname = [get_path '/' fname '.mat'];
if exist(fname,'file')
   eidors_msg('MK_LIBRARY_MODEL: Using stored model');
   if exist('OCTAVE_VERSION');
      load(file_in_loadpath(fname));
   else
      load(fname);
   end
   out = fmdl;
end

function store_model(fmdl,fname)
   fname = [get_path '/' fname '.mat'];
   if exist('OCTAVE_VERSION');
      savver = '-v7';
   else
      savver = '-v7.3';  
   end
   save(fname,savver,'fmdl');

function out = build_if_needed(cmd,str)
out = load_stored_model(str);
if isempty(out)
   if ~iscell(cmd)
      cmd = {cmd};
   end 
   for i = 1:length(cmd)
      if i ==1
         eval(['out = ' cmd{i} ';']);
      else
         eval(cmd{i});
      end
   end
   store_model(out,str);
end

%%%%%
% Lists predefined models (append when adding)
function out = list_predef_model_strings
out = {
    'adult_male_16el';
    'adult_male_32el';
    'adult_male_16el_lungs';
    'adult_male_32el_lungs';
    'adult_male_grychtol2016a_1x32';
    'adult_male_grychtol2016a_2x16';
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
    'lamb_newborn_16el';
    'lamb_newborn_32el';
    'lamb_newborn_16el_organs';
%     'lamb_newborn_32el_organs';
    'beagle_16el';
    'beagle_32el';
    'beagle_16el_lungs';
    'beagle_32el_lungs';
    'beagle_16el_rectelec';
    'beagle_32el_rectelec';
    'beagle_16el_lungs_rectelec';
    'beagle_32el_lungs_rectelec';
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
    case 'adult_male_grychtol2016a_1x32'
        out = mk_thorax_model_grychtol2016a('1x32_ring');
        out = out.fwd_model;
    case 'adult_male_grychtol2016a_2x16'
        out = mk_thorax_model_grychtol2016a('2x16_planar');
        out = out.fwd_model;
        
    case 'cylinder_16x1el_coarse'
       out = build_if_needed(...
          'ng_mk_cyl_models([10,15],[16,5],[0.5,0,0.18])', str);
    case 'cylinder_16x1el_fine' 
       out = build_if_needed(...
          'ng_mk_cyl_models([10,15,1.1],[16,5],[0.5,0,0.15])',str);
    case 'cylinder_16x1el_vfine' 
        out = build_if_needed(...
           'ng_mk_cyl_models([10,15,0.8],[16,5],[0.5,0,0.08])',str);
    case 'cylinder_16x2el_coarse' 
        out = build_if_needed(...
           'ng_mk_cyl_models([30,15],[16,10,20],[0.5,0,0.18])',str);
    case 'cylinder_16x2el_fine'  
        out = build_if_needed(...
           'ng_mk_cyl_models([30,15,1.5],[16,10,20],[0.5,0,0.15])',str);
    case 'cylinder_16x2el_vfine' 
        out = build_if_needed(...
           'ng_mk_cyl_models([30,15,0.8],[16,10,20],[0.5,0,0.08])',str);
        
        
    case 'neonate_16el'
        out = mk_library_model({'neonate','boundary'},[16 1 0.5],[0.1 0 -1 0 60],0.08,49);
    case 'neonate_32el'
        out = mk_library_model({'neonate','boundary'},[32 1 0.5],[0.06 0 -1 0 60],0.08,49);
    case 'neonate_16el_lungs'
        out = mk_library_model({'neonate','boundary','left_lung','right_lung'},[16 1 0.5],[0.1 0 -1 0 60],0.08,49);
    case 'neonate_32el_lungs'
        out = mk_library_model({'neonate','boundary','left_lung','right_lung'},[32 1 0.5],[0.06 0 -1 0 60],0.08,49);
        
    case 'pig_23kg_16el'
        out = mk_library_model({'pig_23kg','boundary'},...
            [16 1 0.5],[0.05 0 -1 0 60],0.08);
    case 'pig_23kg_32el'
        out = mk_library_model({'pig_23kg','boundary'},...
            [32 1 0.5],[0.05 0 -1 0 60],0.08);
    case 'pig_23kg_16el_lungs'
        out = mk_library_model({'pig_23kg','boundary','lungs(1:2:end,:)'},...
            [16 1 0.5],[0.05 0 -1 0 60],0.08);
    case 'pig_23kg_32el_lungs'
        out = mk_library_model({'pig_23kg','boundary','lungs(1:2:end,:)'},...
            [32 1 0.5],[0.05 0 -1 0 60],0.08);    

    case 'lamb_newborn_16el'
%        out = build_if_needed(...
%           {['ng_mk_extruded_model({208,208*',...
%             'shape_library(''get'',''lamb_newborn'',''boundary'')',...
%             ',0,10},[16,1.995,104],[1])'],'out.nodes = out.nodes/204;'}, ...
%           str);
         out = mk_library_model({'lamb_newborn','boundary'},[16,1.995,104],[1],10,0,208);
    case 'lamb_newborn_32el'
         % Very sensitive to the .980 offset. This is the only I can find that works.
         out = mk_library_model({'lamb_newborn','boundary'},[32,1.980,104],[1],15,50,208);
         out.electrode = out.electrode([2:32,1]);
    case 'lamb_newborn_16el_organs'
       out = mk_library_model({'lamb_newborn','boundary','lungs','heart'},[16,1.995,104],[1],10,0,208);
%     case 'lamb_newborn_32el_organs'
    case 'beagle_16el';
      scale = 49;
      out = mk_library_model({'beagle','boundary'}, ...
         [16 1 scale*0.5],[2,0,0.10],10,0,49);
    case 'beagle_32el';
      scale = 49;
      out = mk_library_model({'beagle','boundary'}, ...
         [32 1 scale*0.5],[2,0,0.10],10,0,49);
    case 'beagle_16el_rectelec';
      scale = 49;
      out = mk_library_model({'beagle','boundary'}, ...
         [16 1 scale*0.5],8*[0.25,1,0.05],10,0,49);
    case 'beagle_32el_rectelec';
      scale = 49;
      out = mk_library_model({'beagle','boundary'}, ...
         [16 1 scale*0.5],8*[0.25,1,0.05],10,0,49);

    case 'beagle_16el_lungs';
      scale = 49;
      out = mk_library_model({'beagle','boundary','left_lung','right_lung'}, ...
         [16 1 scale*0.5],[2,0,0.10],10,0,49);

    case 'beagle_16el_lungs_rectelec';
      scale = 49;
      out = mk_library_model({'beagle','boundary','left_lung','right_lung'}, ...
         [16 1 scale*0.5],8*[0.25,1,0.05],10,0,49);

    case 'beagle_32el_lungs';
      scale = 49;
      out = mk_library_model({'beagle','boundary','left_lung','right_lung'}, ...
         [32 1 scale*0.5],[2,0,0.10],10,0,49);

    case 'beagle_32el_lungs_rectelec';
      scale = 49;
      out = mk_library_model({'beagle','boundary','left_lung','right_lung'}, ...
         [32 1 scale*0.5],8*[0.25,1,0.05],10,0,49);
         
    otherwise
        error('No such model');
end
%give the model a name
out.name = str;


function str = make_filename(shape, elec_pos, elec_shape, ...
                             maxsz, nfft, scale);
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
if ~isempty(nfft)
    str = [str '_nfft_' num2str(nfft)];
end
if ~isempty(scale)
    str = [str '_scale' num2str(scale)];
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
n_models = numel(models);
sqrt_n_models = ceil(sqrt(n_models));
for i = 1:numel(models)
    eidors_msg('\n\n\n DOING  MODEL (%s)\n\n\n',models{i},0);
    mdl = mk_library_model(models{i});
    img = mk_image(mdl,1);
    try   
        n = numel(mdl.mat_idx); 
    catch
        n =1; 
    end
    if n >1
        for j = 2:n
            img.elem_data(mdl.mat_idx{j}) = 0.25;
        end
    end
    subplot(sqrt_n_models, sqrt_n_models,i);
    show_fem(img,[0,1,0]); axis off;
    title(models{i},'Interpreter','none');
    drawnow
end


out = mk_library_model({'pig_23kg','boundary','lungs(1:2:end,:)'},[32 1 0.5],[0.05 0 -1 0 60],0.08);

