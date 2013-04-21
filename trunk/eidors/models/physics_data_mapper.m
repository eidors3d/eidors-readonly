function img = physics_data_mapper(img, reverse)
%PHYSICS_DATA_MAPPER maps img.physics data to elem or node data
% img = physics_data_mapper(img) will only work if there is only a single
% physics data field on the img, throwing an error otherwise. 
% For multi-physics images, img.physics_data_mapper must be specified. It
% can be either a name of the appropriate physics data field on the img, or
% a function that will handle the data differently.
% 
% To avoid confusion, the physics currently in use by image data will be
% specified in a descriptive string tag img.current_physics
%
% img = physics_data_mapper(img, TRUE) will reverse the process updating
% the current physics with the data from img.elem_data or img.node_data.
%
% Examples:
%  imdl = mk_common_model('a2c2',8);
%  c = 3*ones(length(imdl.fwd_model.elems),1);
%  img = eidors_obj('image','fwd_model',imdl.fwd_model);
%  img.conductivity.elem_data = c;
%  img = physics_data_mapper(img);
%
%  img.resistivity.elem_data = 1./c;
%  img.physics_data_mapper = 'resistivity';
%  img = physics_data_mapper(img);
%
%  img.physics_data_mapper = 'some_function_that_takes_an_image';
%  img = physics_data_mapper(img);
%
% KNOWN ISSUES:
%   - doesn't fully support image arrays
%
% See also: MK_IMAGE, SUPPORTED_PHYSICS

% (C) 2012 Bartlomiej Grychtol. 
% Licenced under GPL version 2 or 3
% $Id$

if isstr(img) && strcmp(img,'UNIT_TEST'); do_unit_test; return; end

if nargin < 2
   reverse = false;
end

if reverse
   img = map_data_to_physics(img);
else
   img = map_physics_to_data(img);
end

end

function img = map_physics_to_data(img);

knownphysics = supported_physics;
flds = fieldnames(img);
phys = ismember(flds,knownphysics);

switch sum(phys)
   case 0
      % no physics found, check for elem_data or node_data
      if ~( isfield(img,'elem_data') || isfield(img,'node_data'))
         error('No physics, elem_data or node_data found on img');
      else
         eidors_msg('@@@ No physics present on img. Assuming conductivity',3);
% STUPID MATLAB CAN'T KEEP SYNTAX STRAIGHT BETWEEN VERSIONS
%        img = setfield(img,{},'current_physics', 'conductivity');
% NEED TO DO THIS, UGLY INEFFICIENT CODE
         for i=1:length(img)
            img(i).current_physics = 'conductivity';
         end
      end
      return      % returns img
   case 1
      ph = find(phys);
      if isfield(img.(flds{ph}),'func')
         % let the function handle the whole process
         img = feval(img.(flds{ph}).func, img);
         return
      else
         img = copy_physics_to_data(img, flds{ph});
      end
   otherwise
      % multiple physics
      try
         phys_fun = img.physics_data_mapper;
      catch
         error('Multiple physics found on img: img.physics_data_mapper required');
      end
      if ~isa(phys_fun, 'function_handle') && ismember(phys_fun,flds);
         img = copy_physics_to_data(img, phys_fun);
      else
         try
            img = feval(phys_fun,img);
         catch
            error('img.physics_data_mapper not understood');
         end
      end

end
end

function img = copy_physics_to_data(img, phys)

phys_flds = fieldnames(img.(phys));
for i = 1:length(phys_flds)
   if isfield(img, phys_flds{i})
      eidors_msg('@@@ Overwriting img.%s',phys_flds{i},3);
      warning('EIDORS:OverwritingData', 'Overwriting img.%s',phys_flds{i});
   end
   % allow elem_data/node_data to be scalar
   if strcmp(phys_flds{i},'elem_data') && numel(img.(phys).elem_data) == 1
      n_elem = calc_num_elems(img.fwd_model);
      img.elem_data = ones(n_elem,1);
      img.elem_data(:) = img.(phys).elem_data;
   elseif strcmp(phys_flds{i},'node_data') 
      img.node_data = ones(size(img.fwd_model.nodes,1),1);
      img.node_data(:) = img.(phys).node_data;
   else
      img.(phys_flds{i}) = img.(phys).(phys_flds{i});
   end
end
img.current_physics = phys;
end

function n = calc_num_elems(fmdl)
if isfield(fmdl, 'coarse2fine')
   n = size(fmdl.coarse2fine,2);
else
   n = length(fmdl.elems);
end
end

function img = copy_data_to_physics(img, phys)
phys_flds = fieldnames(img.(phys));
try
   for i = 1:length(phys_flds)
      img.(phys).(phys_flds{i}) = img.(phys_flds{i});
      img = rmfield(img, phys_flds{i});
   end
   img.current_physics = [];
catch
   error('Fields specified by img.%s missing from img.',phys);
end
end

function img = map_data_to_physics(img)
try 
   curphys = img.current_physics;
catch
   if isfield(img,'elem_data') || isfield(img, 'node_data')
      return % old type image, nothing to do
   else
      error('img.current_physics required');
   end
end
if ismember(curphys, fieldnames(img))
   img = copy_data_to_physics(img,curphys);
elseif strcmp(curphys,'conductivity')
   % current_physics is conductivity, but there's no img.conductivity
   if ~any(ismember(fieldnames(img),supported_physics))
      % we're dealing with an old physics-oblivious image
      % nothing to do
      img.current_physics = [];
   else
      error('Cannot reverse %s mapping',curphys);
   end
else
   try
      % user-provided physics_data_mapper must provide the reverse
      img = feval(curphys,img,1);
   catch 
      error('physics_data_mapper %s failed to reverse', curphys);
   end
end
end

function do_unit_test
   ll = eidors_msg('log_level');
   eidors_msg('log_level',3);
   
   imdl = mk_common_model('a2c2',8);
   c = 3*ones(length(imdl.fwd_model.elems),1);
   img = eidors_obj('image','test_image','fwd_model',imdl.fwd_model)
   img.conductivity.elem_data = c;
   img = physics_data_mapper(img)
   
   img = physics_data_mapper(img) % give a msg at log_level 3
   
   img = physics_data_mapper(img,1) % reverse
   
   img.resistivity.elem_data = 1./c;
   img.physics_data_mapper = 'resistivity';
   img = physics_data_mapper(img)
   img.elem_data(1) %0.3333
   
   img.physics_data_mapper = @class; % some function that takes an image
   physics_data_mapper(img)
   
   eidors_msg('log_level',ll);
end
