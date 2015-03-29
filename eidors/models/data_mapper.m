function img = data_mapper(img, reverse)
%DATA_MAPPER maps img.params data to elem or node data
% img = data_mapper(img) will only work if there is only a single
% params data field on the img, throwing an error otherwise.
% For multi-parametrization images, img.data_mapper must be specified. It
% can be either a name of the appropriate parametrization data field on the img, or
% a function that will handle the data differently.
%
% To avoid confusion, the parametrization currently in use by image data will be
% specified in a descriptive string tag img.current_params
%
% img = data_mapper(img, TRUE) will reverse the process updating
% the current parametrization with the data from img.elem_data or img.node_data.
%
% Examples:
%  imdl = mk_common_model('a2c2',8);
%  c = 3*ones(length(imdl.fwd_model.elems),1);
%  img = eidors_obj('image','fwd_model',imdl.fwd_model);
%  img.conductivity.elem_data = c;
%  img = data_mapper(img);
%
%  img.resistivity.elem_data = 1./c;
%  img.data_mapper = 'resistivity';
%  img = data_mapper(img);
%
%  img.data_mapper = 'some_function_that_takes_an_image';
%  img = data_mapper(img);
%
% KNOWN ISSUES:
%   - doesn't fully support image arrays
%
% See also: MK_IMAGE, SUPPORTED_PARAMS

% (C) 2012 Bartlomiej Grychtol.
% Licenced under GPL version 2 or 3
% $Id$

if isstr(img) && strcmp(img,'UNIT_TEST'); img = do_unit_test; return; end

if nargin < 2
   reverse = false;
elseif isstr(reverse) && strcmpi(reverse, 'pack')
   reverse = false;
else % 'unpack'
   reverse = true;
end
% need to handle image arrays
for i = 1:length(img)
   if reverse
      tmp(i) = unpack_params(img(i));
   else
      tmp(i) = pack_params(img(i));
   end
end
img = tmp;
end

function img = pack_params(img);

flds = fieldnames(img);
prms = ismember(flds, supported_params);

switch sum(prms)
   case 0
      % no params found, check for elem_data or node_data
      if ~( isfield(img,'elem_data') || isfield(img,'node_data'))
         error('No params, elem_data or node_data found on img');
      else
         if isfield(img,'current_params') && ~isempty(img.current_params)
                eidors_msg('@@@ Careful! Image already mapped. Doing nothing',4);
           return;
         end

         eidors_msg('@@@ No params present on img. Assuming conductivity',4);
% STUPID MATLAB CAN'T KEEP SYNTAX STRAIGHT BETWEEN VERSIONS
%        img = setfield(img,{},'current_params', 'conductivity');
% NEED TO DO THIS, UGLY INEFFICIENT CODE
         for i=1:length(img)
            img(i).current_params = 'conductivity';
         end
      end
      return      % returns img
   case 1
      ph = find(prms);
      if isfield(img.(flds{ph}),'func')
         % let the function handle the whole process
         img = feval(img.(flds{ph}).func, img);
         return
      else
         img = pack_params_using_params_str(img, flds{ph});
      end
   otherwise
      % multiple params
      try
         prms_fun = img.data_mapper;
      catch
         error('Multiple params found on img require: img.data_mapper = function or {"string","array"}');
      end
      if ~isa(prms_fun, 'function_handle') && ismember(prms_fun,flds);
         img = pack_params_using_params_str(img, prms_fun); % assume a string array
      else
         try
            img = feval(prms_fun,img);
         catch
            error('img.data_mapper not understood');
         end
      end

end
end

function img = pack_params_using_params_str(img, prms)

prms_flds = fieldnames(img.(prms));
for i = 1:length(prms_flds)
   if isfield(img, prms_flds{i})
      eidors_msg('@@@ Overwriting img.%s',prms_flds{i},4);
      warning('EIDORS:OverwritingData', 'Overwriting img.%s',prms_flds{i});
   end
   % allow elem_data/node_data to be scalar
   if strcmp(prms_flds{i},'elem_data') && numel(img.(prms).elem_data) == 1
      n_elem = calc_num_elems(img.fwd_model);
      img.elem_data = ones(n_elem,1);
      img.elem_data(:) = img.(prms).elem_data;
   elseif strcmp(prms_flds{i},'node_data')
      img.node_data = ones(size(img.fwd_model.nodes,1),1);
      img.node_data(:) = img.(prms).node_data;
   else
      img.(prms_flds{i}) = img.(prms).(prms_flds{i});
   end
end
img.current_params = prms;
end

function n = calc_num_elems(fmdl)
if isfield(fmdl, 'coarse2fine')
   n = size(fmdl.coarse2fine,2);
else
   n = length(fmdl.elems);
end
end

function img = copy_data_to_params(img, prms)
prms_flds = fieldnames(img.(prms));
try
   for i = 1:length(prms_flds)
      img.(prms).(prms_flds{i}) = img.(prms_flds{i});
      img = rmfield(img, prms_flds{i});
   end
   img = rmfield(img, 'current_params');
catch
   error('Fields specified by img.%s missing from img.',prms);
end
end

function img = unpack_params(img)
try
   curprms = img.current_params;
catch
   if isfield(img,'elem_data') || isfield(img, 'node_data')
      return % old type image, nothing to do
   else
      error('img.current_params required');
   end
end
if ismember(curprms, fieldnames(img))
   img = copy_data_to_params(img,curprms);
elseif strcmp(curprms,'conductivity')
   % current_params is conductivity, but there's no img.conductivity
   if ~any(ismember(fieldnames(img),supported_params))
      % we're dealing with an old params-oblivious image
      % nothing to do
      img = rmfield(img, 'current_params');
   else
      error('Cannot reverse %s mapping',curprms);
   end
else
   try
      % user-provided data_mapper must provide the reverse
      img = feval(curprms,img,1);
   catch
      error('data_mapper %s failed to reverse', curprms);
   end
end
end

function pass = do_unit_test
   pass = 1;
   pass = do_unit_test_legacy() && pass;
   pass = do_unit_test_undo()   && pass;
   disp('');
   if pass
      disp('TEST: overall PASS');
   else
      error('TEST: overall FAIL');
   end
end

function pass = do_unit_test_legacy
   pass = 1;
   ll = eidors_msg('log_level');
   eidors_msg('log_level',3);

   imdl = mk_common_model('a2c2',8);
   c = 3*ones(length(imdl.fwd_model.elems),1);
   img = eidors_obj('image','test_image','fwd_model',imdl.fwd_model);
   img.conductivity.elem_data = c;

   fprintf('TEST: mv img.conductivity.elem_data to img.elem_data\n');
   imgi = data_mapper(img, 'pack');
   if ~isfield(imgi, 'elem_data')
     imgi
     disp('TEST FAIL: missing elem_data');
     pass = 0;
   end

   imgi2 = data_mapper(img, 'pack'); % give a msg at log_level 3

   imgii = data_mapper(imgi,'unpack'); % unpack
   if ~isequaln(imgii, img)
      img
      imgii
      fprintf('TEST FAIL: pack to unpack returned different result\n');
      pass = 0;
   end

   img.resistivity.elem_data = 1./c;
   img.data_mapper = 'resistivity';
   img = data_mapper(img);
   if abs(img.elem_data(1) - 1/3) > 1e-4
      img.elem_data(1)
      pass = 0;
      disp('TEST FAIL: expected img.elem_data to be ~1/3');
   end

   img.data_mapper = @unit_test_passthrough; % some function that takes an image
   imgf = data_mapper(img);
   % if this works without failing, then we are happy
   if ~isequaln(imgf, img)
      img
      imgf
      fprintf('TEST FAIL: pass-through function failed\n');
      pass = 0;
   end

   eidors_msg('log_level',ll);
   if ~pass
      disp('TEST: ---------------------------------------');
   end
end

function img = unit_test_passthrough(img1,rev)
  img = img1;
end

function pass = do_unit_test_undo
   pass = 1;

   imdl = mk_common_model('c2c', 16);
   img = mk_image(imdl);
   d = img.elem_data;
   img0 = rmfield(img, 'elem_data');

   disp('TEST: img.elem_data');
   img = img0; img.elem_data = d;
   pass = test_undo(img, d, 'elem_data') && pass;

   disp('TEST: img.node_data');
   img = img0; img.node_data = d;
   pass = test_undo(img, d, 'node_data') && pass;

   disp('TEST: img.conductivity.elem_data');
   img = img0; img.conductivity.elem_data = d;
   pass = test_undo(img, d, 'elem_data') && pass;
% TODO AB below here blows up
%   disp('TEST: img.conductivity.node_data');
%   img = img0; img.conductivity.node_data = d;
%   pass = test_undo(img, d) && pass;

%   disp('TEST: img.resistivity.node_data');
%   img = img0; img.resistivity.node_data = d;
%   pass = test_undo(img, d) && pass;

%   disp('TEST: img.conductivity');
%   img = img0; img.conductivity = d; img.params_mapper = 'conductivity';
%   pass = test_undo(img, d) && pass;
end

% img, d=data, dl=data location (once packed)
function pass = test_undo(img, d, dl)
   pass = 1;
   imgp = data_mapper(img);
   imgup = data_mapper(imgp,'unpack');
   if ~isequaln(img, imgup)
      img
      imgp
      imgup
      disp('TEST FAIL: the initial image and the image after being packed, then unpacked differ');
      pass = 0;
   end
   if ~isequaln(imgp.(dl), d)
      disp(['TEST FAIL: the packed data imgp.' dl ' seems broken']);
      pass = 0;
   end
   if ~pass
      disp('TEST: ---------------------------------------');
   end
end
