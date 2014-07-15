function img = params_mapper(img, rev)
% img = PARAMS_MAPPER(img, [rev])
%
% img = params_mapper(img); % will return an image with packed img.inv_params
% img = params_mapper(img, 1); % will unpack the previously packed img
%
% A packed img should have the following fields:
%  img.inv_params = []; % a numeric vector
%     .current_inv_params = {'xx', @f}; % the sources of the inv_params
%     .inv_params_sel = []; % [start stop] indices,
%                           % one row per img.current_inv_params{i}
% The packed img should not have any of the source data.
%
% An unpacked img should have the following fields:
%  img.xx = []; % a numeric vector
%     .asdf = [1xN struct]; % data used to calculate inv_params (inverse problem's parameters) for @f
%     .params_mapper = {'asdf', @f}; % the sources of the inv_params
% An unpacked img should have none of the inv_params fields.
%
% WARNING: For supporting legacy operations, this function calls data_mapper(),
% to maintain backwards compatilibity. Legacy operation is detected when there
% is no img.params_mapper field found.
%
% LEGACY: img.conductivity.elem_data or
%         img.conductivity.node_data     <-> img.inv_params
% This function handles
%         img.elem_data or
%         img.node_data     <-> img.inv_params
% The return route is saved as img.current_inv_params as either 'node_data' or
% 'elem_data'.
% The data_mapper() function, called within, deals with moving
%         img.conductivity.elem_data or
%         img.conductivity.node_data     <-> img.elem_data or img.node_data
% The return route is saved in img.current_params.
%
% LIMITATIONS:
%  - img arrays are not handled, as in img(1:3).inv_params
%  - can not handle more than one @f img.params_mapper function_handle
%  - if using a function handle, the data is left in img.inv_params though it
%    should really have a home...
%
% SEE ALSO: data_mapper, convert_img_units

% (C) 2014 Alistair Boyle, Bartlomiej Grychtol. License: GPL version 2 or version 3

if isstr(img) && strcmp(img,'UNIT_TEST'); img = do_unit_test; return; end

if nargin < 2
   img = pack_params(img);
else
   img = unpack_params(img);
end

% PACK
function img = pack_params(img)
   % if no parameter mapper has been provided, handle the munging of data in a
   % sane fashion using the legacy 'data_mapper' function
   if ~isfield(img, 'params_mapper') % LEGACY support
      % move img.conductivity.elem_data -> img.elem_data OR
      % move img.conductivity.node_data -> img.node_data
% TODO why not replace img? what does data_mapper do that is bad?
      jnk = data_mapper(img);
      % move img.elem_data OR img.node_data -> img.inv_params
      % save where we can from in img.current_inv_params
      try
         img.inv_params = jnk.node_data;
         img.current_inv_params = 'node_data';
         if isfield(img, 'node_data') img = rmfield(img, 'node_data'); end
      catch
         img.inv_params = jnk.elem_data;
         img.current_inv_params = 'elem_data';
         if isfield(img, 'elem_data') img = rmfield(img, 'elem_data'); end
      end
      img.inv_params_sel = [1 length(img.inv_params)]; % [start end] == full vector
      img.current_params = jnk.current_params; % needed to reverse
   else % MODERN handling
      % Otherwise, someone has provided an explicit list on how to do the
      % mapping, follow those instructions:
      pm = img.params_mapper;
      if ~iscell(pm) % if its a lonely function handle, wrapper it in a cell array
        pm = {pm};
      end
      inv_params = []; % where we are appending all our data together
      len = []; % lengths of collected inv_params prior to appending, one row per pm{i}
      for i = 1:length(pm)
         % there is some clever-er function to do it: fine, be that way
         if isa(pm{i},'function_handle')
            % get data
            pv = feval(pm{i}, img);
            % we don't know enough about where it came from to delete the data: too bad
         elseif isstr(pm{i})
            % get data
            pn = ['img.' pm{i}]; % parameterization name
            pv = eval(pn); % parameterization values
            % delete empty structures left behind from moving the data
            pnl = strsplit(pn,'.'); % split source name on '.'
            while( (length(pnl) >= 2) && (length(eval(strjoin(pnl,'.'))) > 0) )
              rmfield(eval(strjoin(pnl(1:end-1))), pnl(end)); pnl = pnl(1:end-1);
            end
            % is now clean
         else
            error('do not know how shuffle around your parameterization into img.inv_params: see "help params_mapper" and img.params_mapper to sort out your structure');
         end
         % append new data
         inv_params = [inv_params; pv];
         len = [len; length(pv)];
      end
      % place updated data in the img structure
      img.inv_params = inv_params;
      img.current_inv_params = pm;
      % len=[10 4 2] --> sel=[1 10; 11 14; 15 16]
      start = cumsum(len)-len+1; stop = cumsum(len);
      img.inv_params_sel = [ start' stop'];
   end

% UNPACK
function img = unpack_params(img)
   if ~isfield(img, 'params_mapper') % LEGACY support
      if ~isfield(img, 'current_inv_params')
         error('Need img.current_inv_params');
      end
      switch img.current_inv_params
         case 'node_data'
            img.node_data = img.inv_params;
         case 'elem_data'
            img.elem_data = img.inv_params;
         otherwise
            error('huh?');
      end
      img = data_mapper(img,1); % unpack
   else % MODERN handling
      % Otherwise, someone has provided an explicit list on how to do the
      % mapping, follow those instructions:
      pm = img.params_mapper;
      if isa(pm,'function_handle') % if its a lonely function handle, wrapper it in a cell array
        pm = {pm};
      end
      if iscell(pm) % its a cell array: go to work
         function_seen = 0;
         for i = 1:length(pm)
            ps = img.inv_params_sel(i,:); % the [start stop] indices within img.inv_params
            if isa(pm,'function_handle')
               % TODO FIXME? AB: this looks totally broken -- how can the function put data back where it came from!
               jnk = img;
               jnk.inv_params = jnk.inv_params(ps(1):ps(2),:); % hide the rest of the inv_params from the function_handle
               pvf = feval(pm{i},jnk,1); % 1=REVERSE; we save function_handle output to pvf for storage to img.inv_params again, once we are done
               % TODO HACK to get out before we break something
               if(function_seen)
                 error('This function can not handle cell arrays of img.params_mapper that contain more than a single function_handle since the result leaves the mapped data in img.inv_params where it would be overwritten by a second function_handler if one were allowed. Multiple function_handles in the img.params_mapper cell array is not allowed at this time: try combining your function_handles into a single function.');
               end
               function_seen = 1;
            else
               % put the data back
               pn = ['img.' pm{i}]; % parameterization name
               eval([pn ' = img.inv_params(ps(1):ps(2),:);']); % parameterization values to their original homes
            end
         end
         % TODO HACK to protect function_handle unpacking operation into the old data
         if function_seen
            img.inv_params = pvf;
         end
      else % sorry...
         error('do not know how shuffle around your parameterization into img.inv_params: see "help params_mapper" and img.params_mapper to sort out your structure');
      end
  end
  % remove parameters after unpacking
  img = rmfield(img,'inv_params');
  img = rmfield(img,'inv_params_sel');
  img = rmfield(img,'current_inv_params');

function pass = do_unit_test()
   pass = 1;
   pass = pass && do_unit_test_undo();
   disp('');
   if pass
      disp('TEST: overall PASS');
   else
      disp('TEST: overall FAIL');
   end

function pass = do_unit_test_undo
   pass = 1;

   imdl = mk_common_model('c2c', 16);
   img = mk_image(imdl);
   d = img.elem_data;
   img0 = rmfield(img, 'elem_data');

   disp('TEST: img.elem_data');
   img = img0; img.elem_data = d;
   pass = test_undo(img, d) && pass;

   disp('TEST: img.node_data');
   img = img0; img.node_data = d;
   pass = test_undo(img, d) && pass;

   disp('TEST: img.conductivity.elem_data');
   img = img0; img.conductivity.elem_data = d;
   pass = test_undo(img, d) && pass;
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

function pass = test_undo(img, d)
   pass = 1;
   imgp = params_mapper(img);
   imgup = params_mapper(imgp,'unpack');
   if ~isequaln(img, imgup)
      img
      imgp
      imgup
      fprintf('TEST FAIL: the initial image and the image after being packed, then unpacked differ\n');
      pass = 0;
   end
   if ~isequaln(imgp.inv_params, d)
      fprintf('TEST FAIL: the packed data imgp.inv_params seems broken\n');
      pass = 0;
   end
%   if ~isequaln(imgp.inv_params_sel, d)
%      fprintf('TEST FAIL: the packed data seems broken\n');
%      pass = 0;
%   end
   if ~pass
      disp('TEST: ---------------------------------------');
   end
