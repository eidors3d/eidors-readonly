function img = map_img(img, arg1, arg2);
% function img  = map_img(img,  out)
% function data = map_img(data, in, out)

   %--------------------------
   % UNIT_TEST?
   if isstr(img) && strcmp(img,'UNIT_TEST'); img = do_unit_test; return; end

   if nargin == 3
      in = arg1;
      out = arg2;
      img = map_data(img, in, out); % convert data directly
      return
   else
      out = arg1;
   end

   err_if_inf_or_nan(img.elem_data, 'img-pre');
   try in = img.current_physics;
   catch in = {'conductivity'};
   end
   % make cell array of strings
   if isstr(in)
      in = {in};
      img.current_physics = in;
   end
   if isstr(out)
      out = {out};
   end

   % if we have mixed data, check that we have a selector to differentiate between them
   if ~isfield(img, 'physics_sel')
      if length(in(:)) == 1
         img.physics_sel = {1:size(img.elem_data,1)};
      else
         error('found multiple physics but no physics_sel cell array in img');
      end
   end

   % create data?! we don't know how
   if length(out(:)) > length(in(:))
      % if we are missing movement, construct it (default = 0 movement)
      if strcmp(out{end}, 'movement')
         if ~isfield(img, 'elem_movement_init')
            error('to create initial movement data we need the inv_model.parameters.elem_movement_init vector with initial movement values');
         end
         e = img.elem_movement_init;
         ei = length(img.elem_data);
         ei = [ei+1:ei+length(e)];
         img.elem_data(ei) = e;
         img.physics_sel{end+1} = ei;
         in{end+1} = 'movement';
         img.current_physics = in;
      end
      % are we still broken -- then error out
      if length(out(:)) > length(in(:))
         error('missing data (more out types than in types)');
      end
   elseif length(out(:)) < length(in(:))
      % delete data: we can do that
      % NOTE that if we are doing this, we always assume its the *last* items in the list
      for i = 1:length(in(:))-length(out(:)) % delete the extra
         img.elem_data(img.physics_sel{end}) = []; % rm elem_data
         img.physics_sel(end) = []; % rm physics_sel
         img.current_physics{end} = []; % rm current_physics
      end
   end

   % the sizes now match, we can do the mapping
   for i = 1:length(out(:))
      % map the data
      x = img.elem_data(img.physics_sel{i});
      x = map_data(x, in{i}, out{i});
      img.elem_data(img.physics_sel{i}) = x;
      img.current_physics{i} = out{i};
   end
   err_if_inf_or_nan(img.elem_data, 'img-post');

   % clean up physics_sel/current_physics if we only have one physics
   if length(img.current_physics(:)) == 1
      img.current_physics = img.current_physics{1};
      img = rmfield(img, 'physics_sel'); % unnecessary since we know its all elem_data
   end

function x = map_data(x, in, out)
   % check that in and out are single strings, not lists of strings
   if ~isstr(in)
      if iscell(in) && (length(in(:)) == 1)
         in = in{1};
      else
         error('expecting single string for map_data() "in" type');
      end
   end
   if ~isstr(out)
      if iscell(out) && (length(out(:)) == 1)
         out = out{1};
      else
         error('expecting single string for map_data() "out" type');
      end
   end

   % quit early if there is nothing to do
   if strcmp(in, out) % in == out
      return; % do nothing
   end

   % resistivity to conductivity conversion
   % we can't get here if in == out
   % we've already checked for log convserions on input or output
   if any(strcmp(in,  {'resistivity', 'conductivity'})) && ...
      any(strcmp(out, {'resistivity', 'conductivity'}))
      x = 1./x; % conductivity <-> resistivity
   % log conversion
   elseif any(strcmp({in(1:3), out(1:3)}, 'log'))
      % log_10 x -> x
      if strcmp(in(1:6), 'log10_')
         if any(x >= log10(realmax)-eps) warning('loss of precision -> inf'); end
         x = map_data(10.^x, in(7:end), out);
      % ln x -> x
      elseif strcmp(in(1:4), 'log_')
         if any(x >= log(realmax)-eps) warning('loss of precision -> inf'); end
         x = map_data(exp(x), in(5:end), out);
      % x -> log_10 x
      elseif strcmp(out(1:6), 'log10_')
         if any(x <= 0 + eps) warning('loss of precision -> -inf'); end
         x = log10(map_data(x, in, out(7:end)));
      % x -> ln x
      elseif strcmp(out(1:4), 'log_')
         if any(x <= 0 + eps) warning('loss of precision -> -inf'); end
         x = log(map_data(x, in, out(5:end)));
      else
         error(sprintf('unknown conversion (log conversion?) %s - > %s', in, out));
      end
   else
      error('unknown conversion %s -> %s', in, out);
   end
   x(x == +inf) = +realmax;
   x(x == -inf) = -realmax;
   err_if_inf_or_nan(x, 'map_data-post');

function err_if_inf_or_nan(x, str);
  if any(isnan(x) | isinf(x))
      error(sprintf('bad %s (%d NaN, %d Inf of %d)', ...
                    str, ...
                    length(find(isnan(x))), ...
                    length(find(isinf(x))), ...
                    length(x)));
  end

% test functions
function pass = do_unit_test
pass = 1;
d = 1;
while d ~= 1 & d ~= 0
  d = rand(1);
end
disp('TEST: map_data()');
elem_types = {'conductivity', 'log_conductivity', 'log10_conductivity', ...
              'resistivity',  'log_resistivity',  'log10_resistivity'};
expected = [d         log(d)         log10(d)      1./d      log(1./d)      log10(1./d); ...
            exp(d)    d              log10(exp(d)) 1./exp(d) log(1./exp(d)) log10(1./exp(d)); ...
            10.^d     log(10.^d )    d             1./10.^d  log(1./10.^d ) log10(1./10.^d ); ...
            1./d      log(1./d  )    log10(1./d)   d         log(d)         log10(d); ...
            1./exp(d) log(1./exp(d)) log10(1./exp(d)) exp(d) d              log10(exp(d)); ...
            1./10.^d  log(1./10.^d)  log10(1./10.^d)  10.^d  log(10.^d)     d ];
for i = 1:length(elem_types)
  for j = 1:length(elem_types)
    pass = test_map_data(d, elem_types{i}, elem_types{j}, expected(i,j), pass);
  end
end

if pass
   disp('TEST: overall PASS');
else
   disp('TEST: overall FAIL');
end

function pass = test_map_data(data, in, out, expected, pass)
fprintf('TEST: map_img(%s -> %s)\n', in, out);
if map_img(data, in, out) ~= expected
   pass = 0;
   fprintf('TEST: FAIL map_img(%s -> %s)\n', in, out);
end
