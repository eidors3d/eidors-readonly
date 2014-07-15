function y = convert_img_units(img,arg1,arg2)
%CONVERT_IMG_UNITS change image data units 
%  img = convert_img_units(img,new_unit) converts img.elem_data or img.node_data
%  expressed in img.current_params into different units for the same 
%  physical property
% 
% Examples: 
%  img = mk_image(mdl, 2, 'resisitivity');
%  img = convert_img_units(img, 'conductivity');
%
%  img = mk_image(mdl, 1);
%  img.current_params = 'log_resistivity';
%  img = convert_img_units(img, 'conductivity');

% (C) 2012-2014 Alistair Boyle and Bartlomiej Grychtol
% License: GPL version 2 or 3
% $Id$

if isstr(img) && strcmp(img,'UNIT_TEST'); y = do_unit_test; return; end

if nargin < 2
  error('EIDORS:ConvertUnits:Input', ...
      'Not enough input arguments.');
end

% convert_img_units(x, in, out);
if isnumeric(img)
  % require 3 arguments
  if nargin < 3
    error('EIDORS:ConvertUnits:Input', ...
      'Two units are required to convert numeric values.');
  else
    x = img;
    cur_unit = arg1;
    new_unit = arg2;
    y = map_data(x, cur_unit, new_unit);
    return
  end
end

% convert_img_units(img, out);
if any(isfield(img,{'elem_data','node_data'}))
  if nargin == 3
    cur_unit = arg1;
    new_unit = arg2;
  else
    try
      cur_unit = img.current_params;
      new_unit = arg1;
    catch
      error(['Don''t know what to convert. '...
              'Need either 3 arguments or img.current_params.']);
    end
  end
else
  if nargin ==3
    cur_unit = arg1;
    img.data_mapper = cur_unit;
    new_unit = arg2;
  else
    new_unit = arg1;
  end
  img = data_mapper(img);
  cur_unit = img.current_params;
  
  % what to do about old parametrization??
%   img = rmfield(img, 'cur_unit');
end


[x do_nodes] = get_data(img);

y = map_data(x, cur_unit, new_unit);


if do_nodes
  img.node_data = y;
else
  img.elem_data = y;
end
img.current_params = new_unit;
y = img;

end

function x = map_data(x, cur_unit, new_unit)

if isempty(cur_unit) || isempty(new_unit)
  error('Unit string must not be empty.');
end

if strcmp(cur_unit, new_unit)
    return %nothing to do 
end

f = str2func(sprintf('%s2%s',cur_unit,new_unit));
try
   x = feval(f,x);
   x = fix_and_test(x);
catch err
  if strcmp(err.identifier,'MATLAB:UndefinedFunction')
    cur_pre = regexp(cur_unit,'^(.*?)_','match');
    if ~isempty(cur_pre) && ismember(cur_pre{1}, {'log_', 'log10_'});
      l = length(cur_pre{1});
      x = reverse(cur_pre{1}, x);
      x = map_data(x, cur_unit(l+1:end), new_unit);
      x = fix_and_test(x);
      return
    end
    new_pre = regexp(new_unit,'^(.*?)_','match');
    if ~isempty(new_pre) && ismember(new_pre{1}, {'log_', 'log10_', 'abs_'});
      l = length(new_pre{1});
      x = map_data(x, cur_unit, new_unit(l+1:end));
      x = apply(new_pre{1}, x);
      x = fix_and_test(x);
      return
    end
    
    %no direct function for conversion, see if we can help
    error('EIDORS:ConversionNotSupported', errorstr(cur_unit, new_unit));
  else
    rethrow(err);
  end
end

end

function x = fix_and_test(x)
  x(x == +inf) = +realmax;
  x(x == -inf) = -realmax;
  err_if_inf_or_nan(x, 'map_data-post');
end

function err_if_inf_or_nan(x, str)
  if any(isnan(x) | isinf(x))
      error(sprintf('bad %s (%d NaN, %d Inf of %d)', ...
                    str, ...
                    length(find(isnan(x))), ...
                    length(find(isinf(x))), ...
                    length(x)));
  end

end

function x = conductivity2resistivity(x)
  x = 1./x;
end

function x = resistivity2conductivity(x)
  x = 1./x;
end

function erstr = errorstr(cur_unit, new_unit)
  erstr = sprintf( ['Unit conversion from %s to %s is not supported.\n'...
             'Please write a function %s2%s.'], cur_unit, new_unit, cur_unit, new_unit);
end

function x = reverse(str, x);
  switch str
    case 'log10_'
      if any(x >= log10(realmax)-eps) warning('loss of precision -> inf'); end
      x = 10.^x;
    case 'log_'
      if any(x >= log(realmax)-eps) warning('loss of precision -> inf'); end
      x = exp(x);
    otherwise
      error('Cannot reverse %s', str);
  end
end


function x = apply(str, x)
  switch str
    case 'log10_'
      if any(x <= 0 + eps) warning('loss of precision -> -inf'); end
      x = log10(x);
    case 'log_'
      if any(x <= 0 + eps) warning('loss of precision -> -inf'); end
      x = log(x);
    case 'abs_'
      x = abs(x);
    otherwise
      error('Cannot apply %s', str);
  end
end
  
function [x, do_nodes] = get_data(img)
  do_nodes = 0;
  try
    x = img.elem_data;
  catch
    try
      x = img.node_data;
      do_nodes = 1;
    catch
      error('No elem_data or node_data found');
    end
  end
end

% test functions
function pass = do_unit_test
  pass = 1;
  d = 1;
  while d ~= 1 & d ~= 0
    d = rand(1);
  end
  disp('TEST: convert_img_units()');
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
end


function pass = test_map_data(data, in, out, expected, pass)
  fprintf('TEST: convert_img_units(%s -> %s)\n', in, out);
  if convert_img_units(data, in, out) ~= expected
     pass = 0;
     fprintf('TEST: FAIL convert_img_units(%s -> %s)\n', in, out);
  end
end
