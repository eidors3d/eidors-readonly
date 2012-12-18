function img = convert_units(img,new_unit)
%CONVERT_UNITS change image data units 
%  img = convert_units(img,new_unit) converts img.elem_data or img.node_data
%  expressed in img.current_physics into different units for the same 
%  physical property
% 
% Example: 
%  img = mk_image(mdl, 2, 'resisitivity');
%  img = convert_units(img, 'conductivity');

% (C) 2012 Bartlomiej Grychtol, License: GPL version 2 or 3
% $Id$

try 
    cur_unit = img.current_physics;
catch
    error('Image must contain the current_physics field');
end

if strcmp(cur_unit, new_unit)
    return %nothing to do
end

f = str2func(sprintf('%s2%s',cur_unit,new_unit));
try
    img = feval(f,img);
catch err
    if strcmp(err.identifier,'MATLAB:UndefinedFunction')
        error('EIDORS:ConversionNotSupported', ...
            'Unit conversion from %s to %s is not supported', cur_unit, new_unit);
    else     
        rethrow(err);
    end
end
end

function img = resistivity2conductivity(img)
    img = apply_func(img,@(x)1./x);
    img.current_physics = 'conductivity';
end

function img = log_conductivity2conductivity(img)
    img = apply_func(img,@(x)exp(x));
    img.current_physics = 'conductivity';
end

function img = log_resistivity2conductivity(img)
    img = apply_func(img,@(x)exp(x));
    img = apply_func(img,@(x)1./x);
    img.current_physics = 'conductivity';
end

function img = conductivity2resistivity(img)
    img = apply_func(img,@(x)1./x);
    img.current_physics = 'resistivity';
end

function img = conductivity2log_conductivity(img)
    img = apply_func(img,@(x)log(x));
    img.current_physics = 'log_conductivity';
end

function img = conductivity2log_resistivity(img)
    img = apply_func(img,@(x)1./x);
    img = apply_func(img,@(x)log(x));
    img.current_physics = 'log_resistivity';
end

function img = apply_func(img,func)
    try 
        img.elem_data = feval(func,img.elem_data);
    catch
        try 
            img.node_data = feval(func,img.node_data); 
        catch
            error('No elem_data or node_data found');
        end
    end
end
