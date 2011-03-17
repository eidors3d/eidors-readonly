function out = shape_library(action, shape, varargin)
%SHAPE_LIBRARY Common shapes for models
%  SHAPE_LIBRARY('list') lists available shapes (strings)
%  
%  SHAPE_LIBRARY('list',SHAPE) lists available outlines for the given shape
%  
%  SHAPE_LIBRARY('get',SHAPE) returns a struct with all the information for 
%  the given shape, including at least the fields:
%     .boundary  - a set of [x y] points defining the boundary
%     .pic.img   - an image from which the contours were extracted
%     .pic.X 
%     .pic.Y     - coordinate vectors to use with IMAGESC
%     .copyright - source, authors and license
%  Additional outlines of internal objects (e.g. lungs) are defined in
%  separate fields akin to .boundary above.
%
%  SHAPE_LIBRARY('get', SHAPE, FIELD [,FIELD2,..]) 
%  SHAPE_LIBRARY('get', SHAPE, {FIELD [,FIELD2,...]})
%  returns a requested outline FIELD from SHAPE. If more than one FIELD is 
%  requested, a cell array of [x y] matries is returned
%
%EXAMPLES:
% shape_library('list');
% shape_library('list','pig_32kg');
% shape_library('get','pig_32kg','trunk','lungs')

% (C) Bartlomiej Grychtol, 2011. License: GPL version 2 or version 3
% $Id$

n_IN = nargin;

switch action
    case 'UNIT_TEST'
        do_unit_test;
        return
    case 'list'
        if n_IN>1
            out = list_shapes(shape);
        else
            out = list_shapes('');
        end
    case 'show'
        show_shape;
    case 'get'
        get_shape;
end

function out = list_shapes(shape)
out = '';
if isempty(shape)
    out = who('*','-file','shape_library.mat');
    disp('Available shapes:');
else
    s = load('shape_library.mat',shape);
    try
       out = fieldnames(s.(shape)); 
    catch, not_found(shape);end
    disp(['Shape ' shape ' contains:']);     
end
disp(out);

function show_shape
disp('show_shape');

function get_shape
disp('get_shape');

function not_found(shape)
eidors_msg(['SHAPE_LIBRARY: Didn''t find shape ' shape]);


function do_unit_test
    shape_library('list');
    shape_library('list','a-shape-we-dont-have');% give error
    shape_library('list','pig_23kg');
    shape_library('show','pig_23kg');
    shape_library('get','pig_23kg','trunk','lungs')
