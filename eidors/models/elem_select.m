function memb_frac = elem_select( fmdl, select_fcn )
% ELEM_SELECT: select element fractions inside a function
%   memb_frac = elem_select( fmdl, select_fcn )
%
%  memb_frac = fraction of each element within the fcn
%  fmdl =      fwd_model structure
%  select_fcn = function to describe membership, 
%              can also be a cell array of functions, but all must accept
%              parameters x, y, and z. ELEM_SELECT then finds elements that
%              satisfy ALL the functions at once
%              finally, select_fcn can be a string which will then
%                be turned into a function via inline(str,'x','y','z')
% OR
% select_fcn = string accepting named variables x,y,z.
%
% parameters
%   fwd_model.elem_select.interp_no  - interpolation density
%
% Example:
%   img = mk_image(mk_common_model('b2d1c',8));
%   select_fcn = inline('(x-0.2).^2+(y-0.5).^2<0.2^2','x','y','z');
%   memb_frac = elem_select( img.fwd_model, select_fcn)
%   img.elem_data = 1 + memb_frac*0.1;
%   show_fem(img);
%
% Example
%   img = mk_image(mk_common_model('b2d1c',8));
%   select_fcn = '(x-0.2).^2+(y-0.5).^2<0.2^2';
%   memb_frac = elem_select( img.fwd_model, select_fcn)
%   img.elem_data = 1 + memb_frac*0.1;
%   show_fem(img);
%
% See Also:
%   mk_c2f_circ_mapping

if ischar(fmdl) && strcmp(fmdl,'UNIT_TEST'); do_unit_test; return; end


% 4 for 2D, 3 for 3D
dims = size(fmdl.nodes,2);
interp_no = 6 - dims;
try 
   interp_no = fmdl.elem_select.interp_no;
end

pts = interp_mesh( fmdl, interp_no );
x = squeeze(pts(:,1,:));
y = squeeze(pts(:,2,:));
if dims ==2;
  z = 0*x;
else
  z = squeeze(pts(:,3,:));
end
if ischar(select_fcn)
%  we have a string, create a function
    select_fcn = inline(select_fcn, 'x','y','z');
    memb_frac = mean( feval(select_fcn,x,y,z), 2);
elseif ~iscell(select_fcn) 
    % the normal case  
    memb_frac = mean( feval(select_fcn,x,y,z), 2);
else
    % many functions case
    memb_val = ones(size(x));
    for i = 1:numel(select_fcn)
        memb_val = memb_val .* feval(select_fcn{i},x,y,z);
    end
    memb_frac = mean(memb_val,2);
end


function do_unit_test;
    imdl = mk_common_model('a2c2',8);
    select_fcn = '(x-0.2).^2+(y-0.5).^2<0.2^2';
    memb_frac = elem_select( imdl.fwd_model, select_fcn);
    unit_test_cmp('a2c2 (string)',find(memb_frac), [5, 10,18,26,27]');

    imdl = mk_common_model('a2c2',8);
    select_fcn = inline('(x-0.2).^2+(y-0.5).^2<0.2^2','x','y','z');
    memb_frac = elem_select( imdl.fwd_model, select_fcn);
    unit_test_cmp('a2c2',find(memb_frac), [5, 10,18,26,27]');


    select_fcn2= inline('y<0.4','x','y','z');
    memb_frac = elem_select( imdl.fwd_model, {select_fcn,select_fcn2});
    unit_test_cmp('a2c2 (2fcns)',find(memb_frac), [5, 10,18,27]');

    imdl = mk_common_model('n3r2',[16,2]);
    select_fcn = inline('(x-0.2).^2+(y-0.5).^2 + (z-1).^2<0.1^2','x','y','z');
    memb_frac = elem_select( imdl.fwd_model, select_fcn);
    unit_test_cmp('n3r2',find(memb_frac), [156 159 162 168 431 434 437 503]');
