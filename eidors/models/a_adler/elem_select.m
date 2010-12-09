function memb_frac = elem_select( fmdl, select_fcn );
% ELEM_SELECT: select element fractions inside a function
%   memb_frac = elem_select( fmdl, select_fcn )
%
%  memb_frac = fraction of each element within the fcn
%  fmdl =      fwd_model structure
%  select_fcn = function to describe membership
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
% See Also:
%   mk_c2f_circ_mapping

if isstr(fmdl) && strcmp(fmdl,'UNIT_TEST'); do_unit_test; return; end


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

memb_frac = mean( feval(select_fcn,x,y,z), 2);

function do_unit_test;
    imdl = mk_common_model('a2c2',8);
    select_fcn = inline('(x-0.2).^2+(y-0.5).^2<0.2^2','x','y','z');
    memb_frac = elem_select( imdl.fwd_model, select_fcn);
    do_indiv_test('a2c2',find(memb_frac), [5, 10,18,26,27]);

    imdl = mk_common_model('n3r2',8);
    select_fcn = inline('(x-0.2).^2+(y-0.5).^2 + (z-1).^2<0.1^2','x','y','z');
    memb_frac = elem_select( imdl.fwd_model, select_fcn);
    do_indiv_test('n3r2',find(memb_frac), [156 159 162 168 431 434 437 503]);


function do_indiv_test(txt,a,b, tol)
   if nargin < 4; tol = 0; end
   fprintf('%10s = ',txt);
   ok='fail';
   try; if all(abs(a - b) <= tol);  ok='ok'; end; end
   disp(ok)
