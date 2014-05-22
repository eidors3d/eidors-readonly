function img = apply_circle_parametrization_3d(img,dir)

if nargin == 1;
   dir = 0;
end
if dir==0
   circ = @(r,x,y,z) 1 - 0.5 * elem_select(img.fwd_model,inline(sprintf('(x-%f).^2 + (y-%f).^2 + (z-%f).^2 <= %f^2',x,y,z,r),'x','y','z'));
   p = img.conductivity.params;
   img.elem_data = circ(p(1),p(2),p(3),p(4));
   img.current_params = 'conductivity';
else
   warning('reverse mapping is not possible');
end
