function img = apply_circle_parametrization(img,dir)

if nargin == 1;
   dir = 0;
end
if dir==0
   circ = @(r,x,y,v) 1 - v * elem_select(img.fwd_model,inline(sprintf('(x-%f).^2 + (y-%f).^2 <= %f^2',x,y,r),'x','y','z'));
   p = img.conductivity.params;
   img.elem_data = circ(p(1),p(2),p(3),p(4));
   img.current_physics = 'conductivity';
else
   warning('reverse mapping is not possible');
end
