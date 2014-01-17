function img = apply_bladder_parametrization(img,dir)

if nargin == 1;
   dir = 0;
end
if dir==0
   bladder = @(r,v) 1 - v * elem_select(img.fwd_model,inline(sprintf('(x-%f).^2 + (y-%f).^2 + (z-%f).^2 <= %f^2',-r-0.01,0,r+0.01,r),'x','y','z'));
   p = img.conductivity.params;
   img.elem_data = bladder(p(1),0.5);
   img.current_physics = 'conductivity';
else
   warning('reverse mapping is not possible');
end
