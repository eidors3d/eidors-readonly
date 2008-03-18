function fmdl= crop_model( axis_handle, fcn_handle );
% CROP_MODEL: Crop away parts of a fem model
%
% USAGE #1: crop display to show model internals
%   crop_model( axis_handle, fcn_handle );
% 
%   if axis_handle==[], then use the current axis
%   examples:
%     crop_model([],  inline('z==3','x','y','z'))
%     crop_model(gca, inline('x+y>0','x','y','z'))
%
% USAGE #2: crop fwd_model to create new fwd_model
%   fmdl_new= crop_model( fwd_model, fcn_handle );
% 
%   example:
%   fmdl2= crop_model(fmdl1, inline('x+y>0','x','y','z'))

% (C) 2006-2008 Andy Adler. License: GPL version 2 or version 3
% $Id: crop_model.m,v 1.16 2008-03-18 15:36:26 aadler Exp $

usage_graphics= 1;
try if axis_handle.type == 'fwd_model'
   usage_graphics= 0;
end; end

if usage_graphics
   if isempty(axis_handle)
      axis_handle= gca;
   end
   crop_graphics_model(axis_handle, fcn_handle);
else
   fmdl= crop_fwd_model(axis_handle, fcn_handle);
end


function crop_graphics_model(axis_handle, fcn_handle);
   kk= get(axis_handle,'Children');
   % only the FEM frame
   %k=kk( find( kk== min(kk) ));

   for k= kk(:)'
      try
         x= get(k,'XData');
         y= get(k,'YData');
         z= get(k,'ZData');
         c= get(k,'CData');
         idx= ~all( feval(fcn_handle,x,y,z) );
         if any(size(c)>[1,1])
            set(k,'Xdata', x(:,idx), ...
                  'Ydata', y(:,idx), ...
                  'Zdata', z(:,idx), ...
                  'Cdata', c(:,idx));
         else
            set(k,'Xdata', x(:,idx), ...
                  'Ydata', y(:,idx), ...
                  'Zdata', z(:,idx));
         end
      end
   end
