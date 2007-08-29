function crop_model( axis_handle, fcn_handle );
% crop_model( axis_handle, fcn_handle );
% Crop away parts of a fem model
% 
% if axis_handle==[], then use the current axis
% examples:
%   crop_model([],  inline('z==3','x','y','z'))
%   crop_model(gca, inline('x+y>0','x','y','z'))

% (C) 2006 Andy Adler. License: GPL version 2 or version 3
% $Id: crop_model.m,v 1.6 2007-08-29 09:12:07 aadler Exp $

if exist('OCTAVE_VERSION');
   warning('show_fem does not support octave');
   return
end

if isempty(axis_handle)
   axis_handle= gca;
end

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
