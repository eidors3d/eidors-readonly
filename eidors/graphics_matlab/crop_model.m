function [fmdl,c2f_idx]= crop_model( axis_handle, fcn_handle );
% CROP_MODEL: Crop away parts of a fem model
%
% USAGE #1: crop display to show model internals
%   crop_model( axis_handle, fcn_handle );
%
%   fcn_handle ==1 where model is *kept*
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
%
%  with two parameters output
% [fmdl,c2f_idx]= crop_model( axis_handle, fcn_handle );
%     c2f_idx maps each elemen in fmdl_new to fwd_model

% (C) 2006-2008 Andy Adler. License: GPL version 2 or version 3
% $Id: crop_model.m,v 1.19 2008-06-11 14:47:37 aadler Exp $

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
   [fmdl,c2f_idx]= crop_fwd_model(axis_handle, fcn_handle);
end

% CROP GRAPHICS
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

% CROP fwd_model
function [fmdl1,c2f_idx]= crop_fwd_model(fmdl0, fcn_handle);
   fmdl1= fmdl0;

% Find nodes to remove
   nodes= fmdl0.nodes;
   [n,d]= size(nodes);
   n2xyz= eye(d,3); 
   xyz= nodes*n2xyz;
   idx0= ~all( feval(fcn_handle,xyz(:,1), ...
                                xyz(:,2), ...
                                xyz(:,3)),2);
% remove these nodes
   fmdl1.nodes(idx0,:) = [];

% renumber nodes, set unused ones to 0
   idx1= zeros(n,1);
   idx1(~idx0)= 1:sum(~idx0);

   fmdl1.elems(:) = idx1(fmdl0.elems);
   remove= any( fmdl1.elems == 0, 2);
   fmdl1.elems(remove,:)= [];
% c2f_idx maps each new elem to its original position
   c2f_idx= find(~remove);

   fmdl1.boundary(:) = idx1(fmdl0.boundary);
   remove= any( fmdl1.boundary == 0, 2);
   fmdl1.boundary(remove,:)= [];

% renumber nodes, set unused ones to 0
if isfield(fmdl1,'electrode');
   for i=1:length(fmdl1.electrode)
      el_nodes= fmdl0.electrode(i).nodes;
      el_nodes(:)= idx1(el_nodes);
      if any(el_nodes==0);
         error('crop_model: nodes in electrode are removed');
      end
      fmdl1.electrode(i).nodes= el_nodes;
   end
end

