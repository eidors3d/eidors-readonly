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
%   if the fcn_handle is a string, a function with x,y,z is created
%     crop_model(gca, 'x+y>0') % same as previous
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
% $Id$

if ischar(axis_handle) && strcmp(axis_handle,'UNIT_TEST'); do_unit_test; return; end

% TODO (update 2 apr 2012):
%  - make crop_model work for 2D fems

if isstr(fcn_handle)
  fcn_handle = inline(fcn_handle,'x','y','z');
end


if isempty(axis_handle)
   axis_handle= gca;
   crop_graphics_model(axis_handle, fcn_handle);
elseif ishandle( axis_handle )
   crop_graphics_model(axis_handle, fcn_handle);
elseif axis_handle.type == 'fwd_model';
   [fmdl,c2f_idx]= crop_fwd_model(axis_handle, fcn_handle);
elseif axis_handle.type == 'image';
   [fmdl,c2f_idx]= crop_fwd_model(axis_handle.fwd_model, fcn_handle);
   keyboard
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
% Keep these nodes
   idx0=  all( feval(fcn_handle,xyz(:,1), ...
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
   n_elecs = length(fmdl1.electrode);
   rm_elec_list = zeros(n_elecs,1);
   for i=1:n_elecs;
      el_nodes= fmdl0.electrode(i).nodes;
      el_nodes(:)= idx1(el_nodes);
      if any(el_nodes==0);
         rm_elec_list(i) = 1;
      end
      fmdl1.electrode(i).nodes= el_nodes;
   end
   if any(rm_elec_list)
      str = sprintf('%d,', find(rm_elec_list));
      eidors_msg('crop_model: removing electrodes (%s)',str(1:end-1),1);
      fmdl1.electrode = fmdl1.electrode( find(~rm_elec_list));
   end
end


function do_unit_test
   imdl = mk_common_model('a2c0',8); fmdl= imdl.fwd_model;
   fmdl = crop_model(fmdl,inline('x<0','x','y','z'));
   unit_test_cmp('crop_model-a2c0-01',length(fmdl.electrode),5);
   unit_test_cmp('crop_model-a2c0-02',size(fmdl.elems),[32,3]);
   unit_test_cmp('crop_model-a2c0-03',size(fmdl.nodes),[25,2]);

   fmdl = crop_model(fmdl,'x<0'); % verify it's same
   unit_test_cmp('crop_model-str-a2c0-01',length(fmdl.electrode),5);
   unit_test_cmp('crop_model-str-a2c0-02',size(fmdl.elems),[32,3]);
   unit_test_cmp('crop_model-str-a2c0-03',size(fmdl.nodes),[25,2]);

   imdl = mk_common_model('n3r2',[16,2]); fmdl= imdl.fwd_model;
   fmdl = crop_model(fmdl,inline('x<0','x','y','z'));
   unit_test_cmp('crop_model-n3r2-01',length(fmdl.electrode),16);
   unit_test_cmp('crop_model-n3r2-02',size(fmdl.elems),[342,4]);
   unit_test_cmp('crop_model-n3r2-03',size(fmdl.nodes),[128,3]);
   fmdl = crop_model(fmdl,inline('z<2','x','y','z'));
   unit_test_cmp('crop_model-n3r2-04',length(fmdl.electrode),8);
   unit_test_cmp('crop_model-n3r2-05',size(fmdl.elems),[114,4]);
   unit_test_cmp('crop_model-n3r2-06',size(fmdl.nodes),[64,3]);




   subplot(331)
   imdl = mk_common_model('n3r2',[16,2]); fmdl= imdl.fwd_model;
   show_fem(fmdl);
   subplot(332)
   show_fem(fmdl); hh= gca; 
   subplot(333)
   show_fem(fmdl);
   crop_model([],inline('z<2','x','y','z'));
   crop_model(hh,inline('x<0','x','y','z'));

   subplot(334)
   fmdlc = fmdl;
   fmdlc = crop_model(fmdlc,inline('x<0','x','y','z'));
   show_fem(fmdlc);

   subplot(337)
   imdl = mk_common_model('d2c2'); fmdl= imdl.fwd_model;
   show_fem(fmdl);
   subplot(338)
   show_fem(fmdl); hh= gca; 
   title('expected fail');
   subplot(339)
   show_fem(fmdl);
   crop_model([],inline('y<0','x','y','z'));
   crop_model(hh,inline('x<0','x','y','z'));
   title('expected fail');
   

