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
% $Id$

if isstr(axis_handle) && strcmp(axis_handle,'UNIT_TEST'); do_unit_test; return; end

% TODO (update 2 apr 2012):
%  - make crop_model work for 2D fems

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
         v = get(k, 'Vertices');
         
         sel = feval(fcn_handle, v(:,1), v(:,2) , v(:,3) );
         sel(end+1) = 1;
         f = get(k, 'Faces');
         f(isnan(f)) = length(sel);
         idx= ~any( sel(f') );
         % show faces fully in
         F = f(idx,:);
         set(k, 'Faces', F);
         c = get(k,'CData');
%          if numel(unique(c)) ~= 1
%             return;
%          end
         
         set(k, 'CData', c(idx));
         
         % find intersected edges
         idx = (1:size(f,2))'; idx(:,2) = circshift(idx,-1);
         edges = reshape( f(:,idx), [], 2);
         
         e2f = reshape(repmat((1:length(f))',1,3)',[],1);
         
         % remove any "edges" with nans
         idx = any(edges' == length(sel));
         edges(idx,:) = [];
         e2f  (idx  ) = [];
         ein = sel(edges);
         idx = sum( ein, 2) == 1;
         % only keep the edges that cross
         edges = edges(idx,:); ein = ein(idx,:);
         e2f   = e2f  (idx);
         if ~isempty(edges);
            % sort them from in to out
            [jnk ix] = sort(ein,2,'ascend');
            ix(ix==2) = size(ix,1)+1;
            ix = (0:length(ix)-1)'*ones(1,2) + ix;
            edges = edges(ix);
            % the function can be anything, so we find a numerical
            % approximation rather than the true solution
            N = 1000;
            LS = linspace(0,1,N+2); LS([1 end]) = [];
            dr = v(edges(:,2),:) - v(edges(:,1),:);
            X = v(edges(:,1),1)*ones(1,N) + dr(:,1)*ones(1,N)*diag(LS);
            Y = v(edges(:,1),2)*ones(1,N) + dr(:,2)*ones(1,N)*diag(LS);
            Z = v(edges(:,1),3)*ones(1,N) + dr(:,3)*ones(1,N)*diag(LS);
            val = feval(fcn_handle, X, Y , Z );
            d = diff(val,1,2);
            d( val(:,1) == 1, 1) = 1;
            d( val(:,N) == 0, N) = 1;
            d = logical(d);
            % why does x = X(d) not work??
            x = X; x = x'; x = x(d'); 
            y = Y; y = y'; y = y(d'); 
            z = Z; z = z'; z = z(d'); 
            % add new vertices
            vlen = size(v,1);
            set(k,'Vertices',[v; x y z]);
            % make edges use the new vertices
            edges(:,2) = vlen + (1:size(edges,1));
            edges(:,3:size(F,2)) = NaN;
%             keyboard
            set(k,'Faces',[F;edges]);
%             c = get(k,'FaceVertexCData');
%             switch length(c)
%                case length(F)
%                   c(end:end+length(edges),:) = 0;
%                   set(k,'FaceVertexCData',c);
%                case length(v)
%                   c(end:end+length(x),:) = 0;
%                   set(k,'FaceVertexCData',c);
%             end
         end
         
         
      end
   end

% CROP fwd_model
function [fmdl1,c2f_idx]= crop_fwd_model(fmdl0, fcn_handle)
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
   

