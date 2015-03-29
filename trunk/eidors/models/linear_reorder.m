function [fwd_model] = linear_reorder(fwd_model,ccw)
%[fwd_model] = linear_reorder(fwd_model,ccw)
%Function to reorder local nodes (counter)clockwise per element
%Input:  - fwd_model structure
%        - ccw = -1 (default) - counter clockwise OR 1 - clockwise   
%Output: - fwd_model structure (only .elems changes)
%NOTE:Function only for linear triangles, since in this case, identity:
%         No. of nodes/element = No. spatial dimensions + 1

% (C) 2011 Michael Crabb. License: GPL version 2 or version 3
% $Id$

if isstr(fwd_model) && strcmp(fwd_model,'UNIT_TEST'); do_unit_test; return; end

if (nargin==1) 
    ccw=-1; %Default specify counter-clockwise nodes
end

   fwd_model =     do_reorder(fwd_model, ccw);
%  fwd_model = old_do_reorder(fwd_model, ccw);
end


function fmdl = do_reorder(fmdl, ccw)
   dim=size( fmdl.nodes,2);
   els=num_elems( fmdl );

   if dim==3 && elem_dim( fmdl ) == 2
      % For a surface in 3D, we need to walk across elements
      fmdl = boundary_align( fmdl );
      return
   end

     xx= reshape( fmdl.nodes( fmdl.elems, 1 ), els, dim+1);
     xx= xx - xx(:,1)*ones(1,dim+1);
     yy= reshape( fmdl.nodes( fmdl.elems, 2 ), els, dim+1);
     yy= yy - yy(:,1)*ones(1,dim+1);
   if dim==2;
     vol = xx(:,2).*yy(:,3) - xx(:,3).*yy(:,2);
   elseif dim==3
     zz= reshape( fmdl.nodes( fmdl.elems, 3 ), els, dim+1);
     zz= zz - zz(:,1)*ones(1,dim+1);

     vol = zz(:,4).*( xx(:,2).*yy(:,3) - xx(:,3).*yy(:,2) ) ...
         - zz(:,3).*( xx(:,2).*yy(:,4) - xx(:,4).*yy(:,2) ) ...
         + zz(:,2).*( xx(:,3).*yy(:,4) - xx(:,4).*yy(:,3) );

% Looks more elegant, but slower
%    vol = sum(zz(:,[4,3,2]).*xx(:,[2,4,3]).*yy(:,[3,2,4]) ,2) ...
%        - sum(zz(:,[4,3,2]).*xx(:,[3,2,4]).*yy(:,[2,4,3]) ,2);

   else
     error('reorder for 2 or 3 dimensions');
   end

   ff = find( sign(vol) == ccw );
   % reverse first two nodes
   fmdl.elems(ff,1:2) = fmdl.elems(ff,[2,1]);
end

function mdl = boundary_align( mdl );
   Ne = num_elems(mdl);
   checked = ( 1:Ne )';
   tocheck = zeros(Ne,1); tocheck(1) = 1; maxptr = 2;
   edges = mdl.elem2edge;

   for i=1:num_elems(mdl)
      curr = tocheck(i);
      curredges = edges(curr,:);
      for j=curredges(:)'
         ff = any( j == edges,2 );
         ff(curr) = 0;
         ff = find(ff); % Should only have 1
         if length(ff)>1;
            error('Only one elem should share this edge. unexpected error.');
         end
         if checked(ff) == 0; % we've already done it
            continue;
         end
         nodes_j = mdl.edges(j,:); % Nodes in this edge
         a1= find(mdl.elems(curr,:) == nodes_j(1));
         a2= find(mdl.elems(curr,:) == nodes_j(2));
         curr_o = rem(3+a1-a2,3);
         b1= find(mdl.elems(ff  ,:) == nodes_j(1));
         b2= find(mdl.elems(ff  ,:) == nodes_j(2));
         ff_o = rem(3+b1-b2,3);
         if (ff_o ~= curr_o); 
            mdl.elems(ff,[b1,b2]) = mdl.elems(ff,[b2,b1]);
         end
         % Fix orientation
         tocheck(maxptr) = ff; maxptr= maxptr+1;     
         checked(ff)= 0; % set as done
      end
   end
      
end 

function fwd_model = old_do_reorder(fwd_model, ccw)
   nodecoords = fwd_model.nodes; %Cache coorindates of nodes [nnodesxnodedim]
   elementnodes = fwd_model.elems; %Cache matrix of elements [eletotalxelenode]

   eletotal = size(elementnodes,1); %No. of elements
   elenode = size(elementnodes,2); %No. of nodes per element
   nodedim = size(fwd_model.nodes,2);
   midpoint = mean(fwd_model.nodes(unique(fwd_model.elems),:));

   for e=1:eletotal; %Loop over all elements
       %Row vector of global nodes [1xelenode]
       enodes = elementnodes(e,:); 
       %Matrix of nodal positions [elenodexdim] (Linear dimension==elenode-1) 
       nd = nodecoords(enodes,:);
       
       % surface meshes need tweaking. Use the midpoint to fit the 3D formula.
       % This will not work for non simply-connected surfaces.
       if elenode == 3 && nodedim == 3
          nd = [nd; midpoint]; 
       end
       %Calculate area of triangle/volume defined by the elements nodes
       %In 2D this is area and in 3D this is volume
       area= det([ones(length(nd),1),nd]);
       areasign=sign(area); 
       
       %If sign is (pos) neg swap two nodes (last two will suffice..)
       if(areasign == ccw) %Swap last two entries of enodes 
           temp=enodes(elenode-1);
           enodes(elenode-1)=enodes(elenode);
           enodes(elenode) = temp;
           %elementnodes(e,:)=enodes; %Put back into elementnodes matrix
       end
       elementnodes(e,:)=enodes; %Put enodes back into elementnodes matrix
   end
   fwd_model.elems=elementnodes; %Reassign fwd_model.elems
end

function do_unit_test
   for i = 1:10
     clear imdl;
     switch i,
       case 1; imdl = mk_common_model('n3r2',[16,2]);
       case 2; imdl = mk_common_model('a2c2',8);
       case 3; imdl = mk_common_model('d3cr',[16,2]);
       case 4; imdl = mk_common_model('f3cr',[16,2]);
     end
     if ~exist('imdl'); continue ; end

     fmdl = imdl.fwd_model;
     vol = test_linear_reorder( fmdl ); ok = std(sign(vol))==0; % not all 
     t = cputime;
     fm0 = linear_reorder(fmdl, -1);
     fm1 = linear_reorder(fmdl, 1);
     t = cputime - t;
     vol0= test_linear_reorder( fm0 );
     vol1= test_linear_reorder( fm1 );

     fprintf('test%02d(t=%4.2f): OK=%d=>(%d,%d)\n',i, t, ...
          ok, all(vol0>0), all(vol1<0));
   end
end
   

function vol = test_linear_reorder(fwd_model)

   dim=size(fwd_model.nodes,2); elee=size(fwd_model.elems,1);

   for e=1:elee
       b=fwd_model.elems(e,:);  [v]=fwd_model.nodes(b,:);
           for i=1:dim
               vv1(i,:)=v(i+1,:)-v(1,:);
           end
       vol(e)=det([vv1]);
   end
end

