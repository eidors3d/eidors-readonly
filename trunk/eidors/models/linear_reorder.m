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

fwd_model = old_do_reorder(fwd_model, ccw);
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

