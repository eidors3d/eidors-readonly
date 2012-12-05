function volint=tet_vol_int(v1,v2)
% volint=tet_vol_int(v1,v2)
% Find vertices of intersection volume of two tetrahedra
% the vertices of the tetrahedra are in 4x3 matrices v1 and v2 and
% need not be consistently oriented

% The algorithm is as follows. First calculate a system of linear
% inequalites so that the tetrahedron is x: A*x -b <=0, then use
% this to see if the tetrahedra have an intersection. If they do 
% the intersection is a convex polyhedron, the vertices of which are the
% intersection of planes corresponding to faces of the tetraherda.
% Other than vertices of one tetrahedron that are in the other tetrahedron 
% these vertices are the intersection of two planes from one tetrahedron
% and one from the other. All such points are calculated and tested to see
% if they lie within both tetrahedra. The convex hull of these points is
% the intersection.  


% (C) 2012 Bill Lionheart. License GPL v2 or v3
% $Id$

if isstr(v1) && strcmp(v1,'UNIT_TEST'); do_unit_test; return; end

% List of choices of 2 from 1 and one from the other
choices =[ 1,2,1;1,2,2;1,2,3;1,2,4;
           1,3,1;1,3,2;1,3,3;1,3,4;
           1,4,1;1,4,2;1,4,3;1,4,4;                      
           2,3,1;2,3,2;2,3,3;2,3,4;
           2,4,1;2,4,2;2,4,3;2,4,4;
           3,4,1;3,4,2;3,4,3;3,4,4];
           
           
           
    
epsilon=1e-10;
[A1,b1]=tet_to_inequal(v1);
[A2,b2]=tet_to_inequal(v2);
i12 = (A1*v2' -b1*ones(1,4))< epsilon;
i21 = (A2*v1' -b2*ones(1,4))< epsilon;
if all(all(i12))
    volint=tet_vol(v2);
elseif all(all(i21))
    volint=tet_vol(v1);
elseif all(~all(i12)) & all(~all(i21));
    % disjoint
    volint=0;
else
   % some intersection, 
   vs=[];
   % add the vertices that are already in both
   vs=[vs;v1(find(all(i21)),:)];
   vs=[vs;v2(find(all(i12)),:)];

   %try all two faces from one intersected with one face from the
   %other
   for i = 1:24
       A=[A1(choices(i,1:2),:);A2(choices(i,3),:)];
       b=[b1(choices(i,1:2));b2(choices(i,3))];
       
       if abs(det(A))>1e-5
        vs=[vs;(A\b)'];
       end
         
       A=[A2(choices(i,1:2),:);A1(choices(i,3),:)];
       b=[b2(choices(i,1:2));b1(choices(i,3))];
    
       if abs(det(A))>1e-5
        vs=[vs;(A\b)'];
       end
       
   end
   if isempty(vs)
      volint=0;
   else
     thosein = find( all(A1 * vs'- b1*ones(1,size(vs,1),1) <epsilon) &  all(A2 * vs'- b2*ones(1,size(vs,1)) <epsilon,1));
     vs=vs(thosein,:);
     vs = 1e-5*round(1e5*vs);
     vs=unique(vs,'rows');
%    vs=uniquetol(vs,1e-10,'rows');
     if size(vs,1)<4
         volint=0;
     else   
       try
          [K,volint]=convhulln(vs); %FIXME: get info from convhulln
       catch err 
          if strcmp(err.identifier, 'MATLAB:qhullmx:DegenerateData')
             volint = 0;
          else 
             rethrow( err); 
          end 
       end
     end
   end    
end
    
end

function vol=tet_vol(v)
   %finds the volume of one tetrahedron
   edges= v(2:end,:)-ones(3,1)*v(1,:);
   vol= abs(det(edges))/6;
end   

function do_unit_test 
   % Test inequalities
   v1=[0,0,0;eye(3)];
   v2 = v1;v2(1,:)=v2(1,:)+0.1;
        
   out =  tet_vol_int(v1,v2);
   correct = 7/60;
   unit_test_cmp('Shifted rightangle tetrahedron', out, correct,1e-14)

   % test volume

   out =  tet_vol(v1);
   correct = 1/6;
   unit_test_cmp('Unit tetrahedron volume', out, correct)

   A=[1,4,-1;3,4,1;0,0,2];
   v3 =(A*v1')'+ ones(4,1)*[1,4,-3];
   volu=out;
   out= tet_vol(v3);
   correct = volu * abs(det(A));
   unit_test_cmp('Scaled shifted tetrahedron volume', out, correct)

tic
   unit_test_smaller;
toc

end

function unit_test_smaller
  f_mdl =  mk_circ_tank(3,[0,1],0 );
  c_mdl =  mk_circ_tank(2,[0,1],0 );

   nef = num_elems(f_mdl);
   nec = num_elems(c_mdl);

   mapping = sparse(nef,nec);

   for f = 1:nef
      vf = f_mdl.nodes(f_mdl.elems(f,:),:);
      for c = 1:nec
         vc = c_mdl.nodes(c_mdl.elems(c,:),:);
%disp([f,c]);%if f==10 && c==2 ; keyboard; end
         mapping(f,c) = tet_vol_int(vc,vf);
      end
   end
end

    
function unit_test_big
   electrodes_per_plane = 16;
   number_of_planes = 2;
   tank_radius = 0.2;
   tank_height = 0.5;
   fine_mdl = ng_mk_cyl_models([tank_height,tank_radius],...
       [electrodes_per_plane,0.15,0.35],[0.01]);

   % Create coarse model for inverse problem
   coarse_mdl_maxh = 0.07; % maximum element size 
   coarse_mdl = ng_mk_cyl_models([tank_height,tank_radius,coarse_mdl_maxh],[0],[]);

   disp('Calculating coarse2fine mapping ...');
   inv3d.fwd_model.coarse2fine = ...
          mk_coarse_fine_mapping_analytic( fine_mdl, coarse_mdl);
end
