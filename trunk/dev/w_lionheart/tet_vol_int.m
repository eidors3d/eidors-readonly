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

% NOTES ON BUGS: 17 dec 2012
% 1. FIXED for 
% vs =
%   -0.5541   -0.1847    0.7388
%   -0.2500   -0.2500    0.5000
%   -0.1847   -0.5541    0.7388
%         0   -0.5000    0.5000
%K>> convhulln(vs)
%Error using qhullmx
%QH6114 qhull precision error: initial simplex is not convex. Distance=0
% 2. For 
%  f_mdl =  mk_circ_tank(2,[0,1],0 );
%  c_mdl =  mk_circ_tank(1,[0,1],0 );
%  c_mdl.nodes = c_mdl.nodes*10;
%  c_mdl.nodes(:,3) = c_mdl.nodes(:,3) -7;
%  c_mdl.nodes(:,2) = c_mdl.nodes(:,2) +5;
%  c_mdl.nodes(:,1) = c_mdl.nodes(:,1) +2;
% all elements should be completely contained. WHy
% is convhulln running?
%
% 3. 
% Why are values different to meshing model




% (C) 2012 Bill Lionheart. License GPL v2 or v3
% $Id$

if isstr(v1) && strcmp(v1,'UNIT_TEST'); do_unit_test; return; end

           
if isstruct(v1)
   vol1 = get_elem_volume(v1);
   nn1 = num_nodes(v1);
   ne1 = num_elems(v1);
   e1 = v1.elems;
   v1 = v1.nodes;
else
   vol1 = tet_vol(v1);
   nn1 = 4;
   ne1 = 1;
   e1 = [1 2 3 4];
end
if isstruct(v2)
   vol2 = get_elem_volume(v2);
   nn2 = num_nodes(v2);
   ne2 = num_elems(v2);
   e2 = v2.elems;
   v2 = v2.nodes;
else
   vol2 = tet_vol(v2);
   nn2 = 4;
   ne2 = 1;
   e2 = [1 2 3 4];
end

volint = sparse(ne2,ne1);

[A1,b1]=tet_to_inequal(v1,e1);
[A2,b2]=tet_to_inequal(v2,e2);
i12 = (A1*v2' -b1*ones(1,nn2))< epsilon; % 4*ne1 X nn2
i21 = (A2*v1' -b2*ones(1,nn1))< epsilon; % 4*ne2 X nn1
% i12 = i12(:,reshape(e2',[],1));          % 4*ne1 X 4*ne2
% i21 = i21(:,reshape(e1',[],1));          % 4*ne2 X 4*ne1
% % is there a prettier way of doing this?
% % TRUE if elem of 2 contained in elem of 1
% ei12 = reshape(all(reshape(reshape(all(reshape(...
%           i12,4,[])),[],4*ne1)',[],4)')',[],ne1)'; % ne1 X ne2
% % TRUE if elem of 1 contained in elem of 2
% ei21 = reshape(all(reshape(reshape(all(reshape(...
%           i21,4,[])),4*ne2,[])',[],4)')',[],ne2)'; % ne2 X ne1
% % TRUE if elem of 2 does NOT have nodes in elem of 1       
% ni12 = reshape(all(reshape(reshape(~all(reshape(...
%           i12,4,[])),[],4*ne1)',[],4)')',[],ne1)'; % ne1 X ne2
% % TRUE if elem of 1 does NOT have nodes in elem of 2       
% ni21 = reshape(all(reshape(reshape(~all(reshape(...
%           i21,4,[])),4*ne2,[])',[],4)')',[],ne2)'; % ne2 X ne1
% disjoint = ni12' & ni21;
i12 =  4 == (i12(1:4:end,:) + i12(2:4:end,:) + i12(3:4:end,:) + i12(4:4:end,:));
i21 =  4 == (i21(1:4:end,:) + i21(2:4:end,:) + i21(3:4:end,:) + i21(4:4:end,:));
i12 = i12(:,reshape(e2',[],1));% ne1 X 4*ne2
i21 = i21(:,reshape(e1',[],1));% ne2 X 4*ne1
s12 = i12(:,1:4:end) + i12(:,2:4:end) + i12(:,3:4:end) + i12(:,4:4:end);
s21 = i21(:,1:4:end) + i21(:,2:4:end) + i21(:,3:4:end) + i21(:,4:4:end);
ei12 = s12 ==4;
ei21 = s21 ==4;
ni12 = s12 ==0;
ni21 = s21 ==0;
disjoint = ni12' & ni21;

if any(ei12(:))
   [el1 el2] = find(ei12);
   volint = volint + sparse(el2,el1,vol2(el2),ne2,ne1);
end
if any(ei21(:))
   [el2 el1] = find(ei21);
   volint = volint + sparse(el2,el1,vol1(el1),ne2,ne1);
end
todo = ~(volint~=0 | disjoint);
   
% if all(all(i12))
%     volint=tet_vol(v2);
% elseif all(all(i21))
%     volint=tet_vol(v1);
% elseif all(~all(i12)) & all(~all(i21));
%     % disjoint
%     volint=0;
% if all(~all(i12)) & all(~all(i21))
%    keyboard
% end
% VECTORISE FROM HERE ON
[y,x]=find(todo);
for i=1:length(x)
        id_x = (x(i)-1)*4 + (1:4);
        id_y = (y(i)-1)*4 + (1:4);
   volint(y(i),x(i)) = calc_int_vol( ...
                                     v1(e1(x(i),:),:), ...
                                     v2(e2(y(i),:),:), ...
                                     i12(x(i),id_y), ...
                                     i21(y(i),id_x), ...
                                     A1(id_x,:),...
                                     A2(id_y,:), ...
                                     b1(id_x,:), b2(id_y,:));
end

%v1,v2,i21,i12, A1, A2, b1,b2);
end

function epsi= epsilon; epsi=1e-10; end

function vs = do_choices_old(A1,A2,b1,b2,vs);
   choices =[ 1,2,1;1,2,2;1,2,3;1,2,4;
              1,3,1;1,3,2;1,3,3;1,3,4;
              1,4,1;1,4,2;1,4,3;1,4,4;                      
              2,3,1;2,3,2;2,3,3;2,3,4;
              2,4,1;2,4,2;2,4,3;2,4,4;
              3,4,1;3,4,2;3,4,3;3,4,4];
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
end
function volint = calc_int_vol(v1,v2,i12,i21, A1, A2, b1,b2);
   % List of choices of 2 from 1 and one from the other
   if 0
       clf
       hold on
       m1.type = 'fwd_model';
       m1.nodes = v1;
       m1.elems = [1 2 3 4];
       m2.type = 'fwd_model';
       m2.nodes = v2;
       m2.elems = [1 2 3 4];
       h1 = show_fem(m1)
       set(h1,'edgecolor','b')
       h2 = show_fem(m2)
       set(h1,'edgecolor','r')
       pause
   end
   % some intersection, 
   vs=[];
   % add the vertices that are already in both
   vs=[vs;v1(find(i21),:)];
   vs=[vs;v2(find(i12),:)];
   vs = do_choices_old(A1,A2,b1,b2,vs);

   %try all two faces from one intersected with one face from the
   %other
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

function vol=tet_vol(v)
   %finds the volume of one tetrahedron
   edges= v(2:end,:)-ones(3,1)*v(1,:);
   vol= abs(det(edges))/6;
end   

function do_unit_test 
%  simple_inequalities_test
tic
c= unit_test_smaller(0); % old c2f
toc

tic
a= unit_test_smaller(1); % for loop
toc

tic
b= unit_test_smaller(2); % vectorised
toc
  unit_test_cmp('unit_test_smaller', b,c, 0.2);
% keyboard
  unit_test_big

end

function simple_inequalities_test
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


end

function c2f= unit_test_smaller( select)
  f_mdl =  mk_circ_tank(2,[0,1],0 );
  c_mdl =  mk_circ_tank(1,[0,1],0 );
  c_mdl.nodes = c_mdl.nodes*10;
  c_mdl.nodes(:,3) = c_mdl.nodes(:,3) -7;
  c_mdl.nodes(:,2) = c_mdl.nodes(:,2) +5;
  c_mdl.nodes(:,1) = c_mdl.nodes(:,1) +2;
  

   nef = num_elems(f_mdl);
   nec = num_elems(c_mdl);

   c2f = sparse(nef,nec);

switch select;
case 1;
   for f = 1:nef
      vf = f_mdl.nodes(f_mdl.elems(f,:),:);
      for c = 1:nec
         vc = c_mdl.nodes(c_mdl.elems(c,:),:);
%disp([f,c]);
%   if f==23 && c==7 ; keyboard; end
         c2f(f,c) = tet_vol_int_org(vc,vf);
      end
   end
   c2f = c2f ./ repmat(get_elem_volume(f_mdl),1,num_elems(c_mdl));
case 2;
   c2f = tet_vol_int(c_mdl,f_mdl);
   c2f = c2f ./ repmat(get_elem_volume(f_mdl),1,num_elems(c_mdl));
case 0;
   c2f = mk_coarse_fine_mapping( f_mdl, c_mdl );
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
   disp('Old c2f:');
   tic
   c = mk_coarse_fine_mapping( fine_mdl, coarse_mdl);
   toc
   disp('Analytic:');
   tic
   a = tet_vol_int(coarse_mdl,fine_mdl);
   toc
%    tic
% %    a = a ./ repmat(get_elem_volume(fine_mdl),1,num_elems(coarse_mdl));
%    toc
   
%    unit_test_cmp('unit_test_big', a,c, 0.2);
   keyboard
end
