function [mapping failed] = mk_c2f_circ_mapping( mdl, xyzr );
% MK_C2F_CIRC_MAPPING: create a mapping matrix from circles/spheres to FEM
% mapping= mk_c2f_circ_mapping( mdl, xyzr );
%
% Mapping approximates elem_data_fine from elem_data_coase
%   elem_data_model = Mapping * elem_data_circles
%
% mdl is coarse fwd_model
% xyzr is the 3xN matrix (2D) or 4xN matrix (3D) of
%      circle centres and radii
%
% this function approximates using points interpolated into elements
%   use mdl.interp_mesh.n_points to control interpolation density  
%
% if a 3xN matrix is specified for a 3D model, then cylindrical
%  shapes (circle extruded in z) are selected
%

% (C) 2009 Andy Adler and Bartlomiej Grychtol. 
% License: GPL version 2 or version 3
% $Id$

if ischar(mdl) && strcmp(mdl,'UNIT_TEST'); do_unit_test; return; end

copt.cache_obj = cache_obj(mdl, xyzr);
copt.fstr = 'mk_c2f_circ_mapping';
[mapping, failed] = eidors_cache(@circ_mapping,{mdl,xyzr},copt);

function [mapping, failed] = circ_mapping(mdl,xyzr,copt)

    failed = false;% not all subfunctions set this
    mdl = fix_model(mdl);
    switch size(xyzr,1)
      case 3; % use for 2D or for cylinder mapping in 3D
         mapping = contained_elems_2d( mdl, xyzr );
         if mdl_dim(mdl) == 2;
            correctmap = pi*xyzr(end,:).^2;
         else
            % No analytic way to calculate correct value (for non extruded models)
            correctmap = (get_elem_volume(mdl)'*mapping);
         end
      case 4; % use for 3D and spherical maps
         [mapping failed] = contained_elems_3d( mdl, xyzr );
         correctmap = 4/3*pi*xyzr(end,:).^3;
      otherwise; error('size of xyzr incorrect');
    end

    % Correct
    vol = get_elem_volume(mdl)';

    mapping = bsxfun(@times, mapping, correctmap./(vol*mapping));
    
% Mapping depends only on nodes and elems - remove the other stuff
function c_obj = cache_obj(mdl, xyzr)
   c_obj = {mdl.nodes, mdl.elems, xyzr};

% Redirector during test code dev
function mapping = contained_elems_2d( mdl, xyr );
%mapping = contained_elems_2d_new( mdl, xyr );
 mapping = contained_elems_2d_old( mdl, xyr );

function mapping = contained_elems_2d_new( mdl, xyr );
   % We fill sparse by columns, (ie adding in CCS storage)
   Nc = size(xyr,      2); % Num circs
   too_far = elems_too_far( mdl, xyr );

   mapping = sparse( num_elems(mdl) , Nc );
   for i=1:Nc
     mapping(:,i) = circ_in_elem_2d(mdl, find( ~too_far(:,i)), ...
                xyr(1,i), xyr(2,i), xyr(3,i));
   end

% look only at elements 'look'
function mapping = circ_in_elem_2d( mdl, look, xc, yc, rc);
   Nt = elem_dim(mdl) + 1; % nodes per simplex
   pirc2 = pi*rc^2;
% start assuming no content
   mapping = sparse(num_elems(mdl),1);
% For each element, find nodes inside
   els = mdl.elems(look,:);
   ndx = reshape(mdl.nodes(els,1) - xc, size(els));
   ndy = reshape(mdl.nodes(els,2) - yc, size(els));
   n_in = (ndx.^2 + ndy.^2) < rc^2; % node inside
% triangles with 3 nodes inside are all in, stop looking at them
   all_n_in = sum(n_in,2) == Nt;
   mapping(look(all_n_in)) = 1;
   look(all_n_in)= []; n_in(all_n_in,:)= [];
% find distance inside each face
   f_in = zeros( length(look), Nt); 
   k=1;  for i= look(:)';
      faces = mdl.elem2face(i,:);
      out   =~mdl.inner_normal(i,:);
      f_norm= mdl.normals( faces, :);
      f_norm(out,:) = -f_norm(out,:);

      f_ctr = mdl.face_centre( faces,:);
      v_ctr = repmat([xc,yc],Nt,1) - f_ctr;
      v_ctr = sum(v_ctr .* f_norm,2)/rc; % >1 in, <-1 out, 
      f_in(k,:) = v_ctr';
   k=k+1;end

% triangles with any sides outside are out
   any_s_out= any(f_in<-1,2);
   look(any_s_out)= [];
   n_in(any_s_out,:) = [];
   f_in(any_s_out,:) = [];

% triangles with 3 sides inside are all in.
   all_s_in = sum(f_in>1,2) == Nt;
   mapping(look(all_s_in)) = pirc2 / ...
            mdl.elem_volume(look(all_s_in));
   look(all_s_in)= [];
   n_in(all_s_in) = [];
   f_in(all_s_in) = [];

% Now, all triangles should be partially in
% Calculate area chopped out 
   fin1 = f_in<1;
   a_out = zeros(size(fin1));
   a_out(fin1) = acos(f_in(fin1));
   a_out(fin1) = (a_out(fin1) - cos(a_out(fin1)).*sin(a_out(fin1)))/pi;
   
   % start with the default. This is accurate if there are
   % no contained nodes, otherwise we need to add back
   % those fractions
   mapping(look) = pirc2 / mdl.elem_volume(look);

   %TODO: rewrite loop to avoid case 0.
   k=1; for i= look(:)';
      vol = pi*rc^2 / mdl.elem_volume(i);
      switch sum(n_in(k,:))
         case 0; % already do this
   
         case 1; 
            nd = mdl.elems(k, n_in(k,:));
            vol = vol + pi_slice(p1,p2,[xc,yc],mdl.nodes(nd,:),rc);

         case 2; 
            nd = mdl.elems(k, n_in(k,:));
            vol = vol ...
             + pi_slice(p1,p2,[xc,yc],mdl.nodes(nd(1),:),rc) ...
             + pi_slice(p1,p2,[xc,yc],mdl.nodes(nd(2),:),rc);

         otherwise; error('cant get here'); 
      end
   k=k+1; end
   

% Calculate the area of a slice
function a = pi_slice(p1,p2,c,p,r)
  a_p12 = 0.5*abs(det([1,p1;1,p2;1,p]));

  a_c12 = 0.5*abs(det([1,p1;1,p2;1,c]));
  np1c  = p1-c; np1c = np1c / norm(np1c);
  np2c  = p2-c; np2c = np2c / norm(np2c);
  ang   = acos( dot(np1c,np2c) );
  area  = ang*r^2/2 - a_c12 + a_p12;
     

function mapping = contained_elems_2d_old( mdl, xyr );
   Ne = size(mdl.elems,1); % Num elems
   Nc = size(xyr,      2); % Num circs
   % We fill sparse by columns, due to CCS storage, this is fairly efficient
   mapping = sparse( Ne, Nc );

   % Interpolate
   n_interp =  7-size(mdl.nodes,2);
   m_pts = interp_mesh( mdl, n_interp); 
   for i=1:Nc
     xc = m_pts(:,1,:) - xyr(1,i);
     yc = m_pts(:,2,:) - xyr(2,i);
     inr= xc.^2 + yc.^2 < xyr(3,i)^2;
     frac= mean(inr,3);
     mapping(:,i) = frac;
   end

   % 1. Get furthest node in each element
   % 2. Get the max edge length of each elem
   % 3. Find elems that are too far
function too_far = elems_too_far( mdl, xyr );
   Ne = num_elems(mdl);
   Nc = size(xyr, 2); % Num circs
   Nt = elem_dim(mdl) + 1; % Elements per simplex
   if 0 % for biggish Nc, this is insane
       nodes = repmat(mdl.nodes,[1,1,Nc]);
       targets = repmat(xyr(1:mdl_dim(mdl),:),[1,1,num_nodes(mdl)]);
       targets = shiftdim(targets,2);
       dist = nodes - targets;
       dist = sqrt(sum(dist.^2,2));
       node_target_dist = squeeze(dist);
       furthest_elem_node_dist = node_target_dist(mdl.elems,:);
       furthest_elem_node_dist = reshape(furthest_elem_node_dist,Ne,Nt,Nc);
       [furthest_elem_node_dist, furthest_elem_node]= max(furthest_elem_node_dist,[],2);
       furthest_elem_node_dist = squeeze(furthest_elem_node_dist);
       furthest_elem_node = squeeze(furthest_elem_node);
       max_edge_len =  repmat(mdl.max_edge_len,1,Nc);
       radius = ones(Ne,1)*xyr(Nt,:);
       too_far = (furthest_elem_node_dist - max_edge_len) > radius;
   else
       too_far = false(Ne,Nc);
       progress_msg('mk_c2f_circ_mapping: prepare models',0,Nc);
       for i = 1:Nc
          progress_msg(i,Nc);
          targets = repmat(xyr(1:mdl_dim(mdl),i),[1,num_nodes(mdl)])';
          dist = mdl.nodes - targets;
          dist = sqrt(sum(dist.^2,2));
          furthest_elem_node_dist =max(dist(mdl.elems),[],2);
          too_far(:,i) = (furthest_elem_node_dist - mdl.max_edge_len) > xyr(Nt,i);
       end
       progress_msg(Inf);
   end
  
   
  

function [mapping failed] = contained_elems_3d( mdl, xyr );
   Ne = size(mdl.elems,1); % Num elems
   Nc = size(xyr,      2); % Num circs
   failed(1:Nc) = false;
   % We fill sparse by columns, due to CCS storage, this is fairly efficient
   mapping = sparse( Ne, Nc );

       % 4. Make a tmp model with only the remaining elems
       % 5. Interpolate
       % 6. Merge
       
       too_far = elems_too_far( mdl, xyr );
       
       tmp = eidors_obj('fwd_model','tmp','nodes',mdl.nodes,'elems',mdl.elems);
       %mapping = sparse( Ne, Nc );
       try   
           n_interp_min = mdl.interp_mesh.n_points;
       catch
           n_interp_min = 4;
       end
       n_interp_max = 10;

   if 0
       % INterpolate
       n_interp = 4; % 7-df
       m_pts = interp_mesh( mdl, n_interp); 
       for i=1:Nc
         mapping(:,i) = contained_elem_pts(m_pts, xyr(:,i));
       end
   else

       progress_msg('mk_c2f_circ_mapping: calculate mapping',0,Nc);
       for i=1:Nc
           progress_msg(i,Nc);
           good = ~too_far(:,i);
           if ~any(good), continue, end %point outside the mesh
           tmp.elems = mdl.elems(good,:);
           n_interp = n_interp_min-1;
           log_level = eidors_msg('log_level',1);
           while(sum(mapping(good,i))==0 && n_interp < n_interp_max-1)
               n_interp = n_interp+1;
               m_pts = interp_mesh( tmp, n_interp);
               mapping(good,i) = contained_elem_pts(m_pts, xyr(:,i));
           end
           eidors_msg('log_level', log_level);
           if (sum(mapping(good,i)) == 0)
               failed(i) = true;
               eidors_msg(['mk_c2f_circ_mapping: Interpolation failed for point ' num2str(i)]);
           end
       end
       progress_msg(Inf);
   end
   
   
function frac= contained_elem_pts(m_pts, xyr);
% This is more clear
%    xc = m_pts(:,1,:) - xyr(1);
%    yc = m_pts(:,2,:) - xyr(2);
%    zc = m_pts(:,3,:) - xyr(3);
%    inr= xc.^2 + yc.^2 + zc.^2 < xyr(4)^2;

% But this is how to stop matlab from wasting memory
     inr = (m_pts(:,1,:) - xyr(1)).^2 + ...% xc =
           (m_pts(:,2,:) - xyr(2)).^2 + ...% yc =
           (m_pts(:,3,:) - xyr(3)).^2;     % zc =
     inpts = inr < xyr(4)^2;

     frac= mean( inpts ,3);
%    FIXME: Octave doesn't like to mean on logical
%    frac= mean( int8( inpts ) ,3);
     if sum(inpts(:))==0
         % TODO: This message is outdated
         eidors_msg(['mk_c2f_circ_mapping: Interpolation failed: increase ', ...
                         'fwd_model.interp_mesh.n_interp']);
     end

function do_unit_test
   fmdl = ng_mk_cyl_models([0,1,.1],[16,1],.03);
   vol = get_elem_volume(fmdl)';
   xyr = ones(3,1)*linspace(-.5,.5,7);
   rr = .1; VV = pi*rr^2; xyr(3,:) = rr;
   [c2f,fail] = mk_c2f_circ_mapping(fmdl,xyr);
   unit_test_cmp('2D #1.1:',sum(fail),0);
   unit_test_cmp('2D #1.2(r=.1):',vol*c2f/VV,1,1e-2);

   fmdl = ng_mk_cyl_models([2,1,.1],[16,1],.03);
   fmdl.nodes(:,3) = fmdl.nodes(:,3) - 1;
   vol = get_elem_volume(fmdl)';
   xyzr = ones(4,1)*linspace(-.5,.5,7);

   rr = .1; VV = pi*4/3*rr^3; xyzr(4,:) = rr;
   [c2f,fail] = mk_c2f_circ_mapping(fmdl,xyzr);
   unit_test_cmp('3D #1.1:',sum(fail),0);
   unit_test_cmp('3D #1.2(r=.1):',vol*c2f/VV,1,1e-2);

   rr = .01; VV = pi*4/3*rr^3; xyzr(4,:) = rr;
   [c2f,fail] = mk_c2f_circ_mapping(fmdl,xyzr);
   unit_test_cmp('3D #1.1(r=.01):',sum(fail),0);
   unit_test_cmp('3D #1.2:',vol*c2f/VV,1,1e-2);

   %2D example
   imdl = mk_common_model('a2c2',16); fmdl=imdl.fwd_model;
   xyc = [0,0.27,0.18;0,-0.1,0.03;0,0.1,0.2;0.1,0.37,0.1]';
   th=linspace(0,2*pi,20)';
   xx=[0*th+1]*xyc(1,:)+sin(th)*xyc(3,:);
   yy=[0*th+1]*xyc(2,:)+cos(th)*xyc(3,:);
   show_fem(fmdl,[0,0,1]); set(line(xx,yy),'LineWidth',2);

   % split over four elements
   rr= 0.1;c2f= mk_c2f_circ_mapping( fmdl, [0;0;rr] );
   tt= zeros(size(c2f)); tt(1:4) = pi*rr^2/4; tt= tt./get_elem_volume(fmdl);
   unit_test_cmp('2D ex 1:',c2f,tt,1e-10);

   % all in element #1
   rr= 0.03;c2f= mk_c2f_circ_mapping( fmdl, [.0;.05;rr]); 
   tt= zeros(size(c2f)); tt(1) = pi*rr^2; tt= tt./get_elem_volume(fmdl);
   unit_test_cmp('2D ex 2:',c2f,tt,1e-10);
      
   %3D example - cylinder
   imdl = mk_common_model('a3cr',16); fmdl=imdl.fwd_model;
   fmdl.nodes = 1.1*fmdl.nodes;
   rr=0.1;c2f= mk_c2f_circ_mapping( fmdl, [0;0;rr]); 
   V = pi*rr^2*(max(fmdl.nodes(:,3)) - min(fmdl.nodes(:,3)));
   unit_test_cmp('3D ex 1 (cylinder):',get_elem_volume(fmdl)'*c2f,V,1e-2);

   %3D example - sphere
   imdl = mk_common_model('a3cr',16); fmdl=imdl.fwd_model;
   rr=0.05;c2f= mk_c2f_circ_mapping( fmdl, [0;0;0;rr]); 
   tt = 4/3*pi*rr^3/24./get_elem_volume(fmdl);
   unit_test_cmp('3D ex 2a:',c2f(193:196),tt(193:196),1e-10);
   unit_test_cmp('3D ex 2b:',c2f(1:64),0);

   imdl = mk_common_model('a3cr',16); fmdl=imdl.fwd_model;
   rr=0.05;c2f= mk_c2f_circ_mapping( fmdl, [0 0;0 0;0 0;rr,rr]); 
   unit_test_cmp('3D ex 3a:',c2f(193:196,:),tt(193:196)*[1,1],1e-10);
   unit_test_cmp('3D ex 3b:',c2f(1:64,:),0);
