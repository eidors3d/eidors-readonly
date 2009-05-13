function mapping = mk_c2f_circ_mapping( mdl, xyzr );
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

% (C) 2009 Andy Adler. License: GPL version 2 or version 3
% $Id$

c_obj = cache_obj(mdl, xyzr);

mapping = eidors_obj('get-cache', c_obj, 'circle_mapping');
if ~isempty(mapping)
    eidors_msg('mk_c2f_circ_mapping: using cached value', 3);
else

    switch size(xyzr,1)
      case 3; mapping = contained_elems_2d( mdl, xyzr );
      case 4; mapping = contained_elems_3d( mdl, xyzr );
      case 5: error('size of xyzr incorrect');
    end

    eidors_obj('set-cache', c_obj, 'circle_mapping', mapping);
    eidors_msg('mk_coarse_fine_mapping: setting cached value', 3);
end

% Mapping depends only on nodes and elems - remove the other stuff
function c_obj = cache_obj(mdl, xyzr)
   c_obj = {mdl.nodes, mdl.elems, xyzr};

function mapping = contained_elems_2d( mdl, xyr );
   Ne = size(mdl.elems,1); % Num elems
   Nc = size(xyr,      2); % Num circs
   % We fill sparse by columns, due to CCS storage, this is fairly efficient
   mapping = sparse( Ne, Nc );

   % INterpolate
   n_interp = 4; % 7-df
   m_pts = interp_mesh( mdl, n_interp); 
   for i=1:Nc
     xc = m_pts(:,1,:) - xyr(1,i);
     yc = m_pts(:,2,:) - xyr(2,i);
     inr= xc.^2 + yc.^2 < xyr(3,i)^2;
     frac= mean(inr,3);
     mapping(:,i) = frac;
   end

function mapping = contained_elems_3d( mdl, xyr );
   Ne = size(mdl.elems,1); % Num elems
   Nc = size(xyr,      2); % Num circs
   % We fill sparse by columns, due to CCS storage, this is fairly efficient
   mapping = sparse( Ne, Nc );

   % INterpolate
   n_interp = 3; % 7-df
   m_pts = interp_mesh( mdl, n_interp); 
   for i=1:Nc
     xc = m_pts(:,1,:) - xyr(1,i);
     yc = m_pts(:,2,:) - xyr(2,i);
     zc = m_pts(:,3,:) - xyr(3,i);
     inr= xc.^2 + yc.^2 + zc.^2 < xyr(4,i)^2;
     frac= mean(inr,3);
     mapping(:,i) = frac;
   end

