function cylinder_mdl = create_cylinder_mesh_ng( basename, ...
            cylinder_radius, cylinder_height, maxh )
%USAGE: cylinder_mdl = create_cylinder_mesh_ng( basename ...
%                           cylinder_radius, cylinder_height, maxh )
%            
% Parameters:
%      - basename : basename for filenames,
%      - cylinder_radius, cylinder_height, 
%      - maxh : maximum element size,
%
% Function creates coarse homogenous cylinder mesh without electrodes, 
% suitable for dual mesh reconstruction 3d-3d.
%
% (C) 2009,  Bartosz Sawicki
% $Id$
% Licenced under the GPLv2 or later

% maxh not given
if nargin < 4
   maxh = 10000;
end

geofn= [basename,'.geo'];
meshfn= [basename,'.vol'];

fid=fopen(geofn,'w');
fprintf(fid,'#Automatically generated cylinder mesh\n');
fprintf(fid,'algebraic3d\n');
fprintf(fid,'solid cuts = plane(0,0,0;0,0,-1) ');
fprintf(fid,'   and plane(0,0,%6.2f;0,0,1);\n', cylinder_height);
fprintf(fid,'solid cyl=cylinder (0,0,0;0,0,1;%6.2f) and cuts -maxh = %6.2f; \n', ...
           cylinder_radius, maxh);
fprintf(fid,'tlo cyl;\n');
fclose(fid);

disp('Calling Netgen. Please wait.....');
call_netgen( geofn, meshfn);

disp(['Now reading back data from file: ' meshfn])        
[srf,vtx,fc,bc,simp,edg,mat_ind] = ng_read_mesh(meshfn);

mdl.nodes = vtx;
mdl.elems = simp;
mdl.boundary = srf;
mdl.name = 'Simple cylinder';
cylinder_mdl= eidors_obj('fwd_model', mdl);         

