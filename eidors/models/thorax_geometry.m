function [out1, out2, out3 ] = thorax_geometry(in1,in2);
% THORAX_GEOMETRY: deform mesh to have a human thorax like shape
% Definition of thorax geometry for adult male (from CT)
% at five levels (values in mm)
%
% Calling: 
%   [x_coord, y_coord, z_mag ] = thorax_geometry(level,normalize);
% OR
%   fwd_model_new = thorax_geometry( fwd_model, level)
% where level is 1..5 for T1 .. T5

% (C) 2009 Andy Adler. License: GPL version 2 or version 3
% $Id$

call_thorax_geom = 1;
try if strcmp(in1.type, 'fwd_model')
   call_thorax_geom = 0;
end; end

if call_thorax_geom
   [out1, out2, out3 ] = thorax_geometry_defs(in1,in2);
else
   [x_coord, y_coord, z_mag ] = thorax_geometry_defs;
   reidx= [13:16, 1:12];
   geo.z_mag = z_mag;
   geo.xy= [x_coord(in2,reidx); y_coord(in2,reidx)]';
   
   [out1]              = deform_cylinder(in1, geo );
end


function [x_coord, y_coord, z_mag ] = thorax_geometry_defs(level,normalize);

    x_coord= [ ...
      0,60,123,173,223,202,144,75,0,-75,-144,-202,-223,-173,-123,-60;
      0,50,105,138,144,144,109,50,0,-50,-109,-144,-144,-138,-105,-50;
      0,52, 99,133,148,141,110,61,0,-61,-110,-141,-148,-133,- 99,-52;
      0,51, 92,129,148,136, 96,47,0,-47,- 96,-136,-148,-129,- 92,-51;
      0,49, 92,128,148,141,111,64,0,-64,-111,-141,-148,-128,- 92,-49 ];
    y_coord= [ ...
      123,116, 91,42,-4,-67,-105,-119,-108,-119,-105,-67,-4,42, 91,116;
      129,132,112,62, 3,-57,-101,-110,-107,-110,-101,-57, 3,62,112,132;
      116,112, 92,53, 3,-48,- 88,-106,-105,-106,- 88,-48, 3,53, 92,112;
      143,130, 99,63,14,-35,- 68,- 82,- 82,- 82,- 68,-35,14,63, 99,130;
      136,128,103,68,23,-25,- 62,- 78,- 80,- 78,- 62,-25,23,68,103,128 ];
    z_mag = 150;

    if nargin>=1;
      x_coord= x_coord(level,:);
      y_coord= y_coord(level,:);
    end

    if nargin>=2 && normalize
      maxv= max([ abs(x_coord(:)); ...
                  abs(y_coord(:)) ]);
      x_coord = x_coord / maxv;
      y_coord = y_coord / maxv;
    end
