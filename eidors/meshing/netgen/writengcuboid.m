function writengcuboid(fid,name,c, dirn,h,w,d)
% writes the specification for a netgen cuboid on fid, named name, centerd on c,
% in the direction given by vector dirn, height  h width w and depth d
% direction is in the xy plane
dirnp = [-dirn(2),dirn(1),0];
bl = c - (d/2)* dirn + (w/2)* dirnp -[0,0,h/2];
tr =c + (d/2)* dirn - (w/2)* dirnp +[0,0,h/2];
fprintf(fid,'solid %s  =plane (%6.3f,%6.3f,%6.3f;0, 0, -1 )\n ',name,bl(1),bl(2),bl(3));
fprintf(fid,'        and plane(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f  )\n ', bl(1),bl(2),bl(3),-dirn(1),-dirn(2),0);
fprintf(fid,'        and plane(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f  )\n ', bl(1),bl(2),bl(3),dirnp(1),dirnp(2),0);
fprintf(fid,'        and plane(%6.3f,%6.3f,%6.3f;0, 0, 1  )\n ',tr(1),tr(2),tr(3));
fprintf(fid,'        and plane(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f  )\n ', tr(1),tr(2),tr(3),dirn(1),dirn(2),0);
fprintf(fid,'        and plane(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f  );\n ', tr(1),tr(2),tr(3),-dirnp(1),-dirnp(2),0);

