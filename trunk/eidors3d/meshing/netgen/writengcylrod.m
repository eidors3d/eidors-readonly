function writengcylrod(fid,name,c, d,rd,ln)
% writes the specification for a netgen cylindrical rod on fid, named name, centerd on c,
% in the direction given by vector d, radius rd  lenght ln
% direction is in the xy plane
 % the direction vector
dirn = d.*(ln/(2*norm(d)));   % normalize
inpt = c - dirn.*(ln/2);
outpt =c + dirn.*(ln/2);


fprintf(fid,'solid %s  =plane(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f  )\n ',name , inpt(1),inpt(2),inpt(3),-dirn(1),-dirn(2),-dirn(3));
fprintf(fid,'       and plane(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f  )\n ',outpt(1),outpt(2),outpt(3),dirn(1),dirn(2),dirn(3));
fprintf(fid,'        and cylinder(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f; %6.3f  );\n ', inpt(1),inpt(2),inpt(3),outpt(1),outpt(2),outpt(3), rd);








