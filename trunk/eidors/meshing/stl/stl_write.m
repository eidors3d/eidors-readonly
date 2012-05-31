function stl_write(fv, name)
%STL_WRITE Create a text STL file from a patch struct
% stl_write(fv, name) where:
%  fv is face-vertex structure array (fv.faces, fv.vertices)
%  name is the file name (character string), if no extension give, 'stl'
%  assumed.

% (C) 2006 Eric Carlson
% Adapted by Bartlomiej Grychtol from:
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/126215
% $Id$

if (strcmp(name((end-3):end), '.stl'))
     label = name(1:(end-4));
else
     label = name;
     name = sprintf('%s.stl', name);
end
v1 = fv.vertices(fv.faces(:,2),:)-fv.vertices(fv.faces(:,1),:);
v2 = fv.vertices(fv.faces(:,3),:)-fv.vertices(fv.faces(:,2),:);
Norms = cross3(v1,v2);
clear v1 v2
v1(:,1:3) = fv.vertices(fv.faces(:,1),1:3);
v2(:,1:3) = fv.vertices(fv.faces(:,2),1:3);
v3(:,1:3) = fv.vertices(fv.faces(:,3),1:3);
fid = fopen(name,'w');
fprintf(fid,'solid %s\n',label);
nf = length(fv.faces); %k = (1:nf)';
for k = 1:nf
    fprintf(fid,'facet normal %5.5f %5.5f %5.5f\n outer loop\n vertex %5.5f %5.5f %5.5f\n vertex %5.5f %5.5f %5.5f\n vertex %5.5f %5.5f %5.5f\nendloop\n endfacet\n', ...
        Norms(k,1),Norms(k,2),Norms(k,3), v1(k,1), v1(k,2), v1(k,3),v2(k,1), v2(k,2), v2(k,3),v3(k,1), v3(k,2), v3(k,3) );
end
fprintf(fid,'endsolid %s\n',label);
fclose(fid);

function M=cross3(r,F)
% function to calculate normalized cross product rxF/|rxF|
% handles (same-size) arrays (n by 3) for r and F
%
       M = [(r(:,2).*F(:,3) - r(:,3).*F(:,2)) ...
            (r(:,3).*F(:,1) - r(:,1).*F(:,3)) ...
            (r(:,1).*F(:,2) - r(:,2).*F(:,1))];
       M_mag = sqrt(sum((M.*M)')');
       M(:,1) = M(:,1)./M_mag;
       M(:,2) = M(:,2)./M_mag;
       M(:,3) = M(:,3)./M_mag;
