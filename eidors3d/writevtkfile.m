function writevtkfile(filename,vtx,tet,val)
%
% Copyright 2003 Andrea Borsic, SC-AIP s.a.s.
% Scientific Computing & Applied Inverse Problems
% www.sc-aip.com
% Released under the Gnu GLP as part of the EIDORS project www.eidors.org
%
% Revision 1.1 Original Date 12/08/02, Modified 11/06/03
%
% This function writes 3D meshes and associated scalar/vector fields to a VTK (Visualization Tool Kit)(www.kitware.com) dataset file
%
% Use: writevtkfile('filename',vtx,tet,val)
%
% 'filename' = name of the vtk file, the ".vtk" extension will be addedd automatically 
% vtx = matrix of the verticies, as in EIDORS 3D
% tet = matrix of the tatrahydra, as simp in EIDORS3D
% val = vector holding the data values
% 
% Depending wether the field is piecewise linear or piecewise constant and if it it scalar or vector, val can assume different valid sizes:
%
% val should have a number of rows equal to the number of tetrahydra or to the number of nodes in the mesh (piecewise constant/piecewise linear dataset)
% val should have a number of columns equal to 1 or 3 (scalar field / vector field)
%
% The the generated VTK file can then be displayed with MaiaVi (http://mayavi.sourceforge.net/)
% which has lots of nice features and is fast.

num_vtx=size(vtx,1);
num_tet=size(tet,1);

filename=strcat(filename,'.vtk');

fid=fopen(filename,'wt');
if fid==-1
    error('Cannot open vtk file');
end % if FID

frewind(fid);

fprintf(fid,'# vtk DataFile Version 3.0\n');
fprintf(fid,'%s Exported from MATLAB\n',date);
fprintf(fid,'ASCII\n');

fprintf(fid,'\n');

fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid,'POINTS %i float\n',num_vtx);

for i=1:num_vtx
    fprintf(fid,'%f %f %f\n',vtx(i,1),vtx(i,2),vtx(i,3));
end % for

fprintf(fid,'\n');

fprintf(fid,'CELLS %i %i\n',num_tet,5*num_tet); % per ogni tetraedro 4 dati + 1 che specifica numpint=4 in totale 5 vedi pag 357

for i=1:num_tet
    fprintf(fid,'%i %i %i %i %i\n',4,tet(i,1)-1,tet(i,2)-1,tet(i,3)-1,tet(i,4)-1);
end % for

fprintf(fid,'\n');

fprintf(fid,'CELL_TYPES %i\n',num_tet);

for i=1:num_tet
    fprintf(fid,'10\n');
end % for

fprintf(fid,'\n');

switch size(val,1)
    
    case num_tet
        
        switch size(val,2)
            
            case 1 % scalar field
                
                fprintf(fid,'CELL_DATA %i\n',num_tet);
                fprintf(fid,'SCALARS Scalars float\n');
                fprintf(fid,'LOOKUP_TABLE default\n');
                
                for i=1:num_tet
                    fprintf(fid,'%f\n',val(i));
                end % for
                
            case 3 % vector field
                
                fprintf(fid,'CELL_DATA %i\n',num_tet);
                fprintf(fid,'VECTORS Vectors float\n');
                % No LUT for vectors !
                
                for i=1:num_tet
                    fprintf(fid,'%f %f %f\n',val(i,1),val(i,2),val(i,3));
                end % for
                
            otherwise
                
                fprintf('The number of columns in the val vector should be 1 or 3, depending if it reperesents a scalar or vector field.\n');
                fprintf('A VTK file containing just mesh information was written.\n');     
                
        end % switch size(val,2)
        
    case num_vtx
        
        switch size(val,2)
            
            case 1
                
                fprintf(fid,'POINT_DATA %i\n',num_vtx);
                fprintf(fid,'SCALARS Scalars float\n');
                fprintf(fid,'LOOKUP_TABLE default\n');
                
                for i=1:num_vtx
                    fprintf(fid,'%f\n',val(i));
                end % for
                
            case 3
                
                fprintf(fid,'POINT_DATA %i\n',num_vtx);
                fprintf(fid,'VECTORS Vectors float\n');
                % No LUT for vectors !
                
                for i=1:num_vtx
                    fprintf(fid,'%f %f %f\n',val(i,1),val(i,2),val(i,3));
                end % for
                
            otherwise
                
                fprintf('The number of columns in the val vector should be 1 or 3, depending if it reperesents a scalar or vector field.\n');
                fprintf('A VTK file containing just mesh information was written.\n');
                
        end % switch size(val,2)
        
    otherwise
        
        fprintf('The number of rows of the val vector did not match the number of rows in vtx or in tet.\n');
        fprintf('A VTK file containing just mesh information was written.\n');
        
end % switch size(val,1)

fclose(fid);

% Bye Bye
