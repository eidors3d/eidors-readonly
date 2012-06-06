function writevtkfile(fn,varargin)
% Useage: writevtkfile(filename, eidors_image );
%      or writevtkfile(filename,model,image)
%      or writevtkfile(fn,vtx,simp,vals);
%
% Example:  img= compare_3d_algs(1);
%           writevtkfile('fname',img);
%
% This function writes 3D meshes and associated scalar/vector/tensor fields to a VTK (Visualization Tool Kit)(www.kitware.com) dataset file
%
%
% 'filename' = name of the vtk file, the ".vtk" extension will be addedd automatically 
% vtx = matrix of the verticies
% elem = matrix desribing the simplecies
% val = matrix holding the field values (scalars/vectros/tensors)
% 
% Depending wether the field is piecewise linear or piecewise constant and if it it scalar, vector or tensor, val can assume different valid sizes:
%
% val should have a number of rows equal to the number of simplecies or to the number of nodes in the mesh (piecewise constant/piecewise linear dataset)
% val should have a number of columns equal to 1,3 or 9 (scalar field / vector field / tensor field)
%
% The generated VTK file can then be displayed with MaiaVi (http://mayavi.sourceforge.net/)
% which has lots of nice features and is generally faster than MATALB 3D graphics
%
% CAVEAT: VTK supports only symmetrix tensors. It is left to the user to ensure that the data he is exporting 
% consists of symmetric tensors (as if such a check were performed by writeVTKfile, performance of the filter would suffer)

if nargin==2
   img  = varargin{1};
   fmdl = img.fwd_model;
   writevtkfile_old(fn, fmdl.nodes, fmdl.elems, img.elem_data);
elseif nargin ==3
%vtx = varargin{1}.nodes;
%simp= varargin{1}.elems;
%vals= varargin{2}.elem_data;
   fmdl = varargin{1};
   img  = varargin{2};
   writevtkfile_old(fn,fmdl.nodes, fmdl.elems, img.elem_data);
elseif nargin==4
   writevtkfile_old(fn,varargin{1},varargin{2},varargin{3})
else
 error('Writevtkfile: Wrong number of arguments');
end 

function writevtkfile_old(filename,vtx,elem,val)
%
% Copyright (c) 2002-2003 Andrea Borsic, SC-AIP s.a.s.
% Scientific Computing & Applied Inverse Problems
% www.sc-aip.com
% Revision 1.2 Original Date 12/08/02
% Modified 11/06/03, exports vectors
% Modified 31/08/03, exports tensors, fixes some problems regarding number rapresentation
% Modified 26/03/03 in order to support both triangular and tetrahedral meshes
% 25/06/2005 Changed name, writevtkfile is now overloaded to handle objects as data WRBL 
% This function writes 3D meshes and associated scalar/vector/tensor fields to a VTK (Visualization Tool Kit)(www.kitware.com) dataset file
%
% Use: writevtkfile_old('filename',vtx,elem,val)
%
% 'filename' = name of the vtk file, the ".vtk" extension will be addedd automatically 
% vtx = matrix of the verticies
% elem = matrix desribing the simplecies
% val = matrix holding the field values (scalars/vectros/tensors)
% 
% Depending wether the field is piecewise linear or piecewise constant and if it it scalar, vector or tensor, val can assume different valid sizes:
%
% val should have a number of rows equal to the number of simplecies or to the number of nodes in the mesh (piecewise constant/piecewise linear dataset)
% val should have a number of columns equal to 1,3 or 9 (scalar field / vector field / tensor field)
%
% The generated VTK file can then be displayed with MaiaVi (http://mayavi.sourceforge.net/)
% which has lots of nice features and is generally faster than MATALB 3D graphics
%
% CAVEAT: VTK supports only symmetrix tensors. It is left to the user to ensure that the data he is exporting 
% consists of symmetric tensors (as if such a check were performed by writeVTKfile, performance of the filter would suffer)

% DA FARE: mettere check sul numero di parametri di ingresso (e sul tipo, es controllare che filename sia una stringa)

num_vtx=size(vtx,1);
num_elem=size(elem,1);

filename=strcat(filename,'.vtk');

fid=fopen(filename,'wt');
if fid==-1
    error('Cannot open vtk file');
end % if FID

frewind(fid);

fprintf(fid,'# vtk DataFile Version 3.0\n');
fprintf(fid,'%s %s exported from MATLAB\n',date,filename);
fprintf(fid,'ASCII\n');

fprintf(fid,'\n');

% Here we write the verticies

fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid,'POINTS %i double\n',num_vtx);

for i=1:num_vtx
    fprintf(fid,'%d %d %d\n',vtx(i,1),vtx(i,2),vtx(i,3));
end % for

fprintf(fid,'\n');

% Here we write the simplecies

switch size(elem,2)
    
    case 3 % 2D triangular mesh
        
        fprintf(fid,'CELLS %i %i\n',num_elem,4*num_elem);
        
        for i=1:num_elem
            fprintf(fid,'%i %i %i %i \n',3,elem(i,1)-1,elem(i,2)-1,elem(i,3)-1);
        end % for
        
        fprintf(fid,'\n');
        
        fprintf(fid,'CELL_TYPES %i\n',num_elem);
        
        for i=1:num_elem
            fprintf(fid,'5\n');
        end % for
        
    case 4 % 3D tetrahedral mesh
        
        fprintf(fid,'CELLS %i %i\n',num_elem,5*num_elem); % per ogni tetraedro 4 dati + 1 che specifica numpint=4 in totale 5 vedi pag 357
        
        for i=1:num_elem
            fprintf(fid,'%i %i %i %i %i\n',4,elem(i,1)-1,elem(i,2)-1,elem(i,3)-1,elem(i,4)-1);
        end % for
        
        fprintf(fid,'\n');
        
        fprintf(fid,'CELL_TYPES %i\n',num_elem);
        
        for i=1:num_elem
            fprintf(fid,'10\n');
        end % for
        
    otherwise
        
        error('Simplecies other than triangles and tetrahydra are not currently supported')
        
end

fprintf(fid,'\n');

switch size(val,1)
    
    case num_elem
        
        switch size(val,2)
            
            case 1 % scalar field, cell data
                
                fprintf(fid,'CELL_DATA %i\n',num_elem);
                fprintf(fid,'SCALARS MatlabExportedScalars double\n');
                fprintf(fid,'LOOKUP_TABLE default\n');
                
                for i=1:num_elem
                    fprintf(fid,'%d\n',val(i));
                end % for
                
            case 3 % vector field, cell data
                
                fprintf(fid,'CELL_DATA %i\n',num_elem);
                fprintf(fid,'VECTORS MatlabExportedVectors double\n');
                % No LUT for vectors !
                
                for i=1:num_elem
                    fprintf(fid,'%d %d %d\n',val(i,1),val(i,2),val(i,3));
                end % for
                
            case 9 % tensor field, cell data
                
                fprintf(fid,'CELL_DATA %i\n',num_elem);
                fprintf(fid,'TENSORS MatlabExportedTensors double\n');
                % No LUT for tensors !
                
                for i=1:num_elem
                    for j=1:3
                        fprintf(fid,'%d %d %d\n',val(i,(j-1)*3+1),val(i,(j-1)*3+2),val(i,(j-1)*3+3));
                    end % for j
                end % for
                
            otherwise
                
                fprintf('The number of columns in the val vector should be 1,3 or 9, depending if it reperesents a scalar, a vector field or a tensor field.\n');
                fprintf('A VTK file containing just mesh information was written.\n');     
                
        end % switch size(val,2)
        
    case num_vtx
        
        switch size(val,2)
            
            case 1 % scalar field, point data
                
                fprintf(fid,'POINT_DATA %i\n',num_vtx);
                fprintf(fid,'SCALARS MatlabExportedScalars double\n');
                fprintf(fid,'LOOKUP_TABLE default\n');
                
                for i=1:num_vtx
                    fprintf(fid,'%d\n',val(i));
                end % for
                
            case 3 % vector field, point data
                
                fprintf(fid,'POINT_DATA %i\n',num_vtx);
                fprintf(fid,'VECTORS MatlabExportedVectors double\n');
                % No LUT for vectors !
                
                for i=1:num_vtx
                    fprintf(fid,'%d %d %d\n',val(i,1),val(i,2),val(i,3));
                end % for
                
            case 9 % tensor field, point data
                
                fprintf(fid,'POINT_DATA %i\n',num_vtx);
                fprintf(fid,'TENSORS MatlabExportedTensors double\n');
                % No LUT for tensors !
                
                for i=1:num_vtx
                    for j=1:3
                        fprintf(fid,'%d %d %d\n',val(i,(j-1)*3+1),val(i,(j-1)*3+2),val(i,(j-1)*3+3));
                    end % for j
                end % for
                
            otherwise
                
                fprintf('The number of columns in the val vector should be 1,3 or 9, depending if it reperesents a scalar, a vector field or a tensor field.\n');
                fprintf('A VTK file containing just mesh information was written.\n');
                
        end % switch size(val,2)
        
    otherwise
        
        fprintf('The number of rows of the val vector did not match the number of rows in vtx or in elem.\n');
        fprintf('A VTK file containing just mesh information was written.\n');
        
end % switch size(val,1)

fclose(fid);

% Bye Bye
