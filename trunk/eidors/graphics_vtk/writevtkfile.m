function writevtkfile(fn,varargin)
% Useage: writevtkfile(filename,model,image)
% or writevtkfile(fn,vtx,simp,vals);
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

if nargin ==3
%vtx = varargin{1}.nodes;
%simp= varargin{1}.elems;
%vals= varargin{2}.elem_data;
writevtkfile_old(fn,varargin{1}.nodes,varargin{1}.elems,varargin{2}.elem_data);
elseif nargin==4
 writevtkfile_old(fn,varargin{1},varargin{2},varargin{3})
else
 error('Writevtkfile: Wrong number of arguments');
end 