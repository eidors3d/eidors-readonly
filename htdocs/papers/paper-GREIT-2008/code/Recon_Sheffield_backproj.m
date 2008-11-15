function GREIT_Sheffield_backproj( savename )
% Reconstruct GREIT images using Sheffield Backprojection algorithm
%
% (C) 2008 Andy Adler. Licensed under GPL v2 or v3
% Reconstruction Matrix is Licensed by University of Sheffield
%   for academic use. See http://eidors3d.sf.net/data_contrib/
% $Id$

   [RM,map] = calc_backproj_matrix;
   normalize_flag = 1;
   save(savename, 'RM','normalize_flag');

function [RM,map] = calc_backproj_matrix;
   [x,y]= meshgrid(1:16,1:16); % Take a slice
   ss1 = (y-x)>1 & (y-x)<15;
   sel1 = abs(x-y)>1 & abs(x-y)<15;
   [x,y]= meshgrid(-15.5:15.5,-15.5:15.5);
   ss2 = abs(x-y)<25 & abs(x+y)<25 & x<0 & y<0 & x>=y ;
   sel2 = abs(x-y)<25 & abs(x+y)<25;
 
   load Sheffield_Backproj_Matrix.mat
   BP  = zeros(16^2, 32^2);
   BP(ss1,ss2) = Sheffield_Backproj_Matrix;
   BP  = reshape(BP, 16,16,32,32); % Build up
   BP  = BP + permute(BP, [2,1,3,4]); % Reciprocity
% FLIP LR
   el= 16:-1:1;            BP= BP + BP(el,el,[32:-1:1],:);
% FLIP UD
   el= [8:-1:1,16:-1:9];   BP= BP + BP(el,el,:,[32:-1:1]);
% Transpose
   el= [12:-1:1,16:-1:13]; BP= BP + permute(BP(el,el,:,:), [1,2,4,3]);
% Final flip to match radiological view (upward toward patient)
   BP = permute(BP, [1,2,4,3]);

   RM= reshape(BP, 256, [])';
   RM= RM(:,sel1);

   map = sel2;
