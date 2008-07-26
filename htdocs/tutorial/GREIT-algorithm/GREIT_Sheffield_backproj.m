function [img,map]= GREIT_Sheffield_backproj( ref_meas, reconst_meas )

   [RM,map] = calc_backproj_matrix;
keyboard

   % Expand ref_meas to the full size of reconst_meas
   num_meas = size(reconst_meas,2);
   ref_meas = ref_meas * ones(1,num_meas);
   dv = ( reconst_meas - ref_meas ) ./ ref_meas; % CHANGE IS HERE:

   % reconst image
   ds = RM*dv;

   img= reshape(ds, 32,32,num_meas);

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
% Final UD flip to match radiological view (upward toward patient)
% Here electrodes are connected CW starting from TDC
   BP= BP(:,:,:,[32:-1:1]);

   RM= reshape(BP, 256, [])';
   RM= RM(:,sel1);

   map = sel2;
