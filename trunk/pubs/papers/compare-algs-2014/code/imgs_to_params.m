function imgs_to_params( animals, RMs )
% imgs_to_params: Convert images to params
% animals = vector of animals
% RMs = list of RM strings 

% TODO: fix ROIs to excude edge effects
   if nargin<1, do_unit_test; return; end

   if isstr(animals) && strcmp(animals,'UNIT_TEST'); do_unit_test; return; end

   fid = fopen('params_out.m','w');
   fprintf(fid,'%% Written at: %s\n', datestr(now));
   fclose(fid);

   for animal = animals(:)'
     for RM = RMs
       img_to_params(animal, RM{1});
     end
   end

function do_unit_test
global locn; locn=0;
   imgs_to_params([5:12],{ 'R0', 'R3D', 'Rbkg', 'Rdif', 'RGR', 'RL1', 'Rlpf', 'RTSVD' ,'Rmv', 'Rn', 'RSBP', 'RTV'});


function img_to_params(animal, RM)
   fprintf('PROCESSING animal=%02d RM=%s\n', animal, RM);
   [img, bdy] = load_animal( animal, RM);
   ROI = calc_ROI(img,bdy);
   calc_write_params( img, bdy, ROI, animal, RM);

global locn; locn=locn+1;
   subplot(15,15,locn);
%  subplot(5,5,locn);
   imagesc(ROI); title(sprintf('%02d:%s', animal,RM)); 
   axis off;

function [img, bdy] = load_animal( animal, RM);
  for i=[1:6,12:20]
    fname = sprintf('../recon/S%02d-%03d-%s.mat',animal,i,RM);
    load(fname,'raster');
    bdy = isnan(raster);
    raster(bdy) = 0;
    img(i).raster = - raster;  % convert non-conductive to + ventilation
  end

function ROI = calc_ROI(img,bdy);

  ROI = false;
  for i=1:6 % the 6 VT images
    imgi = img(i).raster;
    maxi = max(imgi(:));
    ROIi = imgi > maxi * 0.20;
    ROI = ROI | ROIi; % accumulate all
  end

% IDEA #1 with imfill
% ROI1 = imfill(~(ROI | bdy), 'holes');
% ROI = ROI & ROI1;

  bdy = conv2(double(bdy), ones(5),'same') > 0;
  ROI = ROI & ~bdy;

% OUTPUT LIKE:
% data( 5,  1).RSBP    =[### ]'; % VT,EIT,ZEEP,21,#1
function  calc_write_params( img, bdy, ROI, animal, RM);
   fid = fopen('params_out.m','a');
   for i = 1:length(img);
      if isempty(img(i).raster); continue;end
      iname =  define_img_no(i);
      imgi = img(i).raster;
      ampli = sum(sum( ROI.*imgi));
      imgi = abs(imgi);
      yaxis = linspace(+1,-1,size(imgi,1))' * ones(1,size(imgi,2));
      CoG =   sum(sum( ROI.*yaxis.*imgi))/sum(sum( ROI.*imgi));
      if isnan(CoG),keyboard,end;
      fprintf(fid,'data(%2d,%3d).%-8s=[ %+10.4f;%+6.4f ]; %% %s\n', ...
                  animal, i, RM, ampli, CoG, iname );

   end
   fclose(fid);
