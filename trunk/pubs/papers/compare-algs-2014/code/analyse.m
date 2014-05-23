function analyse(opts) % PEEP-FIo2
% Analyse PEEP - FIO2 data from Inez Frerichs/ Uni Kiel
% usage: analyse(opts)
%    opts: generate end insp and end exp for each breath
% $Id$

% Create cropped images;
% for i in S*.png ; do convert -crop '63x247+0-184' $i C-$i ; done
% rm C-*_sig.png
% montage C-S12-*.png -geometry +0+0 -tile 4x6 out.png

calc_colours('mapped_colour',125);
calc_colours('backgnd',.9*[1,1,1]);
calc_colours('greylev',-.1);
if nargin==0; opts = 0; end
switch opts;
  case 0; analyse_exps(analyse_dirs);
  case 1; % Taxonomy figure
  case 2; mk_figure2;
  case 3; mk_figure3;

  case 99; mk_drift_figure;

  otherwise; error('huh?');
end

function printncrop(fname);
     fname = ['../paper/',fname];
     print('-dpdf',fname);
     system(['LD_LIBRARY_PATH="" pdfcrop ',fname,' ',fname]);

function mk_figure2;
% Finite element models used. Electrode nodes are indicated in
% green.  A: 2D circular uniform FEM (R0) B: 3D cylindrical
% uniform FEM (R3D) C: 3D cylindrical FEM with lung regions (Rbkg)
% D: 2D circular uniform FEM with electrode movement (Rmv)

   mdl1 = select_inv_model('paper-2D');
   mdl2 = select_inv_model('paper-3D');
   mdl3 = select_inv_model('paper-3D-lungs');
   lungs = mdl3.mat_idx{2};
   mdl3 = mk_image(mdl3,1); mdl3.elem_data(lungs) = 0.9;
   mdl3.calc_colours.ref_level=1;
   mdl3.calc_colours.clim=.3;
   mdl4 = mk_image(mdl1,1);
   mdl4.calc_colours.ref_level=1;
   th= linspace(0,2*pi,17)'; th(end) = [];
   mdl4.elem_data = [mdl4.elem_data;.01*sin(th);-.01*cos(th)];
   LW = 0.2;
subplot(1,4,1);
    hh = show_fem(mdl1);
    set(hh,'LineWidth',LW);
    axis image;
    text(-.9,-.94,'A');
    axis off
subplot(1,4,2);
    hh=show_fem(mdl2);  view(0,50);
    set(gca,'YTickLabel',[]);
    set(hh,'LineWidth',LW,'FaceColor','w');
    axis image;
    text(-.9,-.9,.1,'B');
    axis off
subplot(1,4,3);
    hh = show_fem(mdl3);  view(0,50);
    set(hh,'LineWidth',LW,'FaceColor','w');
    set(gca,'XTick',[-.5,0,.5]);
    set(gca,'YTickLabel',[]);
    text(-.9,-.9,.1,'C');
    axis off
%   axis equal;
subplot(1,4,4);
    hh = show_fem_move(mdl4);
    set(hh,'LineWidth',LW,'FaceColor','w');
    set(gca,'XTick',[-.5,0,.5]);
    set(gca,'YTickLabel',[]);
    axis image; axis([-1,1,-1,1]); 
        text(-.9,-.94,'D');
    axis off
%   axis equal;
   printncrop('fig-FEMs.pdf');

function mk_figure3;
calc_colours('mapped_colour',125);
calc_colours('backgnd',.9*[1,1,1]);
calc_colours('greylev',.01);
   [imdl_v, ROIv] = select_inv_model(1000, 'S12', []);

   FRATE= 13;
   fname= 'S12/S12-006.get';
   [d_insp,d_expi, W] = get_vent_meas( fname, imdl_v, ROIv, FRATE, []);

   %TODO: Comment out line 93 text in find_frc.m
   set(gca,'YTickLabel',[]);
   set(gca,'XTickLabel',{'0','10','20','30','40','50',''});
   printncrop('fig-EndExpiInspi.pdf');

function mk_drift_figure
   [imdl_v, ROIv] = select_inv_model(1, 'S12', []); FRATE= 13;
   for diri = analyse_dirs; dir= diri{1};
      dd1 = eidors_readdata(sprintf('%s/%s-001.get',dir,dir));
      dd2 = eidors_readdata(sprintf('%s/%s-002.get',dir,dir));
      dd3 = eidors_readdata(sprintf('%s/%s-003.get',dir,dir));
      dd4 = eidors_readdata(sprintf('%s/%s-004.get',dir,dir));
      dd5 = eidors_readdata(sprintf('%s/%s-005.get',dir,dir));
      dd6 = eidors_readdata(sprintf('%s/%s-006.get',dir,dir));
      ii = inv_solve( imdl_v, mean(dd1,2), [dd1,dd2,dd3,dd4,dd5,dd6]);   
      [einsp,eexpi] = find_frc( ii, ROIv, FRATE, ['/',dir], 1, []);   
      axis tight
      print_convert(['Drift-',dir,'.png']);
   end


   return
   figure(3)
   [di061,de061, W] = get_vent_meas( 'S06/S06-001.get', imdl_v, ROIv, FRATE, []);
   [di062,de062, W] = get_vent_meas( 'S06/S06-002.get', imdl_v, ROIv, FRATE, []);
   [di071,de071, W] = get_vent_meas( 'S07/S07-001.get', imdl_v, ROIv, FRATE, []);
   [di072,de072, W] = get_vent_meas( 'S07/S07-002.get', imdl_v, ROIv, FRATE, []);
   show_slices(inv_solve(imdl_v, ...
        [de061, de062, de061, de071, de072, de071, di071], ...
        [di061, di062, de062, di071, di072, de072, di072]))



function dirs = analyse_dirs;
 dirs= { 'S05', 'S06','S07', 'S08', 'S09', 'S10', 'S11', 'S12'}

function analyse_exps(analyse_dirs)
calc_colours('mapped_colour',125);
calc_colours('backgnd',.9*[1,1,1]);
calc_colours('greylev',.01);
[fid,fidb,fidi] = openoutfiles;

% ! find -name \*.png -exec convert -crop 64x64+0+0 {} {} \;

for dir = analyse_dirs
   dir= dir{1}; %Why is matlab stupid like this
   fprintf('ANALYSING[ %s ]\n', dir);
   htmlnewdir(fidb,dir); htmlnewdir(fidi,dir);

   [expis, insps, W] = get_insp_expi_data( dir, fidb);
%    continue
   for i=[1:6,12:20]
      [iname, Dt] =  define_img_no(i,expis,insps);

      imfrac= sprintf('%s-0%02d-',dir,i);

      imname={};
    for imdl_no = [1000:1011] % DOMODELS
         [imdl, ROI, mdlname] = select_inv_model(imdl_no,dir,W );
         if imdl_no == 1005 && i == 1
             write_icov(imdl,dir);
%            continue
         end
         [this_im_name, raster] = recon_write(imdl, Dt, imfrac, mdlname);
         imname= {imname{:}, this_im_name};

         write_analyse(raster, imfrac, mdlname, iname);
      end
     
      htmlimage(fidi,imname);
   end
   
end
htmltail(fidb)
htmltail(fidi)
system('find -name "*_sig.png" -exec convert -trim "{}" PNG8:"{}" ";"');

function write_icov(imdl, dir)
fid = fopen('icov_vals.txt','a');
fprintf(fid,'%s', dir);
fprintf(fid,', %f', full(diag(imdl.meas_icov))');
fprintf(fid,'\n');
fclose(fid);

function [expis, insps, W] = get_insp_expi_data( dir, fidb)
   [imdl_v, ROIv] = select_inv_model(1009, dir, []);

   FRATE= 13;
   expis= [];
   insps= [];
   for i= 1:6
         nn= sprintf('../%s/%s/%s-00%d.get', 'raw-data', dir, dir, i);
         [d_insp,d_expi, W] = get_vent_meas( nn, imdl_v, ROIv, FRATE, fidb);
         imname = sprintf('../img/%s-00%d_sig.png',dir,i);
         print_convert(imname);
         htmlimage(fidb,{imname});

         if length(d_expi) == 256;
            d_expi= d_expi(imdl_v.fwd_model.meas_select);
            d_insp= d_insp(imdl_v.fwd_model.meas_select);
         end
         expis= [expis, d_expi]; % keep the exps
         insps= [insps, d_insp];
   end




function fid= htmlheader(fname, extra);
   fid= fopen(fname,'w');
   fprintf(fid,'<HTML><BODY>\n');
   fprintf(fid,'<TABLE><TR><TH> Subject\n');
   if nargin==1;
      for i=1:6
         fprintf(fid,'<TH>Meas 00%d',i);
      end
   else
      for i= extra
         fprintf(fid,'<TH><FONT SIZE="-1">%s</FONT>',i{1});
      end
   end


function htmlnewdir(fid,dir);
   fprintf(fid,'\n<TR>\n<TH><b>%s</b>',dir);

function htmlimage(fid,imname)
%  fprintf(fid,'<TD><img src="%s"><BR>%s-00%d\n',imname,dir,i);
   fprintf(fid,'<TD>');
   for im= imname
      str= im{1};
      fprintf(fid, '<img src="%s"><BR>', str);
   end
   fprintf(fid,'\n');

function htmltail(fid)
   fprintf(fid,'</TABLE>\n');
   fprintf(fid,'</BODY></HTML>\n');
   fclose(fid);

% raster -> 64x64 
% alg    -> '-bp'
% imfrac -> 'S12-001'
% iname  -> 'VT,EIT,ZEEP ...'
% data(12,001).ghnn = [1,2,3,4];

function write_analyse(raster, imfrac, alg, iname);
   ani_no= str2num(imfrac(2:3));
   exp_no= str2num(imfrac(5:7));
   fid= fopen('analyse_out.m','a');
   fprintf(fid,'data(%2d,%3d).%-8s=[',ani_no,exp_no,alg   );

   ss= raster(:)'*[ ...
        calc_ROI_splits( raster, 32,  0), ... % All of it
        calc_ROI_splits( raster, 32, 10), ... % All of it - bottom
        calc_ROI_splits( raster, 25, 10), ...
        calc_ROI_splits( raster, 20, 10)];

   fprintf(fid,' %+8.5f', ss);
   fprintf(fid,']''; %% %s\n', iname); %Stupid matlab concatenates
   fclose(fid);

function All_Top_Bot= calc_ROI_splits( raster, rad, bot)
   xctr= mean([1,size(raster,1)]);
   yctr= mean([1,size(raster,2)]);
   ROI= calc_lung_roi( raster, yctr, xctr, rad, bot);

   sROI= sum(ROI,2)/sum(ROI(:));
   [jnk, smid]= min(abs(cumsum(sROI) - 0.5));

   ROIa = ROI                      ;
   ROIb = ROI; ROIb(1:smid    ,:)=0;
   ROIt = ROI; ROIt(1+smid:end,:)=0;

% We don't want the average in the roi, we want the 
%  ROIa= ROIa(:)/sum(ROIa(:));
%  ROIb= ROIb(:)/sum(ROIb(:));
%  ROIt= ROIt(:)/sum(ROIt(:));
ROIa= ROIa(:);
ROIb= ROIb(:);
ROIt= ROIt(:);

   All_Top_Bot= [ROIa, ROIt, ROIb];

function [imname,raster]= recon_write(imdl, Dt, imfrac, alg) 
   imname= ['../img/',imfrac,alg,'.png'];
   vent_img = inv_solve( imdl, Dt(:,1), Dt(:,2)); 
   if size(Dt,2)>2
      vent_img2 = inv_solve( imdl, Dt(:,3),Dt(:,4));
      vent_img.elem_data = vent_img.elem_data - vent_img2.elem_data;
   end
   clf; axes('position',[0,0,1,1]); axis xy;
   vent_img.calc_colours.ref_level= 0;

   out_img= show_slices(vent_img);
   drawnow
% STUPID CRAP MATLAB DOESN'T LISTEN TO GCF settings!!
   imwrite(out_img, colormap, imname);

   raster = calc_slices(vent_img);

% SAVE EACH IMAGE AS A RASTER HERE
   matname= ['../recon/' imfrac,alg,'.mat'];
   save(matname,'raster');

   raster(isnan(raster))= 0;

function [imname,raster]= recon_write_lr(imdl, Dt, imfrac, alg) 
   imname= [imfrac,alg,'.png'];
   vent_img = inv_solve( imdl, Dt(:,1),Dt(:,2));
   if size(Dt,2)>2
      vent_img2 = inv_solve( imdl, Dt(:,3),Dt(:,4));
      vent_img.elem_data = vent_img.elem_data - vent_img2.elem_data;
   end
   raster = calc_slices(vent_img);
   raster(isnan(raster))= 0;
LIM= 0.25;

   max_v = min(vent_img.elem_data);
   lim_img= vent_img.elem_data<LIM*max_v;
   vent_img.elem_data = [vent_img.elem_data, max_v*lim_img];

   vent_img.calc_colours.ref_level= 0;
   out_img= show_slices(vent_img);
   lung_roi= calc_lung_roi( out_img, 32.5, 64+32.5, 25, 10);
   out_img = out_img + 50*lung_roi.*sign(50 - out_img);

%  out_img= [out_img;126*ones(25+64+40,64)];
   out_img= [out_img(1:64,:); 126*ones(64+25,64); ...
             out_img(65:(128-10),:);126*ones(40,64)];
if 0
   imwrite(out_img, colormap, imname);
else
   if ~exist('OCTAVE_VERSION')
      pos = get(gcf,'paperposition');
      set(gcf,'paperposition',[0,0,size(out_img')/64]);
      ppos = get(gcf,'papersize');
      set(gcf,'papersize',        [size(out_img')/64]);
   end
   clf; axes('position',[0,0,1,1]);

   out_img= first_row( out_img, raster, 64*1+25);
   if 0 % defined from image
      ROI= raster<LIM*min(min(raster));
   else % simple anatomy
      ROI= calc_lung_roi( raster, 32.5, 32.5, 25, 10);
   end
   lungstart = 8;
   out_img= second_row( out_img, raster, lungstart, ROI, 64*3+25-10);


   image(out_img); axis off

%  bottom= 64*3+25; line([0,64],(bottom-0.5)*[1,1],'Color',[0,0,0])
      
   puttext(raster, 64)
   
   if exist('OCTAVE_VERSION');
      print(imname,'-dpng','-S200,75');
   else;
      print('-dpng','-r64',imname);
      set(gcf,'paperposition',pos);
      set(gcf,'papersize',ppos);
   end
end

function puttext(raster,vposn)
% amt in each slice
   q(1)= sum(sum( raster(1:32,1:32)));
   q(2)= sum(sum( raster(1:32,33:64)));
   q(3)= sum(sum( raster(33:64,1:32)));
   q(4)= sum(sum( raster(33:64,33:64)));

   k=1;for i=[3,1:2,4]
     txt= sprintf('%2.0f',100*q(i)/sum(q));
     text(14*(k-1)+3 + 5*(i==2) - 2*(i==1), ...
          vposn+7+(i>2)*10,txt,'FontSize',12);
     k=k+1;
   end

   text(32,vposn+21,sprintf('%3.1f',100*abs(max(raster(:)))), ...
          'FontSize',12,'HorizontalAlignment','center', ...
          'FontWeight','Bold')

function out_img= first_row( out_img, raster, bottom)
%  axes('position',[0,0,1,1/(3+25/64)])
   ri= zeros(8,1);
   le= zeros(8,1);
   for i=1:8
      le(i)= sum(sum(raster( 8*i-(0:7),  1:32)));
      ri(i)= sum(sum(raster( 8*i-(0:7), 33:64)));
   end
   maxv= max([abs(le);abs(ri)]) + 1e-4; % prevent / zeros
   le=le/maxv;
   ri=ri/maxv;
   for i=1:8
      ri_dist = 0:floor(31*abs(ri(i))); 
      out_img(bottom + 8*i-(0:7), 33 + ri_dist) = 128+96*sign(ri(i));
      le_dist = 0:floor(31*abs(le(i))); 
      out_img(bottom + 8*i-(0:7), 32 - le_dist) = 128+96*sign(le(i));
   end

function out_img= second_row( out_img, raster, lungstart, ROI, bottom)
   ri= zeros(40,1);
   le= zeros(40,1);
   for i=(1:40)
      llim = sum(ROI(i + lungstart, 1:32));
      if llim>0
         le(i)= sum(raster( i + lungstart, 32-(0:llim)))/llim;
      end
      rlim = sum(ROI(i + lungstart,33:64));
      if rlim>0
         ri(i)= sum(raster( i + lungstart, 33+(0:rlim)))/rlim;
      end
   end
   maxv= max([abs(le);abs(ri)]) + 1e-4; % prevent / zeros
   le=le/maxv;
   ri=ri/maxv;
   for i=(1:40)
      ri_dist = 0:floor(31*abs(ri(i))); 
      out_img(bottom + i, 33 + ri_dist) = 128+96*sign(ri(i));
      le_dist = 0:floor(31*abs(le(i))); 
      out_img(bottom + i, 32 - le_dist) = 128+96*sign(le(i));
   end

function lung_roi= calc_lung_roi( out_img, xc, yc, rad, ybot );
   lung_roi = logical(zeros(size(out_img)));
   [x,y]= meshgrid(1:size(out_img,2),1:size(out_img,1));
   r= sqrt( (x-xc).^2 + (y-yc).^2 );
   lung_roi( r<rad & y<yc+ybot ) = 1;

function [fid,fidb,fidi] = openoutfiles
fid= fopen('analyse_out.m','w');
   fprintf(fid,'% analyse output: %s',datestr(now)); fclose(fid);
fidb= htmlheader('../html/breathing.html');
fidi= htmlheader('../html/images.html', ...
         {'V<sub>T,<br>ZEEP,21,#1</sub>', 'V<sub>T,<br>PEEP,21,#1</sub>' ...
          'V<sub>T,<br>ZEEP,100</sub>',   'V<sub>T,<br>PEEP,100</sub>' ...
          'V<sub>T,<br>ZEEP,21,#2</sub>', 'V<sub>T,<br>PEEP,21,#2</sub>' ...
          '&Delta;EELV <br><sub>Z&minus;P,21,#1</sub>', ...
          '&Delta;EELV <br><sub>Z&minus;P,21,#2</sub>', ...
          '&Delta;EELV <br><sub>Z&minus;P,100</sub>', ...
          '&Delta;EELV <br><sub>Z&minus;Z,21</sub>', ...
          '&Delta;EELV <br><sub>P&minus;P,21</sub>', ...
          '&Delta;EELV <br><sub>Z21&minus;Z100,#1</sub>', ...
          '&Delta;EELV <br><sub>Z21&minus;Z100,#2</sub>', ...
          '&Delta;EELV <br><sub>P21&minus;P100,#1</sub>', ...
          '&Delta;EELV <br><sub>P21&minus;P100,#2</sub>'});

%        {'V<sub>T,EIT,ZEEP,21,#1</sub>', 'V<sub>T,EIT,PEEP,21,#1</sub>' ...
%         'V<sub>T,EIT,ZEEP,100</sub>',   'V<sub>T,EIT,PEEP,100</sub>' ...
%         'V<sub>T,EIT,ZEEP,21,#2</sub>', 'V<sub>T,EIT,PEEP,21,#2</sub>' ...
%.. %     '(E&minus;I)<sub>M2</sub>&minus;(E&minus;I)<sub>M1</sub>', ...
%.. %     '(E&minus;I)<sub>M4</sub>&minus;(E&minus;I)<sub>M3</sub>', ...
%         '&Delta;EELV<sub>EIT,Z&minus;P,21,#1</sub>', ...
%         '&Delta;EELV<sub>EIT,Z&minus;P,21,#2</sub>', ...
%         '&Delta;EELV<sub>EIT,Z&minus;P,100</sub>', ...
%         '&Delta;EELV<sub>EIT,Z&minus;Z,21</sub>', ...
%         '&Delta;EELV<sub>EIT,P&minus;P,21</sub>', ...
%         '&Delta;EELV<sub>EIT,Z21&minus;Z100,#1</sub>', ...
%         '&Delta;EELV<sub>EIT,Z21&minus;Z100,#2</sub>', ...
%         '&Delta;EELV<sub>EIT,P21&minus;P100,#1</sub>', ...
%         '&Delta;EELV<sub>EIT,P21&minus;P100,#2</sub>'});
