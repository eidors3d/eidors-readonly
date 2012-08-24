function show_pseudosection( fwd_model, data, orientation)
%SHOW_PSEUDOSECTION: show a pseudo-section image of data
%
% OPTIONS:
%   fwd_model - fwd model of the domain
%   data      - data to display
%   orientation - orientation of the pseudo-section
%      orientation = 'profile': show section in +z direction from 
%         the electrode plane
%      orientation = 'gallery': show section in gallery outside electrode plane
%      orientation (other options not yet available)

% (C) 2005-2008 Nolwenn Lesparre. License: GPL version 2 or version 3
% $Id$

if isstr(fwd_model) && strcmp(fwd_model,'UNIT_TEST'); do_unit_test; return; end 

[fmdl,data_tomel,I] = geomFactor(data_tomel,fmdl)
plotPseudoSectionProfile(data,v,normalisation,elec_posn,name)

function plotPseudoSectionProfile(data,v,normalisation,elec_posn,name)
   fs= 20;
   
   
   vi= rem(data(:,3:6),1000);
   vi= vi-min(vi(:))+1;
   
   elecNumber= vi;
   
   xps= (elec_posn(elecNumber(:,1),1)+elec_posn(elecNumber(:,2),1))/2;
   a= abs(elecNumber(:,1)-elecNumber(:,2));
   
   de= elec_posn(1,1)-elec_posn(2,1);
   
   % Identiy reciprocal data(elecNumber(:,1)-elecNumber(:,2)) > abs(elecNumber(:,3)-elecNumber(:,4))
   R= find(abs(elecNumber(:,1)-elecNumber(:,2)) < abs(elecNumber(:,3)-elecNumber(:,4)));
   
   xps(R)= (elec_posn(elecNumber(R,3),4)+elec_posn(elecNumber(R,4),4))/2;
   a(R)= abs(elecNumber(R,3)-elecNumber(R,4));
   
   zps= -a*de/2;
   
   k= normalisation;
   app_resistivity=  k.*v;
   
   P= zps+1i*xps;
   [Pu,iu,ju]= unique(P);
   
   xu= xps(iu);
   zu= zps(iu);
   app_resistivityu= iu*0;
   for i= 1:length(Pu)
       app_resistivityu(i)= mean(app_resistivity(ju==i));
   end
   
   aa= find(app_resistivityu<0);
   figure;scatter(xu,zu,100,(app_resistivityu),'filled');  colorbar
   xlabel('Radius (cm)','fontsize',fs,'fontname','Times');
   ylabel('Depth (m)','fontsize',fs,'fontname','Times')
   axis equal; axis tight; ylim([-130 0]) %caxis([0 12]);
   if isempty(app_resistivityu(app_resistivityu<0));  
       caxis([5 120])
   else
       load('matFiles/BlueWhitePinkMap','mycmap')
       caxis([-10 10]); colormap(mycmap)
   end
   set(gca,'fontsize',fs,'fontname','Times')
   im= ['Figures/PseudoSection' name];
   saveas(gcf, [im '.pdf'], 'pdf');
   
   % figure;scatter(xu,zu,100,log10(app_resistivityu),'filled');  colorbar
   % xlabel('Radius (cm)','fontsize',fs,'fontname','Times');
   % ylabel('Depth (m)','fontsize',fs,'fontname','Times')
   % axis equal; axis tight; %caxis([0 12]);
   % set(gca,'fontsize',fs,'fontname','Times')
   % im= ['Figures/LogPseudoSection' name];
   % saveas(gcf, [im '.pdf'], 'pdf');
   
   % figure;scatter(xpsd,zpsd,150,(app_resistivity),'filled');  colorbar
   % xlabel('Radius (cm)','fontsize',fs,'fontname','Times');
   % ylabel('Depth (m)','fontsize',fs,'fontname','Times')
   % axis tight; caxis([0 10]);
   % set(gca,'fontsize',fs,'fontname','Times','DataAspectRatio',[100 1 1])
   % im= ['Figures/PseudoSection' name];
   % saveas(gcf, [im '.pdf'], 'pdf');
   % 
   % 
   % figure;scatter(xpsd(1:115),zpsd(1:115),150,(app_resistivity(1:115)),'filled');  colorbar
   % xlabel('Radius (cm)','fontsize',fs,'fontname','Times');
   % ylabel('Depth (m)','fontsize',fs,'fontname','Times')
   % axis tight; caxis([0 10]);
   % set(gca,'fontsize',fs,'fontname','Times','DataAspectRatio',[100 1 1])
   % im= ['Figures/PseudoSection' name];
   % saveas(gcf, [im '.pdf'], 'pdf');
   
   return
   x= min(xps):0.5:max(xps);
   z= min(zps):0.5:max(zps);
   
   App= TriScatteredInterp(xu,zu,app_resistivityu);
   [X,Z]= meshgrid(x,z);
   % Appz= App(X,Z);
   
   % figure; pcolor(X,Z,Appz);colorbar; shading flat
   % xlabel('Radius (cm)','fontsize',fs,'fontname','Times');
   % ylabel('Depth (m)','fontsize',fs,'fontname','Times')
   % axis equal; axis tight; %caxis([0 12]);
   % set(gca,'fontsize',fs,'fontname','Times')
   % im= ['Figures/PseudoSectionInterp' name];
   % saveas(gcf, [im '.pdf'], 'pdf');
   % 
   % figure; pcolor(X,Z,log10(Appz));colorbar; shading flat
   % xlabel('Radius (cm)','fontsize',fs,'fontname','Times');
   % ylabel('Depth (m)','fontsize',fs,'fontname','Times')
   % axis equal; axis tight; %caxis([0 12]);
   % set(gca,'fontsize',fs,'fontname','Times')
   % im= ['Figures/PseudoSectionInterp' name];
   % saveas(gcf, [im '.pdf'], 'pdf');
   
   AppL= TriScatteredInterp(xu,zu,log10(app_resistivityu));
   [X,Z]= meshgrid(x,z);
   AppzL= AppL(X,Z);
   AppzLF= medfilt2(AppzL,[20 20]);
   
   figure; pcolor(X,Z,AppzLF);colorbar; shading flat
   xlabel('Radius (cm)','fontsize',fs,'fontname','Times');
   ylabel('Depth (m)','fontsize',fs,'fontname','Times')
   axis equal; axis tight; %caxis([0 12]);
   set(gca,'fontsize',fs,'fontname','Times')
   im= ['Figures/PseudoSectionInterp' name];
   saveas(gcf, [im '.pdf'], 'pdf');

end
   
function [fmdl,data_tomel,I] = geomFactor(data_tomel,fmdl)
   
   % Compute the geometrical factor corresponding to each data set and record
   % it in the 11th column of the data_tomel
   data_tomel(:,8)= data_tomel(:,8)*100./data_tomel(:,7);
   data_tomel(:,7)= data_tomel(:,7)*0+100;
   data_tomel(:,9)= data_tomel(:,8)./data_tomel(:,7);
   
   % if ~isfield(fmdl,'stimulation') || size(fmdl.stimulation,2) < size(data_tomel,1)
   stim = mk_stim_patterns_tomel(fmdl.misc.elec_posn,data_tomel);
   fmdl.stimulation = stim;
   % end
   if size(data_tomel,2) < 11 || sum(data_tomel(:,11)) == size(data_tomel,1) 
       img1= mk_image(fmdl,1);
       vh1= fwd_solve(img1);
       normalisation= 1./vh1.meas;
       data_tomel(:,11)=  normalisation;
   else
       normalisation= data_tomel(:,11);
   end
   I= speye(length(normalisation));
   I(1:size(I,1)+1:size(I,1)*size(I,1))= normalisation;
end

function do_unit_test
end
