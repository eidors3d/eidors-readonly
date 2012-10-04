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

% [fmdl,data_tomel,I] = geomFactor(data,fwd_model)
ifstrcmp(orientation,'HorizontalDownaward')

switch(upper(orientation))
    case 'HORIZONTALDOWNWARD'; plotPseudoSectionProfileDown(fwd_model,data)
    case 'CIRCULAROUTSIDE'; plotPseudoSectionCircularOut(fwd_model,data)
  otherwise;
    error('No orientation of type "%s" available', upper(pat_type));
end

end

function plotPseudoSectionProfileDown(fmdl,data)
   fs= 20;
   A= zeros(length(fmdl.stimulation),1);
   B= zeros(length(fmdl.stimulation),1);
   M= zeros(length(fmdl.stimulation),1);
   N= zeros(length(fmdl.stimulation),1);
   for i=1:length(fmdl.stimulation)
       A(i)= find(fmdl.stimulation(1,i).stim_pattern<0);
       B(i)= find(fmdl.stimulation(1,i).stim_pattern>0);
       M(i)= find(fmdl.stimulation(1,i).meas_pattern<0);
       N(i)= find(fmdl.stimulation(1,i).meas_pattern>0);
   end
   elecNumber= [A B M N];
   
   elec_posn= zeros(length(fmdl.electrode),3);
   for i=1:length(fmdl.electrode)
   elec_posn(i,:)= mean(fmdl.nodes(fmdl.electrode(1,i).nodes,:),1);
   end
 
   
   xps= (elec_posn(elecNumber(:,1),1)+elec_posn(elecNumber(:,2),1))/2;
   a= abs(elecNumber(:,1)-elecNumber(:,2));
   de= elec_posn(1,1)-elec_posn(2,1);
   
   % Identiy reciprocal data(elecNumber(:,1)-elecNumber(:,2)) > abs(elecNumber(:,3)-elecNumber(:,4))
   R= find(abs(elecNumber(:,1)-elecNumber(:,2)) < abs(elecNumber(:,3)-elecNumber(:,4)));
   
   if ~isempty(R)
       xps(R)= (elec_posn(elecNumber(R,3),4)+elec_posn(elecNumber(R,4),4))/2;
       a(R)= abs(elecNumber(R,3)-elecNumber(R,4));
   end
   zps= -abs(a*de/2);
      
   P= zps+1i*xps;
   [Pu,iu,ju]= unique(P);
   
   xu= xps(iu);
   zu= zps(iu);
   du= iu*0;
   for i= 1:length(Pu)
       du(i)= mean(data(ju==i));
   end
   
   figure;scatter(xu,zu,100,(du),'filled');  colorbar
   xlabel('Distance (m)','fontsize',fs,'fontname','Times');
   ylabel('Depth (m)','fontsize',fs,'fontname','Times')
    axis equal; axis tight;
   set(gca,'fontsize',fs,'fontname','Times')
end



function plotPseudoSectionCircularOut(fmdl,data)
   fs= 20;
   A= zeros(length(fmdl.stimulation),1);
   B= zeros(length(fmdl.stimulation),1);
   M= zeros(length(fmdl.stimulation),1);
   N= zeros(length(fmdl.stimulation),1);
   for i=1:length(fmdl.stimulation)
       A(i)= find(fmdl.stimulation(1,i).stim_pattern<0);
       B(i)= find(fmdl.stimulation(1,i).stim_pattern>0);
       M(i)= find(fmdl.stimulation(1,i).meas_pattern<0);
       N(i)= find(fmdl.stimulation(1,i).meas_pattern>0);
   end
   elecNumber= [A B M N];
   
   elec_posn= zeros(length(fmdl.electrode),3);
   for i=1:length(fmdl.electrode)
   elec_posn(i,:)= mean(fmdl.nodes(fmdl.electrode(1,i).nodes,:),1);
   end
 
   
   xps= (elec_posn(elecNumber(:,1),1)+elec_posn(elecNumber(:,2),1))/2;
   a= abs(elecNumber(:,1)-elecNumber(:,2));
   de= elec_posn(1,1)-elec_posn(2,1);
   
   % Identiy reciprocal data(elecNumber(:,1)-elecNumber(:,2)) > abs(elecNumber(:,3)-elecNumber(:,4))
   R= find(abs(elecNumber(:,1)-elecNumber(:,2)) < abs(elecNumber(:,3)-elecNumber(:,4)));
   
   if ~isempty(R)
       xps(R)= (elec_posn(elecNumber(R,3),4)+elec_posn(elecNumber(R,4),4))/2;
       a(R)= abs(elecNumber(R,3)-elecNumber(R,4));
   end
   zps= -abs(a*de/2);
   
%    k= normalisation;
%    app_resistivity=  k.*v;
   app_resistivity= data;
   
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
%    xlabel('Radius (cm)','fontsize',fs,'fontname','Times');
%    ylabel('Depth (m)','fontsize',fs,'fontname','Times')
    axis equal; axis tight;%ylim([-130 0]) %caxis([0 12]);
   set(gca,'fontsize',fs,'fontname','Times')

   
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

function do_unit_test
shape_str = ['solid top    = plane(0,0,0;0,1,0);\n' ...
             'solid mainobj= top and orthobrick(-100,-200,-100;410,10,100) -maxh=20.0;\n'];
e0 = linspace(0,310,64)';
elec_pos = [e0,0*e0,0*e0,1+0*e0,0*e0,0*e0];
elec_shape= [0.1,0.1,1];
elec_obj = 'top';
fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);
%fmdl.nodes = fmdl.nodes(:,[1,3,2]);
% fmdl.stimulation= stim_pattern_geophys( 64, 'Wenner', {'spacings', 1:32} );


% spacing= [1 1 1 2 3 4 6 8 8 11 12 14 17];
% multiples= [1 2 3 2 5/3 6/4 7/6 1 10/8 1 13/12 15/14 1];
spacing= [1 1 1 2 3 3 4 4 5 6 6 7 8 8 9 10 10 11 12 12 13 14 14 15 16 17];
multiples= [1 2 3 2 1 5/3 1 2  1 1 7/6 1 1 10/8 1 1 12/10 1 1 13/12 1 1 15/14 1 1 1];
% fmdl.stimulation= stim_pattern_geophys( 64, 'Schlumberger', {'spacings', spacing,'multiples',multiples} );
fmdl.stimulation= stim_pattern_geophys( 64, 'DipoleDipole', {'spacings', spacing,'multiples',multiples} );


img1= mk_image(fmdl,1);
vh1= fwd_solve(img1);
normalisation= 1./vh1.meas;
I= speye(length(normalisation));
I(1:size(I,1)+1:size(I,1)*size(I,1))= normalisation;


img = mk_image(fmdl,0+ mk_c2f_circ_mapping(fmdl,[100;-50;0;50])*100);
img.elem_data(img.elem_data==0)= 0.1;
dd  = fwd_solve(img);
show_pseudosection( fmdl, I*dd.meas, '')

end
