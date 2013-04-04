function fwd_model= show_pseudosection( fwd_model, data)
%SHOW_PSEUDOSECTION: show a pseudo-section image of data
%
% OPTIONS:
%   fwd_model - fwd model of the domain
%   data      - data to display
%
% OPTIONS:
% fwd_model.show_pseudosection.orientation - of the pseudo-section
%      orientation = 'horizontaldownward': section in +z direction from 
%         the electrode plane (default)
%      orientation = 'circularoutside': gallery outside electrode plane
%      orientation = 'vertical'
%      orientation = 'circularinside'

% (C) 2005-2008 Nolwenn Lesparre. License: GPL version 2 or version 3
% $Id$

if ischar(fwd_model) && strcmp(fwd_model,'UNIT_TEST'); do_unit_test; return; end 

if size(fwd_model.electrode,2) <= 16
    fwd_model.misc.sizepoint= 500;
elseif size(fwd_model.electrode,2) <= 32
    fwd_model.misc.sizepoint= 200;
else
    fwd_model.misc.sizepoint= 50;
end

try
   orientation= fwd_model.show_pseudosection.orientation;
catch
   orientation = 'horizontaldownward';
end

if iscell(orientation) 
    if strcmp(orientation{2},'yz')
    fwd_model.nodes= fwd_model.nodes(:,[2 3 1]);
    elseif strcmp(orientation{2},'xz')
     fwd_model.nodes= fwd_model.nodes(:,[1 3 2]);
    end
    orientation= orientation{1};
end

switch(upper(orientation))
    case 'HORIZONTALDOWNWARD';  [depth,location]= plotPseudoSectionProfileDown(fwd_model,data);
    case 'VERTICAL';            [depth,location]= plotPseudoSectionProfileVert(fwd_model,data);
    case 'CIRCULAROUTSIDE';     [depth,location]= plotPseudoSectionCircularOut(fwd_model,data);
    case 'CIRCULARINSIDE';      [depth,location]= plotPseudoSectionCircularIn(fwd_model,data);
  otherwise;
    error('No orientation of type "%s" available', upper(orientation));
end

fwd_model.show_pseudosection.depth= depth;
fwd_model.show_pseudosection.location= location;

end

function [a,xps]= plotPseudoSectionProfileDown(fmdl,data)
   fs= 20;
   
   [elec_posn,elecNumber]= electrodesPosition(fmdl);
   
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
   
   scatter(xu,zu,fmdl.misc.sizepoint,(du),'filled','MarkerEdgeColor','k');
   xlabel('Distance (m)','fontsize',fs,'fontname','Times');
   ylabel('Pseudo-depth (m)','fontsize',fs,'fontname','Times')
    axis equal; axis tight;
   set(gca,'fontsize',fs,'fontname','Times')
end

function [a,zps]= plotPseudoSectionProfileVert(fmdl,data)
   fs= 20;
   
   [elec_posn,elecNumber]= electrodesPosition(fmdl);
  
   zps= (elec_posn(elecNumber(:,1),3)+elec_posn(elecNumber(:,2),3))/2;
   a= abs(elec_posn(elecNumber(:,1),3)-elec_posn(elecNumber(:,2),3));
%    a= abs(elecNumber(:,1)-elecNumber(:,2));
%    de= abs(elec_posn(1,3)-elec_posn(2,3));
   
   % Identiy reciprocal data(elecNumber(:,1)-elecNumber(:,2)) > abs(elecNumber(:,3)-elecNumber(:,4))
   R= find(abs(elecNumber(:,1)-elecNumber(:,2)) < abs(elecNumber(:,3)-elecNumber(:,4)));
   
   if ~isempty(R)
       zps(R)= (elec_posn(elecNumber(R,3),3)+elec_posn(elecNumber(R,4),3))/2;
       a(R)= abs(elec_posn(elecNumber(R,3),3)-elec_posn(elecNumber(R,4),3));
%        a(R)= abs(elecNumber(R,3)-elecNumber(R,4));
   end
   xps= a/3;
%    xps= abs(a*de/2);
      
   P= zps+1i*xps;
   [Pu,iu,ju]= unique(P);
   
   xu= xps(iu);
   zu= zps(iu);
   du= iu*0;
   for i= 1:length(Pu)
       du(i)= mean(data(ju==i));
   end
   
   scatter(xu,zu,fmdl.misc.sizepoint,(du),'filled','MarkerEdgeColor','k'); 
   xlabel('Pseudo distance (m)','fontsize',fs,'fontname','Times');
   if zps(1)<0
       ylabel('Depth (m)','fontsize',fs,'fontname','Times')
   else
       ylabel('Height (m)','fontsize',fs,'fontname','Times')
   end
   axis equal; axis tight;
   set(gca,'fontsize',fs,'fontname','Times')
end

function [r_point,th_point]= plotPseudoSectionCircularIn(fmdl,data)
fs= 20;
[elec_posn,elecNumber] = electrodesPosition(fmdl);
[r_point,th_point,r] = polarPosition(elecNumber,elec_posn);

r_point= (r-r_point/pi);%*9*pi*r/10+pi*r/10;

[x_point,y_point]= pol2cart(th_point,r_point);

P= x_point+1i*y_point;
[Pu,iu,ju]= unique(P);
   
xu= x_point(iu);
zu= y_point(iu);
du= iu*0;
for i= 1:length(Pu)
    du(i)= mean(data(ju==i));
end
   
scatter(xu,zu,fmdl.misc.sizepoint,(du),'filled','MarkerEdgeColor','k');  
xlabel('X (m)','fontsize',fs,'fontname','Times');
ylabel('Y (m)','fontsize',fs,'fontname','Times')
axis equal; axis tight;
set(gca,'fontsize',fs,'fontname','Times')
end


function [r_point,th_point]= plotPseudoSectionCircularOut(fmdl,data)
fs= 20;
[elec_posn,elecNumber] = electrodesPosition(fmdl);
[r_point,th_point,r] = polarPosition(elecNumber,elec_posn);

[x_point,y_point]= pol2cart(th_point,r_point+r);

P= x_point+1i*y_point;
[Pu,iu,ju]= unique(P);
   
xu= x_point(iu);
zu= y_point(iu);
du= iu*0;
for i= 1:length(Pu)
    du(i)= mean(data(ju==i));
end
   
scatter(xu,zu,fmdl.misc.sizepoint,(du),'filled','MarkerEdgeColor','k');  colorbar
xlabel('X (m)','fontsize',fs,'fontname','Times');
ylabel('Y (m)','fontsize',fs,'fontname','Times')
axis equal; axis tight;
set(gca,'fontsize',fs,'fontname','Times')
end


function [elec_posn,elecNumber] = electrodesPosition(fmdl)

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
   
end

function [r_point,th_point,r] = polarPosition(elecNumber,elec_posn)

TH= elecNumber*0;
R= elecNumber*0;
for i= 1:4
 [TH(:,i),R(:,i)]= cart2pol(elec_posn(elecNumber(:,i),1),elec_posn(elecNumber(:,i),2));
end

th_point= (TH(:,3)+TH(:,4))/2;
idx_discontinuity= find((TH(:,3)>pi/2 & TH(:,4)<-pi/2) | (TH(:,3)<-pi/2 & TH(:,4)>pi/2));
th_point(idx_discontinuity)= (TH(idx_discontinuity,3)+2*pi+TH(idx_discontinuity,4))/2;
th_point(th_point>pi)= th_point(th_point>pi)-2*pi;
th_point(th_point<-pi)= th_point(th_point<-pi)+2*pi;

Arc_length_AM= abs(TH(:,3)-TH(:,1));
idx_discontinuity= find((TH(:,3)>pi/2 & TH(:,1)<-pi/2));
Arc_length_AM(idx_discontinuity)= 2*pi+TH(idx_discontinuity,1)-TH(idx_discontinuity,3);
idx_discontinuity= find((TH(:,1)>pi/2 & TH(:,3)<-pi/2));
Arc_length_AM(idx_discontinuity)= 2*pi+TH(idx_discontinuity,3)-TH(idx_discontinuity,1);
Arc_length_AM= Arc_length_AM.*(R(:,1)+R(:,3))/2;

Arc_length_NB= abs(TH(:,4)-TH(:,2));
idx_discontinuity= find((TH(:,4)>pi/2 & TH(:,2)<-pi/2));
Arc_length_NB(idx_discontinuity)= 2*pi+TH(idx_discontinuity,2)-TH(idx_discontinuity,4);
idx_discontinuity= find((TH(:,2)>pi/2 & TH(:,4)<-pi/2));
Arc_length_NB(idx_discontinuity)= 2*pi+TH(idx_discontinuity,4)-TH(idx_discontinuity,2);
Arc_length_NB= Arc_length_NB.*(R(:,2)+R(:,4))/2;

r= mean(R(:));
r_point= (Arc_length_AM+Arc_length_NB)/2;
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
subplot(221);
show_pseudosection( fmdl, I*dd.meas )

subplot(222);
fmdl.show_pseudosection.orientation = 'horizontaldownward';
show_pseudosection( fmdl, I*dd.meas )

subplot(223);
fmdl.show_pseudosection.orientation = 'vertical';
show_pseudosection( fmdl, I*dd.meas )

end
