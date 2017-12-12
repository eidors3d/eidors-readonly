function fwd_model= show_pseudosection( fwd_model, data, orientation)
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

if ~isfield(fwd_model, 'misc') || ~isfield(fwd_model.misc,'sizepoint')
  if size(fwd_model.electrode,2) <= 16
      fwd_model.misc.sizepoint= 500;
  elseif size(fwd_model.electrode,2) <= 32
      fwd_model.misc.sizepoint= 200;
  else
      fwd_model.misc.sizepoint= 50;
  end
end

if nargin<3 % otherwise orientation is specified
try
   orientation= fwd_model.show_pseudosection.orientation;
catch
   orientation = 'horizontaldownward';
end
end

if ~isfield(fwd_model,'show_pseudosection') || ~isfield(fwd_model.show_pseudosection,'fs')
    fwd_model.show_pseudosection.fs=16;
end

if iscell(orientation) 
    if strcmp(orientation{2},'yz')
    fwd_model.nodes= fwd_model.nodes(:,[2 3 1]);
    rotation= 'yz';
    elseif strcmp(orientation{2},'xz')
     fwd_model.nodes= fwd_model.nodes(:,[1 3 2]);
     rotation= 'xz';
    end
    orientation= orientation{1};
end
location34= [];
electrodesLocation= [];
switch(upper(orientation))
    case 'HORIZONTALDOWNWARD';  [depth,location]= plotPseudoSectionProfileDown(fwd_model,data);
    case 'VERTICAL';            [depth,location]= plotPseudoSectionProfileVert(fwd_model,data);
    case 'SPECIFICHORIZONTALDOWNWARD';  [depth,location,location34,electrodesLocation]= plotPseudoSectionProfileSpecificDown(fwd_model,data); 
    case 'SPECIFICHORIZONTALUPWARD';  [depth,location,L,d]= plotPseudoSectionProfileSpecificUp(fwd_model,data);
    case 'CIRCULAROUTSIDE';     [depth,location]= plotPseudoSectionCircularOut(fwd_model,data);
    case 'CIRCULARINSIDE';      [depth,location]= plotPseudoSectionCircularIn(fwd_model,data);
    case 'CIRCULARVOLCANO';     [depth,location]= plotPseudoSectionCircularInVolcano(fwd_model,data);
    case 'HORIZONTALLAYERS';	[depth,location,electrodesLocation]= plotPseudoSectionHorizontalLayers(fwd_model,data);
  otherwise;
    error('No orientation of type "%s" available', upper(orientation));
end

if exist('rotation','var') 
    if strcmp(rotation,'yz')
    fwd_model.nodes= fwd_model.nodes(:,[3 1 2]);  
    elseif strcmp(rotation,'xz')
     fwd_model.nodes= fwd_model.nodes(:,[1 2 3]);
    end
end

fwd_model.show_pseudosection.depth= depth;
fwd_model.show_pseudosection.location= location;
fwd_model.show_pseudosection.location34= location34;
fwd_model.show_pseudosection.electrodesLocation= electrodesLocation;
% fwd_model.show_pseudosection.L= L;
% fwd_model.show_pseudosection.d= d;
end

function [depth,location,electrodesLocation]= plotPseudoSectionHorizontalLayers(fmdl,data)
   fs= fmdl.show_pseudosection.fs;
   depthPrecision= fmdl.show_pseudosection.depthPrecision;
   depthRatio= fmdl.show_pseudosection.depthRatio;
   
   xprecision= fmdl.show_pseudosection.xprecision;
   yprecision= fmdl.show_pseudosection.yprecision;
   locFun= fmdl.show_pseudosection.locFun;
   locVar= fmdl.show_pseudosection.locVar;
   
   depthFun= fmdl.show_pseudosection.depthFun;
   depthVar= fmdl.show_pseudosection.depthVar;
      
   if ~isfield(fmdl.show_pseudosection,'uniqueDir');
       uniqueDir= 'x';
   else
       uniqueDir= fmdl.show_pseudosection.uniqueDir;
   end
   [elec_posn,elecNumber]= electrodesPosition(fmdl);
   if isfield(fmdl,'misc') && isfield(fmdl.misc,'elec_posn')
       elec_posn= fmdl.misc.elec_posn;
   end
   me= (mod(elec_posn,1));
   me(me==0)= 1;
   precision_position= min(me(:));
   xposition_elec= reshape(elec_posn(elecNumber,1),[],4); 
   yposition_elec= reshape(elec_posn(elecNumber,2),[],4); 
   zposition_elec= reshape(elec_posn(elecNumber,3),[],4); 
   zposition_elec= zposition_elec-min(zposition_elec(:));
   
    electrodesLocation= xposition_elec;
   if isfield(fmdl.show_pseudosection,'elecsUsed')
       elecsUsed= fmdl.show_pseudosection.elecsUsed;
   else
       elecsUsed= 1:4;
   end
   
   location12x= mean(xposition_elec(:,elecsUsed(1:2)),2);
   location34x= mean(xposition_elec(:,elecsUsed(3:4)),2);
   location1234x= mean(xposition_elec(:,elecsUsed(1:4)),2);
   location12y= mean(yposition_elec(:,elecsUsed(1:2)),2);
   location34y= mean(yposition_elec(:,elecsUsed(3:4)),2);
   location1234y= mean(yposition_elec(:,elecsUsed(1:4)),2);
   
   depth12x= diff(xposition_elec(:,elecsUsed(1:2)),1,2);
   depth34x= diff(xposition_elec(:,elecsUsed(3:4)),1,2);
   depth12y= diff(yposition_elec(:,elecsUsed(1:2)),1,2);
   depth34y= diff(yposition_elec(:,elecsUsed(3:4)),1,2);
   
   x= cell(length(locVar),1);
   for i= 1:length(locVar)
        x{i}= eval(locVar{i});
    end
   location= locFun(x);
   
   y= cell(length(depthVar),1);
   for i= 1:length(depthVar)
        y{i}= eval(depthVar{i});
    end
   depth= depthFun(y);
   [LU,is,js]= unique(location);
   dL= min(diff(LU));
   for j= 1:length(is)
       locationis= location(js==j);
       depthis= depth(js==j);
       if length(unique(depthis)) ~=length(depthis)
           [depthisU,isl,jsl]= unique(depthis);
           for k= 1:length(isl)
               if length(find(jsl==k))>1
                   shift= (linspace(-dL/2,dL/2,length(find(jsl==k))))';
                   locationis(jsl==k)= locationis(jsl==k)+shift;
               end
           end
           location(js==j)= locationis;
       end
   end
   
   P= [depth location];
   Pu= unique(P,'rows','stable');
   if size(P,1)~= size(Pu,1)
       disp(' ')
       disp('P not unique !!!')
       disp(' ')
       keyboard
   end
      
   scatter(location,depth,fmdl.misc.sizepoint,data,'filled','MarkerEdgeColor','k');
   hold on; 
   xlabel('Distance (m)','fontsize',fs,'fontname','Times','interpreter','latex');
   ylabel('Distance (m)','fontsize',fs,'fontname','Times','interpreter','latex')
    axis equal; axis tight;
   bx= get(gca,'xlim'); by= get(gca,'ylim'); 
   set(gca,'fontsize',fs,'fontname','Times','dataAspectRatio',[1 (by(2)-by(1))/(bx(2)-bx(1))*3 1])
end


function [depth,location,location34,electrodesLocation]= plotPseudoSectionProfileSpecificDown(fmdl,data)
   fs= fmdl.show_pseudosection.fs;
   [elec_posn,elecNumber]= electrodesPosition(fmdl);
   xposition_elec= reshape(elec_posn(elecNumber,1),[],4); 
   yposition_elec= reshape(elec_posn(elecNumber,2),[],4); 
   zposition_elec= reshape(elec_posn(elecNumber,3),[],4); 
   zposition_elec= zposition_elec-min(zposition_elec(:));
   if isfield(fmdl.show_pseudosection,'minx') && isfield(fmdl.show_pseudosection,'maxx') && ...
       isfield(fmdl.show_pseudosection,'miny') && isfield(fmdl.show_pseudosection,'maxy')
       minx= fmdl.show_pseudosection.minx;
       maxx= fmdl.show_pseudosection.maxx;
       miny= fmdl.show_pseudosection.miny;
       maxy= fmdl.show_pseudosection.maxy;
   else
       xposition_elec= xposition_elec-min(xposition_elec(:));
       yposition_elec= yposition_elec-min(yposition_elec(:));
       minx= min(xposition_elec(:));
       maxx= max(xposition_elec(:));
       miny= min(yposition_elec(:));
       maxy= max(yposition_elec(:));
   end

   if isfield(fmdl.show_pseudosection,'xToDeduce') && isfield(fmdl.show_pseudosection,'yToDeduce')
       if strcmp(fmdl.show_pseudosection.xToDeduce,'minx')
           xToDeduce= minx;
       elseif strcmp(fmdl.show_pseudosection.xToDeduce,'maxx')
           xToDeduce= maxx;
       end
       if strcmp(fmdl.show_pseudosection.yToDeduce,'miny')
           yToDeduce= miny;
       elseif strcmp(fmdl.show_pseudosection.yToDeduce,'maxy')
           yToDeduce= maxy;
       end
   else
       xToDeduce= minx;
       yToDeduce= miny;
   end
   rposition_elec= sqrt((xposition_elec-xToDeduce).^2 + ...
       (yposition_elec-yToDeduce).^2);
   
   if isfield(fmdl.show_pseudosection,'elecsUsed')
       elecsUsed= fmdl.show_pseudosection.elecsUsed;
   else
       elecsUsed= 3:4;
   end
   
   if isfield(fmdl.show_pseudosection,'dirAxis')
       dirAxis= fmdl.show_pseudosection.dirAxis;
   else
       dirAxis= 'X';
   end
   
   if isfield(fmdl.show_pseudosection,'depthRatio')
       depthRatio= fmdl.show_pseudosection.depthRatio;
   else
       depthRatio= 3;
   end
   
    if isfield(fmdl.show_pseudosection,'depthPrecision')
       depthPrecision= fmdl.show_pseudosection.depthPrecision;
   else
       depthPrecision= 1;
    end
    
    if isfield(fmdl.show_pseudosection,'depthLevel')
       depthLevel= fmdl.show_pseudosection.depthLevel;
   else
       depthLevel= 0;
    end
   
   switch dirAxis
       case 'X'
           location= mean(xposition_elec(:,elecsUsed),2);
           depth12= -abs(xposition_elec(:,elecsUsed(1))-xposition_elec(:,elecsUsed(2)))/depthRatio;
           depth34= -abs(xposition_elec(:,elecsUsed(3))-xposition_elec(:,elecsUsed(4)))/depthRatio;
           electrodesLocation= unique(round(xposition_elec/depthPrecision)*depthPrecision);
       case 'Y'
           location= mean(yposition_elec(:,elecsUsed),2);
           depth= -abs(yposition_elec(:,elecsUsed(1))-yposition_elec(:,elecsUsed(2)))/depthRatio;
           electrodesLocation= unique(round(yposition_elec/depthPrecision)*depthPrecision);
       case 'Z'
           location= mean(zposition_elec(:,elecsUsed),2);
           depth= -abs(zposition_elec(:,elecsUsed(1))-zposition_elec(:,elecsUsed(2)))/depthRatio;
           electrodesLocation= unique(round(zposition_elec/depthPrecision)*depthPrecision);
       case 'R'
           location= mean(rposition_elec(:,elecsUsed),2);
           location12= mean(rposition_elec(:,elecsUsed(1:2)),2);
           location34= mean(rposition_elec(:,elecsUsed(3:4)),2);
           depth12= abs(rposition_elec(:,elecsUsed(1))-rposition_elec(:,elecsUsed(2)));
           depth34= abs(rposition_elec(:,elecsUsed(3))-rposition_elec(:,elecsUsed(4)));
           electrodesLocation= rposition_elec; % unique(round(rposition_elec/depthPrecision)*depthPrecision);
   end
   
   precision= floor(min(diff(unique(rposition_elec(:)))));
  if precision==0
      precision= ceil(min(diff(unique(rposition_elec(:)))));
  end
   
   depth12= round(depth12/precision)*precision;
   depth34= round(depth34/precision)*precision;
   location12= round(location12/precision)*precision;
   location34= round(location34/precision)*precision;

   miSdepth12= min(diff(unique(depth12)));
   if ~isempty(miSdepth12)
       depth34reS= (depth34-min(depth34))/max(depth34)*miSdepth12;
   else
       depth34reS= 0;
   end
   depth= depth12 +depth34reS;
   depth= round(-(depth)/depthRatio/precision)*precision;
     
   [depthU,is,js]= unique(depth);
   for j= 1:length(is)
       changeDepth= 0;
       clear locationis depthis
       location34is= location34(js==j);
       location12is= location12(js==j);
       depthis= depth(js==j);
       datais= data(js==j);
       if length(unique(location34is))==length(location34is)
           locationis= location34is;
       else
           [location34isU,isl,jsl]= unique(location34is);
           for k= 1:length(isl)
               shift= (location12is(jsl==k)-mean(location12is(jsl==k)))/(min(depth34));
               locationis(jsl==k)= round((location34is(jsl==k)+shift)/precision)*precision;
           end
           if length(unique(locationis))~=length(locationis)
               depthisShifted= depth(js==j)*0;
               for k= 1:length(isl)
                   shift= ((location12is(jsl==k)-mean(location12is(jsl==k)))/(min(depth34)*2));
                   if length(shift)==2 && round(shift(1)/precision)*precision==0
                       shift= ((location12is(jsl==k)-mean(location12is(jsl==k)))/(min(depth34)*2));
                   end
                   if length(shift)==3
                       xtremeShift= min([abs(min(shift)) abs(max(shift))]);
                       if round(xtremeShift(1)/precision)*precision==0
                           xtremeShift= max([abs(min(shift)) abs(max(shift))]);
                       end
                       shift= [-xtremeShift 0 xtremeShift]';
                   end
                   depthisShifted(jsl==k)= depthis(jsl==k)+shift;
                   changeDepth= 1;
               end
           end
       end
        location(js==j)= locationis;
        if changeDepth
        depth(js==j)= depthisShifted;
        end
   end
   
   depth= depth+depthLevel;
   scatter(location,depth,fmdl.misc.sizepoint,data,'filled','MarkerEdgeColor','k');
   hold on;
   xlabel('Distance (m)','fontsize',fs,'fontname','Times','interpreter','latex');
   ylabel('Pseudo-depth (m)','fontsize',fs,'fontname','Times','interpreter','latex')
    axis equal; axis tight;
   bx= get(gca,'xlim'); by= get(gca,'ylim'); 
   set(gca,'fontsize',fs,'fontname','Times','dataAspectRatio',[1 (by(2)-by(1))/(bx(2)-bx(1))*3 1])
end

function [depth,location,L,d]= plotPseudoSectionProfileSpecificUp(fmdl,data)
   fs= fmdl.show_pseudosection.fs;
   [elec_posn,elecNumber]= electrodesPosition(fmdl);
   xposition_elec= reshape(elec_posn(elecNumber,1),[],4); 
   yposition_elec= reshape(elec_posn(elecNumber,2),[],4); 
   zposition_elec= reshape(elec_posn(elecNumber,3),[],4); zposition_elec= zposition_elec-min(zposition_elec(:));
   if isfield(fmdl.show_pseudosection,'minx') && isfield(fmdl.show_pseudosection,'maxx') && ...
       isfield(fmdl.show_pseudosection,'miny') && isfield(fmdl.show_pseudosection,'maxy')
       minx= fmdl.show_pseudosection.minx;
       maxx= fmdl.show_pseudosection.maxx;
       miny= fmdl.show_pseudosection.miny;
       maxy= fmdl.show_pseudosection.maxy;
   else
       xposition_elec= xposition_elec-min(xposition_elec(:));
       yposition_elec= yposition_elec-min(yposition_elec(:));
       minx= min(xposition_elec(:));
       maxx= max(xposition_elec(:));
       miny= min(yposition_elec(:));
       maxy= max(yposition_elec(:));
   end
   
   if isfield(fmdl.show_pseudosection,'xToDeduce') && isfield(fmdl.show_pseudosection,'yToDeduce')
       if strcmp(fmdl.show_pseudosection.xToDeduce,'minx')
           xToDeduce= minx;
       elseif strcmp(fmdl.show_pseudosection.xToDeduce,'maxx')
           xToDeduce= maxx;
       end
       if strcmp(fmdl.show_pseudosection.yToDeduce,'miny')
           yToDeduce= miny;
       elseif strcmp(fmdl.show_pseudosection.yToDeduce,'maxy')
           yToDeduce= maxy;
       end
   else
       xToDeduce= minx;
       yToDeduce= miny;
   end
   rposition_elec= sqrt((xposition_elec-xToDeduce).^2 + ...
       (yposition_elec-yToDeduce).^2);
     
%    rposition_elec= sqrt((xposition_elec-min(xposition_elec(:))).^2 + ...
%        (yposition_elec-max(yposition_elec(:))).^2);
   
   if isfield(fmdl.show_pseudosection,'elecsUsed')
       elecsUsed= fmdl.show_pseudosection.elecsUsed;
   else
       elecsUsed= 3:4;
   end
   
   if isfield(fmdl.show_pseudosection,'dirAxis')
       dirAxis= fmdl.show_pseudosection.dirAxis;
   else
       dirAxis= 'X';
   end
   
   if isfield(fmdl.show_pseudosection,'depthRatio')
       depthRatio= fmdl.show_pseudosection.depthRatio;
   else
       depthRatio= 3;
   end
   
    if isfield(fmdl.show_pseudosection,'depthPrecision')
       depthPrecision= fmdl.show_pseudosection.depthPrecision;
   else
       depthPrecision= 1;
    end
    
    if isfield(fmdl.show_pseudosection,'depthLevel')
       depthLevel= fmdl.show_pseudosection.depthLevel;
   else
       depthLevel= 0;
    end
   
   switch dirAxis
       case 'X'
           location= mean(xposition_elec(:,elecsUsed),2);
           depth= abs(xposition_elec(:,elecsUsed(1))-xposition_elec(:,elecsUsed(2)))/depthRatio;
           electrodesLocation= unique(round(xposition_elec/depthPrecision)*depthPrecision);
       case 'Y'
           location= mean(yposition_elec(:,elecsUsed),2);
           depth= abs(yposition_elec(:,elecsUsed(1))-yposition_elec(:,elecsUsed(2)))/depthRatio;
           electrodesLocation= unique(round(yposition_elec/depthPrecision)*depthPrecision);
       case 'Z'
           location= mean(zposition_elec(:,elecsUsed),2);
           depth= abs(zposition_elec(:,elecsUsed(1))-zposition_elec(:,elecsUsed(2)))/depthRatio;
           electrodesLocation= unique(round(zposition_elec/depthPrecision)*depthPrecision);
       case 'R'
           location= mean(rposition_elec(:,elecsUsed),2);%keyboard
           depth= abs(rposition_elec(:,elecsUsed(1))-rposition_elec(:,elecsUsed(2)))/depthRatio;
           electrodesLocation= unique(round(rposition_elec/depthPrecision)*depthPrecision);
   end
   L= (rposition_elec(:,2)-rposition_elec(:,1))/2;
   d= (rposition_elec(:,4)-rposition_elec(:,3))/2;
   
   
   depth= depth+depthLevel; 
   P= depth+1i*location;
   [Pu,iu,ju]= unique(round(P/depthPrecision)*depthPrecision);
   
   if length(Pu) < length(P)
       lu= location(iu);
       zu= depth(iu);
       du= iu*0;
       for i= 1:length(Pu)
           du(i)= mean(data(ju==i));
       end
       if length(fmdl.misc.sizepoint)>length(Pu)
           su= iu*0;
%            misu= iu*0;
%            masu= iu*0;
           for i= 1:length(Pu)
           su(i)= min(fmdl.misc.sizepoint(ju==i));
%            misu(i)= min(fmdl.misc.sizepoint(ju==i));
%            masu(i)= max(fmdl.misc.sizepoint(ju==i));
           end
           fmdl.misc.sizepoint= su;
%            UU= [su misu masu masu-misu];
       end
   else
       lu= location;
       zu= depth;
       du= data;
   end
   
   [Puu]= unique(lu+1i*zu);
   if length(Puu) ~= length(lu)
   keyboard    
   end
   scatter(lu,zu,fmdl.misc.sizepoint,(du),'filled','MarkerEdgeColor','k');
%    hold on; plot(electrodesLocation,electrodesLocation*0,'x','Color',[0 0.5 0])
   xlabel('Distance (m)','fontsize',fs,'fontname','Times','interpreter','latex');
   ylabel('Pseudo-depth (m)','fontsize',fs,'fontname','Times','interpreter','latex')
%     axis equal; axis tight;
   set(gca,'fontsize',fs,'fontname','Times')
end



function [zps,xps]= plotPseudoSectionProfileDown(fmdl,data)
   fs= fmdl.show_pseudosection.fs;
   
   [elec_posn,elecNumber]= electrodesPosition(fmdl);
   xposition_elec= reshape(elec_posn(elecNumber,1),[],4);
   xps= mean(xposition_elec,2);
   
   AB= abs(xposition_elec(:,2)-xposition_elec(:,1));
   MN= abs(xposition_elec(:,4)-xposition_elec(:,3));
   AN= abs(xposition_elec(:,4)-xposition_elec(:,1));
   BM= abs(xposition_elec(:,3)-xposition_elec(:,2));
   AM= abs(xposition_elec(:,3)-xposition_elec(:,1));
   BN= abs(xposition_elec(:,4)-xposition_elec(:,2));
   
   a= MN; %max([AB MN AN BM AM BN ],[],2);
   zps= -a/3;
   
%    a= abs(elecNumber(:,1)-elecNumber(:,2));
%    de= elec_posn(1,1)-elec_posn(2,1);
%    
%    % Identiy reciprocal data(elecNumber(:,1)-elecNumber(:,2)) > abs(elecNumber(:,3)-elecNumber(:,4))
%    R= find(abs(elecNumber(:,1)-elecNumber(:,2)) < abs(elecNumber(:,3)-elecNumber(:,4)));
%    
%    if ~isempty(R)
%        xps(R)= (elec_posn(elecNumber(R,3),4)+elec_posn(elecNumber(R,4),4))/2;
%        a(R)= abs(elecNumber(R,3)-elecNumber(R,4));
%    end
     
   P= zps+1i*xps;
   [Pu,iu,ju]= unique(P);
   
   if length(Pu) < length(P)
       xu= xps(iu);
       zu= zps(iu);
       du= iu*0;
       for i= 1:length(Pu)
           du(i)= mean(data(ju==i));
       end
   else
       xu= xps;
       zu= zps;
       du= data;
   end
       
   scatter(xu,zu,fmdl.misc.sizepoint,(du),'filled','MarkerEdgeColor','k');
   xlabel('Distance (m)','fontsize',fs,'fontname','Times','interpreter','latex');
   ylabel('Pseudo-depth (m)','fontsize',fs,'fontname','Times','interpreter','latex');
   axis equal; axis tight;
   bx= get(gca,'xlim'); by= get(gca,'ylim');
   set(gca,'fontsize',fs,'fontname','Times','dataAspectRatio',[1 (by(2)-by(1))/(bx(2)-bx(1))*3 1])

end

function [xps,zps]= plotPseudoSectionProfileVert(fmdl,data)
   fs= fmdl.show_pseudosection.fs;
   
   [elec_posn,elecNumber]= electrodesPosition(fmdl);
  
   zposition_elec= reshape(elec_posn(elecNumber,3),[],4);
   zps= mean(zposition_elec,2);
   
   AB= abs(zposition_elec(:,2)-zposition_elec(:,1));
   MN= abs(zposition_elec(:,4)-zposition_elec(:,3));
   AN= abs(zposition_elec(:,4)-zposition_elec(:,1));
   BM= abs(zposition_elec(:,3)-zposition_elec(:,2));
   AM= abs(zposition_elec(:,3)-zposition_elec(:,1));
   BN= abs(zposition_elec(:,4)-zposition_elec(:,2));
   
   a= max([AB MN AN BM AM BN ],[],2);
   xps= a/3;
      
   P= zps+1i*xps;
   [Pu,iu,ju]= unique(P);
   
   xu= xps(iu);
   zu= zps(iu);
   du= iu*0;
   for i= 1:length(Pu)
       du(i)= mean(data(ju==i));
   end
   
   scatter(xu,zu,fmdl.misc.sizepoint,(du),'filled','MarkerEdgeColor','k'); 
   xlabel('Pseudo distance (m)','fontsize',fs,'fontname','Times','interpreter','latex');
   if zps(1)<0
       ylabel('Depth (m)','fontsize',fs,'fontname','Times','interpreter','latex')
   else
       ylabel('Height (m)','fontsize',fs,'fontname','Times','interpreter','latex')
   end
   axis equal; axis tight;
   set(gca,'fontsize',fs,'fontname','Times')
end

function [r_point,th_point]= plotPseudoSectionCircularIn(fmdl,data)
fs= fmdl.show_pseudosection.fs;
[elec_posn,elecNumber] = electrodesPosition(fmdl);
[A,th_point,r,xc,yc] = polarPosition(elecNumber,elec_posn);
a= max(A,[],2);
r_point= a/2;
r_point= (r-r_point/pi);%*9*pi*r/10+pi*r/10;

[x_point,y_point]= pol2cart(th_point,r_point);
x_point= x_point+xc; y_point= y_point+yc;

P= x_point+1i*y_point;
[Pu,iu,ju]= unique(P);
   
xu= x_point(iu);
zu= y_point(iu);
du= iu*0;
for i= 1:length(Pu)
    du(i)= mean(data(ju==i));
end
   
scatter(xu,zu,fmdl.misc.sizepoint,(du),'filled','MarkerEdgeColor','k');  
xlabel('X (m)','fontsize',fs,'fontname','Times','interpreter','latex');
ylabel('Y (m)','fontsize',fs,'fontname','Times','interpreter','latex')
axis equal; axis tight;
set(gca,'fontsize',fs,'fontname','Times')
end

function [dMN,th_point]= plotPseudoSectionCircularInVolcano(fmdl,data)
fs= fmdl.show_pseudosection.fs;
[elec_posn,elecNumber] = electrodesPosition(fmdl);
[xposition_elec,yposition_elec,zposition_elec] = electrodesPositionABMN(elecNumber,elec_posn);
dMN= sqrt((xposition_elec(:,4)-xposition_elec(:,3)).^2 + ...
    (yposition_elec(:,4)-yposition_elec(:,3)).^2 + (zposition_elec(:,4)-zposition_elec(:,3)).^2);

[A,r_bary,th_point,r,xc,yc,TH,R] = polarPosition(elecNumber,elec_posn);
[th_bary,r_bary,th_bary34,r_bary34,r,xc,yc,barycenter_x34,barycenter_y34] = bary34(elecNumber,elec_posn);
% a= A(:,2);
% r_point= a;
r= 400; 
r_point= (r-r_bary)/(max(r-r_bary)-min(r-r_bary))*(r-5);%*9*pi*r/10+pi*r/10;
r_point= r_point-min(r_point)+5;
[x_point,y_point]= pol2cart(th_bary34,r_point);
x_point= x_point+xc; y_point= y_point+yc;
P= x_point+1i*y_point;
[Pu,iu,ju]= unique(P);
xu= x_point(iu); zu= y_point(iu); du= iu*0;
for i= 1:length(Pu);     du(i)= mean(data(ju==i)); end
figure; scatter(xu,zu,fmdl.misc.sizepoint,(du),'filled','MarkerEdgeColor','k');  
text(xu(iu==1),zu(iu==1),'1')
hold on; plot(xposition_elec,yposition_elec,'kp')
xlabel('X (m)','fontsize',fs,'fontname','Times','interpreter','latex');
ylabel('Y (m)','fontsize',fs,'fontname','Times','interpreter','latex')
axis equal; axis tight;
set(gca,'fontsize',fs,'fontname','Times')
end

function [r_point,th_point]= plotPseudoSectionCircularOut(fmdl,data)
fs= fmdl.show_pseudosection.fs;
[elec_posn,elecNumber] = electrodesPosition(fmdl);
[A,r_bary,th_point,r,xc,yc] = polarPosition(elecNumber,elec_posn);
a= max(A,[],2);
r_point= a/2;
r_point= (r+r_point);
[x_point,y_point]= pol2cart(th_point,r_point+r);
x_point= x_point+xc; y_point= y_point+yc;

P= x_point+1i*y_point;
[Pu,iu,ju]= unique(P);
   
xu= x_point(iu);
zu= y_point(iu);
du= iu*0;
for i= 1:length(Pu)
    du(i)= mean(data(ju==i));
end
   
scatter(xu,zu,fmdl.misc.sizepoint,(du),'filled','MarkerEdgeColor','k');  colorbar
xlabel('X (m)','fontsize',fs,'fontname','Times','interpreter','latex');
ylabel('Y (m)','fontsize',fs,'fontname','Times','interpreter','latex')
axis equal; axis tight;
set(gca,'fontsize',fs,'fontname','Times')
end


function [elec_posn,elecNumber] = electrodesPosition(fmdl)
    stim= fmdl.stimulation;
    stimulationMatrix= [];
    for i = 1:length(stim);
        nmp= size(stim(i).meas_pattern, 1);
        [idxIN,idxJN]= find(stim(i).meas_pattern<0);
        [idxIP,idxJP]= find(stim(i).meas_pattern>0);
        stimulationMatrix= [stimulationMatrix; [ 0*(1:nmp)'+find(stim(i).stim_pattern<0) ...
            0*(1:nmp)'+find(stim(i).stim_pattern>0)] idxJN(idxIN) idxJP(idxIP)];
    end
%    A=  stimulationMatrix(:,1);
%    A=  stimulationMatrix(:,1);
%    A=  stimulationMatrix(:,1);
%    A=  stimulationMatrix(:,1);
%    
%    A= zeros(length(fmdl.stimulation),1);
%    B= zeros(length(fmdl.stimulation),1);
%    M= zeros(length(fmdl.stimulation),1);
%    N= zeros(length(fmdl.stimulation),1);
%    for i=1:length(fmdl.stimulation)
%        A(i)= find(fmdl.stimulation(1,i).stim_pattern<0);
%        B(i)= find(fmdl.stimulation(1,i).stim_pattern>0);
%        M(i)= find(fmdl.stimulation(1,i).meas_pattern<0);
%        N(i)= find(fmdl.stimulation(1,i).meas_pattern>0);
%    end
%    elecNumber= [A B M N];
   
   elecNumber= stimulationMatrix;
   
   elec_posn= zeros(length(fmdl.electrode),size(fmdl.nodes,2));
   for i=1:length(fmdl.electrode)
   elec_posn(i,:)= mean(fmdl.nodes(fmdl.electrode(1,i).nodes,:),1);
   end
   
end

function [xposition_elec,yposition_elec,zposition_elec] = electrodesPositionABMN(elecNumber,elec_posn)
xposition_elec= reshape(elec_posn(elecNumber,1),[],4);
yposition_elec= reshape(elec_posn(elecNumber,2),[],4);
zposition_elec= reshape(elec_posn(elecNumber,3),[],4);
end

function [A,r_bary,th_bary,r,xc,yc,TH,R] = polarPosition(elecNumber,elec_posn)
TH= elecNumber*0;
R= elecNumber*0;

xposition_elec= reshape(elec_posn(elecNumber,1),[],4);
yposition_elec= reshape(elec_posn(elecNumber,2),[],4);
rx= (max(xposition_elec(:))-min(xposition_elec(:)))/2;
ry= (max(yposition_elec(:))-min(yposition_elec(:)))/2;
rmean= (rx+ry)/2;

xc= mean(elec_posn(:,1)); yc= mean(elec_posn(:,2)); r = rmean;

barycenter_x= mean(xposition_elec,2);
barycenter_y= mean(yposition_elec,2);

[th_bary,r_bary]= cart2pol(barycenter_x-xc,barycenter_y-yc);

for i= 1:4
    [TH(:,i),R(:,i)]= cart2pol(elec_posn(elecNumber(:,i),1)-xc,elec_posn(elecNumber(:,i),2)-yc);
end 

Arc_length_AB= (TH(:,2)-TH(:,1));
idx_discontinuity= find(TH(:,2)>TH(:,1) & TH(:,2)<0);
Arc_length_AB(idx_discontinuity)= (-TH(idx_discontinuity,1)+TH(idx_discontinuity,2));
idx_discontinuity= find(TH(:,1)>TH(:,2));
Arc_length_AB(idx_discontinuity)= (TH(idx_discontinuity,1)-TH(idx_discontinuity,2));
idx_discontinuity= find(TH(:,1)>TH(:,2) & TH(:,1)<0);
Arc_length_AB(idx_discontinuity)= (-TH(idx_discontinuity,2)+TH(idx_discontinuity,1));
Arc_length_AB(Arc_length_AB>=pi+0.005)= 2*pi-Arc_length_AB(Arc_length_AB>=pi+0.005);
Arc_length_AB= Arc_length_AB.*(R(:,1)+R(:,2))/2;

Arc_length_MN= (TH(:,4)-TH(:,3));
idx_discontinuity= find((TH(:,4)>TH(:,3) & TH(:,4)<0));
Arc_length_MN(idx_discontinuity)= -TH(idx_discontinuity,3)+TH(idx_discontinuity,4);
idx_discontinuity= find((TH(:,4)<TH(:,3) & TH(:,4)<0));
Arc_length_MN(idx_discontinuity)= +TH(idx_discontinuity,3)-TH(idx_discontinuity,4);
idx_discontinuity= find((TH(:,3)>TH(:,4)) & TH(:,3)>0);
Arc_length_MN(idx_discontinuity)= TH(idx_discontinuity,3)-TH(idx_discontinuity,4);
% idx_discontinuity= find((TH(:,3)>TH(:,4)) &  TH(:,3)<0);
% Arc_length_MN(idx_discontinuity)= -TH(idx_discontinuity,3)+TH(idx_discontinuity,4);
Arc_length_MN(Arc_length_MN>=pi+0.005)= 2*pi-Arc_length_MN(Arc_length_MN>=pi+0.005);
Arc_length_MN= Arc_length_MN.*(R(:,3)+R(:,4))/2;

Arc_length_AM= (TH(:,3)-TH(:,1));
idx_discontinuity= find((TH(:,3)>TH(:,1) & TH(:,3)<0));
Arc_length_AM(idx_discontinuity)= -TH(idx_discontinuity,1)+TH(idx_discontinuity,3);
idx_discontinuity= find((TH(:,1)> TH(:,3)));
Arc_length_AM(idx_discontinuity)= TH(idx_discontinuity,1)-TH(idx_discontinuity,3);
idx_discontinuity= find((TH(:,1)> TH(:,3) & TH(:,1) <0));
Arc_length_AM(idx_discontinuity)= -TH(idx_discontinuity,1)+TH(idx_discontinuity,3);
Arc_length_AM(Arc_length_AM>=pi+0.005)= 2*pi-Arc_length_AM(Arc_length_AM>=pi+0.005);
Arc_length_AM= Arc_length_AM.*(R(:,1)+R(:,3))/2;

Arc_length_AN= (TH(:,4)-TH(:,1));
idx_discontinuity= find((TH(:,4)>TH(:,1)<0 & TH(:,4)<0));
Arc_length_AN(idx_discontinuity)= -TH(idx_discontinuity,1)+TH(idx_discontinuity,4);
idx_discontinuity= find((TH(:,1)>TH(:,4)));
Arc_length_AN(idx_discontinuity)= TH(idx_discontinuity,1)-TH(idx_discontinuity,4);
idx_discontinuity= find((TH(:,1)>TH(:,4) & TH(:,1)<0));
Arc_length_AN(idx_discontinuity)= TH(idx_discontinuity,4)-TH(idx_discontinuity,1);
Arc_length_AN(Arc_length_AN>=pi+0.005)= 2*pi-Arc_length_AN(Arc_length_AN>=pi+0.005);
Arc_length_AN= Arc_length_AN.*(R(:,1)+R(:,4))/2;

Arc_length_NB= (TH(:,4)-TH(:,2));
idx_discontinuity= find((TH(:,4)>TH(:,2)));
Arc_length_NB(idx_discontinuity)= -TH(idx_discontinuity,2)+TH(idx_discontinuity,4);
idx_discontinuity= find((TH(:,2)> TH(:,4)));
Arc_length_NB(idx_discontinuity)= -TH(idx_discontinuity,4)+TH(idx_discontinuity,2);
idx_discontinuity= find((TH(:,2)> TH(:,4) &  TH(:,4)<0));
Arc_length_NB(idx_discontinuity)= -TH(idx_discontinuity,4)+TH(idx_discontinuity,2);
Arc_length_NB(Arc_length_NB>=pi+0.005)= 2*pi-Arc_length_NB(Arc_length_NB>=pi+0.005);
Arc_length_NB= Arc_length_NB.*(R(:,2)+R(:,4))/2;

Arc_length_BM= (TH(:,3)-TH(:,2));
idx_discontinuity= find((TH(:,3)>TH(:,2) & TH(:,2)<0));
Arc_length_BM(idx_discontinuity)= -TH(idx_discontinuity,2)+TH(idx_discontinuity,3);
idx_discontinuity= find((TH(:,2)> TH(:,3)));
Arc_length_BM(idx_discontinuity)= -TH(idx_discontinuity,3)+TH(idx_discontinuity,2);
idx_discontinuity= find((TH(:,2)> TH(:,3)& TH(:,3)<0));
Arc_length_BM(idx_discontinuity)= TH(idx_discontinuity,3)-TH(idx_discontinuity,2);
Arc_length_BM(Arc_length_BM>=pi+0.005)= 2*pi-Arc_length_BM(Arc_length_BM>=pi+0.005);
Arc_length_BM= Arc_length_BM.*(R(:,2)+R(:,3))/2;

A= [Arc_length_AB Arc_length_MN Arc_length_AM Arc_length_AN Arc_length_NB Arc_length_BM];
end



function [th_bary,r_bary,th_bary34,r_bary34,r,xc,yc,barycenter_x34,barycenter_y34] = bary34(elecNumber,elec_posn)
xposition_elec= reshape(elec_posn(elecNumber,1),[],4);
yposition_elec= reshape(elec_posn(elecNumber,2),[],4);
rx= (max(xposition_elec(:))-min(xposition_elec(:)))/2;
ry= (max(yposition_elec(:))-min(yposition_elec(:)))/2;
rmean= (rx+ry)/2;

xc= mean(elec_posn(:,1)); yc= mean(elec_posn(:,2)); r = rmean;

barycenter_x34= mean(xposition_elec(:,3:4),2);
barycenter_y34= mean(yposition_elec(:,3:4),2);
[th_bary34,r_bary34]= cart2pol(barycenter_x34-xc,barycenter_y34-yc);

barycenter_x= mean(xposition_elec,2);
barycenter_y= mean(yposition_elec,2);
[th_bary,r_bary]= cart2pol(barycenter_x-xc,barycenter_y-yc);


end

function do_unit_test
shape_str = ['solid top    = plane(0,0,0;0,1,0);\n' ...
             'solid mainobj= top and orthobrick(-100,-200,-100;410,10,100) -maxh=20.0;\n'];
e0 = linspace(0,310,64)';
elec_pos = [e0,0*e0,0*e0,1+0*e0,0*e0,0*e0];
elec_shape= [0.1,0.1,1];
elec_obj = 'top';
fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

spacing= [1 1 1 2 3 3 4 4 5 6 6 7 8 8 9 10 10 11 12 12 13 14 14 15 16 17];
multiples= [1 2 3 2 1 5/3 1 2  1 1 7/6 1 1 10/8 1 1 12/10 1 1 13/12 1 1 15/14 1 1 1];
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
