function show_fem( mdl, options )
% show_fem( mdl, options )
% SHOW_FEM: show the EIDORS3D finite element model
% mdl is a EIDORS3D 'model' or 'image' structure
%
% options specifies a set of options
%   options(1) => show colourbar
%   options(2) => number electrodes
%   options(3) => colourbar scale
%
% set ref_level for conductivities with
% calc_colours('ref_level', ref_level)

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: show_fem.m,v 1.38 2006-07-27 01:03:23 camilgomez Exp $

if exist('OCTAVE_VERSION');
   warning('show_fem does not support octave');
   return
end

if nargin <= 1
end

do_colourbar=0;
number_electrodes=0;
if nargin >=2
    optionstr= zeros(1,100);
    optionstr(1:length(options)) = options;
    do_colourbar=      optionstr(1);
    number_electrodes= optionstr(2);
    if optionstr(3)~=0
        scl_colourbar= optionstr(3);
    else
        scl_colourbar=[];
    end
end

% if we have an only img input, then define mdl
name= mdl.name;
if strcmp( mdl.type , 'image' )
   img= mdl;
   mdl= img.fwd_model;
end

cla;
set(gcf, 'Name', name(:)');

if size(mdl.nodes,2)==2
   hax= gca;
   pax= get(hax,'position');
   if exist('scl_colourbar')
       cscale = scl_colourbar;
   else
   cscale= [];
   end
   if exist('img');
      colours= calc_colours(img.elem_data, cscale, do_colourbar);
   else
      colours= [1,1,1]; % white elements if no image
   end
   cla;
   show_2d_fem( mdl, colours );
   show_electrodes_2d(mdl);

   set(hax,'position', pax);
   view(0, 90); axis('xy'); grid('off');
elseif size(mdl.nodes,2)==3
   cscale= [];
   show_3d_fem( mdl );

   if exist('img')
       show_inhomogeneities( img.elem_data , mdl);
       if do_colourbar
           calc_colours(img.elem_data, cscale, do_colourbar);
       end
   end

   show_electrodes_3d(mdl, number_electrodes);
else
   error(['model is not 2D or 3D']);
end

function show_electrodes_2d(mdl)
    if ~isfield(mdl,'electrode'); return; end

    ee= get_boundary( mdl );
    ctr_x= mean(mdl.nodes(:,1));
    ctr_y= mean(mdl.nodes(:,2));

% scale away from model
    S= 1.02;

for e=1:length(mdl.electrode)
    elec_nodes= mdl.electrode(e).nodes;

    vx= (mdl.nodes(elec_nodes,1) - ctr_x)*S;
    vy= (mdl.nodes(elec_nodes,2) - ctr_y)*S;
    % sort nodes around the model (to avoid crossed lines)
    [jnk,idx] = sort(atan2( vy, vx ));
    ecolour = electr_colour( e );
    line(vx(idx)+ctr_x,vy(idx)+ctr_y,  ...
         'LineWidth', 2, 'Color', [1 0 0], ...
         'Marker','.','MarkerSize',20,'MarkerEdgeColor',ecolour);
end

function show_electrodes_3d(mdl, number_electrodes);
% show electrode positions on model
if ~isfield(mdl,'electrode'); return; end

ee= get_boundary( mdl );
for e=1:length(mdl.electrode)
    elec_nodes= mdl.electrode(e).nodes;

    colour= electr_colour( e);

    if length(elec_nodes) == 1  % point electrode model
        vtx= mdl.nodes(elec_nodes,:);
        line(vtx(1),vtx(2),vtx(3), ...
            'Marker','h','MarkerSize',24, ...
            'MarkerFaceColor',colour, 'MarkerEdgeColor', colour);
        if number_electrodes
            text(vtx(1),vtx(2),vtx(3), num2str(e));
        end
    else
        % find elems on boundary attached to this electrode
        nn=ones(size(ee,1),1)*mdl.electrode(e).nodes;
        oo=ones(1,size(nn,2));
        ec=zeros(size(ee));
        for i=1:size(ec,2);
           ec(:,i) = any( (ee(:,i)*oo==nn)' )';
        end
        sels= find(all(ec'));

        for u=1:length(sels)
            ee= get_boundary( mdl );
            paint_electrodes(sels(u), ee, ...
                             mdl.nodes, colour, ...
                             number_electrodes);
        end
        if number_electrodes
            el_nodes= mdl.nodes(unique(mdl.boundary(sels,:)),:);
            hh=text( mean(el_nodes(:,1)), ...
                     mean(el_nodes(:,2)), ...
                     mean(el_nodes(:,3)), num2str(e) );
            set(hh,'FontWeight','bold','Color',[.6,0,0]);
        end
    end
end

function show_inhomogeneities( elem_data, mdl)
% show
hold('on');
repaint_inho(elem_data, 0, mdl.nodes, mdl.elems); 
camlight('left');
lighting('none'); % lighting doesn't help much
hold('off');

function paint_electrodes(sel,srf,vtx, colour, show_num);
%function paint_electrodes(sel,srf,vtx);
%
% plots the electrodes red at the boundaries.
%
% sel = The index of the electrode faces in the srf matrix
%       sel can be created by set_electrodes.m 
% srf = the boundary faces (triangles)
% vtx = The vertices matrix.

l = srf(sel,1); m = srf(sel,2); n = srf(sel,3);

Xs = [vtx(l,1);vtx(m,1);vtx(n,1)];
Ys = [vtx(l,2);vtx(m,2);vtx(n,2)];
Zs = [vtx(l,3);vtx(m,3);vtx(n,3)];

h=patch(Xs,Ys,Zs, colour);
% need 'direct' otherwise colourmap is screwed up
set(h, 'FaceLighting','none', 'CDataMapping', 'direct' );

function show_3d_fem( mdl, options )
   ee= get_boundary( mdl );
   hh= trimesh(ee, mdl.nodes(:,1), ...
                   mdl.nodes(:,2), ...
                   mdl.nodes(:,3));
   set(hh, 'EdgeColor', [0,0,0]);
   axis('image');
   set(gcf,'Colormap',[0 0 0]);
   hidden('off');

function show_2d_fem( mdl, colours )
  
  el_pos= avg_electrode_posn( mdl );  
% plot_2d_mesh(mdl.nodes', mdl.elems', el_pos', .95, [0,0,0])

  S= 1; %.95; % shrink factor
  elem= mdl.elems'; e= size(elem,2);
  Xs=zeros(3,e);
  Xs(:)=mdl.nodes(elem(:),1);
  Xs= S*Xs+ (1-S)*ones(3,1)*mean(Xs);
  Ys=zeros(3,e);
  Ys(:)=mdl.nodes(elem(:),2);
  Ys= S*Ys+ (1-S)*ones(3,1)*mean(Ys);
  h= patch(Xs,Ys,zeros(3,e),colours);
  set(h, 'FaceLighting','none', 'CDataMapping', 'direct' );

  max_x= max(mdl.nodes(:,1)); min_x= min(mdl.nodes(:,1));
  max_y= max(mdl.nodes(:,2)); min_y= min(mdl.nodes(:,2));
  axis([ mean([max_x,min_x]) + 0.55*[-1,1]*(max_x-min_x), ...
         mean([max_y,min_y]) + 0.55*[-1,1]*(max_y-min_y) ]);



function old_code_3d_plot
  domesnum= 0;
  donodenum= 0;
  dotext=0;
  if exist('OCTAVE_VERSION')
    gset('nokey');
    dotext=0;
    domesnum= 0;
    donodenum= 0;
  end

  xxx=zeros(4,e); xxx(:)=NODE(1,ELEM(:));
  xxx= S*xxx+ (1-S)*ones(4,1)*mean(xxx);
  yyy=zeros(4,e); yyy(:)=NODE(2,ELEM(:));
  yyy= S*yyy+ (1-S)*ones(4,1)*mean(yyy);
  zzz=zeros(4,e); zzz(:)=NODE(3,ELEM(:));
  zzz= S*zzz+ (1-S)*ones(4,1)*mean(zzz);
  plot3( xxx([1 2 3 1 4 2 3 4],:), ...
         yyy([1 2 3 1 4 2 3 4],:), ...
         zzz([1 2 3 1 4 2 3 4],:),'k' );
  hold on;
  xy=NODE(1:2,MES(1,:));
  plot3(arrow*xy,arrow*[0 1;-1 0]*xy,ones(8,1)*NODE(3,MES(1,:)),'b') 
  hold off;
  axis([ [-1.1 1.1]*max(NODE(1,:)) [-1.1 1.1]*max(NODE(2,:)) ...
         [-1.1 1.1]*max(NODE(3,:))  ])


function plot_2d_mesh(NODE,ELEM,el_pos, S, options)
%  if options(1) -> do mesh num
%  if options(2) -> do node num
%  if options(3) -> do text
%  S=.95; % shrink simplices by this factor

  e=size(ELEM,2);
  d=size(ELEM,1);
  R= zeros(1,e);


  arrow= [1.02 0;1.06 .05;1.06 .02;1.1 .02; ...
          1.1 -.02; 1.06 -.02;1.06 -.05;1.02 0];
  % arrow= [1.06 0;1.18 .15;1.18 .06;1.3 .06; ...
  %          1.3 -.06; 1.18 -.06;1.18 -.15;1.06 0];
  xxx=zeros(3,e); xxx(:)=NODE(1,ELEM(:));
  xxx= S*xxx+ (1-S)*ones(3,1)*mean(xxx);
  xxx= [xxx;xxx(1,:)];

  yyy=zeros(3,e); yyy(:)=NODE(2,ELEM(:));
  yyy= S*yyy+ (1-S)*ones(3,1)*mean(yyy);
  yyy= [yyy;yyy(1,:)];

  if isempty(el_pos)
      plot(xxx,yyy,'b');
  elseif exist('OCTAVE_VERSION')
      plot(arrow*xy,arrow*[0 1;-1 0]*xy,'m');
      hold('on')
      idx= find(R==0);
      if ~isempty(idx); plot(xxx(:,idx),yyy(:,idx),'c'); end
      idx= find(R>0);
      if ~isempty(idx); plot(xxx(:,idx),yyy(:,idx),'r'); end
      idx= find(R<0);
      if ~isempty(idx); plot(xxx(:,idx),yyy(:,idx),'b'); end
      hold('off')
  else
      xy= el_pos;
      hh=plot([xxx;xxx(1,:)],[yyy;yyy(1,:)],'b', ...
              arrow*xy,arrow*[0 1;-1 0]*xy,'r');
      set(hh(find(R>0)),'Color',[1 0 0],'LineWidth',2);
      set(hh(find(R<0)),'Color',[0 0 1],'LineWidth',2);
  end

  if options(1) %domesnum
    mesnum= reshape(sprintf(' %-2d',[0:15]),3,16)';
    text(NODE(1,MES(1,:))*1.08,NODE(2,MES(1,:))*1.08,mesnum, ...
         'Color','green','HorizontalAlignment','center');
  end %if domesnum

  if options(2) %donodenum
    nodenum= reshape(sprintf('%3d',1:length(NODE)),3,length(NODE))';
    text(NODE(1,:),NODE(2,:),nodenum, ...
         'Color','yellow','HorizontalAlignment','center','FontSize',14);
  end %if domesnum

  if options(3) %dotext
    numeros= reshape(sprintf('%4d',[1:e]),4,e)';
%   decal= ( 2-0.2*floor(log10([1:e]')))*sqrt(axis*axis'/4)/100;
    decal= .02;
    xcoor=mean(NODE(2*ELEM-1))';
    ycoor=mean(NODE(2*ELEM))';
    text(xcoor-decal,ycoor,numeros,'FontSize',8, ...
         'HorizontalAlignment','center');

  end  % if nargin~=0
  axis([ [-1.1 1.1]*max(NODE(1,:)) [-1.1 1.1]*max(NODE(2,:)) ])

function  mes= avg_electrode_posn( mdl )
   if ~isfield(mdl,'electrode'); mes=[]; return; end
   n_elec= length( mdl.electrode );
   nodes = mdl.nodes;
   n_dims= size(nodes,2);
   mes= zeros( n_elec, n_dims );
   for i= 1:n_elec
      e_nodes=  mdl.electrode(i).nodes;
      mes(i,:)= mean( nodes(e_nodes,:) , 1 );
   end

function colour= electr_colour( e);
    if e==1;
       colour = [0,.6,0]; % light green
    else
       colour = [0,.3,0]; % dark green
    end

function ee= get_boundary( mdl )
   if isfield(mdl,'boundary')
       ee= mdl.boundary;
   else
       % calc and cache boundary
       pp = np_fwd_parameters( mdl );
       ee = pp.srf;
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003, A. Adler 2004
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
