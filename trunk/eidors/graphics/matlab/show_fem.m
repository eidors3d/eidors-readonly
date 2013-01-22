function hh=show_fem( mdl, options)
% SHOW_FEM: show the EIDORS3D finite element model
% hh=show_fem( mdl, options )
% mdl is a EIDORS3D 'model' or 'image' structure
% hh= handle to the plotted model
%
% options may be specified by a list
%
% options specifies a set of options
%   options(1) => show colourbar
%   options(2) => show numbering on electrodes
%   options(3) => number elements (==1) or nodes (==2);
%
% for detailed control of colours, use
%    img.calc_colours."parameter" = value
% see help for calc_colours.
%
% for control of colourbar, use img.calc_colours.cb_shrink_move
%
% parameters
%     fwd_model.show_fem.linecolour
%
% to change line properties:
%      hh=show_fem(...); set(hh,'EdgeColor',[0,0,1];

% (C) 2005-2011 Andy Adler. License: GPL version 2 or version 3
% $Id$

% TESTS
switch nargin
   case 0;
     error('Insufficient parameters for show_fem');
   case 1;
     if isstr(mdl) && strcmp(mdl,'UNIT_TEST'); do_unit_test; return; end
     options = [];
end
[img,mdl,opts] = proc_params( mdl, options );

if ~ishold
    cla;
    axis auto;
end


switch opts.dims
   case 2;    hh= show_2d(img,mdl,opts);
   case 3;    hh= show_3d(img,mdl,opts);
   otherwise; error('model is not 2D or 3D');
end

if opts.show_elem_numbering
   xyzc= interp_mesh(mdl);
   placenumbers(xyzc, 7, [0,0,0],'none');
end
if opts.show_node_numbering
   xyzc= mdl.nodes;
   placenumbers(xyzc, 7, [0.0,0,0.5],[1,1,1]);
end

if nargout == 0; clear hh; end

if ~ishold
   axis equal;
   axis tight;
end

function placenumbers(xyzc, fontsize, colour, bgcolour)
   xyzc= xyzc * eye(size(xyzc,2),3); %convert to 3D
   for i= 1:size(xyzc,1)
      text(xyzc(i,1),xyzc(i,2), xyzc(i,3), num2str(i), ...
            'HorizontalAlignment','center', ...
            'FontSize',fontsize,'Color',colour,'BackgroundColor',bgcolour);
   end

function [img,mdl,opts] = proc_params( mdl, options );

   opts.do_colourbar=0;
   opts.number_electrodes=0;
   opts.show_numbering  =0;
   opts.show_elem_numbering = 0;
   opts.show_node_numbering = 0;
   if nargin >=2
       % fill in default options
       optionstr= zeros(1,100);
       optionstr(1:length(options)) = options;

       opts.do_colourbar=      optionstr(1);
       opts.number_electrodes= optionstr(2);
       switch optionstr(3)
          case 1; opts.show_elem_numbering = 1;
          case 2; opts.show_node_numbering = 1;
          case 3; opts.show_elem_numbering = 1;
                  opts.show_node_numbering = 1;
       end
   end

   
   % if we have an only img input, then define mdl
   if strcmp( mdl(1).type , 'image' )
      img= mdl;
      mdl= img(1).fwd_model;
   else 
      img = [];
   end
   
   opts.transparency_thresh = calc_colours('transparency_thresh');
   try
       opts.transparency_thresh = img.calc_colours.transparency_thresh;
   end
   
   opts.dims = size(mdl.nodes,2);

% 2D Case
function hh= show_2d(img,mdl,opts)
   hax= gca;
   pax= get(hax,'position');
   if ~isempty(img)
      colours= calc_colours(img, [] );
   else
      colours= [1,1,1]; % white elements if no image
   end
   hh= show_2d_fem( mdl, colours );
   show_electrodes_2d(mdl, opts.number_electrodes);

   set(hax,'position', pax);
   view(0, 90); axis('xy'); grid('off');

% IN MATLAB 7 WE NEED TO RERUN THIS BECAUSE STUPID STUPID
% MATLAB WILL RESET THE COLOURBAR EVERY TIME WE RUN PATCH!!!
   if exist('img','var') && opts.do_colourbar;
      colours= calc_colours(img, [], opts.do_colourbar);
      % Matlab is so weird. It puts the first colorbar in the wrong place
      %   sometimes ...  (tested with 6.5 and with 7.8)
      %   The trick is to never try to move it on the first go
      %   OR we reset it and then replace it. STUPID STUPID

     % Here's the magic trick I found. Force a drawnow,
     %     then delete and recreate

     % Because of a change in matlab colorbar somewhere around
     % 2012, none of this compensation code works any more. We
     % don't really understand what to do anymore, but this
     % seems best, for now ...
    % if ~exist('OCTAVE_VERSION')
    %    drawnow; colorbar('delete');
    %    colours= calc_colours(img, [], opts.do_colourbar);
    % end
   end

   

% 3D Case
function hh= show_3d(img,mdl,opts)
   hh= show_3d_fem( mdl );

   if ~isempty(img)
       elem_data = get_img_data(img);
       show_inhomogeneities( elem_data , mdl, img, opts);
       if opts.do_colourbar
           calc_colours(img, [], opts.do_colourbar);
       end
   end
   if size(mdl.elems,2) == 3
      show_electrodes_surf(mdl, opts.number_electrodes);
   else
      show_electrodes_3d(mdl, opts.number_electrodes);
   end

function show_electrodes_2d(mdl, number_electrodes)
    if ~isfield(mdl,'electrode'); return; end

    ee= get_boundary( mdl );
    ctr_x= mean(mdl.nodes(:,1));
    ctr_y= mean(mdl.nodes(:,2));

% scale away from model

for e=1:length(mdl.electrode)
    if isfield(mdl.electrode(e),'pos') && ~isfield(mdl.electrode(e),'nodes')
        vx = mdl.electrode(e).pos(:,1) - ctr_x;
        vy = mdl.electrode(e).pos(:,2) - ctr_y;
        idx = 1:length(vx);
    else
        elec_nodes= mdl.electrode(e).nodes;
        
        S= 1.00;
        vx= (mdl.nodes(elec_nodes,1) - ctr_x)*S;
        vy= (mdl.nodes(elec_nodes,2) - ctr_y)*S;
        % sort nodes around the model (to avoid crossed lines)
        [jnk,idx] = sort(unwrap(atan2( vy, vx )));
    end
        
    ecolour = electr_colour( e );
    if numel(vx) == 1
       % Point Electrode Models: put a circle around the node
       line(vx(idx)+ctr_x,vy(idx)+ctr_y,  ...
            'LineWidth', 2, 'LineStyle','-','Color', ecolour, ...
            'Marker','o','MarkerSize', 6,'MarkerEdgeColor',ecolour);
    else
       % Complete/Shunt Electrode Models (multiple nodes per electrode)
       %  put a line along the edges that form the electrode
       line(vx(idx)+ctr_x,vy(idx)+ctr_y,  ...
            'LineWidth', 3, 'LineStyle','-','Color', ecolour, ...
            'Marker','none','MarkerSize', 6,'MarkerEdgeColor',ecolour);
    end
    if number_electrodes
       S= 1.05;
       vx= (mdl.nodes(elec_nodes,1) - ctr_x)*S;
       vy= (mdl.nodes(elec_nodes,2) - ctr_y)*S;
       switch number_electrodes
          case {1 true}
             txt = num2str(e);
          case 2
             try, txt = mdl.electrode(e).label; end
       end
       hh= text(mean(vx)+ctr_x, mean(vy)+ctr_y, txt);
       set(hh, 'HorizontalAlignment','center', 'FontWeight','bold');
    end
end

function show_electrodes_surf(mdl, number_electrodes)
    if ~isfield(mdl,'electrode'); return; end

    ee= get_boundary( mdl );
    ctr_x= mean(mdl.nodes(:,1));
    ctr_y= mean(mdl.nodes(:,2));
    ctr_z= mean(mdl.nodes(:,3));
% scale away from model

for e=1:length(mdl.electrode)
    if isfield(mdl.electrode(e),'pos') && ~isfield(mdl.electrode(e),'nodes')
        vx = mdl.electrode(e).pos(:,1) - ctr_x;
        vy = mdl.electrode(e).pos(:,2) - ctr_y;
        vz = mdl.electrode(e).pos(:,3) - ctr_z;
        idx = 1:length(vx);
    else
        elec_nodes= mdl.electrode(e).nodes;
        
        S= 1.00;
        vx= (mdl.nodes(elec_nodes,1) - ctr_x)*S;
        vy= (mdl.nodes(elec_nodes,2) - ctr_y)*S;
        vz= (mdl.nodes(elec_nodes,3) - ctr_z)*S;
        % sort nodes around the model (to avoid crossed lines)
        % TODO: figure out what to do in different directions
        [jnk,idx] = sort(unwrap(atan2( vy, vx )));
    end
    ecolour = electr_colour( e );
    if numel(vx) == 1
       % Point Electrode Models: put a circle around the node
       line(vx(idx)+ctr_x,vy(idx)+ctr_y, vz(idx)+ctr_z,  ...
            'LineWidth', 2, 'LineStyle','-','Color', ecolour, ...
            'Marker','o','MarkerSize', 6,'MarkerEdgeColor',ecolour);
    else
       % Complete/Shunt Electrode Models (multiple nodes per electrode)
       %  put a line along the edges that form the electrode
       line(vx(idx)+ctr_x,vy(idx)+ctr_y, vz(idx)+ctr_z, ...
            'LineWidth', 3, 'LineStyle','-','Color', ecolour, ...
            'Marker','none','MarkerSize', 6,'MarkerEdgeColor',ecolour);
    end
    if number_electrodes
       S= 1.05;
       vx= (mdl.nodes(elec_nodes,1) - ctr_x)*S;
       vy= (mdl.nodes(elec_nodes,2) - ctr_y)*S;
       vz= (mdl.nodes(elec_nodes,3) - ctr_z)*S;
       switch number_electrodes
          case {1 true}
            txt = num2str(e);
          case 2
             try, txt = mdl.electrode(e).label; end
       end
       hh= text(mean(vx)+ctr_x, mean(vy)+ctr_y,mean(vz)+ctr_z,txt);
       set(hh, 'HorizontalAlignment','center', 'FontWeight','bold');
    end
end

function show_electrodes_3d(mdl, number_electrodes)
% show electrode positions on model
if ~isfield(mdl,'electrode'); return; end

ee= get_boundary( mdl );
for e=1:length(mdl.electrode)
    colour= electr_colour( e);
    
    if isfield(mdl.electrode(e),'pos') && ~isfield(mdl.electrode(e),'nodes')
        show_electrodes_surf(mdl, number_electrodes);
        return
    end
    elec_nodes= mdl.electrode(e).nodes;
    
    
    if length(elec_nodes) == 1  % point electrode model
        vtx= mdl.nodes(elec_nodes,:);
        line(vtx(1),vtx(2),vtx(3), ...
            'Marker','o','MarkerSize',12, ...
            'MarkerFaceColor',colour, 'MarkerEdgeColor', colour);
        if number_electrodes
            hh= text(vtx(1),vtx(2),vtx(3), num2str(e));
            set(hh, 'HorizontalAlignment','center', 'FontWeight','bold');
        end
    else
        % find elems on boundary attached to this electrode
        nn=ones(size(ee,1),1)*mdl.electrode(e).nodes(:)';
        oo=ones(1,size(nn,2));
        ec=zeros(size(ee));
        for i=1:size(ec,2);
           ec(:,i) = any( ee(:,i*oo)==nn, 2);
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
            switch number_electrodes
               case {1 true}
                  txt = num2str(e);
               case 2
                  try, txt = mdl.electrode(e).label; end
            end
            hh=text( mean(el_nodes(:,1)), ...
                     mean(el_nodes(:,2)), ...
                     mean(el_nodes(:,3)), txt );
            set(hh,'FontWeight','bold','Color',[.6,0,0]);
        end
    end
end

function show_inhomogeneities( elem_data, mdl, img, opt)
% show
if size(elem_data,2)>1
   eidors_msg('warning: show_fem only shows first image',1);
end
repaint_inho(elem_data(:,1), 'use_global' , mdl.nodes, mdl.elems, ...
    opt.transparency_thresh, img); 
if ~exist('OCTAVE_VERSION');
camlight('left');
lighting('none'); % lighting doesn't help much
end

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

function hh= show_3d_fem( mdl, options )
   ee= get_boundary( mdl );
   hh= trimesh(ee, mdl.nodes(:,1), ...
                   mdl.nodes(:,2), ...
                   mdl.nodes(:,3));
   set(hh, 'EdgeColor', [0,0,0]);
   axis('image');
   hidden('off');

function hh=show_2d_fem( mdl, colours, imgno )
  szcolr = size(colours);
  if nargin<3;
      imgno = 1;
      if szcolr(1:2)>1; 
          eidors_msg('warning: show_fem only shows first image',1);
      end
  end

  if szcolr(1:2) == [1,3]  % no reconstruction  - just use white
     hh= patch('Faces',mdl.elems,'Vertices',mdl.nodes, 'facecolor',colours);
  elseif size(colours,1) == num_elems(mdl);
     colours = squeeze(colours(:,imgno,:));

% THE OCTAVE GNUPLOT INTERPRETER HAS BUGS - recommend FLTK
%    ver= eidors_obj('interpreter_version');
% 
%    if ver.isoctave && size(colours,2)==1 %% Octave bug on CDataMapping
%       colours = colours /(2*calc_colours('mapped_colour')+1);     
%    end

     hh= patch('Faces',mdl.elems,'Vertices',mdl.nodes, 'facecolor','flat', ...
               'facevertexcdata',colours,'CDataMapping','direct'); 
  elseif size(colours,1) == num_nodes(mdl);
     colours = squeeze(colours(:,imgno,:));
     % 'interp' not supported in octave
     hh= patch('Faces',mdl.elems,'Vertices',mdl.nodes, 'facecolor','interp', ...
               'facevertexcdata',colours,'CDataMapping','direct'); 
  else
    eidors_msg('warning: image elements and mesh do not match. Showing grey',1);
    colours= [1,1,1]/2; % Grey to show we're not sure
    hh= patch('Faces',mdl.elems,'Vertices',mdl.nodes, 'facecolor',colours);
  end

  set(hh, 'FaceLighting','none', 'CDataMapping', 'direct' );

function hh=show_2d_fem_oldmatlab( mdl, colours, imgno )
  szcolr = size(colours);
  if nargin<3;
      imgno = 1;
      if szcolr(1:2)>1;  eidors_msg('warning: show_fem only shows first image',1); end
  end
  
% el_pos= avg_electrode_posn( mdl );  
% plot_2d_mesh(mdl.nodes', mdl.elems', el_pos', .95, [0,0,0]);

  S= 1; %.95; % shrink factor
  elem= mdl.elems'; e= size(elem,2);
  node= mdl.nodes'; n= size(node,2);
  Xs=zeros(3,e);
  Xs(:)=mdl.nodes(elem(:),1);
  Xs= S*Xs+ (1-S)*ones(3,1)*mean(Xs);
  Ys=zeros(3,e); Ys(:)=mdl.nodes(elem(:),2);
  Ys= S*Ys+ (1-S)*ones(3,1)*mean(Ys);
  Zs = zeros(size(Xs));
  if szcolr(1:2) == [1,3]  % no reconstruction  - just use white
     hh= patch(Xs,Ys,zeros(3,e),colours);
  elseif size(colours,1) == e % defined on elems
% THE STUPID MATLAB 7 WILL RESET THE COLOURBAR WHENEVER YOU DO A PATCH. DAMN.
     colours = permute(colours(:,imgno,:),[2,1,3]);
  elseif size(colours,1) == n  % multiple images
     colours = colours(elem(:), imgno, :);
     colours = reshape( colours, size(elem,1), size(elem,2), []);
  else
    eidors_msg('warning: image elements and mesh do not match. Showing grey',1);
    colours= [1,1,1]/2; % Grey to show we're not sure
  end

  hh= patch(Xs,Ys,Zs,colours);
  set(hh, 'FaceLighting','none', 'CDataMapping', 'direct' );
  % FOR NODE RGB MATLAB SCREWS UP THE COLOURS FOR US (ONCE AGAIN). THERE MUST
  % BE SOME KIND OF SECRET FLAG@@@

  max_x= max(mdl.nodes(:,1)); min_x= min(mdl.nodes(:,1));
  max_y= max(mdl.nodes(:,2)); min_y= min(mdl.nodes(:,2));
  axis([ mean([max_x,min_x]) + 0.55*[-1,1]*(max_x-min_x), ...
         mean([max_y,min_y]) + 0.55*[-1,1]*(max_y-min_y) ]);



function old_code_3d_plot
  domesnum= 0;
  donodenum= 0;
  dotext=0;

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
   mes= zeros( n_elec, mdl_dim(mdl) )
   for i= 1:n_elec
      e_nodes=  mdl.electrode(i).nodes;
      mes(i,:)= mean( nodes(e_nodes,:) , 1 );
   end

function colour= electr_colour( e);
    if e==1;
       colour = [0,.7,0]; % light green electrode #1
    elseif e==2
       colour = [0,.5,0]; % mid-green electrode #2
    else
       colour = [0,.3,0]; % dark green
    end

function ee= get_boundary( mdl )
   if isfield(mdl,'boundary')
       ee= mdl.boundary;
   else
       % calc and cache boundary
       ee = find_boundary( mdl.elems );
   end


% TESTS:
function do_unit_test
     ver= eidors_obj('interpreter_version');
   clf

   img=calc_jacobian_bkgnd(mk_common_model('a2c0',8)); 
   img.elem_data=rand(size(img.fwd_model.elems,1),1);
   subplot(3,4,1); show_fem(img.fwd_model,[0,0,1]) 
   title('regular mesh numbered');

if ~ver.isoctave 
   imgn = rmfield(img,'elem_data');
   imgn.node_data=rand(size(img.fwd_model.nodes,1),1);
   subplot(3,4,9); show_fem(imgn) 
   title('interpolated node colours');
end

   img2(1) = img; img2(2) = img;
   subplot(3,4,2); show_fem(img,[1]);
   title('colours with legend');
   subplot(3,4,3); show_fem(img2,[0,1]);
   title('colours with legend');
   img.calc_colours.mapped_colour = 0; % USE RGB colours
   subplot(3,4,4); show_fem(img,[0,1,1]);
   title('RGB colours');
   subplot(3,4,4); show_fem(img);
   title('RGB colours');

   img.elem_data = [1:10];
   subplot(3,4,12);show_fem(img); %Should show grey
   title('error -> show grey');

if ~ver.isoctave
   imgn.calc_colours.mapped_colour = 0; % USE RGB colours
   subplot(3,4,10);show_fem(imgn,[0,1]) 
   title('interpolated node colours');


   subplot(3,4,11);hh=show_fem(imgn); set(hh,'EdgeColor',[0,0,1]);
   title('with edge colours');

end

   img3=calc_jacobian_bkgnd(mk_common_model('n3r2',[16,2]));
   img3.elem_data= randn(828,1);                       
   subplot(3,4,5); show_fem(img3.fwd_model) 
   subplot(3,4,6); show_fem(img3,[1])
   subplot(3,4,7); show_fem(img3,[1,1])
   subplot(3,4,8); show_fem(img3,[1,1,1])

