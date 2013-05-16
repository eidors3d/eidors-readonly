function [hf,hh] = show_fem_move( img, move, scale, options )
% SHOW_FEM_MOVE   Plot EIT finite element model (FEM) and movement 
%    vectors.
% Args: img     - eidors_obj of type 'image'
%       move    - FEM node coordinate matrix of movement vectors (optional)
%       scale   - factor to scale movement arrows (optional)
%       options - options array passed on to show_fem()

% (C) 2005 Andy Adler. License GPL v2 or v3
% $Id$

if isstr(img) && strcmp(img,'UNIT_TEST'); do_unit_test; return; end

% Check for single argument call
if nargin == 1
    move = [];
    scale = 20;
    options = [0,0,[]];
elseif nargin == 2
    scale = 20;
    options = [0,0,[]];
elseif nargin == 3
    options = [0,0,[]];
end

% Extract forward model parameters
fwdp = fwd_model_parameters( img.fwd_model );
try
    fwdp.n_elem = size(img.fwd_model.coarse2fine,2);
end
% Verify if img is partitioned by conductivity and move submatrices
if length(img.elem_data) > fwdp.n_elem
    move = reshape( ...
        img.elem_data( fwdp.n_elem+(1:fwdp.n_elec*fwdp.n_dims) ), ...
        fwdp.n_elec, fwdp.n_dims);
    img.elem_data = img.elem_data(1:fwdp.n_elem);    
end

% Plot FEM with conductivity elements with or without colourbars
hf = show_fem(img, options); % Show colourbar

% Plot movement vectors on electrodes
if ~isempty(move)
    % e_nodes is the average position of each electrodes nodes
    e_nodes = [];
    electr= img.fwd_model.electrode;
    nodes = img.fwd_model.nodes;
    switch size(move,1)
      case length(electr)
        for i=1:length(electr)
           e_nodes(i,:) = mean(nodes(electr(i).nodes,:),1);
        end
    
      case num_nodes(img)
        % Keep only electrode node movement coordinates
        e_nodes = [electr(:).nodes];
        if size(move,1) == fwdp.n_node;
            move = move(e_nodes,:);
        end
        e_nodes = nodes(e_nodes,:);

      otherwise;
        error('movement vector doesn''t match model');
    end
    move = move- ones(size(move,1),1)*mean(move);
    
    % Render movement arrows for each electrode
    hold on;
    if nargin < 3
        scale = 20;
    end
    hh = working_quiver(e_nodes, scale*move);
    set(hh,'Color', [0,.3,0], 'LineWidth', 2, 'Clipping', 'off');
    hold off;
end
% Format output figure
axis('off'); 
axis('image'); 
%axis(1.3*[-1,1,-1,1]); % let it take it's own space

function hh= working_quiver( nn,mm )
% WORKING_QUIVER   Matlab has made a completely imcompatible
% quiver function which you can't call properly with different
% versions of matlab.
%
% Last I checked, the V7 version of quiver was horrible, so
%  we use the v6 one.

% TODO: Write a new, fixed quiver function that can do 3D

ver= eidors_obj('interpreter_version');
if ver.isoctave || ver.ver<7;
    hh = quiver( nn(:,1),nn(:,2), mm(:,1),mm(:,2),0);
else
    warning('off','MATLAB:quiver:DeprecatedV6Argument');
    hh = quiver('v6', nn(:,1),nn(:,2), mm(:,1),mm(:,2),0);
end

function do_unit_test;
   subplot(231);
   img = mk_image(mk_common_model('a2c2',8));
   img.elem_data = [img.elem_data;.01*randn(16,1)];
   show_fem_move(img);
   title('move at electrodes');

   subplot(232);
   img = mk_image(mk_common_model('a2c2',8));
   show_fem_move( img, img.fwd_model.nodes*[-1,0;0,1], .2);
   title('move at electrodes (via all move)');

   subplot(233);
   img = mk_image(mk_common_model('a2C2',8));
   img.elem_data = [img.elem_data;.01*randn(16,1)];
   show_fem_move(img);
   title('move at electrodes ctrs');

   subplot(234);
   img = mk_image(mk_common_model('a2C0',8));
   img.elem_data(1:2) = 1.1;
   show_fem_move( img, img.fwd_model.nodes*[-1,0;0,1], .2);
   title('move at electrodes nodes');

   subplot(235);
   img = mk_image(mk_common_model('n3r2',[16,2]));
   img.elem_data(400)= 0.9;
   img.elem_data = [img.elem_data;.01*randn(3*32,1)];
   show_fem_move(img); view([-14,62]);
   title('BUG:move at electrodes ctrs');

   subplot(236);
   img = mk_image(mk_common_model('n3r2',[16,2]));
   move = img.fwd_model.nodes*diag([-1,1,1]);
   show_fem_move( img, move, .2); view([-14,62]);
   title('BUG:move at electrodes ctrs');
