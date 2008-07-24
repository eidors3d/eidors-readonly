function show_fem_move( img, move, scale, options )
% SHOW_FEM_MOVE   Plot EIT finite element model (FEM) and movement 
%    vectors.
% Args: img     - eidors_obj of type 'image'
%       move    - FEM node coordinate matrix of movement vectors (optional)
%       scale   - factor to scale movement arrows (optional)
%       options - options array passed on to show_fem()

% $Id: show_fem_move.m,v 1.6 2007-08-30 03:37:32 aadler Exp $

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
fwdp = aa_fwd_parameters( img.fwd_model );

% Verify if img is partitioned by conductivity and move submatrices
if length(img.elem_data) > fwdp.n_elem
    move = reshape( ...
        img.elem_data( fwdp.n_elem+(1:fwdp.n_elec*fwdp.n_dims) ), ...
        fwdp.n_elec, fwdp.n_dims);
    img.elem_data = img.elem_data(1:fwdp.n_elem);    
end

% Plot FEM with conductivity elements with or without colourbars
show_fem(img, options); % Show colourbar

% Plot movement vectors on electrodes
if ~isempty(move)
    % Group all electrode node indicies into a single vector
    e_nodes = cell2mat({img.fwd_model.electrode(:).nodes});
    
    % Keep only electrode node movement coordinates
    if size(move,1) == fwdp.n_node;
        move = move(e_nodes,:);
    end
    move = move- ones(size(move,1),1)*mean(move);
    
    % Render movement arrows for each electrode
    hold on;
    if nargin < 3
        scale = 20;
    end
    nodes = img.fwd_model.nodes;
    hh = working_quiver(nodes(e_nodes,1), nodes(e_nodes,2), ...
        scale*move(:,1), scale*move(:,2), 0 );
    set(hh,'Color', [0,.3,0], 'LineWidth', 2, 'Clipping', 'off');
    hold off;
end
% Format output figure
axis('off'); 
axis('image'); 
axis(1.3*[-1,1,-1,1]);

function hh= working_quiver( varargin )
% WORKING_QUIVER   Matlab has made a completely imcompatible
% quiver function which you can't call properly with different
% versions of matlab.

v = version;
octave = exist('OCTAVE_VERSION') | str2num(v(1)) < 7;
if octave
    hh = quiver( varargin{:} );
else
    hh = quiver('v6', varargin{:} );
end
