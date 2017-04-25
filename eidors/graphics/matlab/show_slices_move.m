function show_slices_move( img, move, move_scale )
% SHOW_SLICES_MOVE   Shows planar slices of a 3D FEM with movement vectors
% if electrodes are visible on the slice.
% Args:     img  - eidors_obj type image
%           move - change in  position vectors for nodes [x,y,z] after movement 

% $Id$

if ischar(img) && strcmp(img,'UNIT_TEST'); do_unit_test; return; end


num_levs = 3;

mdl = img.fwd_model;
elecs = length(mdl.electrode);
pos = zeros(elecs,3);
if exist('move') && size(move,1) == size(mdl.nodes,1)
    move_nodes = move;
    move= zeros(elecs,3);
end
for e = 1:elecs;
    if isfield(mdl.electrode(e),'pos')
       pos(e,:) = mean(mdl.electrode(e).pos);
    else
       elec_nodes = mdl.electrode(e).nodes;
       pos(e,:) = mean(mdl.nodes(elec_nodes,:),1);
    end
    if exist('move_nodes');
        move(e,:) = mean( move_nodes(elec_nodes,:), 1);
    end
end

elec_zmax = max(pos(:,3));
elec_zmin = min(pos(:,3));

levels = inf*ones(num_levs,3);
levels(:,3) = [elec_zmax; elec_zmax-(elec_zmax - elec_zmin)/2; elec_zmin];
levels(:,4) = ones(num_levs,1);
levels(:,5) = (1:num_levs)';

show_slices( img, levels );

isize = calc_colours('npoints');
toplayer = pos(:,3) < levels(2,3); % origin is top-left corner of image
xofs = .5 + isize*(.5);
xscale = (100/104)*isize/2;
yofs = .5 + isize*(.5 + (num_levs-1)*toplayer);
yscale = -(100/104)*isize/2;
vx = pos(:,1) * xscale + xofs;
vy = pos(:,2) * yscale + yofs;
ecolour = [0,.3,0];
hh= line(vx, vy);
set(hh, 'LineStyle','none','Marker','.', ...
    'MarkerSize',10,'MarkerEdgeColor',ecolour);

pp = fwd_model_parameters( mdl );
if nargin == 1;
    move = [];
end
if length(img.elem_data) > pp.n_elem
    move = reshape( ...
        img.elem_data( pp.n_elem+(1:pp.n_elec*pp.n_dims) ), ...
        pp.n_elec, pp.n_dims);
end

if ~isempty( move )
    nodes = img.fwd_model.nodes;

    % zero mean movement for each electrode row; this is not quite legit
    idx = 1:16;
    move(idx,:) = move(idx,:)- ones(16,1)*mean(move(idx,:));
    if length(img.fwd_model.electrode) == 32
        idx = 17:32;
        move(idx,:) = move(idx,:)- ones(16,1)*mean(move(idx,:));
    end
    hold on;
    if nargin < 3
        move_scale = 20*isize;
    end
    hh = working_quiver( vx, vy, move_scale*move(:,1), ...
        - move_scale*move(:,2), 0 );
    set(hh,'Color', [0,.3,0], 'LineWidth', 2, 'Clipping', 'off');
    hold off;
end

function hh= working_quiver( varargin )
% WORKING_QUIVER   Matlab has made a completely imcompatible
% quiver function which you can't call properly with different
% versions of matlab.

% Unfortunately, the idea of using 'v6' is now a warning.
% We can't do too many version checks, giving up!

hh = quiver( varargin{:} );

%   v = version;
%   octave = exist('OCTAVE_VERSION')
%   if octave
%       hh = quiver( varargin{:} );
%   else
%       hh = quiver('v6', varargin{:} );
%   end


function do_unit_test
   subplot(221);
   img = mk_image( mk_common_model('n3r2',[16,2]) );
   show_slices_move(img);

   subplot(222);
   for i=1:length(img.fwd_model.electrode)
     e_node = img.fwd_model.electrode(i).nodes;
     move(i,:) = 0.01*mean( img.fwd_model.nodes(e_node,:),1);
   end
   show_slices_move(img,move);

   subplot(223);
   show_slices_move(img,move, 1000);
