function mdl = ng_mk_2d_model(varargin)
%NG_MG_2D_MODELS create a 2D mesh with Netgen via the in2d interface
% mdl = ng_mk_2d_model(shape)
% mdl = ng_mk_2d_model(shape, elec_pos)
% mdl = ng_mk_2d_model(shape, elec_pos, elec_shape)
%
% SHAPE can be:
%  - xy (Nx2)             : a counter- clockwise list of points in 2D 
%                           defining the outer contour
%  - {xy, xy1, xy2, ...}  : allows specifying additional counter-clockwise 
%                           loops  xy1, xy2, etc, which represent holes in  
%                           the bigger contour xy contour
%  - {..., maxsz}         : specifies maximum element size of the mesh.
%                           If absent, mesh paremeters are controlled by
%                           the ng.opt file in the current directory.
%
% WARNING: Specifying maxsz overwrites the ng.opt in the current directory.
%
% ELEC_POS (optional) defines electrodes:
%  - ep (Nx2)              : a list of points in 2D (will be projected on
%                            closest edge of the first contour specified in
%                            SHAPE
%  - ep (1x1) = N          : the number of equidistant electrodes to create
%                            with first electrode on the first point in XY
%                            and counter-clockwise ordering. Specify a
%                            negative number for clockwise ordering.
%  - ep (1x2) = [N offset] : specify offset of the first electrode with
%                            respect to the first point of XY
%                            (clockwise if negative, counter-clockwise
%                            otherwise)
%  - {ep1, ep2, ...}       : allows specifying electrodes on the internal 
%                            contours specified in SHAPE. Use an empty 
%                            array [] if a contour has no electrodes
%
% ELEC_SHAPE (optional) defines the electrode shape
%  - es (1x2) = [wd rfnum] : WD defines width of the electrode (default: 0 
%                            i.e. point electrode
%                            RFNUM controls amount of local refinement
%                            around the electrode.
%  - es (NEx2)             : specifies the above for each electrode
%                            individually
% 
% NOTE: smaller MAXSZ generally requires a lower RFNUM than a coarser mesh
% would.
%
% Examples:
%
% xy = [0 0;  1 0; 1 1; 0 1];
% mdl = ng_mk_2d_model({xy, 0.25 + 0.5*xy});
% mdl = ng_mk_2d_model({xy, 0.25 + 0.5*xy, 0.1}, [0.5 1; 0.5 0; 0 0.5]);
% mdl = ng_mk_2d_model({xy, 0.25 + 0.5*xy, 0.1}, [5, 0.25]);
% mdl = ng_mk_2d_model({xy, 0.25 + 0.5*xy, 0.1}, [-5, 0.25]);
% mdl = ng_mk_2d_model({xy, 0.25 + 0.5*xy, 0.1}, [5, -0.25]);
% mdl = ng_mk_2d_model({xy, 0.25 + 0.5*xy, 0.1}, {[5, -0.25], [4 0.1]});
% mdl = ng_mk_2d_model({xy, 0.1 + 0.25*xy, 0.4 + 0.5*xy, 0.1}, {[5, -0.25], [4 0.1]});
% mdl = ng_mk_2d_model({xy, 0.1 + 0.25*xy, 0.4 + 0.5*xy, 0.1}, {[5, -0.25], [4 0.1], [4]});
% mdl = ng_mk_2d_model({xy, 0.1 + 0.25*xy, 0.4 + 0.5*xy, 0.1}, {[5, -0.25], [], [4]});
% mdl = ng_mk_2d_model({xy, 0.1 + 0.25*xy, 0.4 + 0.5*xy, 0.1}, {[5, -0.25], [], [4]},[0 30]);
% mdl = ng_mk_2d_model({xy, 0.25 + 0.5*xy, 0.1}, [5, 0.25],[0.2,10;0 20; 0 30; 0 20; 0 10]);
% mdl = ng_mk_2d_model({xy, 0.25 + 0.5*xy, 0.1}, [5, 0],[0.2,10;0 20; 0 20; 0 20; 0 20]);
% mdl = ng_mk_2d_model({xy, 0.1 + 0.25*xy, 0.4 + 0.5*xy, 0.1}, {[5, -0.25], [], [4]},...
%     [0.2,10;0 20; 0 20; 0 20; 0 20; 0 20; 0 20; 0.2 20; 0 20]);


% (C) 2012-2013 Bartlomiej Grychtol, License: GPL version 2 or version 3
% $Id$


if ischar(varargin{1}) && strcmp(varargin{1}, 'UNIT_TEST'), mdl = do_unit_test; return, end 

[shape, elec_pos, elec_shape] = process_input(varargin{:});

mdl = eidors_cache(@ng_mk_2d_model_do,{shape, elec_pos, elec_shape});



function [shape, elec_pos, elec_shape] = process_input(shape, elec_pos, elec_shape)

if ~iscell(shape)
   shape = {shape};
end

if nargin < 2
    elec_pos = [];
end
if ~iscell(elec_pos)
    elec_pos = {elec_pos};
end

if nargin < 3
    elec_shape = [0 10]; % point electrode
end
if ~iscell(elec_shape)
    elec_shape = {elec_shape};
end
if numel(elec_shape) == 1 && numel(elec_pos) > 1
    elec_shape(2:numel(elec_pos)) = elec_shape(1);
end


function mdl = ng_mk_2d_model_do(shape, elec_pos, elec_shape)

[shape,i_wrote_ng_opt] = process_maxsz(shape);

points = [];
eidx = [];
eref = [];
for i = 1:length(shape)
   lp = length(points);
   ls = length(shape{i});
   if i <= numel(elec_pos) && ~isempty(elec_pos{i})
       [pp e_idx elec_pos{i} e_ref] = integrate_elecs(shape{i},elec_pos{i},elec_shape{i});
       ls = length(pp);
       points  = [points; pp];
   else
       e_idx = zeros(1,length(shape{i}));
       e_ref = [];
       points = [points; shape{i}];
   end
   if ~isempty(eidx)
       eidx = [eidx max(double(eidx))*(e_idx>0)+e_idx];
   else 
       eidx = e_idx;
   end
   eref = [eref; e_ref];
   seg{i} = repmat([0 1],ls,1) + lp + repmat((1:ls)',1,2);
   seg{i}(end,2) = lp + 1;
end

fnamebase = tempname;
fnamein2d = [fnamebase, '.in2d'];
fnamevol =  [fnamebase, '.vol'];

write_in2d_file(fnamein2d, points, seg, eidx, eref);

call_netgen( fnamein2d, fnamevol);
if i_wrote_ng_opt; delete('ng.opt'); end

[mdl,mat_idx] = ng_mk_fwd_model( fnamevol, [], 'ng', []);

delete(fnamein2d); 
delete(fnamevol); 

mdl.nodes(:,3) = [];
if ~isempty(elec_pos{1})
    mdl = find_electrodes(mdl, points(find(eidx),:), nonzeros(eidx));
end
mdl.boundary = find_boundary(mdl);
if isfield(mdl, 'electrode')
    for i = 1:length(mdl.electrode)
        mdl.electrode(i).z_contact = 0.01;
    end
end

function [shape,i_wrote_ng_opt] = process_maxsz(shape)
maxsz = [];
if numel(shape{end})==1
    maxsz = shape{end};
    shape(end)=[];
end
if ~isempty(maxsz)
    ng_write_opt('meshoptions.fineness',6,'options.meshsize',maxsz);
    i_wrote_ng_opt = true;
else
    i_wrote_ng_opt = false;
end

function mdl = find_electrodes(mdl, elec_pts, e_idx)

opt.boundary_face = 1;
mdl = fix_model(mdl, opt); % in case there are multi-point electrodes

nel = max(e_idx);
npts = length(elec_pts);
nn  = length(mdl.nodes);
e = elec_pts';
d = repmat(e(:)',nn,1) - repmat(mdl.nodes,1,npts);
d = sqrt(d(:,1:2:end).^2 + d(:,2:2:end).^2);
for j = 1:nel
    epts = find(e_idx==j);
    for k = 1:length(epts)
        [val mdl.electrode(j).nodes(k)] = min(d(:,epts(k)));
    end
    if numel(mdl.electrode(j).nodes) > 1
        mdl.electrode(j).nodes = fill_in_elec_nodes(mdl, mdl.electrode(j).nodes);
    end
end

function nds = fill_in_elec_nodes(mdl,enodes)
fcs = mdl.faces(mdl.boundary_face,:);
% fcs are ordered such that all(fcs(:,1) < fcs(:,2))
% we assume that enodes are also sorted
nds(1) = enodes(1);
for i = 1:length(enodes)-1
    startnode  = enodes(i);
    targetnode = enodes(i+1);
    nextnode   = startnode;
    while nextnode ~= targetnode
        % find the two faces the startnode is on
        % consider which of the two nodes at their other ends is closer to
        % targetnode
        idx = find(fcs(:,1) == nextnode);
        switch numel(idx)
            case 2
                c1 = fcs(idx(1),2);
                c2 = fcs(idx(2),2);
            case 1
                c1 = fcs(idx(1),2);
                idx(2) = find(fcs(:,2) == nextnode);
                c2 = fcs(idx(2),1);
            case 0
                idx = find(fcs(:,2) == nextnode);
                c1 = fcs(idx(1),1);
                c2 = fcs(idx(2),1);
            otherwise
                error('huh?');
        end
        d1 = sqrt(sum((mdl.nodes(c1,:) - mdl.nodes(targetnode,:)).^2,2));
        d2 = sqrt(sum((mdl.nodes(c2,:) - mdl.nodes(targetnode,:)).^2,2));
        if d1 < d2
            nextnode = c1;
        else
            nextnode = c2;
        end
        nds(end+1) = nextnode;
    end
end
    

    



function [newpoints eidx elec_pos e_ref] = integrate_elecs(points, elec_pos, elec_shape)


n_elecs = size(elec_pos,1);
if n_elecs == 1
    % the number of electrodes was specified, positions need to be found
    n_elecs = elec_pos(1);
    start = 0;
    try start = elec_pos(2); end
    elec_pos = equidistant_elec_pos(points, n_elecs, start);
    n_elecs = size(elec_pos,1);
end

if size(elec_shape,1) == 1;
    elec_shape = repmat(elec_shape,n_elecs,1);
end

newpoints = points;
eidx = zeros(1, length(points));
eref = zeros(1, length(points));

for i = 1:n_elecs
    A = newpoints;
    B = circshift(newpoints,-1);
    AB = B-A;    L = sqrt(sum((AB).^2,2));

    % 1. find the closest edge
    % 2. add between the endpoints
    E = elec_pos(i,:);
    AE = repmat(E,size(A,1),1) - A;
    r = sum(AE .* AB,2)./L.^2;
    P = A + r*[1 1].*AB; % E projected on each edge
    D = sqrt(sum((repmat(E, size(A,1),1)-P).^2,2));
    D(r>1) = Inf; D(r<0) = Inf;
    [jnk e] = min(D); % index of closest edge
    
    if elec_shape(i,1) == 0 % point electrode
        if r(e) == 0
            eidx(e) = i;
            eref(e) = elec_shape(i,2);
        elseif r(e) == 1
            if e==length(A);
                eidx(1) = i;
                eref(1) = elec_shape(i,2);
            else
                eidx(e+1) = i;
                eref(e+1) = elec_shape(i,2);
            end
        else
            newpoints = [newpoints(1:e,:); P(e,:); newpoints(e+1:end,:)];
            eidx      = [eidx(1:e) i eidx(e+1:end)];
            eref      = [eref(1:e) elec_shape(i,2) eref(e+1:end)];
        end
    else % multi-point electrode
        % e is the first node of the edge the centre lies on
        
        % 1. Need the perimeter coordinate of the centre
        p = sqrt(sum((circshift(newpoints,-1) - newpoints).^2,2));
        vec = [0; cumsum(p)];
        L = vec(end); % total length
        ctr = vec(e) + r(e)*(vec(e+1) - vec(e));
        e_fr = linspace(ctr-elec_shape(i,1)/2 , ctr+elec_shape(i,1)/2,2);
        e_fr = rem(e_fr, L); % wrap around
        e_fr(e_fr<0) = L + e_fr(e_fr<0); % wrap around
        for j = 1:numel(e_fr)
            k = find(vec <= e_fr(j), 1, 'last');
            if k == length(vec)
                % handle the case where electrode falls on the last point
                eidx(1) = i;
                eref(1) = elec_shape(i,2);
            else
                r = (e_fr(j) - vec(k)) / (vec(k+1) - vec(k));
                if r == 0
                    eidx(k) = i;
                    eref(k) = elec_shape(i,2);
                end
                jnkpts = newpoints; jnkpts(end+1,:) = jnkpts(1,:);
                p = newpoints(k,:) + r * (jnkpts(k+1,:) - newpoints(k,:));
                if k < length(eidx)
                    newpoints = [newpoints(1:k,:); p; newpoints(k+1:end,:)];
                    eidx = [eidx(1:k) i eidx(k+1:end)];
                    eref = [eref(1:k) elec_shape(i,2) eref(k+1:end)];
                else
                    eidx = [eidx i ];
                    eref = [eref elec_shape(i,2)];
                    newpoints = [newpoints; p];
                end
                vec = [vec(1:k); e_fr(j); vec(k+1:end)];
            end
        end
        
    end
        
end
e_ref = nonzeros(eref);


function elec_pos = equidistant_elec_pos(points, n_elecs, start)
% 1. Calculate the perimeter
p = sqrt(sum((circshift(points,-1) - points).^2,2));
vec = [0; cumsum(p)];
L = vec(end); % total length

if n_elecs > 0
    e_fr = linspace(start, L+start, n_elecs+1); e_fr(end) = [];
else
    e_fr = linspace(L+start, start, -n_elecs+1); e_fr(end) = [];
    n_elecs = -n_elecs;
end
e_fr = rem(e_fr, L); % wrap around
e_fr(e_fr<0) = L + e_fr(e_fr<0); % wrap around
elec_pos = NaN(n_elecs,2);
points(end+1,:) = points(1,:);
for i = 1:n_elecs
    j = find(vec <= e_fr(i), 1, 'last');
    if j == length(vec)
        % handle the case where electrode falls on the last point
        elec_pos(i,:) = points(1);
    else
        r = (e_fr(i) - vec(j)) / (vec(j+1) - vec(j));
        elec_pos(i,:) = points(j,:) + r * (points(j+1,:) - points(j,:));
    end
end



function write_in2d_file(fname,points, seg, e_idx, e_ref)

if length(e_idx) < length(points);
    e_idx(length(points)) = 0;
end

refine = ones(length(points),1);
if any(e_idx)
    refine(logical(e_idx)) = e_ref;  % refinement factor (somehow)
end
fid = fopen(fname,'w');
fprintf(fid, '%s\n','splinecurves2dv2');
fprintf(fid, '%d\n',6); % global grading factor, 6 should force use of ng.opt
fprintf(fid, '%s\n','points');
for i = 1:length(points)
   fprintf(fid, '%d   %f   %f   %f\n',i,points(i,:), refine(i));
end
fprintf(fid,'%s\n','segments');
% here we assume the first loop is the boundary, all the others are holes
domains = [ 1 0];
for i = 1:length(seg)
   if i > 1
      domains = [0 1];
   end
   for j = 1:length(seg{i})
      fprintf(fid,'%d   %d   %d   %d   %d -bc=%d\n',domains, 2, seg{i}(j,:),i);
   end
end
fclose(fid);

function mdl = do_unit_test
xy = [0 0;  1 0; 1 1; 0 1];
for i = 1:15
    switch i
        case 1
            mdl = ng_mk_2d_model({xy, 0.25 + 0.5*xy});
        case 2
            mdl = ng_mk_2d_model({xy, 0.25 + 0.5*xy, 0.1}, [0.5 1; 0.5 0; 0 0.5]);
        case 3
            mdl = ng_mk_2d_model({xy, 0.25 + 0.5*xy, 0.1}, [5, 0.3]);
        case 4
            mdl = ng_mk_2d_model({xy, 0.25 + 0.5*xy, 0.1}, [5, 0.25]);
        case 5
            mdl = ng_mk_2d_model({xy, 0.25 + 0.5*xy, 0.1}, [-5, 0.25]);
        case 6 
            mdl = ng_mk_2d_model({xy, 0.25 + 0.5*xy, 0.1}, [5, -0.25]);
        case 6
            mdl = ng_mk_2d_model({xy, 0.25 + 0.5*xy, 0.1}, {[5, -0.25], [4 0.1]});
        case 7
            mdl = ng_mk_2d_model({xy, 0.1 + 0.25*xy, 0.1}, {[5, -0.25], [4 0.1]});
        case 8
            mdl = ng_mk_2d_model({xy, 0.1 + 0.25*xy, 0.4 + 0.5*xy, 0.1}, {[5, -0.25], [4 0.1]});
        case 9
            mdl = ng_mk_2d_model({xy, 0.1 + 0.25*xy, 0.4 + 0.5*xy, 0.1}, {[5, -0.25], [4 0.1], [4]});
        case 10
            mdl = ng_mk_2d_model({xy, 0.1 + 0.25*xy, 0.4 + 0.5*xy, 0.1}, {[5, -0.25], [], [4]});
        case 11 
            mdl = ng_mk_2d_model({xy, 0.1 + 0.25*xy, 0.4 + 0.5*xy, 0.1}, {[5, -0.25], [], [4]},[0 30]);
        case 12
            mdl = ng_mk_2d_model({xy, 0.25 + 0.5*xy, 0.1}, [5, 0.25],[0.2,10;0 20; 0 30; 0 20; 0 10]);
        case 13
            mdl = ng_mk_2d_model({xy, 0.25 + 0.5*xy, 0.1}, [5, 0],[0.2,10;0 20; 0 20; 0 20; 0 20]);
        case 14
            mdl = ng_mk_2d_model({xy, 0.1 + 0.25*xy, 0.4 + 0.5*xy, 0.1}, {[5, -0.25], [], [4]},...
                [0.2,10;0 20; 0 20; 0 20; 0 20; 0 20; 0 20; 0.2 20; 0 20]);
        case 15
            mdl = ng_mk_2d_model({xy, 0.1 + 0.25*xy, 0.4 + 0.5*xy, 0.05}, {[5, -0.25], [], [4]},...
                [0.2,10;0 20; 0 20; 0 20; 0 20; 0 20; 0 20; 0.2 20; 0 20]);
    end
    show_fem(mdl,[0 1 0]);
    drawnow
end
