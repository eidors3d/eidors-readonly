function mdl = ng_mk_2d_model(shape, elec_pos)
%NG_MG_2D_MODELS create a 2D mesh with Netgen via the in2d interface
% mdl = ng_mk_2d_model(xy) creates a 2D model from a clockwise list of
% points
%
% mdl = ng_mk_2d_model({xy, xy1, xy2, ...}) allows specifying additional
% clockwise loops xy1, xy2, etc, which represent holes in the bigger xy
% contour
%
% mdl = ng_mk_2d_model({..., maxsz}) specifies maximum element size
%
% Example:
%
%  xy = [0 0;  1 0; 1 1; 0 1];
%  ng_mk_2d_model({xy, 0.25 + 0.5*xy});
%

% (C) 2012 Bartlomiej Grychtol, License: GPL version 2 or version 3
% $Id$

%TODO: Add support for electrodes

if nargin < 2
    elec_pos = [];
end

mdl = [];
if ischar(shape) && strcmp(shape, 'UNIT_TEST'), do_unit_test; return, end 
if ~iscell(shape)
   shape = {shape};
end

shape = process_maxsz(shape);

points = [];
e_idx = [];

for i = 1:length(shape)
   lp = length(points);
   ls = length(shape{i});
   if i == 1 && ~isempty(elec_pos)
       [points e_idx elec_pos e_ref] = integrate_elecs(shape{i},elec_pos);
       ls = length(points);
   else
       points = [points; shape{i}];
   end
   seg{i} = repmat([0 1],ls,1) + lp + repmat((1:ls)',1,2);
   seg{i}(end,2) = lp + 1;
end

% plot(points(:,1),points(:,2),'bo');
% hold on
% plot(elec_pos(:,1),elec_pos(:,2),'rx');
% hold off
% return

points = [points; elec_pos];
write_in2d_file('tmp2.in2d',points, seg, e_idx, e_ref);

call_netgen( 'tmp2.in2d', 'tmp2.vol');
[mdl,mat_idx] = ng_mk_fwd_model( 'tmp2.vol', [], 'ng', []);
mdl.nodes(:,3) = [];
if ~isempty(elec_pos)
    mdl = find_electrodes(mdl, elec_pos);
end

function shape = process_maxsz(shape)
maxsz = [];
if numel(shape{end})==1
    maxsz = shape{end};
    shape(end)=[];
end
if ~isempty(maxsz)
    ng_write_opt('meshoptions.fineness',6,'options.meshsize',maxsz);
end

function mdl = find_electrodes(mdl, elec_pos)
% threshold for accepting node as correct position
THRESH = 1e-6 * max(max(mdl.nodes) - min(mdl.nodes));

nel = size(elec_pos,1);
nn  = length(mdl.nodes);
e = elec_pos';
d = repmat(e(:)',nn,1) - repmat(mdl.nodes,1,nel);
d = sqrt(d(:,1:2:end).^2 + d(:,2:2:end).^2);
for i = 1:nel
    [val mdl.electrode(i).nodes] = min(d(:,i));
    if 0 %val > THRESH
        mdl.electrode(i).nodes = [];
        warning('No nodes for electrode %d found',i);
    end
end
   


function [newpoints eidx elec_pos e_ref] = integrate_elecs(points, elec_pos)
e_ref = 10;
if iscell(elec_pos)
    e_ref = elec_pos{2};
    elec_pos = elec_pos{1};
end
n_elecs = size(elec_pos,1);
if n_elecs == 1
    % the number of electrodes was specified, positions need to be found
    n_elecs = elec_pos(1);
    start = 0;
    try start = elec_pos(2); end
    elec_pos = equidistant_elec_pos(points, n_elecs, start);
    n_elecs = size(elec_pos,1);
end

newpoints = points;
eidx = false(1, length(points));
for i = 1:n_elecs
    A = newpoints;
    B = circshift(newpoints,-1);
    AB = B-A;
    L = sqrt(sum((AB).^2,2));

    % 1. find the closest edge
    % 2. add between the endpoints
    E = elec_pos(i,:);
    AE = repmat(E,size(A,1),1) - A;
    r = sum(AE .* AB,2)./L.^2;
    P = A + r*[1 1].*AB; % E projected on each edge
    D = sqrt(sum((repmat(E, size(A,1),1)-P).^2,2));
    [jnk e] = min(D); % index of closest edge
    if r(e) == 0
        eidx(e) = true;
    elseif r(e) == 1
        if e==length(A);
            eidx(1) = true;
        else
            eidx(e+1) = true;
        end
    else
        newpoints = [newpoints(1:e,:); P(e,:); newpoints(e+1:end,:)];
        eidx      = [eidx(1:e) true eidx(e+1:end)];
    end
end

function elec_pos = equidistant_elec_pos(points, n_elecs, start)
% 1. Calculate the perimeter
p = sqrt(sum((circshift(points,-1) - points).^2,2));
vec = [0; cumsum(p)];
L = vec(end); % total lenght
vec = vec;

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
    if j == length(points)
        % handle the case where electrode falls on the last point
        elec_pos(i,:) = points(j);
    else
        r = (e_fr(i) - vec(j)) / (vec(j+1) - vec(j));
        elec_pos(i,:) = points(j,:) + r * (points(j+1,:) - points(j,:));
    end
end



function write_in2d_file(fname,points, seg, e_idx, e_ref)
if length(e_idx) < length(points);
    e_idx(length(points)) = false;
end

refine = ones(length(points),1);
if any(e_idx)
    refine(e_idx) = e_ref;  % refinement factor (somehow)
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

function do_unit_test
xy = [0 0;  1 0; 1 1; 0 1];
% ng_mk_2d_model({xy, 0.25 + 0.5*xy});
% mdl = ng_mk_2d_model({xy, 0.25 + 0.5*xy, 0.1}, [0.5 1; 0.5 0; 0 0.5]);
% this one fails for no good reason:
% mdl = ng_mk_2d_model({xy, 0.25 + 0.5*xy, 0.1}, [5, 0.3]);
% mdl = ng_mk_2d_model({xy, 0.25 + 0.5*xy, 0.1}, [5, 0.25]);
% mdl = ng_mk_2d_model({xy, 0.25 + 0.5*xy, 0.1}, [-5, 0.25]);
mdl = ng_mk_2d_model({xy, 0.25 + 0.5*xy, 0.1}, [5, -0.25]);

show_fem(mdl,[0 1 0]);