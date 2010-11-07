% fwd_mdl = deform_cylinder( fwd_mdl, niv);
% Deform the boundary of the cylinder to make it like a torso
% niv= 1.. 5 => Torso shape from T5 - T12
% xyz_expand - rescale xyz - default should be [1];
function fwd_mdl = deform_cylinder( fwd_mdl, geo);
    NODE= fwd_mdl.nodes';
    a_max= size(geo.xy,1);
    ab_geo=sqrt(sum(([ geo.xy; geo.xy(1,:) ]').^2)');
    nn= zeros(size(NODE));
    for i=1:size(NODE,2);
      angle = rem(a_max*atan2( NODE(2,i), ...
            NODE(1,i) )/2/pi+a_max,a_max)+1;
      fl_angl = floor(angle + eps );
      fac=(  (fl_angl + 1 - angle) * ab_geo(fl_angl) + ...
             (angle - fl_angl)     * ab_geo(fl_angl + 1)  );
      nn(1:2,i)= NODE(1:2,i)* fac;
    end  %for i=1:size
    if size(nn,1) == 3;
       nn(3,:) = NODE(3,:)*geo.z_mag;
    end

    xyz_expand = 1;
    fwd_mdl.nodes = nn'*eye(xyz_expand);

function fmdl_out = deform_cylinder2( fmdl, deform_pararms)
% fmdl_out = deform_cylinder( fmdl_in, deform_pararms)
% DEFORM_CYLINDER: deform cylinder from circle to
%   nearly circular shape. Need to be a bit careful,
%   because too much deformation will make for poorly
%   shaped FEM elements. 
% 
% fmdl_out = output (deformed) fmdl
% fmdl_in  = input (non-deformed) fmdl
%
% deform_params values
%    .axes (rotate axes before and after deform)
%       deformation is in x,y axis
%       either as vector [1,3,2] or 2x2 (2D) or 3x3 mat
%       Default = [1,2,3]
%    .mode    = 'add' = add deformation 
%               'mul' = multiply by deformation 
%    .xycentre  = [xctr, yctr] (DEFAULT = 0,0)
%
%    .xyspecpos = [xelec, yelec]
%         specified new electrode positions
%
%    .zerorad = radius at which deform goes to zero
%        NOT IMPLEMENTED YET

if isstr(fmdl) && strcmp(fmdl,'UNIT_TEST'); do_unit_test;return;end

pp = proc_params( deform_pararms, fmdl )

function pp = proc_params( deform_pararms, fmdl )
   d = size(fmdl.nodes,2);
   if isfield(deform_pararms,'axes')
     dpax = deform_pararms.axes;
     [r,c] = size(dpax);
     if r==1;
        pp.axes = sparse(dpax,1:c,1,c,c); 
     else 
        pp.axes = sparse(dpax);

     end
   else 
     pp.axes = speye(d);
   end

   if isfield(deform_pararms,'mode')
     switch deform_pararms.mode
        case 'add'; pp.mode = 1;
        case 'mul'; pp.mode = 2;
        otherwise; error(['deform mode not understood']);
     end
   else
     pp.mode = 2;
   end

   if isfield(deform_params,'xycentre')
     pp.xctr = deform_pararms.xycentre(1);
     pp.yctr = deform_pararms.xycentre(2);
   else
     pp.xctr = 0;
     pp.yctr = 0;
   end

   % deform_pararms.xyspecpos must be given
   xy= deform_pararms.xyspecpos;
  [pp.th_s,pp.r_s] = cart2pol(xy(:,1)-pp.xctr, ...
                              xy(:,2)-pp.yctr);
   
   % Find rad, theta of electrodes
   for i=1:length(fmdl.electrode);
      enodes = fmdl.electrode(i).nodes;
      fepos(i,:) = mean(fmdl.nodes(enodes,:),1);
   end

  [pp.th_e,pp.r_e] = cart2pol(fepos(:,1)-pp.xctr, ...
                              fepos(:,2)-pp.yctr);

function deform(fmdl, pp);
th = mean([the,thm],2);
dr = re-rm;
dr = [dr;dr;dr];
th = [th-2*pi;th;th+2*pi];
[th,idx]=sort(th); dr= dr(idx);
tpf = polyfit(th,idx,10);
plot(th,dr,'*');

[thn,rn] = cart2pol( fmdl.nodes(:,2)-yc, fmdl.nodes(:,3)-zc ); 
drn = interp1(th,dr,thn);
rn = rn + drn;
[nodesy, nodesz] = pol2cart(thn,rn);
fmdl.nodes(:,2:3) = [nodesy+yc, nodesz+zc];

hh=plot(elec_pos(:,2),elec_pos(:,3),'*-',fepos(:,2),fepos(:,3),'o-');
set(hh,'LineWidth',2);
hold on;
f1ml = fmdl; f1ml.nodes = f1ml.nodes(:,[2,3,1]);
show_fem(f1ml);
hold off;

function do_unit_test
