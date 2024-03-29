function [elec,sels, electrodes] = ng_tank_find_elec(srf,vtx,fc,centres)
%[elec,sels, electrodes] = ng_tank_find_elec(srf,vtx,fc,centres);
%
% This function Tries to find the electrdes given the x y x coords of their centres.
%
% Version 5.0
% B.D.Grieve - 13/02/2002 + modyfication by lmazurk
% WRBL added default as prevous choice 20/1/2004
% WRBL deleted ground plane 05/12/2005
% WRBL derived automatic version ditto
% AA   speedup and fix to not output zeros
% AA  add ability to define electrodes with a function
%
% srf      = The boundary surfaces
% vtx      = The vertices matrix
% fc       = A one column matrix containing the face numbers
% elsrf    = Cell array of indices matrices mapping into vtx each electrode face
% sels     = The indices into the srf matrix of the selected electrode faces
% elec     = The EIDORS-3D electrode matrix of dimensions NxM, where 
%            where N: no. of electrodes, M: 3 * max no. of faces per electrode
%            [ kept for backward compatibility. Use electrodes output instead]
% centres(k,:)=[ x,y,z ] of kth electrode
%  OR 
% centres = struct where centres(k).centre is electrode centre and
%             centres(k).fcn = fcn of vtx and ctr which is
%             zero outside and one inside electrode
%    example: ctr_fcn = inline( 'sum((vtx-ones(size(vtx,1),1)*ctr).^2,2)<1', ...
%                               'vtx','ctr');
%            [ctr_param(1:nn).fcn]     = deal( ctr_fcn );
%             centres=  mat2cell( [x(:),y(:),z(:)], ones(nn,1),3);
%            [ctr_param(1:nn).centre] = deal( centres{:} );
%    this form allows for more complicated electrode shapes
%          
% electrodes = EIDORS V3.x electrodes structure

% WARNING! the 'fc' variable in this function is what everywhere else is
% called 'bc', and that's how it should be used.

% (C) 2002-2006. Licenced under the GPL
% $Id$

if  isstruct(centres)
   % Calc centre of each surface
   for d=1:size(vtx,2)
      srfctr(:,d) = mean(reshape( vtx(srf,d), size(srf) ),2);
   end
   for e=1:length(centres)
      inside = feval(centres(e).fcn,srfctr,centres(e).centre);
      this_el = srf(inside,:);
      electrodes(e).nodes     = unique( this_el(:) )';
      electrodes(e).z_contact = 0.1; % set placeholder value
   end
   elec= NaN; sels= NaN; % No backward compatible for this mode
elseif size(centres,2)==3;
   [elec,sels, electrodes] = find_elec_centres(srf,vtx,fc,centres);
else
   error('don`t understand format of centres');
end


function [elec,sels, electrodes] = find_elec_centres(srf,vtx,fc,centres);

for loop1 = 1:max(fc)
    % Create a logical array (lgelfc) to determine which faces are electrodes
    lgelfc(loop1) = logical(0);
    
%   [fcsrf,fci] = ng_extract_face(srf,vtx,fc,loop1);
    fci  = find( fc == loop1 );
    fcsrf= srf(fci,:); % should be vertex numbers for this face
    fcsrf= unique(fcsrf); fcsrf= fcsrf(fcsrf>0);
    coordsforthisface= vtx(fcsrf,:);
    face_coords{loop1}= coordsforthisface;
    ttlfcsrf(loop1) = {fcsrf};
    
end


[sels,lgelfc] = find_selected_face(centres, face_coords, lgelfc);



% Extract from the total face indices matrix (ttlfcsrf) the
% faces which are electrodes and store them in the cell
% array (elsrf)
elsrf = ttlfcsrf(lgelfc);

if 0
   % Display each electrode in turn as a wire mesh
   figure
   set(gcf,'Name','Wire Mesh Electrode Faces')
   for loop1 = 1:size(elsrf,2)
       trimesh(elsrf{loop1},vtx(:,1),vtx(:,2),vtx(:,3),'EdgeColor','red')
       title(['Electrode ' num2str(loop1) ': red'])
       axis equal
       axis(mshaxs)
       view(45,10)
       hidden off
       hold on
       pause(0.75)
   end
   title('Electrodes: red,')
   hidden off
   pause(2)
end
% Convert elsrf into the EIDORS-3D matrix electrode matrix format

nmel = size(elsrf,2); 
for loop1 = 1:nmel
    nmfc(loop1) = size(elsrf{loop1},1);
end
% Initiate electrode matrix (elec) & pad with zeros
elec = zeros(nmel,3*max(nmfc));
% Put electrode surface information into elec
for loop1 = 1:nmel
    el_idx= sels(loop1);
    this_el= ttlfcsrf{el_idx}';
    l_this_el= prod(size(this_el));
    elec(loop1, 1:l_this_el) = this_el(:)';

    electrodes(loop1).nodes     = unique( this_el(:) )';
    electrodes(loop1).z_contact = 0.1; % set placeholder value
end


% Find the electrode node which is closest to the specified point
function [sels,lgelfc] = find_selected_face(centres, face_coords, lgelfc) 
   sels = [];
   elecn_idx= [];   
   elecnodes= [];   
   for i=1:length(face_coords)
       elecn_idx = [elecn_idx; i*ones(length(face_coords{i}),1)];
       elecnodes = [elecnodes; face_coords{i}];
   end
   elecsep = sum(min(abs(diff(centres))));
   for ielec = 1:size(centres,1)
   % Find the distance from the centre of faces to this electrode
       dists =  (elecnodes(:,1) - centres(ielec,1)).^2 + ...
                (elecnodes(:,2) - centres(ielec,2)).^2 + ...
                (elecnodes(:,3) - centres(ielec,3)).^2;
       dmin = min(dists); %iface is closest face to this electrode.
       % take the first, closest ones
       iface = find(dmin + elecsep >= dists);
       if length(iface)>1 % found electrode and background. Take smallest
          tryfaces = unique( elecn_idx(iface));
          ff=[]; for i=1:length(tryfaces)
            ff(i) = max( dists( elecn_idx == tryfaces(i)));
          end
          [~, i] = min(ff);
          iface = tryfaces(i);
       else
          iface = elecn_idx(iface);
       end
       lgelfc(iface) = true;
       if sum(lgelfc) ~= ielec
%           disp(ielec);
          error('Electrode #%d not found', ielec);
       end
       sels(ielec)= iface;
   %   disp([ielec, iface, d]);
   %   now remove that face so we dont use it again
       ff = find(elecn_idx == iface);
       elecn_idx(ff,:) = [];
       elecnodes(ff,:) = [];
   end

function [sels,lgelfc] = find_selected_face_old(centres, face_coords, lgelfc) 
   sels = [];
   for i=1:length(face_coords)
       centreofface(i,:)= mean(face_coords{i});
   end
   for ielec = 1:size(centres,1)
   % Find the distance from the centre of faces to this electrode
       dists =  (centreofface(:,1) - centres(ielec,1)).^2 + ...
                (centreofface(:,2) - centres(ielec,2)).^2 + ...
                (centreofface(:,3) - centres(ielec,3)).^2;
       [d,iface] = min(dists); %iface is closest face to this electrode.
       lgelfc(iface) = logical(1);
       if sum(lgelfc) ~= ielec
          disp(ielec);
          error('Electrode #%d not found', ielec);
       end
       sels(ielec)= iface;
   %   disp([ielec, iface, d]);
   end
