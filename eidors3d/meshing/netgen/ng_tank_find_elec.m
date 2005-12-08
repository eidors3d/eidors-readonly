function[elec,sels] = ng_tank_find_elec(srf,vtx,fc,centres);
%function[elec,sels] = ng_tank_find_elec(srf,vtx,fc,centres);
%
% This function Tries to find the electrdes given the x y x coords of their centres.
%
% Version 5.0
% B.D.Grieve - 13/02/2002 + modyfication by lmazurk
% WRBL added default as prevous choice 20/1/2004
% WRBL deleted ground plane 05/12/2005
% WRBL derived automatic version ditto
%
% srf      = The boundary surfaces
% vtx      = The vertices matrix
% fc       = A one column matrix containing the face numbers
% elsrf    = Cell array of indices matrices mapping into vtx each electrode face
% sels     = The indices into the srf matrix of the selected electrode faces
% elec  = The EIDORS-3D electrode matrix of dimensions NxM, where 
%            where N: no. of electrodes, M: 3 * max no. of faces per electrode
% centres(k,:)=[ x,y,z ] of kth electrode


sels = [];



for loop1 = 1:max(fc)
    % Create a logical array (lgelfc) to determine which faces are electrodes
    lgelfc(loop1) = logical(0);
    
    [fcsrf,fci] = ng_extract_face(srf,vtx,fc,loop1);
%   size(fcsrf) % should be vertex numbers for this face
    coordsforthisface= vtx(fcsrf,:);
    centreofface(loop1,:)= mean(coordsforthisface);
    ttlfcsrf(loop1) = {fcsrf};
    
end

for ielec = 1:size(centres,1)
% Find the distance from the centre of faces to this electrode
    dists =  (centreofface(:,1) - centres(ielec,1)).^2 + ...
             (centreofface(:,2) - centres(ielec,2)).^2 + ...
             (centreofface(:,3) - centres(ielec,3)).^2;
    [d,iface] = min(dists); %iface is closest face to this electrode.
    lgelfc(iface) = logical(1);
    sels(ielec)= iface;
end


% Extract from the total face indices matrix (ttlfcsrf) the
% faces which are electrodes and store them in the cell
% array (elsrf)
elsrf = ttlfcsrf(lgelfc);

close(gcf)
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
    this_el= ttlfcsrf{el_idx};
    for loop2 = 1:size(this_el,1)
        elec(loop1,loop2*3 + (-2:0))=this_el(loop2,:);
    end
end
