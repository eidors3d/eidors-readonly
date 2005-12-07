function[elec,sels] = ng_tank_select_elec(srf,vtx,fc,mshaxs);
%function[elec,sels] = ng_tank_select_elec(srf,vtx,fc,mshaxs);
%
% This function takes the wire frame model and the face
% data and sequentially requests if each face is an
% electrode or earth. From this the indices matrices (sels & sgnd)
% are concatenated so as to map all the selected electrodes and 
% grounded planes into the srf matrix. In addition the cell array
% (elsrf) is created containing the electrode surfaces and ground
% plane indices matrices which map directly into vtx.
%
% Version 5.0
% B.D.Grieve - 13/02/2002 + modyfication by lmazurk
% WRBL added default as prevous choice 20/1/2004
% WRBL deleted ground plane 05/12/2005
%
% srf      = The boundary surfaces
% vtx      = The vertices matrix
% fc       = A one column matrix containing the face numbers
% mshaxs   = Axes details for plotting wire frame
% elsrf    = Cell array of indices matrices mapping into vtx each electrode face
% sels     = The indices into the srf matrix of the selected electrode faces
% elec  = The EIDORS-3D electrode matrix of dimensions NxM, where 
%            where N: no. of electrodes, M: 3 * max no. of faces per electrode
%


sels = [];


figure
set(gcf,'Name','Object Faces','Colormap',[0 0 0])
elfc_prev = 'N';
max(fc)

for loop1 = 1:max(fc)
    % Create a logical array (lgelfc) to determine which faces are electrodes
    lgelfc(loop1) = logical(0);
    
    [fcsrf,fci] = ng_extract_face(srf,vtx,fc,loop1);
    if ~isempty(fcsrf)
        % Add this face's vtx indices matrix to the cell array ttlfcsrf
        ttlfcsrf(loop1) = {fcsrf};
        % Plot this face
        trimesh(fcsrf,vtx(:,1),vtx(:,2),vtx(:,3))
        title(['Face number: ' num2str(loop1) ' of ' num2str(max(fc))])
        axis equal
        axis(mshaxs)
        view(45,10)
	qstr =sprintf('Is this face an electrode? Y/N [%s] ',elfc_prev);
        elfc = input(qstr,'s');
	if isempty(elfc)
	   elfc = elfc_prev;
	else
	   elfc_prev=elfc;
	end
        if  (elfc=='y' | elfc=='Y')
            lgelfc(loop1) = logical(1);
            sels = [sels; fci]; % Concatenate indices into sels for this face
        end
    end
end


% Extract from the total face indices matrix (ttlfcsrf) the
% faces which are electrodes and store them in the cell
% array (elsrf)
elsrf = ttlfcsrf(lgelfc);

close(gcf)
% Display each electrode in turn as a wire mesh
figure
set(gcf,'Name','Wire Mesh Electrode Faces')
for loop1 = 1:size(elsrf,2)-1
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

% Convert elsrf into the EIDORS-3D matrix electrode matrix format

nmel = size(elsrf,2);
for loop1 = 1:nmel
    nmfc(loop1) = size(elsrf{loop1},1);
end
% Initiate electrode matrix (elec) & pad with zeros
elec = zeros(nmel,3*max(nmfc));
% Put electrode surface information into elec
for loop1 = 1:nmel
    for loop2 = 1:size(elsrf{loop1},1)
        elec(loop1,loop2*3-2:loop2*3)=elsrf{loop1}(loop2,:);
    end
end
