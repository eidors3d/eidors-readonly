function[srf,vtx,fc,bc,simp,edg,mat_ind] = ng_read_mesh(filename)
%[srf,vtx,fc,bc,simp,edg,mat_ind] = ng_read_mesh(filename)
% Function to read in a mesh model from NetGen and saves it in
% five arrays; surface (srf), veritices (vtx), face no. (fc)
% volume (simp) and edges (edg)
% 
% Version 4.0
% B.D.Grieve - 27/01/2002 + modifications by lmazurk
%
% EIDORS's srf array is a subset of NetGen's surface element data
% (columns 6:8). The first column of the surface element data also
% ascribes a face number to each surface which is saved as the fc
% array. Each line of the srf array contains 3 indices to define
% a triangle mapped on to the three dimensional vtx array.
% EIDORS's vtx array is a direct equivalent to NetGen's pointer data.
%
%
% srf      = The surfaces indices into vtx
% simp     = The volume indices into vtx
% vtx      = The vertices matrix
% fc       = A one column matrix containing the face numbers
% edg      = Edge segment information
% filename = Name of file containing NetGen .vol information
% mat_ind  = Material index


%filename = input('Please enter the mesh filename [e.g. demo.vol]: ','s');


mdca = textread(filename,'%s');

% Retrieve pointers and data length information from .vol file
% Nested loops extract relevant information into arrays
for loop3 = 1:size(mdca,1)
    if strcmp(mdca(loop3),'surfaceelementsgi')
        lngse = str2num(cat(2,mdca{loop3+1}));
        startse = loop3+2;
        % Put Surface Element data into the array srf
        for loop1 = 0:lngse-1
            for loop2 = 0:10
                se(loop1+1,loop2+1) = str2num(cat(2,mdca{startse+loop1*11+loop2}));
            end
        end
        disp([num2str(lngse) ' surface elements retrieved'])
    end
    if strcmp(mdca(loop3),'volumeelements')
        lngve = str2num(cat(2,mdca{loop3+1}));
        startve = loop3+2;
        % Put Volume Element data into the array ve
        for loop1 = 0:lngve-1
            for loop2 = 0:5
                ve(loop1+1,loop2+1) = str2num(cat(2,mdca{startve+loop1*6+loop2}));
            end
        end
        disp([num2str(lngve) ' volume elements retrieved'])
    end
    if strcmp(mdca(loop3),'edgesegmentsgi2')
        lnges = str2num(cat(2,mdca{loop3+1}));
        startes = loop3+2;
        % Put Edge Segment data into the array edg
        for loop1 = 0:lnges-1
            for loop2 = 0:11
                es(loop1+1,loop2+1) = str2num(cat(2,mdca{startes+loop1*12+loop2}));
            end
        end
        disp([num2str(lnges) ' edge segments retrieved'])
    end
    if strcmp(mdca(loop3),'points')
        lngp = str2num(cat(2,mdca{loop3+1}));
        startp = loop3+2;
        % Put Points data into the array vtx
        for loop1 = 0:lngp-1
            for loop2 = 0:2
                vtx(loop1+1,loop2+1) = str2num(cat(2,mdca{startp+loop1*3+loop2}));
            end
        end
        disp([num2str(lngp) ' vertices retrieved'])
    end
end

srf = se(:,6:8);
fc = se(:,1);
simp = ve(:,3:6);
edg = es;
mat_ind=ve(:,1);
bc = se(:,2);
