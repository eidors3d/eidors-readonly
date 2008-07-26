function[srf,vtx,fc,bc,simp,edg,mat_ind] = ng_read_mesh(filename)
%[srf,vtx,fc,bc,simp,edg,mat_ind] = ng_read_mesh(filename)
% Function to read in a mesh model from NetGen and saves it in
% five arrays; surface (srf), veritices (vtx), face no. (fc)
% volume (simp) and edges (edg)
% 
% Version 4.0
% B.D.Grieve - 27/01/2002 + modifications by lmazurk
% A.Adler - 2006 mods to run quicker
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

% $Id$
% (C) 2002-2006 (C) Licenced under the GPL

eidors_msg('ng_read_mesh',3);


fid = fopen(filename,'r');
while 1
    tline = fgetl(fid);
    if ~ischar(tline); fclose(fid); break; end

    if     strcmp(tline,'surfaceelementsgi')
       se= get_lines_with_numbers( fid, 11);
    elseif strcmp(tline,'volumeelements')
       ve= get_lines_with_numbers( fid, 6);
    elseif strcmp(tline,'edgesegmentsgi2')
       es= get_lines_with_numbers( fid, 12);
    elseif strcmp(tline,'points')
       vtx= get_lines_with_numbers( fid, 3);
    end
end


srf = se(:,6:8);
fc = se(:,1);
simp = ve(:,3:6);
edg = es;
mat_ind=ve(:,1);
bc = se(:,2);

function mat= get_lines_with_numbers( fid, n_cols);
   tline = fgetl(fid);
   n_rows = sscanf(tline,'%d');
   mat= fscanf(fid,'%f',[n_cols,n_rows])';


