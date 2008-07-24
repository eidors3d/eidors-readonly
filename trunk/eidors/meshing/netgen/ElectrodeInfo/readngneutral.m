function [vtx,simp,surf] = readngneutral(filename);
% [vtx,simp,surf] = readngneutral(filename)
% Reads a Netgen neutral volume file and extracts 
% vertex coordinates and simplices
% together with surf which is a list of 'boundary condition index' followed 
% by three indecies for a triangle

% WRB Lionheart 23/07/2004
% Part of the EIDRORS 3D project
% Todo: 
%          Make more robust
%          Include surface data
%          Read specification of vol file rather than reverse 
%              engineering
fid = fopen(filename,'r');
if fid==-1
    error(['Cannot open netgen  neutral volume file ',filename])
end
% Vertices
tline = fgetl(fid); 
nvtx = sscanf(tline,'%d');
for ivtx = 1:nvtx
    tline = fgetl(fid);
    row = sscanf(tline,'%f')';
    vtx(ivtx,:)=row;
end
% elements
tline = fgetl(fid); 
nsimp = sscanf(tline,'%d');
for isimp = 1:nsimp
    tline = fgetl(fid);
    row = sscanf(tline,'%d')';
    simp(isimp,:)=row(2:end);
end
% Boundary surfacess
tline = fgetl(fid); 
nbd= sscanf(tline,'%d');
for ibd = 1:nbd
    tline = fgetl(fid);
    row = sscanf(tline,'%d')';
    surf(ibd,:)=row;
end
fclose(fid);
