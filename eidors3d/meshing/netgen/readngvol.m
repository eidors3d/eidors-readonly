function [vtx,simp,surf] = readngvol(filename);
% [vtx,simp,surf] = readngvol(filename)
% Reads a Netgen native format volume mesh file and extracts 
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
    error(['Cannot open netgen  mesh file ',filename])
end

tline = fgetl(fid);
while  isempty(strfind(tline,'surfaceelement'))
    tline = fgetl(fid);
    if ~ischar(tline),error('Unexpected end of file'),end
end
% Boundary surfacess
tline = fgetl(fid);
nbd= sscanf(tline,'%d')
for ibd = 1:nbd
    tline = fgetl(fid);
    row = sscanf(tline,'%d')';
    surf(ibd,:)=[row(1),row(6:8)];
end

%volumeelements
tline = fgetl(fid);
while  isempty(strfind(tline,'volumeelement'))
    tline = fgetl(fid);
    if ~ischar(tline),error('Unexpected end of file'),end
end

% elements
tline = fgetl(fid); 
nsimp = sscanf(tline,'%d')
for isimp = 1:nsimp
    tline = fgetl(fid);
    row = sscanf(tline,'%d')';
    simp(isimp,:)=row(3:6);
end

%volumeelements
tline = fgetl(fid);
while  isempty(strfind(tline,'points'))
    tline = fgetl(fid);
end

% Vertices
tline = fgetl(fid); 
nvtx = sscanf(tline,'%d')
for ivtx = 1:nvtx
    tline = fgetl(fid);
    row = sscanf(tline,'%f')';
    vtx(ivtx,:)=row;
end
fclose(fid);

    
