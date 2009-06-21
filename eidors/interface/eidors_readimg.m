function img = eidors_readimg( fname, format )
% EIDORS readimg - read reconstructed image files from various EIT equipment
%    manufacturers
%
% Currently the list of supported file formats is:
%    - NATIVE "e3d" file format
%        format = "e3d" or "NATIVE"
%    - MCEIT (Goettingen / Viasys) "igt" file format 
%        format = "IGT" or "MCEIT"
%
% Usage
% img = eidors_readimg( fname, format )
%     img   = eidors image structure
%     img.elem_data = reconstructed image matrix NumPixels x NumFrames
%     fname = file name
%
%  If format is unspecified, we attempt to autodetect.

% (C) 2009 by Bartlomiej Grychtol. Licensed under GPL v2 or v3
% $Id$ 

if ~exist(fname,'file')
   error([fname,' does not exist']);
end

if nargin < 2
% unspecified file format, autodetect
   dotpos = find(fname == '.');
   if isempty( dotpos ) 
      error('file format unspecified, can`t autodetect');
   else
      dotpos= dotpos(end);
      format= fname( dotpos+1:end );
   end
end
fmt= lower(format);

switch fmt
    case {'igt','mceit'}
        img = mceit_readimg( fname );
    case {'e3d','native'}
        img = native_readimg( fname );
    otherwise
        error('eidors_readdata: file "%s" format unknown', fmt);
end




%%
function img = mceit_readimg( fname )
% mceit_readimg - reads in IGT files. 
img.name = ['Read from ' fname];
img.type = 'image';
tempmdl = mk_common_gridmdl('backproj');
img.fwd_model = tempmdl.fwd_model;

fid = fopen(fname,'r');
igt = fread(fid, inf,'4*float');
fclose(fid);

igt = reshape(igt, [], 912);

img.elem_data = igt';


%%
function img = native_readimg( fname )
% native_readimg - reads in native E3D files.
% E3D file is a zipped matlab v6 compatible .mat file called "e3d.temp"
% containing one eidors image struct variable named "img".

files = unzip(fname);
if numel(files) > 1
    error(['File %s is not a proper E3D file. '...
        'The archive contains more than one file'],fname);
end

tempfile = files{1};
S = load(tempfile,'-mat');
delete(tempfile);
if numel(fieldnames(S)) > 1
     warning(['File %s is not a proper E3D file. '...
        'The mat file contains more than one variable. Strugling on.'],fname);
end

if isfield(S,'img')
    img = S.img;
else
    error(['File %s is not a proper E3D file. '...
        'The mat file does not contain an "img" variable.'],fname);
end



