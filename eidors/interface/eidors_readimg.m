function img = eidors_readimg( fname, format )
% EIDORS readimg - read reconstructed image files from various EIT equipment
%    manufacturers
%
% Currently the list of supported file formats is:
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
    otherwise
        error('eidors_readdata: file "%s" format unknown', fmt);
end




%%
function img = mceit_readimg( fname );
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
