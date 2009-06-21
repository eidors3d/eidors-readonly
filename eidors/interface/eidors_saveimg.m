function eidors_saveimg( img, fname, format )
% EIDORS saveimg - save reconstructed image files in formats
%    of various EIT equipment manufacturers
% eidors_saveimg( img, fname, format )
%
% Currently the list of supported file formats is:
%    - NATIVE "e3d" file format
%        format = "E3D"
%    - MCEIT (Goettingen / Viasys) "igt" file format 
%        format = "IGT"
%
% Usage
% eidors_saveimg( img,fname,format )
%     img   = eidors image structure
%     fname = file name
%
%  If format is unspecified, we attempt to autodetect.

% (C) 2009 by Bartlomiej Grychtol. Licensed under GPL v2 or v3
% $Id: eidors_readimg.m 1842 2009-06-20 21:28:12Z nightflight84 $ 


switch nargin
    case 2
        fmt = detect_format(fname);
        if isempty( fmt )
            error('file format unspecified, can`t autodetect');
        end
    case 3
        fmt1 = detect_format(fname);
        fmt = lower(format);
        if isempty(fmt1);
            fname = [fname '.' fmt];
        else
            if ~strcmp(fmt1, fmt)
                error('The extension specified in file name doesn''t match the file format');
            end
        end
    otherwise
       error('Usage: eidors_saveimg( img , fname, format )');       
end


switch fmt
    case 'igt'
        mceit_saveimg( img, fname );
    case 'e3d'
        native_saveimg( img, fname);
    otherwise
        error('eidors_readdata: file "%s" format unknown', fmt);
end




%%
function fmt = detect_format( fname ) 

dotpos = find(fname == '.');
if isempty( dotpos )
    fmt = [];
else
    dotpos= dotpos(end);
    format= fname( dotpos+1:end );
    fmt= lower(format);
end



%%
function fid = open_file( fname );

if exist(fname,'file')
    disp('File already exists.');
    reply = input('Overwrite? Y/N [Y]: ', 's');
    if isempty(reply), reply = 'Y'; end
    reply = lower(reply);
    
    if ~strcmp(reply,'y');
        fid = -1;
        return;
    end
end
fid = fopen( fname ,'w');






%%
function mceit_saveimg( img, fname );
% mceit_readimg - saves IGT files. 

fid = open_file( fname );
if fid < 0
    error('Cannot open file.');
end

n = size(img.elem_data,1);
if n == 912
   %already the right format
   fwrite(fid,img.elem_data','4*float');
else
   data = img2igt(img);
   fwrite(fid, data , '4*float');
end

fclose(fid);

%%
function native_saveimg( img, fname )
% native_saveimg - saves E3D file.
% E3D file is a zipped matlab v6 compatible .mat file called "e3d.temp"
% containing one eidors image struct variable named "img".

% save temporary mat file
if ~exist('OCTAVE_VERSION') && str2double(version('-release')) < 14
    save('e3d.temp', 'img');
else
    save('e3d.temp', 'img', '-v6');
end
zip('temp.zip','e3d.temp');
movefile('temp.zip',fname);
delete e3d.temp 
