function eidors_saveimg( img, fname, format, pp )
% EIDORS saveimg - save reconstructed image files in formats
%    of various EIT equipment manufacturers
% eidors_saveimg( img, fname, format, params )
%
% Currently the list of supported file formats is:
%    - MCEIT (Goettingen / Viasys) "igt" file format 
%        format = "IGT"
%    - SenTec / iBex - the "mat" format for SenTec's ibeX software
%        format = "ZRI.MAT"
%
% Usage
% eidors_saveimg( img,fname,format )
%     img   = eidors image structure
%     fname = file name
%
% For the "zri.mat" format, use (minumum)
%  p.imageRate = FR;
%  p.patient.ROI.Inside = thorax_ROI*100;
%  eidors_saveimg( img, 'filename.zri.mat','zri.mat', p)
%
%  If format is unspecified, we attempt to autodetect.

% (C) 2009 by Bartlomiej Grychtol. Licensed under GPL v2 or v3
% $Id$ 

if ischar(img) && strcmp(img,'UNIT_TEST'); do_unit_test; return; end

switch nargin
    case 2
        fmt = detect_format(fname);
        if isempty( fmt )
            error('file format unspecified, can`t autodetect');
        end
    case {3,4}
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
       error('Usage: eidors_saveimg( img , fname, format,pp )');       
end


switch fmt
    case 'igt'
        mceit_saveimg( img, fname );
    case 'e3d' % Will be deprecated
        native_saveimg( img, fname);
    case 'zri.mat' % Will be deprecated
        zriDmat_saveimg( img, fname,pp);
    otherwise
        error('eidors_readdata: file "%s" format unknown', fmt);
end




%%
function fmt = detect_format( fname ) 

dotpos = find(fname == '.');
if isempty( dotpos )
    fmt = [];
else
    dotpos= dotpos(1);
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

function zriDmat_saveimg( img, fname,pp);
   imgs = calc_slices(img);
   data = fill_in_params_zriDmat(pp, imgs);
   save(fname, 'data');

function po = fill_in_params_zriDmat(pi, imgs);
   sz_imgs = size(imgs);
   po.patient.ROI.RightLung =zeros(sz_imgs([1,2]));
   po.patient.ROI.LeftLung = zeros(sz_imgs([1,2]));
   po.patient.ROI.Heart =    zeros(sz_imgs([1,2]));

   % put to dummy because they are missing
   po.patient.halfChest = 'NaN';
   po.patient.height = 'NaN';
   po.patient.weight = 'NaN';
   po.patient.gender = 'NaN';

   po.measurement.Position.transversal = zeros(1, sz_imgs(3));
   po.measurement.Position.longitudinal= zeros(1, sz_imgs(3));
   po.measurement.ImageQuality =      100*ones(1,sz_imgs(3));
   po.measurement.ElectrodeQuality =     zeros(sz_imgs(3),32);
   po.measurement.ZeroRef = imgs;

   po.injctionPattern= 'NaN';
   po.SensorBelt.NumEl= 'NaN';

   po.measurement.CompositValue=squeeze(sum(sum(imgs,2),1));

   po = replacefields(po, pi, 1 );

function po = replacefields(po, pi, d )
% fprintf('%d:',d);
  if ~isstruct(pi)
    po = pi;
%   disp('=');
    return
  end
  for ff= fieldnames(pi)'; fn= ff{1};
%   disp(fn);
    if isfield(po,fn);
       po.(fn) = replacefields(po.(fn), pi.(fn), d+1 );
    else 
       po.(fn) = pi.(fn);
    end
  end

% The 'native' image format is deprecated - May 2019 - aa
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

function do_unit_test
   load montreal_data_1995;
   imdl = mk_common_model('c2t3',16);
   imdl.hyperparameter.value = 0.1;
   imgr = inv_solve(imdl,zc_resp(:,1),zc_resp);
   imgr.calc_colours.npoints = 32;
   p.imageRate = 7;
   p.patient.ROI.Inside = 100*ones(32);
   show_slices(imgr);
   eidors_saveimg(imgr,'mtldata.zri.mat','zri.mat',p);

   load 'mtldata.zri.mat'
   unit_test_cmp('zri.mat',eidors_var_id(data.patient), ...
         'id_62ECCCE89770EB9BB62D825054E4CE475C243C2A');
   unit_test_cmp('zri.mat',eidors_var_id(data.SensorBelt), ...
         'id_E57BF6A94F6263AAA469AD2B3A70631FD0FCE2A6');
   unit_test_cmp('zri.mat',eidors_var_id(data), ...
         'id_EEBA193E428A3A27E21447155A139BBD1E33CB7B');
