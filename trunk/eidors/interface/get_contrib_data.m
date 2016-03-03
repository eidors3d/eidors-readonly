function p = get_contrib_data( contrib , file )
%GET_CONTRIB_DATA Get files from the EIDORS data_contrib repository
%   GET_CONTRIB_DATA(CONTRIB, FILE) returns the path to the FILE from
%   directory CONTRIB.
%
%   GOT_CONTRIB_DATA('list') lists the valid contributions
%
%   If the file in question is not found locally, GET_CONTRIB_DATA offers
%   the option to download it to the current directory.
%
% See also: WEBSAVE

% (C) 2016 Bartlomiej Grychtol. License: GPL version 2 or version 3
% $Id$

% check for data_contrib
have_local_data = exist(get_local_path,'dir');

if nargin>0 && ischar(contrib) && strcmp(contrib,'list')
   disp(list_directories);
end


ls = list_directories(have_local_data);

if ~ismember(contrib,ls)
   error('Unknown contribution %s.',contrib);
end

if have_local_data
   p = [get_local_path filesep contrib filesep file];
   if ~exist(p, 'file')
      fprintf('File %s not found.\n',p);
      p = download_file(contrib,file,[get_local_path filesep contrib]);
   end
else % attempt to get from the web
   fprintf('The requested file %s is absent and needs to be downloaded from the web\n',file);
   p = download_file(contrib,file,cd);
end
end

function p=download_file(contrib,file,path)
   p = [get_remote_address '/' contrib '/' file];
   while 1
      s = input('Try to download now? Y/N [Y]','s');
      if isempty(s), s='Y'; end
      switch s
         case {'n','N'}
            fprintf('The file can be obtained from <a href="%s">%s</a>\n',p,p);
            break
         case {'y','Y'}
            try
               p = websave([path filesep file],p);
            catch
               fprintf('Download of <a href="%s">%s</a> failed!\n',p,p);
               try delete([path filesep file '.html']); end % cleanup
               p = [];
            end
            break
         otherwise
            fprintf('Response not understood\n');
      end
   end
end

function p = get_local_path
   % where am i
   p = fileparts(mfilename('fullpath'));
   p = strrep(p,['eidors' filesep 'interface'],'');
   p = [p 'htdocs' filesep 'data_contrib'];
end

function p = get_remote_address
   p = 'http://eidors.org/data_contrib';
end

function ls = list_directories(have_local_data)
   if nargin==0, have_local_data=false; end
   if have_local_data
      d = dir(get_local_path);
      d(1:2) = []; % . and ..
      isdir = cell2mat({d.isdir});
      ls = {d.name};
      ls = ls(isdir)';
   else % static list
      ls = {'Dargaville-ICM2010'
         'ab_2d_thorax_model'
         'at-head-mesh'
         'at-thorax-mesh'
         'bb-human-arm'
         'cg-2012-ards-recruitment'
         'cg_deforming_tank_phantom'
         'cg_normal_breathing_2006'
         'chb_child_acute_lung_injury'
         'chb_hfov_porcine'
         'db_backproj_matrix'
         'dg_geophysical_EIT'
         'ds_tank_phantom'
         'evaluation-phantom'
         'gh_pleural_cavity'
         'if-experimental-lung-injury'
         'if-neonate-spontaneous'
         'if-peep-acute-lung-injury'
         'jn_chest_phantom'
         'yg-ventilation-colourmap'
         };
   end

end