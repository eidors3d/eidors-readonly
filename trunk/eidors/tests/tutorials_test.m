function [errors warnings] = tutorials_test(directory,tcount)
% TUTORIAL_TEST: Run all tutorials
% [errors, warnings] = tutorial_test(directory) will run all *.m files in 
% the specified directory (../../htdocs/tutorial/ if omitted). The status
% of each file will be printed to the screen as well as in a
% tutorial_status.txt file in the current directory (using the diary
% functionality). Delete tutorial_status.txt between succesive calls to
% avoid appending. 
% OUTPUTS: 
%  errors -- a cell array containing the files that caused errors
%            (including the directory and the first error message)
% warnings -- a cell array containing the files that caused warnings 
%            (including directories)

% (C) 2011 Bartlomiej Grychtol. License: GPL version 2 or version 3
% $Id$

if nargin<1 || isempty(directory)
    directory = '../../htdocs/tutorial/';
    %directory = cd;
end
diary tutorial_status.txt

global my_dir; my_dir = cd;
cd(directory);
global tut_dir; tut_dir = cd;
global D; D = genpath('.');
global d;

warning off all
eidors_msg('log_level',0);
pause off

% a little consistency would have been nice!
global tut_dlm;
if strcmp(computer,'PCWIN') || strcmp(computer,'PCWIN64')
    tut_dlm = ';';
else
    tut_dlm = ':';
end

global tut_count; tut_count = 0;
global last_tut;
if nargin >1,
    last_tut = tcount;
else last_tut = 0;
end
[d,D]= strtok(D,tut_dlm);
while ~isempty(d)
    if ~isempty(strfind(d,'.svn')); 
        [d,D]= strtok(D,tut_dlm); 
        continue; 
    end;
    cd(d);
    global F;F = dir('*.m');
    F = struct2cell(F); F = sortrows(F(1,:));% assume sorted by name
    global errors; errors={};
    global e_count; e_count = 1;
    global w_count; w_count = 1;
    global warnings; warnings={};


    while length(F) > 0
        tutname = F{1};
        if any(strcmp( tutname, skiplist)); break; end
        %tutname(end-3) = [];
        tutname(end-2) = '*'; % assume tutorials differ by one char
        T = dir(tutname);
        while length(T) > 0
            calc_colours('defaults');
            if length(T(1).name)>length(tutname)+1 % allow some leeway
                T(1) = [];
                continue
            end
            name = T(1).name(1:end-2);
            skip = 0;
            if isfunction(name)
                skip = 1;
            else
                tut_count = tut_count +1;
            end
            if skip || tut_count < last_tut 
                T(1) = [];
                F(1) = [];
                continue
            end
            fprintf([num2str(tut_count,'%04d') '  ' strrep(d,'\','\\') '\\' name]);
            lastwarn('');
            save tmp
            try
                evalin('base',evalc(name));
                if ~isempty(lastwarn)
                    fprintf(' WARNING(S)\n');
                    warnings{w_count,1} = T(1).name;
                    warnings{w_count,1} = d;
                    w_count = w_count+1;
                else
                    fprintf(' OK\n');
                end
                load tmp
                delete('tmp.mat');
                global warnings w_count;
            catch
                load tmp
                delete('tmp.mat');
                errors{e_count,1} = T(1).name;
                L = lasterror();
                errors{e_count,2} = d;
                errors{e_count,3} = L.message;
                e_count = e_count + 1;
                fprintf(' ERROR = (%s)\n',L.message);
            end

            T(1) = [];
            F(1) = [];

        end
        clear; clf; %close all; %close all crashes vnc :(
        global F tut_count e_count errors warnings w_count d D my_dir tut_dir tut_dlm last_tut
    end
    cd(tut_dir);
    [d,D]= strtok(D,tut_dlm);
end
% if numel(errors)>0
%     errors
% end
warning on all
eidors_msg('log_level',2);
pause on
clear global F tut_count e_count errors warnings w_count
cd(my_dir)
diary off

function sl = skiplist
  sl = {'netgen_accuracy01.m', ...
        'common_models01.m', ...
        'Script_For_Tutorial.m', ...
        'cg_ards_recruitment_01.m'};

function flag = isfunction(fname)
fid = fopen([fname '.m']);
l = fgetl(fid);
fclose(fid);
t = strtok(l);
flag = strcmp(t,'function');

