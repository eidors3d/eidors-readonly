function [errors warnings] = tutorials_test(directory)
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

if nargin<1
    directory = '../../htdocs/tutorial/';
    %directory = cd;
end
diary tutorial_status.txt

global my_dir; my_dir = cd;
cd(directory);
global tut_dir; tut_dir = cd;
global D; D = genpath('.');

warning off all
eidors_msg('log_level',0);
pause off

    


[d,D]= strtok(D,';');
while ~isempty(d)
    if ~isempty(strfind(d,'.svn')); [d,D]= strtok(D,';'); continue; end;
    cd(d);
    global F;F = dir('*.m');
    F = struct2cell(F); F = sortrows(F(1,:));% assume sorted by name
    global tut_count; tut_count = 1;
    global errors; errors={};
    global e_count; e_count = 1;
    global w_count; w_count = 1;
    global warnings; warnings={};


    while length(F) > 0
        tutname = F{1};
        %tutname(end-3) = [];
        tutname(end-2) = '*'; % assume tutorials differ by one char
        T = dir(tutname);
        while length(T) > 0
            if length(T(1).name)>length(tutname)+1 % allow some leeway
                T(1) = [];
                continue
            end
            name = T(1).name(1:end-2);
            fprintf([strrep(d,'\','\\') '\\' name]);
            lastwarn('');
            try
               % evalc(name);
                if ~isempty(lastwarn)
                    fprintf(' WARNING(S)\n');
                    warnings{w_count,1} = T(1).name;
                    warnings{w_count,1} = d;
                    w_count = w_count+1;
                else
                    fprintf(' OK\n');
                end
            catch
                errors{e_count,1} = T(1).name;
                L = lasterror();
                errors{e_count,2} = d;
                errors{e_count,3} = L.message;
                e_count = e_count + 1;
                fprintf(' ERROR\n');
            end
            T(1) = [];
            F(1) = [];
            tut_count = tut_count +1;
        end
        clear; clf; close all;
        global F tut_count e_count errors warnings w_count d D my_dir tut_dir
    end
    cd(tut_dir);
    [d,D]= strtok(D,';');
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