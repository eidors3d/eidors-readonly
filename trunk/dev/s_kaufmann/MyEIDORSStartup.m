%% Save Path and jump to EIDORS path
curdir = cd;
cd C:\Repos\eidors\eidors

%% Disable Warnings
warning off MATLAB:dispatcher:nameConflict

%% Start EIDORS and increase cash
eidors_startup({'s_kaufmann'})
eidors_cache('cache_size', 8*1024*1024*1024);
eidors_cache

eidors_msg('log_level',3);

%% Restore old path
cd(curdir)
clear curdir
