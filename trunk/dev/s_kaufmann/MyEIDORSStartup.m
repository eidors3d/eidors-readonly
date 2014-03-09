try 
    % Check if EIDORS runs already
    mk_stim_patterns(1, 1, 1, 1, '', 1);
    
    % When we are here EIDORS runs and we can return
    return;
catch
    %% Save Path and jump to EIDORS path
    curdir = cd;
    cd C:\Repos\eidors\eidors

    %% Disable Warnings
    warning off MATLAB:dispatcher:nameConflict

    %% Start EIDORS and increase cash
    eidors_startup({'s_kaufmann'})
    eidors_cache('cache_size', 4*1024*1024*1024);
    eidors_cache

    %% Restore old path
    cd(curdir)
    clear curdir
end;