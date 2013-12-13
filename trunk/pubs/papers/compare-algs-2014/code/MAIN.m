% To reproduce the analysis in the associated paper, follow these steps.

% 1. Obtain a copy of EIDORS from eidors.org
% 2. In Matlab, type "run /path/to/eidors/startup.m"
% 3. Run the following code:
warning off EIDORS:Deprecated % we renamed some files in EIDORS, it's OK
analyse; % reconstructs all data with all algorithms and stores the images
imgs_to_params; % extract measures from each image
params_stats;   % run statistical tests for expected physiological changes
% results appear as ../paper/stat_tests_out.tex
