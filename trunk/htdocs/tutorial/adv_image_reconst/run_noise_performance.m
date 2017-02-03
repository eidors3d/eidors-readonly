%% EIDORS tutorial about noise performance / hyperparameter selection
%
% Fabian Braun, December 2016
%


%% part 1 - create forward models of the human thorax
run noise_performance_01.m


%% part 2 - generate noisy EIT voltage measurements (fwd_solve)
run noise_performance_02.m


%% part 3 - create reconstruction models and select hyperparameter with different approaches
% part 3a - Gauss-Newton
clear USE_GREIT_NOT_GN;
run noise_performance_03.m

% part 3b - Gauss-Newton
USE_GREIT_NOT_GN = true;
run noise_performance_03.m

