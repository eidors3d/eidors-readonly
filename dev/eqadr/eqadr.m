function [imdl_comp, removed, scores] = eqadr(dataIn, imdl, opt)
% -------------------------------------------------------------------------
% DESCRIPTION:
%   [imdl_comp, removed, scores] = eqadr(dataIn, imdl, opt)
% -------------------------------------------------------------------------
% PARAMETERS:
%   dataIn (array, cell, or struct):
%       array:
%       	first output argument of eidors_readdata.
%       cell:
%           contains multiple EIT data arrays
%       struct:
%           each field of the struct is a file identifier that contains
%           .eit.data or .data subfields. The loading function will search
%           for (ID).eit.data subfield first. If this is not found, it will
%           attemp to load data from (ID).data subfield, for each field in
%           dataIn
%   imdl (EIDORS inverse model object):
%       Inverse model object generated with keep_intermediate_results =
%       true
%   opt (struct):
%       recognized fields:
%           thresh (float):
%               score threshold for data rejection between 0 and 1. Default
%               is 0.25
%           type (string):   
%               'elec' to find the worst n electrodes, or 'meas' to find
%               the worst n measurements. Default is 'elec'
%           n (int):      
%               Number of electrodes to be excluded. The default is the
%               maximum number of electrodes or measurements (excludes all
%               above threshold)
% -------------------------------------------------------------------------   
% RETURNS:
%   imdl_comp (EIDORS inverse model object):
%       inverse model whose reconstruction matrix has been modified with a
%       series of rank one updates that compensate for the removal of the
%       removed electrodes.
%   worst (array):
%       Contains the worst n electrodes or measurements (specifiable) by
%       number, sorted from worst to best.
%   scores (array):
%       The score associated with each electrode or measurement in worst.
%       Score are the proportion of measurements that were identified as
%       faulty when the electrode was measuring, or the proportion of
%       faulty observations for each individual measurement.
% -------------------------------------------------------------------------   
% AUTHOR:
%   Mark Campbell
%   Carleton University
%   markacampbell@sce.carleton.ca
% -------------------------------------------------------------------------
% VERSION:
%   1.0.0
% -------------------------------------------------------------------------
% (C) 2019-2021 Mark Campbell. License: GPL version 2 or version 3
% -------------------------------------------------------------------------

% check options and set parameters
params = parse_options(opt, imdl);

% load and format data
data = eqadr_load_data(dataIn);

% perform data quality assessment
[worst, scores] = worst_n_elecs(data, imdl, params);
idx = scores >= params.thresh;
removed = worst(idx);
scores = scores(idx);

% perform measurement rejection
imdl_comp = comp_RM_bad_elec(imdl, removed, params.type);

end % end function

% ======================================================================= %

function params = parse_options(optIn, imdl)
    assert(isfield(imdl, 'solve_use_matrix'), "The inverse model must contain the 'solve_use_matrix' field. Re-generate the inverse model after setting the 'keep_intermediate_results' field in the inverse model options struct to true.");
    params = optIn;
    params.nElecs = length(imdl.fwd_model.electrode);
    params.mm = find(imdl.fwd_model.meas_select);
    params.nMeas = sum(imdl.fwd_model.meas_select);
    if ~isfield(optIn, 'thresh')
        params.thresh = 0.25;
    end % end if
    if ~isfield(optIn, 'type')
        params.type = 'elec';
    end % end if
    if ~isfield(optIn, 'n')
        params.n = inf;
    end % end if
    if strcmp(params.type, 'elec')
        params.n = min( [params.nElecs, params.n] );
    elseif strcmp(params.type, 'meas')
        params.n = min( [params.nMeas, params.n] );
    end % end if
end

% ======================================================================= %

function dataOut = eqadr_load_data(dataIn)
    dataOut = {};
    if isstruct(dataIn)
        fn = fieldnames(dataIn);
        nFiles = 0;
        idx = 1;
        for i = 1:length(fn)
            try
                dataOut{idx} = dataIn.(fn{i}).eit.data;
            catch
                try
                    dataOut{idx} = dataIn.(fn{i}).data;
                catch
                    continue
                end
            end
            idx = idx + 1;
            nFiles = nFiles + 1;
        end % end for i
    elseif iscell(dataIn)
        dataOut = dataIn;
    else
        dataOut = {dataIn};
    end % end if
end
