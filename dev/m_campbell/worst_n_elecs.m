function [worst, scores, mmScores, stimScores] = worst_n_elecs(D, imdl, opt)
% -------------------------------------------------------------------------
% DESCRIPTION:
%   [worst, scores, mmScores, stimScores] = worst_n_elecs(D, imdl, opt)
%
% Finds the noisiest n electrodes across all data files in D based on
% number of clipped measurements occuring while that electrode was
% measuring. The scores represent the proportion of total measurements made
% using that electrode that were clipped or faulty. Scores range from 0 to
% 1.
% -------------------------------------------------------------------------
% PARAMETERS:
%   D:
%       An array of EIT data in the form given by eidors_readdata. The
%       data can be a single file, multiple files within a cell array, or
%       within the .data or .eit.data fields of a struct.
%   imdl:
%       Inverse model object
%   opt:
%       Struct of optional parameters
%       n:      Number of electrodes to be excluded
%       type:   Denotes whether to find the worst n electrodes or
%               measurements.
% -------------------------------------------------------------------------   
% RETURNS:
%   worst:
%       An array containing the worst n electrodes or measurements
%       (specifiable) by number, sorted from worst to best.
%   scores:
%       The score associated with each electrode or measurement in worst,
%       in the same order as worst. Proportion of measurements that were
%       faulty when the electrode was measuring, or proportion of faulty
%       measurements for each individual measurement.
%   mmScores:
%       The score associated with each measurement.
%   stimScores:
%       The proportion of measurements that were faulty when each electrode
%       was stimulating
% -------------------------------------------------------------------------   
% AUTHOR:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
% -------------------------------------------------------------------------
% VERSION:
%   1.2.0
% -------------------------------------------------------------------------
% (C) 2019-2020 Mark Campbell. License: GPL version 2 or version 3
% -------------------------------------------------------------------------

% find worst n electrodes across all sequences
if nargin == 2
    opt = struct;
end
data = {};
if isstruct(D)
    fn      = fieldnames(D);
    n_files = length(fn);
    for i = 1:n_files
        try
            data{i} = D.(fn{i}).eit.data;
        catch
            try
                data{i} = D.(fn{i}).data;
            catch
                keyboard;
            end
        end
    end % end for i
elseif iscell(D)
    n_files = length(D);
    data    = D;
else
    n_files = 1;
    data    = D;
end % end if

nElecs = length(imdl.fwd_model.electrode);
try
    n = opt.n;
catch
    n = 32;
end
try
    type = opt.type;
catch
    type = 'elec';
end

elecScores  = zeros(nElecs, n_files);
mmScores    = zeros( length(find(imdl.fwd_model.meas_select)), n_files );
stimScores  = zeros(nElecs, n_files);

for i = 1: n_files
    if n_files == 1
        [elecScores, mmScores(:,i), stimScores(:,i)] = elec_clip_scores(data, imdl);
    else
        [elecScores, mmScores(:,i), stimScores(:,i)] = elec_clip_scores(data{i}, imdl);
    end % end if
    elecScores(:,i) = resolve_scores(elecScores, imdl);
%         plot(elecScores); hold on; 
%         plot(stimScores(:,i));
%         plot(elec_scores(:,i));
%         keyboard;
end % end for

if n_files > 1
    elecScores  = mean(elecScores, 2);
    mmScores    = mean(mmScores, 2);
    stimScores  = mean(stimScores, 2);
end % end if

if strcmp(type, 'elec')
    [hi_lo_scores, order] = sort(elecScores, 'descend');    % elecs represents electrodes
elseif strcmp(type, 'meas')
    [hi_lo_scores, order] = sort(mmScores, 'descend');      % elecs represents measurements
end % end if

worst = order(hi_lo_scores > 0);

if length(worst) > n
    worst = worst(1:n);
end % end if

if nargout > 1
    scores = hi_lo_scores(1:n);
end % end if

end % end function


function resolved = resolve_scores(elecScores, imdl)
% If one electrode is fauly, its partners will appear faulty as well. This
% function will detect if one electrode is making its parterns look bad,
% and resolve its partners' scores
resolved    = elecScores;
nElec       = size(imdl.fwd_model.electrode, 2);
stimPairs   = zeros(nElec, 2);

for i = 1:nElec
    stimPairs(i, 1) = find(imdl.fwd_model.stimulation(i).stim_pattern == -1);
    stimPairs(i, 2) = find(imdl.fwd_model.stimulation(i).stim_pattern == 1);    
end % end for

for i = 1:nElec
    [r,c]       = find(stimPairs == i);
    partners    = unique( stimPairs(r,c) );
    partners    = partners(partners ~= i);
    for j = 1:length(partners)
        p           = partners(j);
        [r,c]       = find(stimPairs == p);
        jPartners   = unique( stimPairs(r,c) );
        jPartners   = jPartners(jPartners ~= p);
        jPartner    = jPartners(jPartners ~= i);
        if elecScores(i) >= ( elecScores(p) + elecScores(jPartner) )
            resolved(p) = elecScores(jPartner);
        end % end if
    end % end for
end % end for
% clf();plot(resolved);xlabel('Electrode Number');ylabel('Resolved Electrode Measurement Scores');xlim([0,32]);ylim([0,.8]);hold on;plot([0,32],[.25,.25],'--');
end % end function


% ======================================================================= %

function [clip_scores, mmScores, stimScores] = elec_clip_scores(data, imdl)
% clip_scores:
%     give electrodes a score based on proportion of clipped measurements
%     occuring while that electrode was measuring. Score ranges from 0 to
%     1.
% mmScores:
%     For each selected measurement, score is given based on portion of total
%     measurements that were clipped (scores range from 0 to 1).

    msel            = imdl.fwd_model.meas_select;
    nElec           = size(imdl.fwd_model.electrode, 2);
    mm              = find(msel);
    elecUsed        = zeros(length(mm),2);
    count           = 1;
    nMeasWithElec   = zeros(nElec, 1);
    stimPairs       = zeros(nElec, 2);
    
    for i = 1:nElec
        [r,c]           = find(imdl.fwd_model.stimulation(i).meas_pattern);
        stimPairs(i, 1) = find(imdl.fwd_model.stimulation(i).stim_pattern == -1);
        stimPairs(i, 2) = find(imdl.fwd_model.stimulation(i).stim_pattern == 1);
        for j = 1:max(r)
           elecUsed(count,:)    = c(r == j)';
           count                = count+ 1;
        end % end for
    end % end for

    for i = 1:nElec
        nMeasWithElec(i) = sum(elecUsed == i, 'all');
    end % end for

    data_       = data(mm,:);
    clipped     = find_clipped_agresive(data_); % mark clipped measurements.
    rDataNeg    = real(data_) < 0;              % mark measurements with negative real component.
    allMarked   = sum( (clipped*1 + rDataNeg*1) > 0, 2);
    mmScores    = allMarked ./ size(data_, 2);  % output
    measPerPair = length(mm) / nElec;

    score       = 1;
    clip_scores = zeros(nElec,1);
    stimScores  = zeros(nElec,1);               % number of clipped measurements when this electrode is stimulating
    
    for i = 1:nElec
        stop    = i * measPerPair;
        start   = stop - measPerPair + 1;
        idx     = stimPairs(i, :);
        stimScores( idx ) = stimScores( idx ) + sum( allMarked(start:stop) );
    end % end for
    stimScores = stimScores ./ (measPerPair * size(data, 2));
%     clf();plot(stimScores);xlabel('Electrode Number');ylabel('Proportion of Measurements Erroneous While Stimulating');xlim([0,32]);ylim([0,.8]);hold on;plot([0,32],[.25,.25],'--');
    while score > 0
        score = max(allMarked);
        worst = find(allMarked == score);
        allMarked(worst) = 0;

        for i=1:length(worst)
            elecs = elecUsed(worst(i),:);

            if ~isempty(elecs)
                clip_scores(elecs) = clip_scores(elecs) + score;
            end % end if

        end % end for

    end % end while

    clip_scores = clip_scores ./ (size(data,2) * nMeasWithElec);

end % end function

% ======================================================================= %

function clipped = find_clipped_agresive(data)
% find total clipped data points for each selected measurement Draw a
% circle and slowly decrease the radius. Tally how many points are outside
% of that circle. When the number of points outside the circle drops
% sharply, we can be confident that all clipped data points are outside of
% this circle. 

ALPHA = 0;
DELTA = 0.05;
ITERS = 21;

insideCirc  = zeros(ITERS, 1);
distFromC   = sqrt( real(data).^2 + imag(data).^2 );
maxDist     = max(distFromC, [], 'all');

for i = 1:ITERS
    ALPHA           = ALPHA + DELTA;
    insideCirc(i)   = sum(distFromC < maxDist * ALPHA, 'all');
end % end for

foundInIter     = diff(insideCirc);
stoppingIter    = find(diff(foundInIter) > 0, 1, 'last');
clipDist        = maxDist * (stoppingIter * DELTA);
clipped         = distFromC >= clipDist;

% clf();subplot(3,1,1);plot(insideCirc,'LineWidth',2);xlim([0,length(insideCirc)]);ylabel('Points Within Circle');
% subplot(3,1,2);semilogy(foundInIter,'LineWidth',2);ylabel('First Derivative');xlim([0,length(insideCirc)]);hold on; semilogy([stoppingIter,stoppingIter],[400,max(foundInIter)],'--k');
% subplot(3,1,3);plot(diff(foundInIter),'LineWidth',2);xlim([0,length(insideCirc)]);hold on;plot([1,length(insideCirc)],[0,0],'--r');ylabel('Second Derivative');xlabel('Iteration');plot([stoppingIter,stoppingIter],[min(diff(foundInIter)),max(diff(foundInIter))],'--k');
end % end function
