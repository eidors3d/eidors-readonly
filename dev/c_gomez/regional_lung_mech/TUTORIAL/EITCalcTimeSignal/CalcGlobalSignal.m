function [eta,s] = CalcGlobalSignal(eitdata,graph)
%CALCGLOBALSIGNAL   Calculates and graphs the mean voltage signal over
%time.
%
% Signature:
% function [eta,s] = CalcGlobalSignal(eitdata,graph)
%
% Input:
% eitdata   struct      scalar      .subject_ID
%                                   .subject_type
%                                   .case_date
%                                   .file_name
%                                   .eit_device
%                                   .number_of_frames
%                                   .frame_rate
%                                   .voltage_data
%
% (graph)   boolean     scalar      graph nu (default true)
%
% Output:
% eta       double      1XN         mean voltage signal over time (V)
% s         boolean     scalar      errors present: true
%
% Copyright C. Gomez-Laberge, November 2010.
% $Id$

% Set error status to 'no errors present'
s = false;
% Save calling path
callpath = pwd;

% Check arguments
switch nargin
    case 1
        % Set (graph) to default value
        graph = true;
    case 2
        % Do nothing
    otherwise
        % Error
        s = true;
        help('EITCalcGlobalSignal')
        display('Error EITCalcGlobalSignal: invalid arguments');
        error('Error EITCalcGlobalSignal: aborting execution');
end

K1 = 1000; % Convert to millivolts for plot

eta = mean(eitdata.voltage_data,1);
if graph
    signallength = length(eta);
    maxtime = (signallength-1)/eitdata.frame_rate;
    time = linspace(0,maxtime,signallength);
    hdl=figure; 
    plot(time,K1*eta,'-k');
    figname = [eitdata.subject_ID ' | '...
        eitdata.case_date ' | '...
        eitdata.file_name];
    set(hdl,'Name',figname);
    ylabel('\eta(t) (mV)')
    xlabel('time (s)')
    % Enable special data tip function to see seconds and indices
    dcm_obj = datacursormode(hdl);
    set(dcm_obj,'UpdateFcn',@mydatatipgraph)
    if isfield(eitdata,'tidalindices')
        maxima = eitdata.tidalindices(1,:);
        minima = eitdata.tidalindices(2,:);
        hold on
        plot(time(maxima),K1*eta(maxima),'ks','MarkerFaceColor','b')
        plot(time(minima),K1*eta(minima),'ko','MarkerFaceColor','r')
        title(sprintf('Mean = %0.2f mV, Range = %1.2f mV.\n', ...
            mean(K1*eta), K1*mean(eta(maxima)-eta(minima))));        
    else
        title(sprintf('Mean = %0.2f mV, Range = %1.2f mV.\n', ...
            mean(K1*eta), K1*(max(eta)-min(eta))));
    end
end

% End of function
end %function