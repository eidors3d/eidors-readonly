function [signal,s] = EITCalcTimeSignal(eitstruct,roi,graph)
%EITCALCTIMESIGNAL   Calculates and graphs the time signal or either the
%voltage data in mV or the image data in fractional change from the mean.
%
%If eitstruct = eitdata, then roi is ignored and the global voltage is
%calculated.
%If eitstruct = eitimages, then the mean time signal on pixels in the roi
%are plotted as percentage change from the mean. If roi is empty, a global
%image time signal is calculated. Currently, only the difference image
%series will be plotted.
%
% Signature:
% function [signal,s] = EITCalcTimeSignal(eitstruct,roi,graph)
%
% Input:
% eitstruct struct      scalar      either eitdata or eitimages
% (roi)     double      Mx2         [x1,y1;x2,y2;...;xM,yM] for now M=1
% (graph)   boolean     scalar      graph nu (default true)
%
% Output:
% signal    double      1XN         time signal in V or % of mean
% s         boolean     scalar      errors present: true
%
% Copyright C. Gomez-Laberge, December 2010.
% $Id: $

% Set error status to 'no errors present'
s = false;
% Save calling path
callpath = pwd;

% Check arguments
switch nargin
    case 1
        % Set (roi) to default value
        roi = []; % global signal
        % Set (graph) to default value
        graph = true;
    case 2
        % Set (graph) to default value
        graph = true;
    case 3
        % Do nothing
    otherwise
        % Error
        s = true;
        help('EITCalcTimeSignal')
        display('Error EITCalcTimeSignal: invalid arguments');
        error('Error EITCalcTimeSignal: aborting execution');
end

K1 = 1000; % Convert to millivolts for plot for eitdata
K2 = 100;  % Convert to percentile units for eitimages

% Check data type
switch eitstruct.structure
    case 'eitdata'
        signal = mean(eitstruct.voltage_data,1);
        if graph
            if strcmp(eitstruct.eit_device,'Draeger')
                K1 = 10; % Draeger operates on mV
            end
            signallength = length(signal);
            maxtime = (signallength-1)/eitstruct.frame_rate;
            time = linspace(0,maxtime,signallength);
            hdl=figure;
            plot(time,K1*signal,'-k');
            figname = [eitstruct.subject_id ' | '...
                eitstruct.case_date ' | '...
                eitstruct.file_name];
            set(hdl,'Name',figname);
            ylabel('\eta(t) (mV)')
            xlabel('time (s)')
            % Enable special data tip function to see seconds and indices
            dcm_obj = datacursormode(hdl);
            set(dcm_obj,'UpdateFcn',@mydatatipgraph)
            if isfield(eitstruct,'tidalindices')
                maxima = eitstruct.tidalindices(1,:);
                minima = eitstruct.tidalindices(2,:);
                hold on
                plot(time(maxima),K1*signal(maxima),'ks','MarkerFaceColor','b')
                plot(time(minima),K1*signal(minima),'ko','MarkerFaceColor','r')
                title(sprintf('Global Voltage Signal (Mean = %0.2f mV, Range = %1.2f mV)\n', ...
                    mean(K1*signal), K1*mean(signal(maxima)-signal(minima))));
            else
                title(sprintf('Global Voltage Signal (Mean = %0.2f mV, Range = %1.2f mV)\n', ...
                    mean(K1*signal), K1*(max(signal)-min(signal))));
            end
        end
    case 'eitimages'
        [xtab,index] = TabulateImageData(eitstruct.difference_image_series,...
            eitstruct.image_mask);
        if isempty(roi)
            roiindex = [1:nnz(eitstruct.image_mask)]'; % the global signal
            globalsignal = true;
        else
            % Convert roi pixel to column index
            roitemp = Coord2Index(roi,size(eitstruct.image_mask));
            maskindex = find(index(:,3));
            roiindex = find(roitemp(1) == maskindex);
            if isempty(roiindex)
                s = true;
                help('EITCalcTimeSignal')
                display('Error EITCalcTimeSignal: pixel coordinates out of bounds');
                error('Error EITCalcTimeSignal: aborting execution');
            end
            globalsignal = false;
        end
        signal = mean(xtab(roiindex,:),1);
        if graph
            signallength = length(signal);
            maxtime = (signallength-1)/eitstruct.frame_rate;
            time = linspace(0,maxtime,signallength);
            hdl=figure;
            plot(time,signal,'-b');
            figname = [eitstruct.subject_id ' | '...
                eitstruct.case_date ' | '...
                eitstruct.file_name];
            set(hdl,'Name',figname);
            ylabel('z(t) (%\Delta from reference impedance)')
            xlabel('time (s)')
            % Enable special data tip function to see seconds and indices
            dcm_obj = datacursormode(hdl);
            set(dcm_obj,'UpdateFcn',@mydatatipgraph)
            if isfield(eitstruct,'tidalindices')
                maxima = eitstruct.tidalindices(1,:);
                minima = eitstruct.tidalindices(2,:);
                hold on
                plot(time(maxima),signal(maxima),'ks','MarkerFaceColor','b')
                plot(time(minima),signal(minima),'ko','MarkerFaceColor','r')
                if globalsignal
                    title(sprintf('Global Image Impedance Change (Mean = %0.2f %%, Range = %0.2f %%)\n',...
                        mean(signal), max(signal)-min(signal)));
                else
                    title(sprintf('Pixel (%d,%d) Impedance Change (Mean = %0.2f %%, Range = %0.2f %%)\n',...
                        roi(1),roi(2),mean(signal), max(signal)-min(signal)));
                end
            else
                if globalsignal
                    title(sprintf('Global Image Impedance Change (Mean = %0.2f %%, Range = %0.2f %%)\n',...
                        mean(signal), max(signal)-min(signal)));
                else
                    title(sprintf('Pixel (%d,%d) Impedance Change (Mean = %0.2f %%, Range = %0.2f %%)\n',...
                        roi(1),roi(2),mean(signal), max(signal)-min(signal)));
                end
            end
        end
    otherwise
        s = true;
        help('EITCalcTimeSignal')
        display('Error EITCalcTimeSignal: unrecognized eitstruct object');
        error('Error EITCalcTimeSignal: aborting execution');
end


% End of function
end %function