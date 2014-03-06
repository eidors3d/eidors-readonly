% Function: ShowPattern
%
% Takes the EIDORS structure stim and meas_sel together with other settings
% to gernerate a CommandFile to accquire Data at once.
%
% Input:    stim - EIDORS stim structure
%           meas_sel - EIDORS meas_sel vector
%           CMDFileName - File name and path of the CMD file
%           SystemSettings - System settings struct with the following members
%               Gain_DAC_PGA - 
%               Gain_Diff_PGA - 
%               UseFFT - When one FFT acquisition is used instead of RAW ADC acquisition
%               FramesToRead - Frames to read per data acquisition
%               NumberOfFFTsToAverage 
%               Delay
%
% Return:   CMD - Command Line Arguments which need to be passed to the
%           embedded system. Organized as Cell arrays (one array per
%           current injection)
%
%           FramesToRead - Vector of FramesToRead per Current injection
%
%
% Todo:     meas_sel is not evaluated!
%
% Written by Steffen Kaufmann <sk@steffen-kaufmann.com>
% March 2014

function  ShowPattern(stim, meas_sel)
    % Constants
    MeasurementNumber = 0;
    NumberOfStimulations = length(stim);
        
    %% Measurements
    for i=1:NumberOfStimulations
        fprintf('\nNo. C+ C- V+ V-\t\tData set - Current Cycle %i\n', i);
        
        % Find current mux setting
        [~, CurrentElectrodePlus] = max(full(stim(i).stim_pattern));
        [~, CurrentElectrodeMinus] = min(full(stim(i).stim_pattern));
                
        VoltageMeasurements = size(full(stim(i).meas_pattern), 1);
                
        for k=1:VoltageMeasurements
            % Check if any stimulation is there
            if(sum(abs(full(stim(i).meas_pattern(k, :)))) == 0)
                continue;
            end;
            
            % Find voltage mux setting
            [~, VoltageElectrodePlus] = max(full(stim(i).meas_pattern(k, :)));
            [~, VoltageElectrodeMinus] = min(full(stim(i).meas_pattern(k, :)));
            
            % Take Measurments
            MeasurementNumber = MeasurementNumber + 1;

            % Current Electrode Plus are the odd Channels (1, 3, ...)
            % Current Electrode Minus are the even Channels (2, 4, ...)
            if (mod(CurrentElectrodePlus, 2)) % Even  => Change order
                MUX.CP = (CurrentElectrodePlus + 1) / 2;
                MUX.CM = (CurrentElectrodeMinus / 2);
            else % Odd
                MUX.CP = (CurrentElectrodeMinus + 1) / 2;
                MUX.CM = (CurrentElectrodePlus / 2);
            end;

            % Prepare Mux variables
            % Voltage Electrode Plus are the even Channels (2, 4, ...)
            % Voltage Electrode Minus are the odd Channels (1, 3, ...)
            if (~mod(VoltageElectrodePlus, 2)) % Even
                MUX.VP = (VoltageElectrodePlus / 2);
                MUX.VM = (VoltageElectrodeMinus + 1) / 2;
            else % Odd => Change order
                MUX.VP = (VoltageElectrodeMinus / 2);
                MUX.VM = (VoltageElectrodePlus + 1) / 2;
            end;           
                          
            fprintf(['%2.u: %2.2i %2.2i %2.2i %2.2i\t\t'...
                '%3.3u\n'],...
                k, MUX.CP*2-1, MUX.CM*2, MUX.VP*2, MUX.VM*2-1, ...
                MeasurementNumber);                   
        end;
    end;
end