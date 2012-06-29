function [stim,meas,s] = BuildStimulationPattern(eit_device)
%BUILDSTIMULATIONPATTERN   Defines electrode simulation and measurement
%patterns based on EIT device used.
%   
% Signature:
% function [stim,meas] = BuildStimulationPattern(eit_device)
%
% Input:
% eit_device        double      char        member of eitdata struct
% 
% Output:   
% stim      double      scalar      description here
% meas      boolean     scalar      description here
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
        % Do nothing
    otherwise
        % Error
        s = true;
        help('BuildStimulationPattern')
        display('Error BuildStimulationPattern: invalid arguments');
        error('Error BuildStimulationPattern: aborting execution');
end

% Define function constants
% Viasys parameters
VELEC   = 16;
VRING   = 1;
VINJ    = '{ad}'; 
VMEAS   = '{ad}';
VOPT1   = 'no_meas_current'; 
VOPT2   = 'no_rotate_meas';
VAMPS   = 10; %mA

% Draeger parameters
DELEC   = 16;
DRING   = 1;
DINJ    = '{ad}';
DMEAS   = '{ad}';
DOPT1   = 'no_meas_current'; 
DOPT2   = 'rotate_meas';
DAMPS   = 5; %mA

switch eit_device
    case 'Viasys'
        [stim,meas] = mk_stim_patterns(VELEC,VRING,VINJ,VMEAS,...
            {VOPT1,VOPT2},VAMPS);
    case 'Draeger'
        [stim,meas] = mk_stim_patterns(DELEC,DRING,DINJ,DMEAS,...
            {DOPT1,DOPT2},DAMPS);
    otherwise
        % Report error with abort
        display('Error BuildStimulationPattern: unknown EIT device');
        s = true;
        error('Error BuildStimulationPattern: aborting execution');
end

% End of function
end %function