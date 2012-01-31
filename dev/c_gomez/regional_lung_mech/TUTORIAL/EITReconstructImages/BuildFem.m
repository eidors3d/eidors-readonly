function [im,s] = BuildFem(subject_type,eit_device)
%BUILDFEM   Builds the FEM for image reconstruction.
%   
% Signature:
% function [im,s] = BuildFem(subject_type,eit_device)
%
% Input:
% subject_type      double      char        member of eitdata struct
% eit_device        double      char        member of eitdata struct
% 
% Output:   
% im                double      struct      EIDORS inv_model
% s                 boolean     scalar      errors present: true
% 
% Copyright C. Gomez-Laberge, November 2010.
% $Id: $

% Set error status to 'no errors present'
s = false;
% Save calling path
callpath = pwd;

% Check arguments
switch nargin
    case 2
        % Do nothing
    otherwise
        % Error
        s = true;
        help('BuildFem')
        display('Error BuildFem: invalid arguments');
        error('Error BuildFem: aborting execution');
end

% Define function constants
GENERATOR = 'Distmesh'; % alternative is 'Netgen'
% Distmesh generator parameters
    DENSITY = 'd2';     % Parameters for mk_common_model
    HUMAN = 't2';
    SWINE = 'c2';
    EVIASYS = 16;
    EDRAEGER = 32;
    % SIM01 =           % Simulation parameters can be added as needed
    % SIM02 =
    % SIM03 =

% Netgen generator parameters 
    % Currently unavailable (see EIDORS for details)
    
switch GENERATOR
    case 'Distmesh'
        % Select FEM shape based on subject type
        switch subject_type
            case 'Human'
                fem_string = [DENSITY HUMAN];
            case 'Swine'
                fem_string = [DENSITY SWINE];
            case 'Test'
                fem_string = [DENSITY SWINE];
            case 'Simulation'
                % Special routine will be needed to determine simulation
                display('Error BuildFem: simulation FEM undefined');
                s = true;
                error('Error BuildFem: aborting execution');
            otherwise
                % Report error with abort
                display('Error BuildFem: unknown subject_type');
                s = true;
                error('Error BuildFem: aborting execution');
        end
        % Build FEM and model electrodes according to device
        switch eit_device
            case 'Viasys'
                im = mk_common_model(fem_string,EVIASYS);
            case 'Draeger'
                im = mk_common_model(fem_string,EDRAEGER);
                im.fwd_model.electrode = im.fwd_model.electrode(2:2:32);
            otherwise
                % Report error with abort
                display('Error BuildFem: unknown EIT device');
                s = true;
                error('Error BuildFem: aborting execution');
        end
    case 'Netgen'
        % Currenty unavailable (see EIDORS for details)
    otherwise
        % Report error with abort
        display('Error BuildFem: unknown FEM generator');
        s = true;
        error('Error BuildFem: aborting execution');
end

% End of function
end %function