function [eitdata,s] = EITReadData(basename,filename)
%READEITDATA   Wrapper function to read any EIT data and return a
%structured output containing the data and meta information.
%   
% Signature:
% function [eitdata,s] = EITReadData(filename)
%
% Input:
% filename  char        1xN         
% 
% Output:   
% eitdata   struct      scalar      .structure
%                                   .subject_ID
%                                   .subject_type
%                                   .case_date
%                                   .file_name
%                                   .eit_device
%                                   .number_of_frames
%                                   .frame_rate
%                                   .voltage_data
%
% s         boolean     scalar      errors present: true
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
        help('EITReadData')
        display('Error EITReadData: invalid arguments');
        error('Error EITReadData: aborting execution');
end

% Define function constants
K1=2;           % file extension length
K2=basename;    % official data analysis directory
K3=13;          % Default Viasys frame rate (issue warning)
filename = [basename,'/',filename]; 

% Open file for reading in little-endian form
fid = fopen(filename,'r','l');
if fid == -1
    s = true;
    display('Error EITReadData: file could not be opened');
    error('Error EITReadData: aborting execution');
end
fclose(fid);

% Determine subject ID
lengthK2 = length(K2);
if strncmp(K2,filename,lengthK2) == 1
    idstring = filename(lengthK2+1:end);
    [studycode,idrem1] = strtok(idstring,'/');
    origin = GetSubjectType(studycode);
    [subjectid,idrem2] = strtok(idrem1(2:end),'/');
    [casedate,idrem3] = strtok(idrem2(2:end),'/');
    file = fliplr(strtok(fliplr(idrem3),'/'));
end

% Determine file type and get data
switch filename(end-K1:end)
    case 'get'
        device = 'Viasys';
        vdtemp = eidors_readdata(filename);
        [vd,s] = FormatViasysData(vdtemp);
        fr = K3;
        % Report warning
        display('Warning EITReadData: Viasys frame rate set to default (13 Hz)'); 
    case 'eit'
        device = 'Draeger';
        [fr,s] = ReadDraegerHeader(filename);
        [vdtemp,s] = ReadDraegerFile(filename);
        [vi,vq,s] = GetQuadratureSignals(vdtemp);
        vd = ApplyQuadratureSignals(vi,vq,'vi');
    otherwise
        % Report error with abort
        display('Error EITReadData: unknown file extension');
        s = true;
        error('Error EITReadData: aborting execution');
end

eitdata = struct('structure','eitdata',...
    'subject_id',subjectid,...
    'subject_type',origin,...
    'case_date',casedate,...
    'file_name',file,...
    'eit_device',device,...
    'number_of_frames',size(vd,2),...
    'frame_rate',fr,... % Will require an ASCII file header reader
    'voltage_data',vd);

% End of function
end %function
