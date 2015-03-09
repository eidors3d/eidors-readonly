% CUSTOM EIDORS STARTUP FILE
% this must be a script so that octave can source it


try
    eidors_startup
catch
    [message id] = lasterr;
    if strcmp(id,'MATLAB:infoXmlValidationError')
        % An old version of MATLAB that we don't know how to write
        % info.xml for
    else
        error.message = message;
        error.identifier = id;
        rethrow(error);
    end
end