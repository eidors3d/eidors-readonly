function datafile_utility(inputFile, outputFile, opt) 
% DATAFILE_UTILITY(INPUTFILE, OUTPUTFILE, OPT) 
% utility to operate on EIT data files
%
% INPUTFILE:  Name of input file
% OUTPUTFILE: Name of output file
% opt.action:
%      LQ4_splitEITfile
%      opt.range = [startFrame, endFrame]
%      opt.transform = fcn_to_call on each frame
% 
% (C) 2018 Andy Adler and Beat Mueller
% License: GPL version 2 or version 3
% $Id$  

function LQ4_splitEITfile(filename, newFilename, opt)
% split EIT file
%
% filename:     name of file
% newFilename:  store file here
% opt.range:            [startFrame endFrame]
% opt.transform:        function to call 
%
% header is same as filename
% if endFrame is bigger than number of frames, only available frames are 
% copied
% only consider fileformat 4

    fid = fopen(filename, 'r', 'ieee-le', 'UTF-8');
    
    format_version = fread(fid, 1, 'int32', 'ieee-le');
    
    if format_version ~= 4
        error('wrong file version')
    end
    
    header_size = fread(fid, 1, 'uint32', 'ieee-le');
    
    tsStart = fread(fid, 1, 'uint64', 'ieee-le');
    tsEnd = fread(fid, 1, 'uint64', 'ieee-le');
    nFrames = fread(fid, 1, 'uint32', 'ieee-le');
    
    % check range:
    if length(opt.range) ~= 2
        error('range has wrong size');
    end
    
    if opt.range(1) > nFrames || opt.range(1) > opt.range(2)
        error('range not well defined');
    end
    
    % check header_size
    if header_size < 2664
        disp('assume header size of 2664');
        header_size = 2664;
    end
    
    % create new file and copy input
    fNewId = fopen(newFilename, 'w', 'ieee-le', 'UTF-8');
    fseek(fid, 0, 'bof');
    
    tmp = fread(fid, header_size, 'int8', 'ieee-le');
    fwrite(fNewId, tmp, 'int8', 'ieee-le');
    
     i = 0;
     storeFrames = 0;
     copiedFrames = 0;
     while fseek(fid, 1,'cof') ~= -1
        fseek(fid, -1,'cof');
        % need to have payload
        tStampsAbs = fread(fid,1,'int64','ieee-le'); 
        ft = fread(fid,1,'int32','ieee-le'); 
        pl = fread(fid,1,'int32','ieee-le');
        
        if ft == 0
            i = i + 1;
            
            if i >= opt.range(1) && i <= opt.range(2)
                storeFrames = 1;
                copiedFrames = copiedFrames + 1;
            else
                storeFrames = 0;
            end
        end
        
        if storeFrames == 1
            % put it into new Files
            fwrite(fNewId, tStampsAbs, 'int64','ieee-le'); 
            fwrite(fNewId, ft, 'int32', 'ieee-le');
            fwrite(fNewId, pl, 'int32', 'ieee-le');
            if pl > 0
                payLoad = fread(fid, pl, 'int8', 'ieee-le');
                fwrite(fNewId, payLoad, 'int8', 'ieee-le');
            end
        elseif pl > 0
            % go to the next frame
            fseek(fid,pl,'cof'); 
        elseif pl < 0
            disp('payloadsize < 0. try to continue');
        end
    end      
    % readjust number of Frame
    fseek(fNewId, 4 + 4 + 8 + 8, 'bof');
    fwrite(fNewId, copiedFrames, 'uint32', 'ieee-le');


    fclose(fid);
    fclose(fNewId);
    
