function datafile_utility(inputFile, outputFile, opt)
% DATAFILE_UTILITY(INPUTFILE, OUTPUTFILE, OPT)
% utility to operate on EIT data files
%
% INPUTFILE: Name of input file
% OUTPUTFILE: Name of output file
% opt.action:
%     LQ4_splitEITfile
%     opt.range = [startFrame, endFrame]
%     opt.transform = fcn_to_call on each frame
%          if opt.transform is numeric, multiply by it.
%     opt.addition  = add opt.addition(:,i) to frame i
%
% To replace EIT files with other content, do
%opt.range  = [1,2];
%opt.transform  = sparse(2048,2048);
%opt.addition = zeros(2048,2);
%opt.addition(1:2:2048,:) = [vh.meas,vi.meas]*1e6;
%
% In order to "shrink" EIT files, do
% PAT = [0,5];
% skip5 = {32,1,PAT,PAT,{'meas_current'},1};
% [tst.stimulation,~] = mk_stim_patterns(skip5{:});
% idx = reciprocity_idx( tst );
% mvr =  idx == max([idx,idx(idx)],[],2);
% mv2 = idx; mv2(~mvr) = [];
% mvr = find(mvr);
% mult = sparse(mvr,mvr,1,1024,1024) +  ...
%        sparse(mvr,mv2,1,1024,1024);
% skip5{4} = 'no_meas_current_next2';
% [~,msel] = mk_stim_patterns(skip5{:});
% mult(~msel,:) = 0;
% 
% opt.transform  = sparse(2048,2048);
% opt.transform(1:2:2048,1:2:2048)  = real(mult);
% %opt.transform(2:2:2048,2:2:2048)  = imaginary
% datafile_utility(infile,outfile,opt);
%
% (C) 2018 Andy Adler and Beat Mueller
% License: GPL version 2 or version 3
% $Id$


switch opt.action
   case 'LQ4_splitEITfile';
      if ~isfield(opt,'transform'); 
         opt.transform = 1; %
      end
      if ~isfield(opt,'addition'); 
         opt.addition = 0; %
      end
      splitEITfile(inputFile, outputFile, opt.range, ...
                   opt.transform, opt.addition);
   otherwise;
      error('opt.action = "%s" not recognized', opt.action);
end

function splitEITfile(filename, newFilename, range, fcall, addvec)
% split EIT file
%
% filename:     name of file
% newFilename:  store file here
% range:        [startFrame endFrame]
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
    if length(range) ~= 2
        error('range has wrong size');
    end
    
    if range(1) > nFrames || range(1) > range(2)
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
            
            if i >= range(1) && i <= range(2)
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
                vi_payload = 64;
                iq_payload = 2048;
                p_startframe = ftell(fid);
                payLoad = fread(fid, pl, 'int8', 'ieee-le');
                p_endframe   = ftell(fid);

                fseek(fid,p_startframe + 15*4 + 12 + 4*vi_payload,'bof');
                iqPayload = fread(fid,iq_payload,'int32','ieee-le');
                fseek(fid,p_endframe,'bof');

                p_startframe = ftell(fNewId);
                fwrite(fNewId, payLoad, 'int8', 'ieee-le');
                p_endframe   = ftell(fNewId);

                fseek(fNewId,p_startframe + 15*4 + 12 + 4*vi_payload,'bof');
                if isnumeric(fcall)
                   iqPayload = full(fcall*iqPayload); % transform it
                else
                   iqPayload = fcall(iqPayload); % transform it
                end
                if size(addvec,2) == 1
                   iqPayload = iqPayload + addvec;
                else
                   iqPayload = iqPayload + addvec(:,i);
                end
                fwrite(fNewId, iqPayload, 'int32', 'ieee-le');
                fseek(fNewId,p_endframe,'bof');
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
    
