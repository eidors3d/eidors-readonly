function [vv, auxdata, stim ]= eidors_readdata( fname, format, frame_range, extra )
% EIDORS readdata - read data files from various EIT equipment
%    manufacturers
%
% [vv, auxdata, stim ]= eidors_readdata( fname, format, frame_range, extra )
%    frame_range is the list of frames to load from the file (ie. 1:100).
%    if the frames are beyond the number in the file, no data is returned
%    to get all frames, use frame_range= []
%
% Currently the list of supported file formats is:
%    - "EIT" Files, may be one of 
%          - Draeger Pulmovista (2008+)
%          - GoeIIMF/Carefusion (2010+)
%          - Swisstom BB2/Pioneer (2010+)
%       The function will attempt to autodetect the format
%       
%    - MCEIT (GoeIIMF / Viasys) "get" file format 
%        format = "GET" or "MCEIT"
%        Note that the data is "untwisted" to correspond to a "no_rotate_meas" stim pattern
%    - MCEIT (GoeIIMF) "get" file format
%        format = "GET-RAW"
%        Data in original order, corresponds to "rotate_meas" stim pattern
%
%    - Draeger (Pulmovista) "EIT" file format (2008+)
%        format = "DRAEGER-EIT"
%
%    - Draeger "get" file format (- 2007 - format for Draeger equipment)
%        format = "DRAEGER-GET"
%
%    - Swisstom EIT equipment "EIT" (Pioneer set and BB2):
%           'LQ1' (2010 - 2011)
%           'LQ2' (2013 - 2014)
%           'LQ4' (2015+)
%
%    - Dixtal file format, from Dixtal inc, Brazil
%        format = 'DIXTAL_encode' extract encoder from provided Dll
%           - Output: [encodepage] = eidors_readdata( path,'DIXTAL_encode');
%              where path= '/path/to/Criptografa_New.dll' (provided with system)
%        format = 'DIXTAL'
%           - output: [vv] = eidors_readdata( fname, 'DIXTAL', [], encodepage );
%    - New Carefusion "EIT" file format
%        format = "GOEIIMF-EIT" or "carefusion"

%    - Sheffield MK I "RAW" file format
%        format = "RAW" or "sheffield"
%    - ITS (International Tomography Systems)
%        format = "ITS" or "p2k"
%    - IIRC (Impedance Imaging Research Center, Korea)
%        format = "txt" or "IIRC"
%    - University of Cape Town formats
%        format = "UCT_SEQ"  UCT sequence file
%           - Output: [ stimulation, meas_select]= eidors_readdata(fname, 'UCT_SEQ')
%        format = "UCT_CAL"  UCT calibration file
%           - Output: [vv, no_cur_caldata_raw ]= eidors_readdata( fname, 'UCT_CAL' )
%                 where no_cur_caldata_raw is data captured with no current
%        format = "UCT_DATA"  UCT data frame file
%           - Output: [vv]= eidors_readdata( fname, 'UCT_DATA' )
%
% Usage
% [vv, auxdata, stim ]= eidors_readdata( fname, format )
%     vv      = measurements - data frames in each column
%     auxdata = auxillary data - if provided by system 
%     stim    = stimulation structure, to be used with
%                fwd_model.stimulation. 
%     fname = file name
%     stim.framerate = acquisition rate (frames/sec) if available
%
%  if format is unspecified, an attempt to autodetect is made

% (C) 2005-09 Andy Adler. License: GPL version 2 or version 3
% $Id$

% TODO:
%   - output an eidors data object
%   - test whether file format matches specified stimulation pattern
%   - todo MCEIT provides curr and volt on driven electrodes.
%       how can this be provided to system?

if ~exist(fname,'file')
   error([fname,' does not exist']);
end

if nargin < 2
% unspecified file format, autodetect
   dotpos = find(fname == '.');
   if isempty( dotpos ) 
      error('file format unspecified, can`t autodetect');
   else
      dotpos= dotpos(end);
      format= fname( dotpos+1:end );
   end
end


auxdata = []; % default, can be overriden if the format has it
fmt = pre_proc_spec_fmt( format, fname );
switch fmt
   case 'mceit';
      [vv,curr,volt,auxdata_out,auxtime,rawdata] = mceit_readdata( fname );
      auxdata.auxdata = auxdata_out;
      auxdata.auxtime = auxtime;
      auxdata.curr    = curr;
      auxdata.volt    = volt;


      if strcmp(lower(format),'get-raw')
          vv= rawdata(1:208,:);
          stim = mk_stim_patterns(16,1,[0,1],[0,1], {'no_meas_current','rotate_meas'}, 1);
      else
          stim = basic_stim(16);
      end
   case 'draeger-get'
      vv = draeger_get_readdata( fname );

      stim = basic_stim(16);

   case 'draeger-eit'
     [fr] = read_draeger_header( fname );
     % Currently Draeger equipment uses this pattern with 5mA injection
     stim = mk_stim_patterns(16,1,[0,1],[0,1],{'rotate_meas','no_meas_current'},.005);
     [stim(:).framerate] = deal(fr);
     [vv] = read_draeger_file( fname );
     auxdata = vv; % the actual voltages
     % Based on the draeger *get file converter
     ft = [.00098242, .00019607];% estimated AA: 2016-04-07
     vv = ft(1)*vv(1:208,:) - ft(2)*vv(322+(1:208),:);

   case {'raw', 'sheffield'}
      vv = sheffield_readdata( fname );

      stim = basic_stim(16);
   case {'p2k', 'its'}
      vv = its_readdata( fname );

      stim = 'UNKNOWN';
   case {'txt','iirc'}
      vv = iirc_readdata( fname );

      stim = basic_stim(16);
   case 'uct_seq'
      [vv,auxdata] = UCT_sequence_file( fname );

      stim = 'UNKNOWN';
   case 'uct_cal'
      [vv,auxdata] = UCT_calibration_file( fname );

      stim = 'UNKNOWN';
   case 'uct_data'
      [vv] = UCT_LoadDataFrame( fname );

      stim = 'UNKNOWN';
   case {'carefusion','goeiimf-eit'}
      [vv, auxdata_out, auxtime] = carefusion_eit_readdata( fname );
      auxdata.auxdata = auxdata_out;
      auxdata.auxtime = auxtime;
  
      stim = mk_stim_patterns(16,1,[0,1],[0,1], {'no_meas_current','rotate_meas'}, .005);

   case 'lq1'
      [vv] = landquart1_readdata( fname );

      stim = mk_stim_patterns(32,1,[0,5],[0,5],{'no_rotate_meas','no_meas_current'},.005);
      
   case {'lq2','lq3'}
      [vv, elecImps, tStampsAbs, tStampsRel] = landquart2_readdata( fname );
      auxdata.elec_impedance = elecImps;
      auxdata.t_abs = tStampsAbs;
      auxdata.t_rel = tStampsRel;

      stim = mk_stim_patterns(32,1,[0,5],[0,5],{'no_rotate_meas','no_meas_current'},.005);
     
   case 'lq4pre'
      [vv] = landquart4pre_readdata( fname );

      stim = mk_stim_patterns(32,1,[0,5],[0,5],{'no_rotate_meas','no_meas_current'},.005);

   case 'lq4'
      [vv, evtlist, elecImps, tStampsAbs, tStampsRel] = landquart4_readdata( fname );
      auxdata.event = evtlist;
      auxdata.elec_impedance = elecImps;
      auxdata.t_abs = tStampsAbs;
      auxdata.t_rel = tStampsRel;

      stim = mk_stim_patterns(32,1,[0,5],[0,5],{'no_rotate_meas','no_meas_current'},.005);
      
   case 'dixtal_encode'
      [vv] = dixtal_read_codepage( fname );
      stim = 'N/A';

   case 'dixtal'
      [vv] = dixtal_read_data( fname, frame_range, extra );
      auxdata = vv(1025:end,:);
      vv      = vv([1:1024],:);
      stim= mk_stim_patterns(32,1,[0,4],[0,4], ... 
              {'meas_current','no_rotate_meas'}, 1);

   otherwise
      error('eidors_readdata: file "%s" format unknown', format);
end

function stim = basic_stim(N_el);
   stim= mk_stim_patterns(16,1,[0,1],[0,1], ... \
          {'no_meas_current','no_rotate_meas'}, 1);

function fmt = pre_proc_spec_fmt( format, fname );
   fmt= lower(format);
   if strcmp(fmt,'get')
      if is_get_file_a_draeger_file( fname)
         fmt= 'draeger-get';
      else
         fmt= 'mceit';
      end
   end

   if strcmp(fmt,'get-raw')
      fmt= 'mceit';
   end
   

   if strcmp(fmt,'eit')
      draeger =   is_eit_file_a_draeger_file( fname );
      swisstom=   is_eit_file_a_swisstom_file( fname );
      carefusion= is_eit_file_a_carefusion_file( fname );
      if carefusion
         eidors_msg('"%s" appears to be in GOEIIMF/Carefusion format',fname,3);
         fmt= 'carefusion';
      elseif draeger
         eidors_msg('"%s" appears to be in Draeger format',fname,3);
         fmt= 'draeger-eit';
      elseif swisstom
         fmt= sprintf('lq%d',swisstom);
         if swisstom == 3.5
            fmt= 'lq4pre';
         end
         eidors_msg('"%s" appears to be in %s format',fname,upper(fmt),3);
      else
         error('EIT file specified, but it doesn''t seem to be a Carefusion file')
      end
   end

function df= is_get_file_a_draeger_file( fname)
   fid= fopen(fname,'rb');
   d= fread(fid,[1 26],'uchar');
   fclose(fid);
   df = all(d == '---Draeger EIT-Software---');

function df= is_eit_file_a_draeger_file( fname );
   fid= fopen(fname,'rb');
   d= fread(fid,[1 80],'uchar');
   fclose(fid);
   ff = findstr(d, '---Draeger EIT-Software');
   if ff;
      df = 1;
      eidors_msg('Draeger format: %s', d(ff(1)+(0:30)),4);
   else 
      df = 0;
   end

function df= is_eit_file_a_swisstom_file( fname );
   fid= fopen(fname,'rb');
   d= fread(fid,[1 4],'uchar');
   fclose(fid);
   
% Note that there are two possible LQ4 formats, either with
%     d=[0,0,0,4] or with d = [4,0,0,0]
%  we call the first one version 3.5
   if any(d(4)==[2,3,4]) && all(d([1,2,3]) == 0);
      df = d(4);
      if d(4)==4; df = 3.5; end
%     eidors_msg('Swisstom format: LQ%d', df, 4);
   elseif any(d(1:4) == [4,0,0,0])
      df = 4;
   else 
      df = 0;
   end

function df = is_eit_file_a_carefusion_file( fname )
   fid= fopen(fname,'rb');
   d= fread(fid,[1 80],'uchar');
   fclose(fid);
   df = 0;
   if d(1:2) ~= [255,254]; return; end
   if d(4:2:end) ~= 0; return; end
   str= setstr(d(3:2:end));
   tests1= '<?xml version="1.0" ?>';
   tests2= '<EIT_File>';
   if ~any(findstr(str,tests1)); return; end
   if ~any(findstr(str,tests2)); return; end
   
   df= 1;
   return

function [vv, auxdata, auxtime] = carefusion_eit_readdata( fname );
   fid= fopen(fname,'rb');
   d= fread(fid,[1 180],'uchar');
   str= setstr(d(3:2:end));
   outv = regexp(str,'<Header>(\d+) bytes</Header>','tokens');
   if length(outv) ~=1; 
      error('format problem reading carefusion eit files');
   end
   header = str2num(outv{1}{1});

% NOTE: we're throwing away the header completely. This isn't
% the 'right' way to read a file.

   fseek(fid, header, -1); % From the beginning
   d_EIT= fread(fid,[inf],'float32');

   fseek(fid, header, -1); % From the beginning
   d_struct= fread(fid,[inf],'int32');  % In the struct, meaning of every int: 1 Type, 2 Waveform channel No., 3 sequence No., 4 size in byte, 5 count
   fclose(fid);
   pt_EIT=1;               % pointer of the EIT data
   vv=[];
   auxdata=[];
   
   while (pt_EIT<=length(d_struct))
      switch d_struct(pt_EIT)
         case 3, % type=3,  EIT transimpedance measurement
           % d_struct(pt_EIT+2), Sequence No., start from 0
           vv(1:256,1+d_struct(pt_EIT+2))= d_EIT(pt_EIT+6:pt_EIT+6+255);
           pt_EIT=pt_EIT+6+256;
         case 7, %type=7, EIT image vector
           pt_EIT=pt_EIT+912+6;        % drop the image           
         case 8, %type=8, Waveform data           
            aux_seg_len = d_struct(pt_EIT+3)/4*d_struct(pt_EIT+4);
           if (d_struct(pt_EIT+1) == 0)     % waveform channel 0 is the analog input
               aux_segment = d_EIT(pt_EIT+6 : pt_EIT+6+aux_seg_len-1);
               auxdata = [auxdata; aux_segment];
           else
               % drop all other waveform channels (see Waves/Waveform in header file for what they are)
           end  
           pt_EIT=pt_EIT+6+aux_seg_len;           
         case 10,% type = 10, unknown type
           %drop
           pt_EIT=pt_EIT+6+d_struct(pt_EIT+3)/4*d_struct(pt_EIT+4);
         otherwise,
           eidors_msg('WARNING: unknown type in carefusion file type');
           pt_EIT=pt_EIT+6+d_struct(pt_EIT+3)/4*d_struct(pt_EIT+4);
       end
   end
   vv=vv(47+[1:208],:);
   auxtime = (cumsum([1 13*ones(1,15)])-1)/208;             % relative time as fractions of the EIT frame: assuming uniform sampling - not sure this is 100% correct(!)
   auxtime = reshape(repmat(1:size(vv,2), length(auxtime), 1),[],1) + repmat(auxtime, 1, size(vv,2))';
   

% Read the old <2006 Draeger "get" file format
function [vv,curr,volt,auxdata] = draeger_get_readdata( fname );
   fid= fopen(fname,'rb');
   emptyctr=0;
   while emptyctr<2
     line= fgetl(fid);
     if isempty(line)
        emptyctr= emptyctr+1;
     else
        eidors_msg('data not understood',0);
        emptyctr=0;
     end
   end
   d= fread(fid,inf,'float',0,'ieee-le');
   fclose(fid);
pause

   if rem( length(d), 256) ~=0
      eidors_msg('File length strange - cropping file',0);
      d=d(1:floor(length(d)/256)*256);
   end

   dd= reshape( d, 256, length(d)/256);
   vv= untwist(dd);

   curr=0.00512*dd(209:224,:);  % Amps
   volt=12*dd(225:240,:); %Vrms
   auxdata= dd(241:255,:);
   auxdata= auxdata(:);
    
function jnk
%[adler@adler01 sept07]$ xxd Sch_Pneumoperitoneum_01_001.get | head -30
%0000000: 2d2d 2d44 7261 6567 6572 2045 4954 2d53  ---Draeger EIT-S
%0000010: 6f66 7477 6172 652d 2d2d 0d0a 2d2d 2d50  oftware---..---P
%0000020: 726f 746f 636f 6c20 4461 7461 2d2d 2d2d  rotocol Data----
%0000030: 2d0d 0a0d 0a44 6174 653a 2020 3138 2d30  -....Date:  18-0
%0000040: 322d 3230 3034 0d0a 5469 6d65 3a20 2031  2-2004..Time:  1
%0000050: 333a 3138 2050 4d0d 0a0d 0a46 696c 656e  3:18 PM....Filen
%0000060: 616d 653a 2020 2020 2020 2020 2053 6368  ame:         Sch
%0000070: 5f50 6e65 756d 6f70 6572 6974 6f6e 6575  _Pneumoperitoneu
%0000080: 6d5f 3031 5f30 3031 2e67 6574 0d0a 4453  m_01_001.get..DS
%0000090: 5020 5365 7269 616c 204e 722e 3a20 2020  P Serial Nr.:
%00000a0: 4549 5430 322f 3035 2f30 3030 360d 0a0d  EIT02/05/0006...
%00000b0: 0a46 7265 7175 656e 6379 205b 487a 5d3a  .Frequency [Hz]:
%00000c0: 2020 2020 3937 3635 362e 330d 0a47 6169      97656.3..Gai
%00000d0: 6e3a 2020 2020 2020 2020 2020 2020 2020  n:
%00000e0: 2020 2031 3134 350d 0a41 6d70 6c69 7475     1145..Amplitu
%00000f0: 6465 3a20 2020 2020 2020 2020 2020 2031  de:            1
%0000100: 3030 300d 0a53 616d 706c 6520 5261 7465  000..Sample Rate
%0000110: 205b 6b48 7a5d 3a20 2020 2035 3030 300d   [kHz]:    5000.
%0000120: 0a50 6572 696f 6473 3a20 2020 2020 2020  .Periods:
%0000130: 2020 2020 2020 2020 2032 300d 0a46 7261           20..Fra
%0000140: 6d65 733a 2020 2020 2020 2020 2020 2020  mes:
%0000150: 2020 2020 2031 300d 0a0d 0a0d 0a
%                                       8bc33e       10........>
%0000160: 3f 0548633e bf20933d e192393d a568ea  ?.Hc>. .=..9=.h.
%0000170: 3c 2530f63c 27e6043d 2043ad3c 25ce93  <%0.<'..= C.<%..
%0000180: 3c aebcce3c 8714643d e3533d3e 65b6e1  <...<..d=.S=>e..
%0000190: 3e 7a62103f 81c4143e 8c99813d 35921d  >zb.?...>...=5..
%00001a0: 3d 8e6b0d3d 690cf93c 9910713c 3c9289  =.k.=i..<..q<<..
%00001b0: 3c f6736f3c 2291453d ad1cab3d 386f15  <.so<".E=...=8o.
%00001c0: 3e e82a143f 2a952e3f f568493e f08a8c  >.*.?*..?.hI>...
%00001d0: 3d e43e0e3d 7040253d 19f4af3c 67fd93  =.>.=p@%=...<g..

function vv = sheffield_readdata( fname );
   fid=fopen(fname,'rb','ieee-le');
   draw=fread(fid,[104,inf],'float32');
   fclose(fid);
   ldat = size(draw,2);

   [x,y]= meshgrid( 1:16, 1:16);
   idxm= y-x;
 % HW gain table
   gtbl = [0,849,213,87,45,28,21,19,21,28,45,87,213,849,0];
   idxm(idxm<=0)= 1;
   idxm= gtbl(idxm);

 % Multiply by gains
   draw = draw .* (idxm(idxm>0) * ones(1,ldat)); 

   vv= zeros(16*16, size(draw,2));
   vv(find(idxm),:) = draw;

 % Add reciprocal measurements
   vv= reshape(vv,[16,16,ldat]);
   vv= vv + permute( vv, [2,1,3]);
   vv= reshape(vv,[16*16,ldat]);

   

function [vv,curr,volt,auxdata,auxtime,rawdata] = mceit_readdata( fname );

   fid= fopen(fname,'rb');
   d= fread(fid,inf,'float');
   fclose(fid);

   if rem( length(d), 256) ~=0
      eidors_msg('File length strange - cropping file',0);
      d=d(1:floor(length(d)/256)*256);
   end

   dd= reshape( d, 256, length(d)/256);
   rawdata = dd;
   no_reciprocity = (dd(39,1)==0);      %104 measurements per frame
   if no_reciprocity
       dd=transfer104_208(dd);
   end
   vv= untwist(dd);

   curr=0.00512*dd(209:224,:);  % Amps
   volt=12*dd(225:240,:); %Vrms
   auxdata= dd(241:256,:);
   auxdata= auxdata(:);
   %input impedance=voltage./current-440;        Ohm
   
   if no_reciprocity
     % remove every 14th, 15th and 16th sample as they are always zero
     auxdata(14:16:end) = []; auxdata(14:15:end) = []; auxdata(14:14:end) = [];
     % only 104 measurements were performed with NON-UNIFORM sampling of AUX in between EIT samples
     % Sampling procedure was thought to be: 13 AUX1 13 AUX2 12 AUX3 11 AUX4 10 AUX5 9 ... 2 AUX13 1 [END OF FRAME]
%      auxtime = (cumsum([1 13 12:-1:2]) - 1) / 104;
     % However, this seems to be a little different, finally the sampling points were determined 
     % empirically by measuring a sawtooth signal (linear ramp) and fitting the timings 
     auxtime = [0.0000,0.1390,0.2535,0.3617,0.4642,0.5486,0.6367,0.7110,0.7780,0.8354,0.8865,0.9304,0.9711];  % mean fit at 25 Hz
%      auxtime = [0.0000,0.1460,0.2789,0.3895,0.4875,0.5788,0.6608,0.7356,0.7986,0.8569,0.9045,0.9454,0.9824,]; % mean fit at 44 Hz
   else
     % all 208 measurements were performed 
     auxtime = (cumsum([1 13*ones(1,15)])-1)/208;             % relative time as fractions of the EIT frame
   end
   % time of AUX signal relative to frame number (starting with 1 - Matlab style)
   auxtime = reshape(repmat(1:size(dd,2), length(auxtime), 1),[],1) + repmat(auxtime, 1, size(dd,2))';
   

function array208=transfer104_208(array104),
% The order in the 208-vector is
% Index   Inject    Measure Pair
% 1         1-2        3-4        (alternate pair is 3-4, 1-2, at location 39)
% 2         1-2        4-5
% 3         1-2        5-6
% ¡­
% 13        1-2        15-16
% 14        2-3        4-5
% 15        2-3        5-6
% ¡­
% 26        2-3        16-1
% 27        3-4        5-6
% ¡­
% 39        3-4        1-2
% 40        4-5        6-7
% ¡­

ind=[39,51,63,75,87,99,111,123,135,147,159,171,183,52,64,76, ...
     88,100,112,124,136,148,160,172,184,196,65,77,89,101,113, ...
     125,137,149,161,173,185,197,78,90,102,114,126,138,150,162, ...
     174,186,198,91,103,115,127,139,151,163,175,187,199,104,116, ...
     128,140,152,164,176,188,200,117,129,141,153,165,177,189, ...
     201,130,142,154,166,178,190,202,143,155,167,179,191,203, ...
     156,168,180,192,204,169,181,193,205,182,194,206,195,207,208];
ro=[1:13, 14:26, 27:38, 40:50, 53:62, 66:74, 79:86, 92:98, ...
    105:110,118:122,131:134,144:146,157:158,170:170];

[x,y]=size(array104);
if x~=256 && y~=256,
    eidors_msg(['eidors_readdata: expectingin an input array ', ...
                'of size 208*n']);
    return;
elseif y==256,
    array104=array104';
    y=x;
end
array208=array104;
for i=1:y,
    array208(ind,i)=array104(ro,i);    
end
   

% measurements rotate with stimulation, we untwist them
function vv= untwist(dd);
   elec=16;
   pos_i= [0,1];
   ELS= rem(rem(0:elec^2-1,elec) - ...
        floor((0:elec^2-1)/elec)+elec,elec)';
   ELS=~any(rem( elec+[-1 0 [-1 0]+pos_i*[-1;1] ] ,elec)' ...
        *ones(1,elec^2)==ones(4,1)*ELS')';
   twist= [               0+(1:13), ...
                         13+(1:13), ...
           39-(0:-1:0),  26+(1:12), ...
           52-(1:-1:0),  39+(1:11), ...
           65-(2:-1:0),  52+(1:10), ...
           78-(3:-1:0),  65+(1: 9), ...
           91-(4:-1:0),  78+(1: 8), ...
          104-(5:-1:0),  91+(1: 7), ...
          117-(6:-1:0), 104+(1: 6), ...
          130-(7:-1:0), 117+(1: 5), ...
          143-(8:-1:0), 130+(1: 4), ...
          156-(9:-1:0), 143+(1: 3), ...
          169-(10:-1:0),156+(1: 2), ...
          182-(11:-1:0),169+(1: 1), ...
          195-(12:-1:0), ...
          208-(12:-1:0) ];
    vv= zeros(256,size(dd,2));
    vv(ELS,:)= dd(twist,:);
   %vv= dd(1:208,:);

% Read data from p2k files (I T S system)
% FIXME: this code is very rough, it works for
%   only eight ring data records
function  vv = its_readdata( fname ) 
   fid= fopen( fname, 'rb', 'ieee-le');
   vv=[];

   % don't know how to interpret header
   header= fread(fid, 880, 'uchar');
   frameno= 0;
   rings= 8;
   while( ~feof(fid) )
       frameno= frameno+1;
       % don't know how to interpret frame header
       framehdr= fread(fid, 40);
       data= fread(fid, 104*rings, 'double');
       vv= [vv, data];
   end

   if 0 % convert a ring to 208
      ringno= 1;
      ld= size(vv,2);
      vx= [zeros(1,ld);vv( ringno*104 + (-103:0) ,: )];
      idx= ptr104_208;
      vv= vx(idx+1,:);
   end

% pointer to convert from 104 to 208 meas patterns
function idx= ptr104_208;
    idx= zeros(16);
    idx(1,:)= [0,0,1:13,0];
    ofs= 13;
    for i= 2:14
        mm= 15-i;
        idx(i,:) = [zeros(1,i+1), (1:mm)+ofs ];
        ofs= ofs+ mm;
    end

    idx= idx + idx';
    
function vv = iirc_readdata( fname );
    fid= fopen( fname, 'r');
    while ~feof(fid)
       line = fgetl(fid);
       if isempty(line)
           continue;
       end
       
       num= regexp(line,'Channel : (\d+)');
       if ~isempty(num)
           channels= str2num( line(num(1):end ) );
           continue;
       end
       
       num= regexp(line,'Frequency : (\d+)kHz');
       if ~isempty(num)
           freqency= str2num( line(num(1):end ) );
           continue;
       end

       num= regexp(line,'Scan Method : (\w+)');
       if ~isempty(num)
           scan_method=  line(num(1):end );
           continue;
       end

       num= regexp(line,'Study : (\w+)');
       if ~isempty(num)
           study=  line(num(1):end);
           continue;
       end
           
       if strcmp(line,'Data');
           data= fscanf(fid,'%f',[4,inf])';
           continue;
       end
    end
    vv= data(:,1) + 1i*data(:,2);
    if length(vv) ~= channels^2
        error('eidors_readdata: data length wrong')
    end

% stimulation is the fwd_model stimulation data structure
% meas_select indicates if data is NOT measures on current electrode
function [stimulations,meas_select] = UCT_sequence_file( fname );
   % (c) Tim Long
   % 21 January 2005
   % University of Cape Town


   % open the file
   fid = fopen(fname, 'rt');

   % check to see if file opened ok
   if fid == -1
         errordlg('File not found','ERROR')  
         return;
   end


   tline = fgetl(fid);             % get the spacer at top of text file

   % the measurement and injection pattern is stored as follows:
   % I1V1:  db #$00,#$0F,#$00,#$00,#$10,#$00,#$00,#$21,#$00 ...
   % I2V2:  db #$11,#$0F,#$11,#$11,#$10,#$11,#$11,#$21,#$11 ...
   % etc
   % need to put all the bytes in a vector

   % tokenlist will store the list of bytes as strings
   tokenlist = [];


   tline = fgetl(fid);             % get first line of data

   while length(tline) ~= 0
       
       % the first few characters in the line are junk
       rem = tline(11:end);            % extract only useful data
       
       % extract each byte
       while length(rem) ~=0
           [token, rem] = strtok(rem, ',');
           tokenlist = [tokenlist; token];
       end
       
       % get the next line in sequence file
       tline = fgetl(fid);
   end

   fclose(fid);

   % got everything in string form... need to covert to number format

   drive_lay = [];
   drive_elec = [];
   sense_lay = [];

   injection_no = 1;
   % for each injection
   for i=1:3:length(tokenlist)
       
       % get injection layer
       tsource_layer = tokenlist(i,3);
       tsink_layer = tokenlist(i,4);
       source_layer = sscanf(tsource_layer, '%x');
       sink_layer = sscanf(tsink_layer, '%x');
       
       drive_lay = [drive_lay; [source_layer sink_layer]];
         
       
       % get drive pair
       tsource_elec = tokenlist(i+1,3);
       tsink_elec = tokenlist(i+1,4);
       source_elec = sscanf(tsource_elec, '%x');
       sink_elec = sscanf(tsink_elec, '%x');
       
       drive_elec = [drive_elec; [source_elec sink_elec]];
       
       
       % get sense layer pair
       tpos_layer = tokenlist(i+2,3);
       tneg_layer = tokenlist(i+2,4);
       pos_sense_layer = sscanf(tpos_layer, '%x');
       neg_sense_layer = sscanf(tneg_layer, '%x');
       
       sense_lay = [sense_lay; [pos_sense_layer neg_sense_layer]];
   end

   n_elec = size(sense_lay,1);
   n_inj  = size(drive_lay,1);       % every injection
   elecs_per_plane = 16; % FIXED FOR THE UCT DEVICE
   raw_index = 0;
   meas_select = [];

   for i=1:n_inj      % for every injection
       stimulations(i).stimulation= 'mA';
       
       % find the injection electrodes
       e_inj_p = drive_lay(i, 1) * elecs_per_plane + drive_elec(i,1) + 1;
       e_inj_n = drive_lay(i, 2) * elecs_per_plane + drive_elec(i,2) + 1;

       % create the stimulation pattern for this injection
       inj = zeros(n_elec, 1);
       inj(e_inj_p) = 1;
       inj(e_inj_n) = -1;
       stimulations(i).stim_pattern = inj;
     
       % the UCT instrument always makes 16 measurements per injection.
       % the +ve and -ve electrodes are always adjacent, but might be on
       % different planes
       meas_pat = [];
       for e = 0:15
           raw_index = raw_index + 1;
           meas = zeros(1, n_elec);   % the measurement electrodes for this sample

           % find the measurement electrodes for this measurement (+ve elec is
           % next to -ve electrode)
           e_meas_p = sense_lay(i,1) * elecs_per_plane + mod(e+1,elecs_per_plane) + 1;
           e_meas_n = sense_lay(i,2) * elecs_per_plane + e + 1;

           % if either of the drive electrodes are equal to any of the sense
           % electrodes, we must not include this sample
           if any( e_meas_p == [e_inj_p, e_inj_n] ) | ...
              any( e_meas_n == [e_inj_p, e_inj_n] ) 
               continue;
           end

           meas(e_meas_p) = -1;
           meas(e_meas_n) = 1;
           meas_select = [meas_select;raw_index];
           
           % add this measurement to the measurement pattern
           meas_pat = [meas_pat; meas];

       end     % for each injection there are actually only 13-16 measurements
       stimulations(i).meas_pattern = sparse(meas_pat);
   end

function [cur_data,no_cur_data] = UCT_calibration_file( fname );
   fid = fopen(fname, 'rb');
   mag_num = fread(fid, 1, 'int32');
   version = fread(fid, 1, 'int32');
%  comments = fread(fid, 2048, 'char'); MATLAB CHANGED CHAR to 2 bytes
   comments = fread(fid, 2048, 'uint8');
   comments = setstr(comments');
   no_of_layers = fread(fid, 1, 'int32');
   uppa_dac = fread(fid, 1, 'float64');
   lowa_dac = fread(fid, 1, 'float64');
   raw_frame_size = fread(fid, 1, 'int32');

   no_cur_data = [];
   cur_data = [];

   for i=1:no_of_layers
       no_cur_data = [no_cur_data;fread(fid, raw_frame_size, 'float64')];
       cur_data = [cur_data;fread(fid, raw_frame_size, 'float64')];
   end
   fclose(fid);

%  no_cur_data = UCT_ShuffleData( no_cur_data);
%  cur_data =    UCT_ShuffleData( cur_data);

function [v_raw] = UCT_LoadDataFrame(infilename)
% [v_raw] = LoadDataFrame(infilename)
%
% Loads the data from the tomography file
%
% (c) Tim Long
% 21 January 2005
% University of Cape Town


% Open the file
fid = fopen(infilename, 'rb');

if fid == -1
      errordlg('File not found','ERROR')  
      return;
end

%%% new file changes
magic_number = fread(fid, 1, 'uint32');

if magic_number ~= 2290649224
   disp('UCT File: wrong file type'); 
end
version = fread(fid, 1, 'uint32');
foffset = fread(fid, 1, 'uint32');
no_of_layers = fread(fid, 1, 'uint32');
frame_size = fread(fid, 1, 'uint32');
fseek(fid, foffset-8, 'cof');

%%% end of new file changes
%%% old file stuff
% no_of_layers = fread(fid, 1, 'uint32');
% frame_size = fread(fid, 1, 'uint32');

frame_no = fread(fid, 1, 'uint32');

v_raw = [];
while feof(fid)==0
    v_raw = [v_raw, fread(fid, frame_size*no_of_layers, 'float64')];
    frame_no = fread(fid, 1, 'uint32');
end

fclose(fid);
% UCT_ShuffleData??

% Read data from the file format develped by Pascal
% Gaggero at CSEM Landquart in Switzerland.
function [vv] = landquart1_readdata( fname );
   FrameSize = 1024;

   enableCounter=1;
   counter=0;

   fid=fopen(fname);

   codec=fread(fid,1,'int'); %read codec
   if codec~=1; error('Codec unexpected value'); end

   nbFileInHeader=fread(fid,1,'int'); %read codec
   
   for i=1:nbFileInHeader
       lenghtFile = fread(fid,1,'int64'); 
       jnk= fread(fid,lenghtFile,'int8');% discard the file 
   end
   
   
   vv= []; 
   while 1; 
       type=fread(fid,1,'int'); %check if not End of file
       if type == -1; break ; end
       
       nel= fread(fid,1,'int');
       sz= fread(fid,1,'int');
       iq=fread(fid,[2,sz/8],'int')';

       vv = [vv, iq*[1;1j]]; %I+ 1j*Q vectors
       
       for j=1:nel-1
           sz= fread(fid,1,'int');
           jnk = fread(fid,sz/4,'int');
       end
       
   end
   
% Read data from the file format develped by Swisstom, Landquart, Switzerland.
function [vv, elecImps, tStampsAbs, tStampsRel] = landquart2_readdata( fname )
   [fid msg]= fopen(fname,'r','ieee-be','UTF-8');
   try
      format_version = fread(fid,1,'int32','ieee-be');
      if format_version ~= 3
         error('unsupported file format version');
      else
         % get expected number of frames and preallocate variables
         fseek(fid,16,'cof');
         nFrames = fread(fid, 1, 'int32', 'ieee-be');
         tStampsAbs = nan(1, nFrames);
         tStampsRel = nan(1, nFrames);
         viPayload = nan(64, nFrames);
         iqPayload = nan(2048, nFrames);
         
         header_size = 2264; 
         fseek(fid,header_size + 8,'bof');
         
         %%% get frame size and length of payload at end of frame
         frame_length = fread(fid, 1, 'int32', 'ieee-be') + 12;
         %%% move back to start of 1. frame
         fseek(fid,header_size,'bof');
         
         %%% Read frames
         i = 1;
         while fseek(fid, 1,'cof') ~= -1
            fseek(fid, -1,'cof');
            % read absolute timestamp (milliseconds resolution)
            tStampsAbs(i) = fread(fid,1,'int64','ieee-le');
            pl = fread(fid,1,'int32','ieee-le');    % payload length (i.e. frame length)
            frame_header = fread(fid,12,'int32','ieee-le'); 
            % read relative timestamp (microseconds resolution)
            tStampsRel(i) = frame_header(5);    
            
            % drop some unintersting parts  
            fseek(fid,12, 'cof');
            
            % get electrode impedance measurement
            viPayload(:,i) = fread(fid,64,'int32','ieee-le');
            % get effective payload (voltage measurements)
            iqPayload(:,i) = fread(fid,2048,'int32','ieee-le');
            fseek(fid,header_size + i*frame_length,'bof');
            i = i +1;
         end

      end
   catch err
      fclose(fid);
      rethrow(err);
   end
   fclose(fid);
   
   i = i-1;
   if i ~= nFrames       
      % remove data which were preallocated but not read
      if i < nFrames       
         tStampsAbs(i+1:end) = [];
         tStampsRel(i+1:end) = [];
         viPayload(:,i+1:end) = [];
         iqPayload(:,i+1:end) = [];
      end
      eidors_msg('"%s": expected %.0f frames but read %.0f',fname, nFrames, i ,3);
   end

   % this is just a simple guess
   amplitudeFactor = 2.048 / (2^20 * 360 * 1000);
   vv = amplitudeFactor * (iqPayload(1:2:end,:) + 1i*iqPayload(2:2:end,:));
   
   elecImps = viPayload(1:2:end,:) + 1i*viPayload(2:2:end,:);

% Read data from the file format develped by Swisstom, Landquart, Switzerland.
function [vv, evtlist, elecImps, tStampsAbs, tStampsRel] = landquart4_readdata( fname )
   evtlist = [];
   [fid msg]= fopen(fname,'r','ieee-le','UTF-8');
   try
      format_version = fread(fid,1,'int32','ieee-le');
      if format_version ~= 4
         error('unsupported file format version');
      else          
         header_size = fread(fid,1,'int32', 'ieee-le');   
         eit_frame_offset = 328; % 60 + 12 + 256 = 328
         iq_payload = 2048;
         vi_payload = 64;
         
         % get expected number of frames and preallocate variables
         fseek(fid,16,'cof');
         nFrames = fread(fid, 1, 'int32', 'ieee-le');
         tStampsAbs = nan(1, nFrames);
         tStampsRel = nan(1, nFrames);
         viPayload = nan(vi_payload, nFrames);
         iqPayload = nan(iq_payload, nFrames);

         fseek(fid,header_size,'bof');
         
         %%% Read frames
         i = 1;
         evti = 1;
         while fseek(fid, 1,'cof') ~= -1
            fseek(fid, -1,'cof');
            % drop frame header and some payload:
            % read absolute timestamp (milliseconds resolution)
            tStampsAbs(i) = fread(fid,1,'int64','ieee-le'); 
            ft = fread(fid,1,'int32','ieee-le'); 
            pl = fread(fid,1,'int32','ieee-le');
            if ft == 1
               % event
               evtlist(evti).timestamp = tStampsAbs(i) ;
               evtlist(evti).frame = i ;
               evti = evti + 1;
               if pl > 0
                  evtlist(evti).eventId = fread(fid, 1, 'int32', 'ieee-le');
                end
            elseif ft == 0
                frame_header = fread(fid,15,'int32','ieee-le'); 
                % read relative timestamp (microseconds resolution)
                tStampsRel(i) = frame_header(5);   
                
                % drop some unintersting parts  
                fseek(fid,12, 'cof');
               
                % get electrode impedance measurement
                viPayload(:,i) = fread(fid,vi_payload,'int32','ieee-le');
                % get effective payload (voltage measurements)
                iqPayload(:,i) = fread(fid,iq_payload,'int32','ieee-le');
                fseek(fid,pl-4*iq_payload-eit_frame_offset,'cof');
                i = i+1;
            elseif pl > 0
               fseek(fid,pl,'cof'); 
            else
               % nothing to do
            end
         end         
      end
   catch err
      fclose(fid);
      rethrow(err);
   end
   fclose(fid);
   
   i = i-1;
   if i ~= nFrames       
      % remove data which were preallocated but not read
      if i < nFrames       
         tStampsAbs(i+1:end) = [];
         tStampsRel(i+1:end) = [];
         viPayload(:,i+1:end) = [];
         iqPayload(:,i+1:end) = [];
      end
      eidors_msg('"%s": expected %.0f frames but read %.0f',fname, nFrames, i ,3);
   end

   % this is just a simple guess
   amplitudeFactor = 2.048 / (2^20 * 360 * 1000);
   vv = amplitudeFactor * (iqPayload(1:2:end,:) + 1i*iqPayload(2:2:end,:));
   
   elecImps = viPayload(1:2:end,:) + 1i*viPayload(2:2:end,:);

function [vv] = landquart4pre_readdata( fname )
   [fid msg]= fopen(fname,'r','ieee-le','UTF-8');
   try
      format_version = fread(fid,1,'int32','ieee-le');
      if format_version ~= 4
         error('unsupported file format version');
      else
         header_size = fread(fid,1,'int32', 'ieee-le'); 
         
         fseek(fid,header_size + 12,'bof');
         %%% get frame size and length of payload at end of frame
         frame_length = fread(fid, 1, 'int32', 'ieee-le') + 16;
         %%% move back to start of 1. frame
         fseek(fid,header_size,'bof');
         
         %%% Read frames
         i = 1;
         while fseek(fid, 1,'cof') ~= -1
            fseek(fid, -1,'cof');
            % drop frame header and some payload: 
            % (16 + 60) + 12 + 256 = 344
            fseek(fid,344, 'cof');
            iqPayload(:,i) = fread(fid,2048,'int32','ieee-le');
            fseek(fid,header_size + i*frame_length,'bof');
            i = i +1;
         end

      end
   catch err
      fclose(fid);
      rethrow(err);
   end
   fclose(fid);

   % this is just a simple guess
   amplitudeFactor = 2.048 / (2^20 * 360 * 1000);
   vv = amplitudeFactor * (iqPayload(1:2:end,:) + 1i*iqPayload(2:2:end,:));
   
%  Output: [encodepage] = eidors_readdata( path,'DIX_encode');
%   where path= '/path/to/Criptografa_New.dll' (provided with system)
function [vv] = dixtal_read_codepage( fname );
   %Fname must be the Criptografa_New dll
   fid = fopen(fname,'rb');
   b1234= fread(fid, [1,4], 'uint8');
   if b1234~= [77, 90, 144, 0];
      error('This does not appear to be the correct Dll');
   end
   fseek(fid, hex2dec('e00'), 'bof');
   encodepage1 = fread(fid, hex2dec('1000'), 'uint8');
   fclose(fid);
   encodepage1 = flipud(reshape(encodepage1,4,[]));
   vv = encodepage1(:);

%        format = 'DIXTAL'
%           - output: [vv] = eidors_readdata( fname, 'DIXTAL', encodepage );
function [vv] = dixtal_read_data( file_name, frame_range, encodepage1 );

   % Get encodepage from the file
   fid=fopen(file_name,'rb');
   n1 = fgetl(fid);
   n2 = fgetl(fid);
   n3 = fgetl(fid); 
   nframes = str2num( n3 );
   
   fseek(fid, hex2dec('800'), 'bof');
   encodepage2 = fread(fid, hex2dec('1000'), 'uint8');
   fclose(fid);

   encodepage = bitxor(encodepage1, encodepage2);

   % Read the file
   dplen = hex2dec('1090');
   start = hex2dec('2800');
   fid=fopen(file_name,'rb');
   if isempty(frame_range)
      frame_range = 0:nframes;
   else
      frame_range = frame_range - 1;
      frame_range(frame_range>nframes) = [];
   end	  
   vv= zeros(dplen/4, length(frame_range));
   k= 1;
   for i = frame_range
      status= fseek(fid, start + i*dplen, 'bof');
      if status~=0; break; end
      datapage = fread(fid, dplen, 'uint8');
      if length(datapage)== 0; 
         vv= vv(:,1:k-1); break; % frame number inside file is wrong
      end
      vv(:,k) = proc_dixtal_data( datapage, encodepage );
	  k=k+1;
   end
   fclose(fid);


function  data = proc_dixtal_data( b, encodepage );

   cryptolen = 1024*4;
   b(1:cryptolen) = bitxor(b(1:cryptolen), encodepage);

   b1 = b(1:4:end,:); %b1 = bitand(b1,127);
   b2 = b(2:4:end,:);
   b3 = b(3:4:end,:);
   b4 = b(4:4:end,:);
   bb  = [b1,b2,b3,b4]*(256.^[3;2;1;0]);

   sgnbit = bitand( bb,  2^31);
   sgnbit = +1*(sgnbit == 0) - 1*(sgnbit == 2^31);

   expbit = bitand( bb, 255*2^23 ) /2^23; 
   expbit = 2.^(expbit - 127);

   fracbt = bitand( bb, 2^23-1)/2^23 + 1;
   data = sgnbit .* expbit .* fracbt;

function [fr] = read_draeger_header( filename );
%READ_DRAEGER_HEADER   Opens and reads the header portion of the Draeger .eit
%file. Current parameter returned is frame rate.
%
% function [fr,s] = ReadDraegerHeader(filename)
% fr        double      scalar      frame rate in hertz

% Determine file version
K0 = 2048;   % bytes to read in char class
K1='Framerate [Hz]:'; % Draeger frame rate line in header
K2=15;                % Length of K1 string
K3=13;                % Following F1 Field for delimiter (look for ^M)

% Open file for reading in little-endian form
fid = fopen(filename,'r','l');
if fid == -1
    error('Error read_draeger_header: can''t open file');
end
header = fread(fid,K0,'*char')';
index = strfind(header,K1);
if ~isempty(index)
    [tok,rem]= strtok(header(index+K2:end),K3);
    fr = str2num(tok); 
else
    error('Error read_draeger_header: frame rate unspecified');
end

if isempty(fr)
    error('Error read_draeger_header: Frame rate could not be read');
end

fclose(fid);


function vd = read_draeger_file(filename)
%READDRAEGERFILE   Opens and reads a Draeger .eit file. Returns an array
%containing the voltage data arranged per frame.
%   
% Input:
% filename  char        1xN         
% 
% Output:   
% vd        double      MxN         M = volt measurements per frame (mV)
%                                   N = number of frames        


   ff = dir(filename);
   filelen = ff.bytes;

   % Open file for reading in little-endian form
   fid = fopen(filename,'r','l');
   if fid == -1
       error('Error read_draeger_file: file could not be opened');
   end

   % Determine file version
   K0 = 128;   % bytes to read in char class

   header = fread(fid,K0,'*char')';
if 0 % THis analysis doesn't seem useful
   if ~isempty(strfind(header,'V3.2'))
       version = '3.2';
   elseif ~isempty(strfind(header,'V4.01'))
       version = '4.01'; 
   elseif ~isempty(strfind(header,'V1.10')) % 2014!!
       version = '1.10'; 
   else
       error('Error read_draeger_file: unknown file version');
   end
end

   % Read Format version
   fseek(fid,0,-1); % beginning of file
   hdr = fread(fid, 8, 'uint8');
   offset1 =  (256.^(0:3))*hdr(5:8) + 16; % "good" data starts at offset #2
   switch sprintf('%02X-',header(1:4));
      case '1F-00-00-00-';
         type='Draeger v31';  SPC = 4112;
      case '20-00-00-00-';
         type='Draeger v32';  SPC = 5200;
      case '33-00-00-00-';
         type='Draeger v51';  SPC = 5495;
      otherwise;
         error('File "%s" is format version %d, which we don''t know how to read', ...
               hdr(1))
   end

   len = floor( (filelen-offset1)/SPC );
   vd= zeros(600,len);
   ss= offset1 + (0:len)*SPC;

   for k= 1:length(ss);
      if fseek(fid,ss(k),-1)<0; break; end
      vd(:,k)=fread(fid,600,'double', 0, 'ieee-le');
   end

   % End of function
   fclose(fid);

