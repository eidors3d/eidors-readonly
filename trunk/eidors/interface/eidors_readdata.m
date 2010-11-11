function [vv, auxdata, stim ]= eidors_readdata( fname, format )
% EIDORS readdata - read data files from various EIT equipment
%    manufacturers
%
% Currently the list of supported file formats is:
%    - MCEIT (Goettingen / Viasys) "get" file format 
%        format = "GET" or "MCEIT"
%    - Draeger "get" file format (older format for Draeger equipment)
%        format = "GET" or "draeger-get"
%    - New Carefusion "EIT" file format
%        format = "EIT" or "carefusion"
%    - Sheffield MK I "RAW" file format
%        format = "RAW" or "sheffield"
%    - ITS (International Tomography Systems)
%        format = "ITS" or "p2k"
%    - IIRC (Impedance Imaging Research Center, Korea)
%        format = "txt" or "IIRC"
%    - University of Cape Town formats
%        format = "UCT_SEQ"  UCT sequence file
%           - Output: [ stimulation, meas_select]= eit_readdata(fname, 'UCT_SEQ')
%        format = "UCT_CAL"  UCT calibration file
%           - Output: [vv, no_cur_caldata_raw ]= eit_readdata( fname, 'UCT_CAL' )
%                 where no_cur_caldata_raw is data captured with no current
%        format = "UCT_DATA"  UCT data frame file
%           - Output: [vv]= eit_readdata( fname, 'UCT_DATA' )
%
% Usage
% [vv, auxdata, stim ]= eit_readdata( fname, format )
%     vv      = measurements - data frames in each column
%     auxdata = auxillary data - if provided by system 
%     stim    = stimulation structure, to be used with
%                fwd_model.stimulation. 
%     fname = file name
%
%  if format is unspecified, we attempt to autodetect

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


switch pre_proc_spec_fmt( format, fname );
   case 'mceit';
      [vv,curr,volt,auxdata_out] = mceit_readdata( fname );
      auxdata.auxdata = auxdata_out;
      auxdata.curr    = curr;
      auxdata.volt    = volt;

      stim = basic_stim(N_el);
   case 'draeger-get'
      vv = draeger_get_readdata( fname );

      stim = basic_stim(N_el);
   case {'raw', 'sheffield'}
      vv = sheffield_readdata( fname );

      stim = basic_stim(N_el);
   case {'p2k', 'its'}
      vv = its_readdata( fname );

      stim = 'UNKNOWN';
   case {'txt','iirc'}
      vv = iirc_readdata( fname );

      stim = basic_stim(N_el);
   case 'uct_seq'
      [vv,auxdata] = UCT_sequence_file( fname );

      stim = 'UNKNOWN';
   case 'uct_cal'
      [vv,auxdata] = UCT_calibration_file( fname );

      stim = 'UNKNOWN';
   case 'uct_data'
      [vv] = UCT_LoadDataFrame( fname );

      stim = 'UNKNOWN';
   case 'carefusion'
      [vv] = carefusion_eit_readdata( fname );
  
      stim = basic_stim(N_el);
   otherwise
      error('eidors_readdata: file "%s" format unknown', fmt);
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

   if strcmp(fmt,'eit')
      if is_eit_file_a_carefusion_file( fname )
         fmt= 'carefusion';
      else
         error('EIT file specified, but it doesn''t seem to be a Carefusion file')
      end
   end

function df= is_get_file_a_draeger_file( fname)
   fid= fopen(fname,'rb');
   d= fread(fid,[1 26],'uchar');
   fclose(fid);
   df = all(d == '---Draeger EIT-Software---');

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

function [vv] = carefusion_eit_readdata( fname );
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
   fseek(fid, header, -1); % From origin
   d= fread(fid,[inf],'float32');
   fclose(fid);

%  Data records are 262 long
   endpt = floor(length(d)/262)*262;
   vv= reshape(d(1:endpt),262,[]);
   vv= vv(36+[1:208],:);
   

% Read the old <2006 Draeger "get" file format
function [vv,curr,volt,auxdata] = draeger_get_readdata( fname );
   fid= fopen(fname,'rb');
   emptyctr=0;
   while emptyctr<2
     line= fgetl(fid);
     if isempty(line)
        emptyctr= emptyctr+1;
     else
        eidors_msg(line,0);
        emptyctr=0;
     end
   end
   d= fread(fid,inf,'float',0,'ieee-le');
   fclose(fid);

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

   

function [vv,curr,volt,auxdata] = mceit_readdata( fname );

   fid= fopen(fname,'rb');
   d= fread(fid,inf,'float');
   fclose(fid);

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
   %input impedance=voltage./current-440;        Ohm

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
