function [vv, auxdata ]= eidors_readdata( fname, format )
% EIDORS readdata - read data files from various EIT equipment
%    manufacturers
%
% Currently the list of supported file formats is:
%    - MCEIT (Goettingen / Viasys) "get" file format 
%        format = "GET" or "MCEIT"
%    - Draeger "get" file format
%        format = "GET" or "draeger"
%    - Sheffield MK I "RAW" file format
%        format = "RAW" or "sheffield"
%    - ITS (International Tomography Systems)
%        format = "ITS" or "p2k"
%    - IIRC (Impedance Imaging Research Center, Korea)
%        format = "txt" or "IIRC"
%
% Usage
% [vv, auxdata ]= eit_readdata( fname, format )
%     vv      = measurements - data frames in each column
%     auxdata = auxillary data - if provided by system 
%     fname = file name
%
%  if format is unspecified, we attempt to autodetect

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
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

fmt= lower(format);
if strcmp(fmt,'get')
    if is_get_file_a_draeger_file( fname)
       fmt= 'draeger';
    else
       fmt= 'mceit';
    end
end

if     strcmp(fmt, 'mceit')  
   [vv,curr,volt,auxdata] = mceit_readdata( fname );
elseif strcmp(fmt, 'draeger')
   vv = draeger_readdata( fname );
elseif strcmp(fmt, 'raw') | strcmp(fmt, 'sheffield')
   vv = sheffield_readdata( fname );
elseif strcmp(fmt, 'p2k') | strcmp(fmt, 'its')
   vv = its_readdata( fname );
elseif strcmp(fmt,'txt') | strcmp(fmt, 'iirc')
   vv = iirc_readdata( fname );
else
   error('eidors_readdata: file "%s" format unknown', fmt);
end
   
function df= is_get_file_a_draeger_file( fname)
   fid= fopen(fname,'rb');
   d= fread(fid,[1 26],'char');
   fclose(fid);
   df = all(d == '---Draeger EIT-Software---');

function [vv,curr,volt,auxdata] = draeger_readdata( fname );
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
