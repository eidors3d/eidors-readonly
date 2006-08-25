function [vv,curr,volt]= eit_readdata( fname, format )
% EIDORS readdata - read data files from various EIT equipment
%    manufacturers
%
% Currently the list of supported file formats is:
%    1. MCEIT (Goettingen / Viasys) "get" file format 
%        format = "GET" or "MCEIT"
%    2. ITS (International Tomography Systems)
%        format = "ITS" or "p2k"
%    3. IIRC (Impedance Imaging Research Center, Korea)
%        format = "txt" or "IIRC"
%
% Usage
%  [vv,curr,volt]= eidors_readdata( fname, format )
%     vv = measurements 208xmeas
%     fname = file name
%
%  if format is unspecified, we attempt to autodetect

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: eidors_readdata.m,v 1.6 2006-08-25 07:02:50 tonginoh Exp $

% TODO:
%   - output an eidors data object
%   - test whether file format matches specified stimulation pattern

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
if     strcmp(fmt, 'get') | strcmp(fmt, 'mceit')  
   [vv,curr,volt] = mceit_readdata( fname );
elseif strcmp(fmt, 'p2k') | strcmp(fmt, 'its')
   vv = its_readdata( fname );
elseif strcmp(fmt,'txt') | strcmp(fmt, 'iirc')
   vv = iirc_readdata( fname );
else
   error('eidors_readdata: file "%s" format unknown', fmt);
end
   

function [vv,curr,volt] = mceit_readdata( fname );
   elec=16;
   pos_i= [0,1];
   ELS= rem(rem(0:elec^2-1,elec) - ...
        floor((0:elec^2-1)/elec)+elec,elec)';
   ELS=~any(rem( elec+[-1 0 [-1 0]+pos_i*[-1;1] ] ,elec)' ...
		     *ones(1,elec^2)==ones(4,1)*ELS')';

   fid= fopen(fname,'rb');
   d= fread(fid,inf,'float');
   fclose(fid);

   if rem( length(d), 256) ~=0
      error('File length strange');
   end

   dd= reshape( d, 256, length(d)/256);

% measurements rotate with stimulation, we untwist them
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

   curr=0.00512*dd(209:224,:);  % Amps
   volt=12*dd(225:240,:); %Vrms
   %input impedance=voltage./current-440;	Ohm

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
       
       num= regexp(line,'Channel : (\d+)','tokens');
       if ~isempty(num)
           channels= str2num( num{1}{1} );
           continue;
       end
       
       num= regexp(line,'Frequency : (\d+)kHz','tokens');
       if ~isempty(num)
           freqency= str2num( num{1}{1} );
           continue;
       end

       num= regexp(line,'Scan Method : (\w+)','tokens');
       if ~isempty(num)
           scan_method=  num{1}{1};
           continue;
       end

       num= regexp(line,'Study : (\w+)','tokens');
       if ~isempty(num)
           study=  num{1}{1};
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
