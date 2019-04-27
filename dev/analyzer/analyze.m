function analyze
%ANALYZE: analyze EIT files for lung parameters

% (C) Andy Adler & Symon Stowe 2018-2019. License: GPL v2 or v3.
% $Id$

% TODO:
%  unify filtering -> Maybe into the "FILE" section

  parse_config;
  write_header;
  write_table_hdr;
  iterate_over_files;
  write_footer;

function iterate_over_files
  global pp; % parameters
  files= dir('*.mat');
  for f= 1:length(files)
     fn = files(f).name;
     disp(fn);
     dd = loadfile(fn);
     fid = outfile;
     fprintf(fid,'<TR><TD>%s',fn(1:end-4));
     for i=1:length(pp.callfns)
        clf;
        reqbreaths = feval(pp.callfns{i},'REQBREATHS?');
    	reqbeats   = feval(pp.callfns{i},'REQBEATS?');
        if reqbreaths && dd.n_breaths==0
           out = '<font size="+2"><center><b>No breaths detected</b></center></font>';
        elseif reqbeats && dd.n_beats==0
    	   out = '<font size="+2"><center><b>No beats detected</b></center></font>';
        else
           out = feval(pp.callfns{i},dd,f);
        end
        fprintf(fid,'<TD>%s',out);
     end
        fprintf(fid,'<TD>%s',fn(1:end-4));
     fclose(fid);
  end

function parse_config 
  global pp; pp=struct();
  pp.callfns = cell(0,0); 
  pp.rotate  = 0;
  pp.slices  = 4;
  pp.min_insp_length  = 0.5; % 1 seconds for horses
  pp.min_insp_length  = 1.0; % 1 seconds for horses
  pp.FRC_search_window = 0.2; % search 100ms for FRC
  pp.FRC_relative_match = 0.2; % match start/end FRC
  pp.min_heart_peak_separation = 0.2; % (seconds) horse min heart beat separation

  pp.flow_window = 10:50;
  pp.colourbar = 'colourbar.png';
  pp.color_map = 'ocean';

  fid = fopen('config.m');
  count = 0;
  while true;
     tline = fgetl(fid);
     if ~ischar(tline); break; end
     count = count+1;
     % Remove leading (and trailing) whitespace
     tline = strtrim(tline); 
     do_locs = regexpi(tline,'DO ');
     config_locs = regexpi(tline,'CONFIG ');
     comment_locs = regexpi(tline,'%');
     if do_locs == 1
        tline = tline(4:end);
        pp.callfns{end+1,1} = tline;
     elseif config_locs == 1 
        eval(tline);
     elseif (comment_locs ~= 1) && (size(tline) > 0)
        warning(['Unsupported input in file: ',fid,' on line: ',num2str(count),'.']) 
     end
  end
  fclose(fid);

  subplot(611);
  mycolormap;
  image(linspace(0,256,32));
  set(gca,'YTick',[]);
  set(gca,'XTick',linspace(0.5,32.5,3));
  set(gca,'XTickLabel',{'0','50%','Max'});
  print_convert(pp.colourbar,struct('resolution',30));

function DO(varargin);
  global pp;
  fname = varargin{1};
  if nargin>1; error('DO expects one argument'); end
  try outstr = feval(fname,'TITLE');
     pp.callfns{end+1} = fname;
  catch
     error('Function %s not available',fname);
  end

function CONFIG(varargin)
  global pp;
  cmd = varargin{1};
  v2str = varargin{2};
  v2num = str2num(v2str);
  switch cmd
    case 'rotate';             pp.rotate= v2num;
    case 'min_insp_length';    pp.min_insp_length=   v2num;
    case 'min_expi_length';    pp.min_expi_length=   v2num;
    case 'FRC_search_window';  pp.FRC_search_window= v2num;
    case 'FRC_relative_match'; pp.FRC_relative_match= v2num;
    case 'slices'           ;  pp.slices= v2num;
    case 'color_map'        ;  pp.color_map = v2str;
    case 'LP_filter'        ;  pp.LP_filter = v2num; 
    case 'model_breaths_ncos'; pp.model_breaths_ncos = v2num;
    case 'FILE';
       fh = hsh(v2str);
       vstr = varargin{4};
       vnum = str2num(vstr);
       switch varargin{3}
          case 'endtime'; pp.( fh ).endtime = vnum;
          case 'LPF_fc';  pp.( fh ).LPF_fc = vnum;
          otherwise error(['CONFIG FILE parameter ', ...
              ' %s not understood'], varargin{3});
       end 
    otherwise;
      error('CONFIG parameter %s not understood', cmd);
  end

% Make assoc arrays work
function hh= hsh( str )
   if exist('OCTAVE_VERSION')==5
      hh= ['H',hash('sha1',str)];
   else
      hh= ['H',rptgen.hash(str)];
   end
   
% Frequency axis for filtering
function fax = freq_axis( dd );
  fax = linspace(0,dd.FR,dd.lD+1);
  fax(end)=[];
  fax(fax>dd.FR/2) = fax(fax>dd.FR/2) - dd.FR;

% Filter in the frequency direction
function s = freq_filt(s,fresp,dim);
  f = fft(s,[],dim);
  if size(s,dim)~= length(fresp);
     error('Incompatible sizes');
  end
  fshape = [1,1,1]; fshape(dim) = size(s,dim);
  f = f .* reshape(fresp,fshape);
  s= ifft(f,[],dim);
  if norm(imag(s(:))) > 1e-13
     error('FFT filter has imag output');
  end
  s = real(s);
  

function dd = filtfile(dd,fname);
  global pp;
  fh = hsh(fname); 
  if ~isfield(pp,fh); return; end % no specific info
  cfg = pp.( fh );
  % make sure these are done in order
  %   so that time cuts occur first
  for fn = {'endtime','starttime','LPF_fc'}
     if ~isfield(cfg,fn{1}); continue; end
     fval = cfg.( fn{1} );
     switch fn{1}
        case 'endtime'; 
           lim = round(fval * dd.FR);
           dd.CV(lim:end) = [];
           dd.tt(lim:end) = [];
           dd.ZR(:,:,lim:end) = [];
           dd.lD = length(dd.CV);
        case 'starttime'; 
           lim = round(fval * dd.FR);
           dd.CV(1:lim) = [];
           dd.tt(1:lim) = [];
           dd.ZR(:,:,1:lim) = [];
           dd.lD = length(dd.CV);
        case 'LPF_fc'; 
% TODO Add windowing and other controls
           fax = freq_axis( dd );
           fresp = abs(fax)<fval;
           dd.ZR = freq_filt(dd.ZR,fresp,3);
           dd.CV = freq_filt(dd.CV,fresp,2);

        otherwise; 
           error 'CONFIG FILE option not understood';
     end
  end


function dd = loadfile(fname);
  global pp;
  in = load(fname);
  dd.ZR = in.data.measurement.ZeroRef;
  dd.CV = in.data.measurement.CompositValue(:)';
  dd.FR = in.data.imageRate;
  dd.lD = length(dd.CV);
  dd = filtfile(dd,fname); % filter if required

  dd.tt = (0:dd.lD-1)/dd.FR;
  dd.breaths = find_frc(dd);
  dd.n_breaths = size(dd.breaths,1);
  dd.beats = find_beats(dd); 
  dd.n_beats = size(dd.beats,1);

  ls = linspace(0,1,10);
  ls = [ls,-fliplr(ls)];
  dd.flow= -conv2(dd.CV,ls,'same');
  ls = reshape(ls,1,1,[]);
  dd.ZF = -convn(dd.ZR,ls,'same');

function out= show_perfusion(dd,ii)
  global pp;
  if ischar(dd) && strcmp(dd,'TITLE');
     out = 'Perfusion Image'; return
  end
  if ischar(dd) && strcmp(dd,'REQBREATHS?');
     out = false; return
  end
  if ischar(dd) && strcmp(dd,'REQBEATS?');
     out = true; return
  end
  clf;
  perfusion = perfusion_calc(dd);
  perfusion(dd.ZR(:,:,1)==0) = NaN;
  mycolormap;
  my_image(perfusion*150/max(perfusion(:))+100);
  fout = sprintf('perfusion_image%03d.png',ii);
  out = sprintf(['<center>max pixel=%1.3f<br>' ...
     '<a href="%s"><img width="200" src="%s">' ...
     '</a><p><img src="%s"></center>'], ...
     max(perfusion(:)), fout, fout, pp.colourbar);
  print_convert(fout);
 
function perfusion = perfusion_calc(dd);
  perfusion = [];
  for i=1:dd.n_beats
    eie = dd.beats(i,[1,2,3]);
    perfusion(:,:,i) = dd.ZR(:,:,eie(2)) - mean(dd.ZR(:,:,eie([1,3])),3);
  end
  perfusion= mean(perfusion,3);

function out= show_apnoea(dd,ii) 
  if ischar(dd) && strcmp(dd,'TITLE');
     out = 'Apnoea Segment'; return
  end
  if ischar(dd) && strcmp(dd,'REQBREATHS?');
     out = false; return
  end
  if ischar(dd) && strcmp(dd,'REQBEATS?');
     out = false; return
  end
  fout = sprintf('apnoea_segment%03d.png',ii);
  out = sprintf( ...
     '<a href="%s"><img width="300" src="%s"></a>',...
     fout, fout);
  clf; subplot(211);
  seg = dd.CV;
  plot(dd.tt,seg,'LineWidth',2); box off; axis tight
  title('All pixels')
  xlabel('seconds')
  subplot(212);
  seg = max_pixel(dd);
  plot(dd.tt,seg,'LineWidth',2); box off; axis tight
  title('Maximum pixels')
  xlabel('seconds')
  print_convert(fout);

function out= show_max_pixel(dd,ii) 
  if ischar(dd) && strcmp(dd,'TITLE');
     out = 'Maximum Pixels'; return
  end
  if ischar(dd) && strcmp(dd,'REQBREATHS?');
     out = false; return
  end
  if ischar(dd) && strcmp(dd,'REQBEATS?');
     out = false; return
  end
  fout = sprintf('max_pixel_segment%03d.png',ii);
  out = sprintf( ...
     '<a href="%s"><img width="300" src="%s"></a>',...
     fout, fout);
  clf; subplot(211);
  pix_wave = max_pixel(dd);
  plot(dd.tt,pix_wave,'LineWidth',2); box off;
  axis tight
  print_convert(fout);

% TODO: Let N be the N maximum pixels (separated watershed-like)
function pix_wave = max_pixel(dd, N)
  if nargin == 1; N=1; end
  stdi = std(dd.ZR,[],3);
  [~,idx] = max(stdi(:));
  ZRr = reshape(dd.ZR,[], size(dd.ZR,3));
  pix_wave = ZRr(idx,:);
  

function out= show_beats(dd,ii)
  global pp;
  if ischar(dd) && strcmp(dd,'TITLE');
     out = 'Beat Detection'; return
  end
  if ischar(dd) && strcmp(dd,'REQBREATHS?');
     out = false; return
  end
  if ischar(dd) && strcmp(dd,'REQBEATS?');
     out = true; return
  end
  % Calculate the heart rate
  last_first = (dd.beats(end,2)-dd.beats(1,2))/dd.FR
  HR = (dd.n_beats-1)/last_first*60; 
  fout = sprintf('beat_detection%03d.png',ii);
  out = sprintf([...
       '<center>Average heart rate=%1.0f bpm<br>' ...
       '<a href="%s"><img width="300" src="%s">' ...
       '</a></center>'], ...
       HR, fout, fout);
  clf; subplot(211); % All pixels
  plot(dd.tt,dd.CV,'LineWidth',2); box off;
  H = (max(dd.CV) - min(dd.CV))/10;
  for i=1:dd.n_beats
     eie = dd.beats(i,[1,3,2,1]);
     line(dd.tt(eie), dd.CV(eie), 'color',[0,0,0],'LineWidth',2);
  end
  axis tight
  title('All pixels')
  xlabel('seconds')
  subplot(212); % Brightest pixels
  pix_wave = max_pixel(dd);
  plot(dd.tt,pix_wave,'LineWidth',2); box off;
  H = (max(pix_wave) - min(pix_wave))/10;
  for i=1:dd.n_beats
     eie = dd.beats(i,[1,3,2,1]);
     line(dd.tt(eie), pix_wave(eie), 'color',[0,0,0],'LineWidth',2);
  end
  axis tight
  title('Maximum pixels')
  xlabel('seconds')
  print_convert(fout);

function beats = find_beats(data)
% Identify the heart beats in the Apnoea signal
% May work on ventilation segments
  global pp;
  % seq = data.CV;
  seq = max_pixel(data);
  % Take an FFT of the data...

  lseq = length(seq);
  Fseq= fft(seq);
  % Cut off freqs
  % TODO - check this and potentially remove for non-frequency domain
  fc_low = 0.35; % Low frequency cutoff (Hz)
  fc_high= 0.65; % High frequency cutoff (Hz)
  fc_low = round(fc_low*lseq/data.FR);
  fc_high = round(fc_high*lseq/data.FR);
  fc_center = floor((fc_high-fc_low)/2+fc_low);
  dwn = fc_center - fc_low+5;
  up  = fc_high - fc_center+5;
  L = up+dwn;% window Length
  mask = zeros(size(Fseq));
  window = blackman(L);
  % Added in some harmonics to keep more of the shape
  adjust = 1;
  shift = 3;
  mask([shift+(fc_center)-dwn:shift+(fc_center)+up-adjust]) = window;
  mask([end-fc_center-up:end-fc_center+dwn-adjust]) = window; % Other end
  % Harmonic 1 - 0.5 blackman
  mask([shift+2*(fc_center)-dwn:shift+2*(fc_center)+up-adjust]) = 0.5*window;
  mask([end-2*fc_center-up:end-2*fc_center+dwn-adjust]) = 0.5*window; % Other end
  % Harmonic 2 - 0.25 blackman
  mask([shift+3*fc_center-dwn:shift+3*fc_center+up-adjust]) = 0.25*window;
  mask([end-3*fc_center-up:end-3*fc_center+dwn-adjust]) = 0.25*window; % Other end
  % harmonic 3 - 0.15 blackman
  mask([shift+4*fc_center-dwn:shift+4*fc_center+up-adjust]) = 0.15*window;
  mask([end-4*fc_center-up:end-4*fc_center+dwn-adjust]) = 0.15*window; % Other end
  % Harmonic 4 - 0.1 blackman
  mask([shift+5*fc_center-dwn:shift+5*fc_center+up-adjust]) = 0.1*window;
  mask([end-5*fc_center-up:end-5*fc_center+dwn-adjust]) = 0.1*window; % Other end
  Fseq = Fseq.*mask;
  seq1= ifft(Fseq);
  if std(imag(seq1))>1e-10; error('Heart Frequency FFT code'); end
  seq1= real(seq1);
  % Flow calc
  flow = diff(seq1);% first differences
  thresh = 0.5*median( abs(flow)); % /2 to allow for the detection of all of the heart beats...
  inout= zeros(lseq,1);
  i = 1;
  inout(i) = sign( flow(i)+eps ); % eps to force not zero
  for i = 2:lseq-1
    cval= inout(i-1);
    if cval*flow(i) > -thresh
       inout(i) = cval;
    else
       inout(i) = -cval;
    end
  end
  inout(lseq) = inout(lseq-1);

% Shorter forward window because detection is lagging
  beat_window_fwd= round(pp.min_heart_peak_separation/5*data.FR);
  beat_window_bak= round(pp.min_heart_peak_separation*data.FR);

  dinout= diff( inout );
  fdiff = find( diff(inout) ); % This has all peaks and troughs from filtered signal
  fdiff(fdiff<=beat_window_bak      )= []; % too close
  fdiff(fdiff>lseq - beat_window_fwd)= []; % too close

  min_beat = fdiff( (dinout(fdiff)>0) ); % Minimim point of abs HR signal
  max_beat = fdiff( (dinout(fdiff)<0) ); % Max point of abs heart signal

  % Basic beat detection done here...

  % Find the best point
  ww= -beat_window_bak:beat_window_fwd;
  for i=1:length(min_beat);
    wind= seq( min_beat(i)+ ww );
    ff = find( wind== min(wind) );
    ff= ff(1)+min(ww)-1;
    min_beat(i)= min_beat(i) + ff;
  end
  for i=1:length(max_beat);
    wind= seq( max_beat(i)+ ww );
    ff = find( wind== max(wind) );
    ff= ff(1)+min(ww)-1;
    max_beat(i)= max_beat(i) + ff;
  end
  beat_sep = pp.min_heart_peak_separation;
  beats = [];
  % cycle through all of the beats
  i=1;e=1; while true;
     if i>=length(max_beat) && e>=length(min_beat)-1; break; end % no beats
     if max_beat(i) < min_beat(e); 
        i=i+1;
     else
        rstr =  sprintf('rejecting beat (%d) [%i-%i-%i]: ', ...
              i, min_beat(e), max_beat(i), min_beat(e+1));
        heart_trough = seq1(min_beat(e+[0:1])); % select the troughs between peaks
    	heart_peak  = seq1(max_beat(i)) - mean(heart_trough); % select the heart peaks
        if max_beat(i) - min_beat(e) < beat_sep/1.2;
           disp([rstr,'too soon after previous beat']);
        elseif seq1(max_beat(i)) < 0; 
          disp([rstr,'The beat is too small']);
         elseif abs(diff(heart_trough))/heart_peak > 0.8;
          disp([rstr,'The baseline is too inconsistant']);
        else %accept beat
           beats(end+1,:) = [min_beat(e), max_beat(i), min_beat(e+1)]; 
        end
        i=i+1; e=e+1;
     end
  end

function out= show_breaths(dd,ii)
  if ischar(dd) && strcmp(dd,'TITLE');
     out = 'Breaths'; return
  end
  if ischar(dd) && strcmp(dd,'REQBREATHS?');
     out = true; return
  end
  if ischar(dd) && strcmp(dd,'REQBEATS?');
     out = false; return
  end
  fout = sprintf('breaths%03d.png',ii);
  out = sprintf( ...
   '<a href="%s"><img width="300" src="%s"></a>',...
   fout, fout);
  
  clf; subplot(211);
  plot(dd.tt,dd.CV,'LineWidth',4); box off;
  H = (max(dd.CV) - min(dd.CV))/10;
  for i=1:dd.n_breaths
     eie = dd.breaths(i,[1,3,2,1]);
     line(dd.tt(eie), dd.CV(eie), 'color',[0,0,0],'LineWidth',2);
  end
  axis tight
  print_convert(fout);

function out= model_breaths(dd,ii)
  if ischar(dd) && strcmp(dd,'TITLE');
     out = 'Model fit breaths'; return
  end
  if ischar(dd) && strcmp(dd,'REQBREATHS?');
     out = true; return
  end
  if ischar(dd) && strcmp(dd,'REQBEATS?');
     out = false; return
  end
  fout = sprintf('model_breaths%03d.png',ii);
  out = sprintf( ...
   '<a href="%s"><img width="300" src="%s"></a>',...
   fout, fout);
  
  xlimits = [0,dd.lD/dd.FR];
  clf; subplot(211);
  plot(dd.tt,dd.CV,'LineWidth',4); box off;
  for i=1:dd.n_breaths
    eie = dd.breaths(i,[1,2,3]);
    Lie = diff(eie)+1;
    Teie = dd.tt(eie);
    Ceie = dd.CV(eie);
    inspL= linspace(Ceie(1),Ceie(2),Lie(1));
    insp = dd.CV(eie(1):eie(2)) - inspL;
    insp = low_order_fourier(insp,4) + inspL;
    line(dd.tt(eie(1):eie(2)), insp, 'color',[0,.5,.1],'LineWidth',2);
  
    expiL= linspace(Ceie(2),Ceie(3),Lie(2));
    expi = dd.CV(eie(2):eie(3)) - expiL;
    expi = low_order_fourier(expi,4) + expiL;
    line(dd.tt(eie(2):eie(3)), expi, 'color',[0,.4,.2],'LineWidth',2);
  end

  H = (max(dd.CV) - min(dd.CV))/10;
  for i=1:dd.n_breaths
     eie = dd.breaths(i,[1,3,2,1]);
     line(dd.tt(eie), dd.CV(eie), 'color',[0,0,0],'LineWidth',2);
  end
  axis tight
  xlim(xlimits);

  subplot(212);
  ls = linspace(0,1,10); ls = ls/sum(ls)/10;
  ls = [ls,-fliplr(ls)];
  for i=1:dd.n_breaths
    eie = dd.breaths(i,[1,2,3]);
    Lie = diff(eie)+1;
    Teie = dd.tt(eie);
    Ceie = dd.CV(eie);
    inspL= linspace(Ceie(1),Ceie(2),Lie(1));
    insp = dd.CV(eie(1):eie(2)) - inspL;
    insp = low_order_fourier(insp,7);% + inspL;
    iflow = -conv2(insp,ls,'same') - diff(Ceie(1:2))/Lie(1);
    line(dd.tt(eie(1):eie(2)), iflow, 'color',[0,.5,.1],'LineWidth',2);
  
    expiL= linspace(Ceie(2),Ceie(3),Lie(2));
    expi = dd.CV(eie(2):eie(3)) - expiL;
    expi = low_order_fourier(expi,7);% + expiL;
    eflow = -conv2(expi,ls,'same') - diff(Ceie(2:3))/Lie(2); ;
    line(dd.tt(eie(2):eie(3)), eflow, 'color',[0,.4,.2],'LineWidth',2);
  end
  xlim(xlimits);
  print_convert(fout);
 

function s = low_order_fourier(s,N);
  f = fft(s);
  f([N+2:end-N]) = 0;
  s= ifft(f);
  if norm(imag(s)) > 1e-13; error('FFT'); end
  s= real(s);


function out= show_volume_n_flow(dd,ii)
  if ischar(dd) && strcmp(dd,'TITLE');
     out = 'Volume/Flow vs t'; return
  end
  if ischar(dd) && strcmp(dd,'REQBREATHS?');
     out = true; return
  end
  if ischar(dd) && strcmp(dd,'REQBEATS?');
     out = false; return
  end
  fout = sprintf('vol_flow%03d.png',ii);
  out = sprintf( ...
   '<a href="%s"><img width="300" src="%s"></a>',...
   fout, fout);
  
  clf; subplot(211);
  dd.CV = detrend(dd.CV);
  idx = min(dd.breaths(:)):max(dd.breaths);
  plot(dd.tt(idx),dd.CV(idx),  ...
       'Color',[0,0,1],'LineWidth',3); box off;
  line(dd.tt(idx),dd.flow(idx), ...
       'Color',[0.5,0.2,1],'LineWidth',3);
  axis tight
  print_convert(fout);

function out= stats_volume_n_flow(dd,ii)
  global pp;
  if ischar(dd) && strcmp(dd,'TITLE');
     out = 'Volume/Flow stats';
     return
  end
  if ischar(dd) && strcmp(dd,'REQBREATHS?');
     out = true; return
  end
  if ischar(dd) && strcmp(dd,'REQBEATS?');
     out = false; return
  end
  TV = []; Pif = []; Pef = [];
  for i=1:dd.n_breaths
     eie = dd.breaths(i,[1,2,3]);
     TV(i) = dd.CV(eie(2)) - mean(dd.CV(eie([1,3])));
     Pif(i) = max(-dd.flow(eie(1):eie(2)));
     Pef(i) = max( dd.flow(eie(2):eie(3)));
  end
  out= sprintf(['<table border=1>' ...
    '<TR><TH>Param<TH>#' ...
    '<TR><TD>&Delta;Z<sub>total</sub><TD>%1.3f' ...
    '<TR><TD>Peak<sub>insp flow</sub><TD>%1.3f' ...
    '<TR><TD>Peak<sub>expi flow</sub><TD>%1.3f' ...
    '</table>'], ...
     mean(TV), mean(Pif), mean(Pef));


function out= flow_volume_global(dd,ii)
  global pp; % parameters
  if ischar(dd) && strcmp(dd,'TITLE');
     out = 'Flow-Volume';
     return
  end
  if ischar(dd) && strcmp(dd,'REQBREATHS?');
     out = true; return
  end
  if ischar(dd) && strcmp(dd,'REQBEATS?');
     out = false; return
  end
  clf;plot(NaN); hold on;
  for i=1:dd.n_breaths
     eie = dd.breaths(i,[1,2,3]);
     eflo= eie(2) + pp.flow_window;
     vol = dd.CV(eie(1):eie(3));
     volD= dd.CV(eflo);
     volD= volD-vol(1);
     vol = vol-vol(1);

     flow=    dd.flow(eie(1):eie(3));
     flowD=   dd.flow(eflo);
     plot(vol,flow,'LineWidth',2,'Color',[0,0,1]); 
     plot(volD,flowD,'LineWidth',4,'Color',[0,0,0]); 
  end
  hold off;
  fout = sprintf('g_flow_vol%03d.png',ii);
  out = sprintf( ...
   '<a href="%s"><img width="300" src="%s"></a>',...
   fout, fout);
  print_convert(fout);

function out= flow_volume_components(dd,ii)
  global pp; % parameters
  if ischar(dd) && strcmp(dd,'TITLE');
     out = 'Flow-Volume Components';
     return
  end
  if ischar(dd) && strcmp(dd,'REQBREATHS?');
     out = true; return
  end
  if ischar(dd) && strcmp(dd,'REQBEATS?');
     out = false; return
  end
  clf;plot(NaN); hold on;
  intersect=[];
  for i=1:dd.n_breaths
     eie = dd.breaths(i,[1,2,3]);
     eflo= eie(2) + pp.flow_window;
     efla= eie(3) - (40:130);
     vol = dd.CV(eie(1):eie(3));
     volD= dd.CV(eflo);
     volA= dd.CV(efla);
     volA= volA-vol(1);   volD= volD-vol(1);   vol = vol-vol(1);
     volA= volA/max(vol); volD= volD/max(vol); vol = vol/max(vol);
     volA= 100*volA;      volD= 100*volD;      vol = 100*vol;

     flow=    dd.flow(eie(1):eie(3));
     flowD=   dd.flow(eflo);
     flowA=   dd.flow(efla);
     plot(vol,flow,'LineWidth',2,'Color',[0,0,1]); 
     plot(volD,flowD,'LineWidth',4,'Color',[0,0.5,0]); 
     plot(volA,flowA,'LineWidth',4,'Color',[0.5,0,0]); 

     pD = polyfit(volD,flowD,1); idx = [30,90];
     plot(idx,polyval(pD,idx),'LineWidth',1,'Color',[0,0.5,0]); 
     pA = mean(flowA); idx = [5,50];
     plot(idx,polyval(pA,idx),'LineWidth',1,'Color',[0.5,0,0]); 
     intersect(i) = (pA-pD(2))/pD(1); 
     plot(intersect(i),pA,'k*');
  end
  hold off; xlim([0,100]);
  fout = sprintf('g_fv_comp%03d.png',ii);
  out = sprintf(['<center>Intercept = %3.1f%%<br>' ...
     '<a href="%s"><img width="300" src="%s">', ...
     '</a></center>'],...
     median(intersect), fout, fout);
  print_convert(fout);

function out= flow_volume_global_slope(dd,ii)
  global pp; % parameters
  if ischar(dd) && strcmp(dd,'TITLE');
     out = 'F-V Slope';
     return
  end
  if ischar(dd) && strcmp(dd,'REQBREATHS?');
     out = true; return
  end
  if ischar(dd) && strcmp(dd,'REQBEATS?');
     out = false; return
  end
  out= '<table border=1><TR><TH>#<TH>Slope';
  for i=1:dd.n_breaths
     eie = dd.breaths(i,[1,2,3]);
     eflo= eie(2) + pp.flow_window;
     vol = dd.CV(eie(1):eie(3));
     volD= dd.CV(eflo);
     volD= volD-vol(1);
     vol = vol-vol(1); 
     flow=    dd.flow(eie(1):eie(3));
     flowD=   dd.flow(eflo);

     pf = polyfit(volD,flowD,1); slope = pf(1);
     out = [out, sprintf( ...
        '<TR><TD>%d<TD>%1.3f', i, slope)];
  end
  out = [out,'</table>'];

function TV = TVcalc(dd);
  TV = [];
  for i=1:dd.n_breaths
    eie = dd.breaths(i,[1,2,3]);
    TV(:,:,i) = dd.ZR(:,:,eie(2)) - mean(dd.ZR(:,:,eie([1,3])),3);
  end
  TV= mean(TV,3);

function IM = my_image(IM)
  global pp;
  switch pp.rotate
    case 0; % nothing
    case 180; IM= flipud(fliplr(IM));
    otherwise; error 'rotation value not understood';
  end
  if nargout == 1; return; end
  image(IM); axis off;

function out= TV_image(dd,ii)
  global pp;
  if ischar(dd) && strcmp(dd,'TITLE');
     out = 'TV Image';
     return
  end
  if ischar(dd) && strcmp(dd,'REQBREATHS?');
     out = true; return
  end
  if ischar(dd) && strcmp(dd,'REQBEATS?');
     out = false; return
  end
  TV = TVcalc(dd);
  TV(dd.ZR(:,:,1)==0) = NaN;
  mycolormap;
  my_image(TV*250/max(TV(:))+5);
  
  fout = sprintf('TV_image%03d.png',ii);
  out = sprintf(['<center>max pixel=%1.3f<br>' ...
     '<a href="%s"><img width="200" src="%s">' ...
     '</a><p><img src="%s"></center>'], ...
     max(TV(:)), fout, fout, pp.colourbar);
  print_convert(fout);

%DO flow_volume_slices
function out= TV_slices(dd,ii)
  global pp;
  if ischar(dd) && strcmp(dd,'TITLE');
     out = 'TV Slices';
     return
  end
  if ischar(dd) && strcmp(dd,'REQBREATHS?');
     out = true; return
  end
  if ischar(dd) && strcmp(dd,'REQBEATS?');
     out = false; return
  end
  TV = TVcalc(dd);
  TV = my_image(TV); %rotate if req'd
  
  out = '<table border="1"><TR><TH>#<TH>Left %<TH>Right %';
  cuts = round(linspace(0,32,pp.slices+1));
  fac  = (32/pp.slices)/sum(TV(:))*100;
  for i=1:pp.slices
     ss = fac*mean(TV(cuts(i)+1:cuts(i+1),:),1);
     out = [out, sprintf(['<TR><TD align="right">%d' ...
             '<TD align="right">%3.1f<TD align="right">%3.1f'], ...
         i, sum(ss(1:16)), sum(ss(17:32)))];
  end
  out = [out,'</table>'];

% set color_map in config file
function mycolormap(cmap)
  if nargin==0;
     global pp;
     colormap( feval(pp.color_map,256) ); 
   % colormap(gray(256));
   % colormap(viridis(256));
   % colormap(ocean(256));
  else
     colormap( cmap );
  end

function [RGB] = blue_red_colours;

   ofs= 1;
   glev= 0.1;

   D= (2*ofs - 1);
   ofs= ofs - 2*(ofs==0);
   F= 3*0.1;
   DF= D*F; D_F= D/F;
   scale_data = [linspace(-3,0,100),linspace(0,3,156)];

   red= DF*abs(scale_data+D_F) - ofs;
   red= red.*(red>0).*(red<1) + (red>=1);

   grn= DF*abs(scale_data    ) - ofs;
   grn= grn.*(grn>0).*(grn<1) + (grn>=1);

   blu= DF*abs(scale_data-D_F) - ofs;
   blu= blu.*(blu>0).*(blu<1) + (blu>=1);
   RGB= [red;grn;blu]';
   
   vd = viridis(2);
   RGB(1,:) = vd(1,:);

function [FV,ROI] = FV_calc(dd);
  global pp;
  FV = TVcalc(dd);
  mFV = max(FV(:));
  FV(dd.ZR(:,:,1)==0) = NaN;
  ROI = FV(:,:,1); ROI(~isnan(ROI))= 0;
   
  for i=1:32; for j=1:32
     if isnan(FV(i,j));
        1; % do nothing
     elseif FV(i,j) < 0.1*mFV
        FV(i,j) = 0;
        ROI(i,j) = 0;
     else
        ROI(i,j) = 1;
        slope = [];
        for l=1:dd.n_breaths
           eie = dd.breaths(l,[1,2,3]);
           eflo= eie(2) + pp.flow_window;
           volD= dd.ZR(i,j,eflo);
           flowD=dd.ZF(i,j,eflo);
           pf = polyfit(squeeze(volD),squeeze(flowD),1);
           slope(l) = pf(1);
        end
        FV(i,j) = mean(slope);
     end
  end; end

function out= flow_volume_image(dd,ii)
  global pp;
  if ischar(dd) && strcmp(dd,'TITLE');
     out = 'Flow-Volume Image';
     return
  end
  if ischar(dd) && strcmp(dd,'REQBREATHS?');
     out = true; return
  end
  if ischar(dd) && strcmp(dd,'REQBEATS?');
     out = false; return
  end
  [FV,ROI] = FV_calc(dd);
  mycolormap;
    
  outimg = FV*200/max(FV(:));
% Issue ... what if some components are negative
  outimg(ROI==1) = outimg(ROI==1) + 50;
  my_image(outimg); axis image;
  str= ['<center>max pixel=%1.3f<br>' ...
   '<a href="%s"><img width="200" src="%s">' ...
   '</a><p><img src="%s"></center>'];
  
  fout = sprintf('FV_image%03d.png',ii);
  out = sprintf(str, max(FV(:)), fout, fout, pp.colourbar);
  print_convert(fout);

function out= flow_volume_image_change(dd,ii)
  global pp;
  if ischar(dd) && strcmp(dd,'TITLE');
     out = 'Flow-Volume Image &Delta;';
     return
  end
  if ischar(dd) && strcmp(dd,'REQBREATHS?');
     out = true; return
  end
  if ischar(dd) && strcmp(dd,'REQBEATS?');
     out = false; return
  end
  [FV,ROI] = FV_calc(dd);
  mycolormap(blue_red_colours);

% TODO Set parameter to choose BL image
  if ii==1; pp.flow_volume_image_bl= FV; end
    
  FVd= FV-pp.flow_volume_image_bl;
  outimg = FVd*100;
  outimg(ROI==1) = outimg(ROI==1) + 100;
  outimg( ROI==1 & outimg<2  ) = 2;
  outimg( ROI==1 & outimg>255) = 255;
  
   my_image(outimg); axis image
  str= ['<center>max &Delta;=%1.3f from BL<br>' ...
  '<a href="%s"><img width="200" src="%s">' ...
  '</a><p><img src="%s"></center>'];
  
  fout = sprintf('FVd_image%03d.png',ii);
  out = sprintf(str, max(FVd(:)), fout, fout, pp.colourbar);
  print_convert(fout);

function out= flow_volume_slices(dd,ii)
  global pp;
  if ischar(dd) && strcmp(dd,'TITLE');
     out = 'Flow-volume Slices';
     return
  end
  if ischar(dd) && strcmp(dd,'REQBREATHS?');
     out = true; return
  end
  if ischar(dd) && strcmp(dd,'REQBEATS?');
     out = false; return
  end
  FV = FV_calc(dd);
  FV = my_image(FV); %rotate if req'd
  FV(isnan(FV)) = 0;
  
  out = '<table border="1"><TR><TH>#<TH>Left %<TH>Right %';
  cuts = round(linspace(0,32,pp.slices+1));
  fac  = (32/pp.slices)/sum(FV(:))*100;
  for i=1:pp.slices
     ss = fac*mean(FV(cuts(i)+1:cuts(i+1),:),1);
     out = [out, sprintf(['<TR><TD align="right">%d' ...
             '<TD align="right">%3.1f<TD align="right">%3.1f'], ...
         i, sum(ss(1:16)), sum(ss(17:32)))];
  end
  out = [out,'</table>'];

function breaths= find_frc( data );
   % FIND_FRC: find candidates for FRC
   %   [einsp,eexpi]= find_frc( data);
   % OUTPUTS: 
   %  einsp - end inspiration
   %  eexpi - end expiration
   global pp;

   %seq= ROI*imgs.elem_data;
   seq = data.CV;
   lseq = length(seq);
   % LPF the signal
   Fseq= fft(seq);

   % Cut off freq
   % each point is frate/2/len Hz
   % want to cut a 0.25Hz = L *frate/2/len; L=CUTOFF *2*len/frate
   Fcutoff = 0.25;
   L = round( Fcutoff * 2*2*lseq/data.FR ); % TODO why is there 2*2? 
   Fseq([1+L+1:end-L])=0; %HPF
   Fseq([1,2,end])=0;     %LPF
   seq1= ifft(Fseq);
   % Test and take real part
   if std(imag(seq1))>1e-10; error('Breathing FFT code'); end
   seq1= real(seq1);
   % Flow calc
   flow = diff(seq1);
   thresh = median( abs(flow));

   inout= zeros(lseq,1);
   i = 1;
   inout(i) = sign( flow(i)+eps ); % eps to force not zero
   for i = 2:lseq-1
     cval= inout(i-1);
     if cval*flow(i) > -thresh
        inout(i) = cval;
     else
        inout(i) = -cval;
     end
   end
   inout(lseq) = inout(lseq-1);

% Shorter forward window because detection is lagging
   breath_window_fwd= round(pp.FRC_search_window/4*data.FR);
   breath_window_bak= round(pp.FRC_search_window*data.FR);

   dinout= diff( inout );
   fdiff = find( diff(inout) );
   fdiff(fdiff<=breath_window_bak      )= []; % too close
   fdiff(fdiff>lseq - breath_window_fwd)= []; % too close

   eexpi= fdiff( (dinout(fdiff)>0) );
   einsp= fdiff( (dinout(fdiff)<0) );

   % Find the best point
   ww= -breath_window_bak:breath_window_fwd;
   for i=1:length(eexpi);
     wind= seq( eexpi(i)+ ww );
     ff = find( wind== min(wind) );
     ff= ff(1)+min(ww)-1;
     eexpi(i)= eexpi(i) + ff;
   end
   for i=1:length(einsp);
     wind= seq( einsp(i)+ ww );
     ff = find( wind== max(wind) );
     ff= ff(1)+min(ww)-1;
     einsp(i)= einsp(i) + ff;
   end

   min_insp_length = round(pp.min_insp_length*data.FR);
   min_expi_length = round(pp.min_expi_length*data.FR);
   breaths = [];
   i=1;e=1; while true;
      if i>=length(einsp) && e>=length(eexpi)-1; break; end
      if einsp(i) < eexpi(e);
         i=i+1;
      else
         rstr =  sprintf('rejecting breath (%d) [%i-%i-%i]: ', ...
               i, eexpi(e), einsp(i), eexpi(e+1));
         FRCs= seq1(eexpi(e+[0:1]));
         TV  = seq1(einsp(i)) - mean(FRCs);
         if einsp(i) - eexpi(e) < min_insp_length 
            disp([rstr,'min_insp_length']);
         elseif eexpi(e+1)-einsp(i) < min_expi_length;
            disp([rstr,'min_expi_length']);
         elseif abs(diff(FRCs))/TV > pp.FRC_relative_match
            disp([rstr,'FRC_relative_match']);
         else %accept breath
            breaths(end+1,:) = [eexpi(e), einsp(i), eexpi(e+1)]; 
         end
         i=i+1; e=e+1;
      end
   end
   %[einsp,eexpi] = remove_some_points( einsp, eexpi, remove_pts );

     
function write_table_hdr;
  global pp;
  fid = outfile;
  fprintf(fid,'<TR><TH>Filename');
  for i=1:length(pp.callfns);
    out = feval(pp.callfns{i}, 'TITLE');
    fprintf(fid,'<TH>%s',out);
  end
  fprintf(fid,'<TH>Filename');
  fclose(fid);

  


function fid=outfile(w);
  if nargin==1; w='w'; else; w='a'; end
  fid = fopen('index.html',w);
  if fid<0; error('file problem'); end
  

function write_header;
  fid = outfile('w');
  fprintf(fid,[ ...
   '<HTML><TITLE>Analysis: %s</TITLE></HTML>' ...
   '<BODY><Table>' ...
   '<TR><TD>Date:<TD>%s' ...
   '<TR><TD>Path:<TD>%s' ...
   '<TR><TD>Software:<TD>%s' ...
   '</Table>' ...
   '<H1>Results</H1>' ...
   '<Table border="1">'], date, date, pwd, '$Date::            $');
  fclose(fid);

function write_footer;
  fid = outfile;
  fprintf(fid,[ ...
   '</Table>' ...
   '</BODY></HTML>']);
  fclose(fid);


function print_convert(filename, varargin)
%PRINT_CONVERT: print figures with anti-aliasing and trim them
%  PRINT_CONVERT(FILENAME,OPT) prints current figure to FILENAME, which 
%  must include extenstion.
%
%  OPT is either a string specifying dpi in the format '-r150' or a struct
%  with (some of) these fields:
%   opt.resolution   = 150;              % 150 dpi (default: 125)
%   opt.pagesize     = [width, height];  % whatever matlab's current units
%   opt.jpeg_quality = 85;               % jpeg quality (default: 85)
%   opt.imwrite_opts = {'BitDepth',8};   % options to IMWRITE (default: '')
%   opt.horz_cut     = 50;               % remove horizontal gaps larger
%                                        % than 50 pixels (default: 0)
%   opt.horz_space   = 10;               % replace the removed horizontal
%                                        % gaps with 10 px of background
%   opt.vert_cut     = 50;               % remove vertical gaps larger
%                                        % than 50 pixels (default: 0)
%   opt.vert_space   = 10;               % replace the removed vertical
%                                        % gaps with 10 px of background
%   opt.supersampling_factor = 2;        % anti-aliasing (default: 1 for
%                                        % images, 2 for graphs). Higher is
%                                        % smoother.
%   opt.crop_slack   = [top,bot,left,right] %don't crop right to boundary
%   opt.figno        = number            % print figure, otherwise gca
%
%  Note that opt.imwrite_opts takes precedence over opt.jpeq_quality.
%
%  print_convert has pre-version 3.7 options, which are deprecated
%
%  Examples
%   print_convert('outname.png')         % uses default options
%   print_convert outname.png
%   print_convert outname.png -r150       % set dpi=150
%   print_convert('outname.png',opts);   % use options
%
% See also IMWRITE
 
% (C) Andy Adler and Bartlomiej Grychtol 2010-2013. License: GPL v2 or v3.
% $Id$
%
%  Compatibility with pre-3.7 features:
%  ------------------------------------
%  PRINT_CONVERT(FILENAME, OPTIONS, PAGEHWR)
%
%  FILENAME : filename to print to file type is the extension
%  OPTIONS  : specify dpi as '-density N'
%  PAGEHWR  : set page hight/width to control shape of figures
%  
%   print_convert('outname.png','-density 150') % resolution to 150 dpi
%   print_convert('outname.png','-density 150', 0.5)
 

pp = parse_options(filename, varargin{:});

tmpnam = [tempname,'.png'];

old_ihc = get(gcf, 'InvertHardcopy');
old_col = get(gcf, 'Color');
set(gcf,'InvertHardCopy','off'); 
set(gcf,'Color','w');

set(gcf,'PaperPosition',pp.posn); % I wish matlab had unwind protect - like octave does!
%%% The -dpng driver is broken in R2014a on linux (and maybe others);
print(pp.figno,'-dpng',pp.resolution,tmpnam);
set(gcf,'PaperPosition',pp.page);
 
set(gcf,'InvertHardCopy',old_ihc); % 
set(gcf,'Color',old_col);
 
im = imread(tmpnam,'png');
delete(tmpnam);

im = bitmap_downsize(im, pp.factor);
im = crop_image(im,pp);
try
   imwrite(im,filename,pp.imwrite_opts{:});
catch e
   eidors_msg(['Call to IMWRITE failed.'...
               'Probably opt.imwrite_opts is incorrect for %s files.'],...
               upper(pp.fmt), 1);
   disp('opt.imwrite_opts:');
   disp(pp.imwrite_opts);
   rethrow(e);
end

function im = crop_image(im,pp)
   tu = pp.crop_slack(1);
   bu = pp.crop_slack(2) + 1; %remove starting one more
   lu = pp.crop_slack(3);
   ru = pp.crop_slack(4) + 1;

   szim = size(im);
   bdr = squeeze(mean(double(im(1,:,:)),2));
 
   isbdr = true(szim(1),szim(2));
   for i=1:szim(3);
     isbdr = isbdr & (im(:,:,i) == bdr(i));
   end
 
   horz = [true,all(isbdr,1),true];
   horzpt = find(diff(horz)) - 1; % first 'true'
   if isempty(horzpt)
      eidors_msg('Image is blank. Cropping aborted.',1);
      return
   end
   im(:,horzpt(end)+ru:end,:)= []; % remove higher first
   if pp.horz_cut >0;
      horz_e_pt = find(diff(horz)==-1) -1; horz_e_pt(1) = [];
      horz_s_pt = find(diff(horz)==+1)   ; horz_s_pt(end) = [];
      idx = find(horz_e_pt - horz_s_pt > pp.horz_cut);
      for i=fliplr(idx) % remove higher first
        im(:,horz_s_pt(i)+pp.horz_space:horz_e_pt(i),:)= [];
      end
   end
   im(:,1:horzpt(1)-lu ,:)= [];
 
   vert = [true,all(isbdr,2)',true];
   vertpt = find(diff(vert)) - 1; % first 'true'
   im(vertpt(end)+bu:end,:,:)= [];
   if pp.vert_cut >0;
      vert_e_pt = find(diff(vert)==-1) -1; vert_e_pt(1) = [];
      vert_s_pt = find(diff(vert)==+1)   ; vert_s_pt(end) = [];
      idx = find(vert_e_pt - vert_s_pt > pp.vert_cut);
      for i=fliplr(idx) % remove higher first
        im(vert_s_pt(i)+pp.vert_space:vert_e_pt(i),:,:)= [];
      end
   end
   im(1:vertpt(1)-tu ,:,:)= [];
 
    
% factor = 1 if all plots only contain images, 2 otherwise
function f = default_factor
   f = 1; % default for images
   sp = get(gcf,'Children'); % subplots
   for i = 1:length(sp)
      obj = get(sp(i),'Children');
      tp  = get(obj,'Type');
      if ~all(strcmp(tp,'image'))
         f = 2;
         return;
      end
   end 
   
function fmt = parse_format(filename)   
    ext = lower(regexp(filename,'(?<=\.).+$','match'));
    if isempty(ext); error('no filename extension detected (%s)',filename); end
    switch ext{1}
       case {'jpg', 'jpeg'}
          fmt = 'jpg';
       case {'j2c', 'j2k', 'jp2'}
          fmt = 'jp2';
       case {'tif','tiff'}
          fmt = 'tif';
       otherwise
          fmt = ext{1};
    end
       
function pp = parse_options(filename,varargin)
   
   pp.fmt = parse_format(filename);

   pp.page = get(gcf,'PaperPosition');
   pp.posn = pp.page;
   pp.jpeg_quality = 85; % default jpeg quality
   pp.imwrite_opts = {}; % initial
   pp.horz_cut = 50;
   pp.horz_space = 10;
   pp.vert_cut = 50;
   pp.vert_space = 10;
   pp.factor = default_factor;
   pp.resolution = sprintf('-r%d',125 * pp.factor);
   pp.crop_slack = [0,0,0,0];
   pp.figno = gcf; % default
   pp.LP_filter = inf; % Don't LP filter
   pp.model_breaths_ncos = 4;
   
   
 
% Old options
   if nargin< 2;   
      return; 
   end
   if nargin>=3
      pp.posn(4) = pp.posn(3)*varargin{2};  
   end
 
   opt = varargin{1};
   if ischar(opt)
      val =regexp(opt,'-density (\d+)','tokens');
      if ~isempty(val);
         pp.resolution = sprintf('-r%d', str2double(val{1}{1}) * pp.factor);
      end
      val =regexp(opt,'-r(\d+)','tokens');
      if ~isempty(val);
         pp.resolution = sprintf('-r%d', str2double(val{1}{1}) * pp.factor);
      end
   elseif isstruct(opt)
     if isfield(opt,'supersampling_factor')
        pp.factor = opt.supersampling_factor;
     end
     if isfield(opt,'resolution');
         pp.resolution = sprintf('-r%d', opt.resolution * pp.factor);
     else
         pp.resolution = sprintf('-r%d',125 * pp.factor);
     end
     if isfield(opt,'pagesize');
         pp.posn(3:4) = opt.pagesize;
     end
     % TODO, this code can copy from opt to pp
     if isfield(opt,'jpeg_quality')
        pp.jpeg_quality = opt.jpeg_quality;
     end
     if strcmp(pp.fmt,'jpg')
        pp.imwrite_opts = {'quality',pp.jpeg_quality}; 
     end
     if isfield(opt,'imwrite_opts');
        pp.imwrite_opts = opt.imwrite_opts;
        if strcmp(pp.fmt,'jpg') && ~any(strcmpi(pp.imwrite_opts,'quality'))
           pp.imwrite_opts(end+1:end+2) = {'quality',pp.jpeg_quality};
        end
     end
    
     if isfield(opt,'horz_cut');
         pp.horz_cut = opt.horz_cut;
     end
     if isfield(opt,'vert_cut');
         pp.vert_cut = opt.vert_cut;
     end
     if isfield(opt,'vert_space');
        if opt.vert_space >= pp.vert_cut;
           warrning('Option vert_space must be smaller than vert_cut. Ingoring');
        else
           pp.vert_space = opt.vert_space;
        end
     end
     if isfield(opt,'horz_space');
        if opt.horz_space >= pp.horz_cut;
           warrning('Option vert_space must be smaller than vert_cut. Ingoring');
        else
           pp.horz_space = opt.horz_space;
        end
     end
     if isfield(opt,'crop_slack');
        pp.crop_slack = opt.crop_slack;
     end
     if isfield(opt,'figno');
        pp.figno = opt.figno;
     end
   else
      error('Can''t parse options');
   end
 

function A = bitmap_downsize(A, factor)
%BITMAP_DONWSIZE scale down an image (anti-aliasing)
% A = downsize(A, N) downsizes image array A by a factor N (natural number)
%
% This function is part of the EXPORT_FIG suite by Oliver Woodford
% http://www.mathworks.com/matlabcentral/fileexchange/23629

% Copyright (C) Oliver Woodford 2008-2012
% License: BSD, see accompanying license.txt
% $Id$

if factor == 1
    % Nothing to do
    return
end
try
    % Faster, but requires image processing toolbox
    A = imresize(A, 1/factor, 'bilinear');
catch
    % No image processing toolbox - resize manually
    % Lowpass filter - use Gaussian as is separable, so faster
    % Compute the 1d Gaussian filter
    filt = (-factor-1:factor+1) / (factor * 0.6);
    filt = exp(-filt .* filt);
    % Normalize the filter
    filt = single(filt / sum(filt));
    % Filter the image
    padding = floor(numel(filt) / 2);
    for a = 1:size(A, 3)
        A(:,:,a) = conv2(filt, filt', single(A([ones(1, padding) 1:end repmat(end, 1, padding)],[ones(1, padding) 1:end repmat(end, 1, padding)],a)), 'valid');
    end
    % Subsample
    A = A(1+floor(mod(end-1, factor)/2):factor:end,1+floor(mod(end-1, factor)/2):factor:end,:);
end



% butter from octave signal package (License: GPL)
function [a, b] = butter (n, w)


  stop = false;
  digital = true;


  %# Prewarp to the band edges to s plane
  if (digital)
    T = 2;       # sampling frequency of 2 Hz
    w = 2 / T * tan (pi * w / T);
  end

  %# Generate splane poles for the prototype Butterworth filter
  %# source: Kuc
  C = 1;  ## default cutoff frequency
  pole = C * exp (1i * pi * (2 * [1:n] + n - 1) / (2 * n));
  if (mod (n, 2) == 1)
    pole((n + 1) / 2) = -1;  ## pure real value at exp(i*pi)
  end
  zero = [];
  gain = C^n;

  %# splane frequency transform
% [zero, pole, gain] = sftrans (zero, pole, gain, w, stop);
% Low-pass code from sftrans
  gain = gain * (C/w)^(length(zero)-length(pole));
  pole = w * pole / C;
  zero = w * zero / C;

  %# Use bilinear transform to convert poles to the z plane
  if (digital)
    [zero, pole, gain] = bilinear (zero, pole, gain, T);
  end

  %# convert to the correct output form
    a = real (gain * poly (zero));
    b = real (poly (pole));

% filtfilt from octave signal package (License: GPL)
function y = filtfilt(b, a, x)

  rotate = (size(x,1)==1);
  if rotate,                    # a row vector
    x = x(:);                   # make it a column vector
  end

  lx = size(x,1);
  a = a(:).';
  b = b(:).';
  lb = length(b);
  la = length(a);
  n = max(lb, la);
  lrefl = 3 * (n - 1);
  if la < n, a(n) = 0; end
  if lb < n, b(n) = 0; end

  %# Compute a the initial state taking inspiration from
  %# Likhterov & Kopeika, 2003. "Hardware-efficient technique for
  %#     minimizing startup transients in Direct Form II digital filters"
  kdc = sum(b) / sum(a);
  if (abs(kdc) < inf) # neither NaN nor +/- Inf
    si = fliplr(cumsum(fliplr(b - kdc * a)));
  else
    si = zeros(size(a)); # fall back to zero initialization
  end
  si(1) = [];

  for (c = 1:size(x,2)) # filter all columns, one by one
    v = [2*x(1,c)-x((lrefl+1):-1:2,c); x(:,c);
         2*x(end,c)-x((end-1):-1:end-lrefl,c)]; # a column vector

    %# Do forward and reverse filtering
    v = filter(b,a,v,si*v(1));                   # forward filter
    v = flipud(filter(b,a,flipud(v),si*v(end))); # reverse filter
    y(:,c) = v((lrefl+1):(lx+lrefl));
  end

  if (rotate)                   # x was a row vector
    y = rot90(y);               # rotate it back
  end

% biliear from signal toolbox (license: GPL)
function [Zz, Zp, Zg] = bilinear(Sz, Sp, Sg, T)

  if nargin==3
    T = Sg;
    [Sz, Sp, Sg] = tf2zp(Sz, Sp);
  elseif nargin!=4
    print_usage;
  end

  p = length(Sp);
  z = length(Sz);
  if z > p || p==0
    error("bilinear: must have at least as many poles as zeros in s-plane");
  end

%# ----------------  -------------------------  ------------------------
%# Bilinear          zero: (2+xT)/(2-xT)        pole: (2+xT)/(2-xT)
%#      2 z-1        pole: -1                   zero: -1
%# S -> - ---        gain: (2-xT)/T             gain: (2-xT)/T
%#      T z+1
%# ----------------  -------------------------  ------------------------
  Zg = real(Sg * prod((2-Sz*T)/T) / prod((2-Sp*T)/T));
  Zp = (2+Sp*T)./(2-Sp*T);
  if isempty(Sz)
    Zz = -ones(size(Zp));
  else
    Zz = [(2+Sz*T)./(2-Sz*T)];
    Zz = postpad(Zz, p, -1);
  end

  if nargout==2, [Zz, Zp] = zp2tf(Zz, Zp, Zg); end

