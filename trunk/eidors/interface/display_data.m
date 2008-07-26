function display_data( data, fwd_model )
% DISPLAY_DATA show measurement structure in web browser
%
% Usage:
%  display_data( data )
%
%   where data is an eidors data object
%
% Usage:
%  display_data( data, fwd_model )
%
%   where data is a matrix M x F (Measurements x Frames)
%   and fwd_model is an eidors fwd model structure
 
% (C) 2005 Andy Adler.  License: GPL version 2 or version 3
% $Id$

if nargin==1
    fwd_model= data.fwd_model;
end

d_name= 'Unknown';
if isfield(data,'name');
    dname= data.name;
end
fm_name= fwd_model.name;

fname= sprintf('%s.html', tempname );
fid= fopen(fname,'wt');
if fid==0; 
    error('can''t open file %s', fname );
end

%
% Write out header
%
fprintf(fid, [ ...
'<HTML><HEAD><TITLE>' ...
'Data display %s. Model %s' ...
'</TITLE></HEAD><BODY>' ...
'<H1>Data display. <tt>data.name</tt>= %s' ...
                  '<tt>fwd_model.name</tt>= %s</H1>'], ...
  d_name, fm_name, d_name, fm_name );

%
% Write out electrodes
%
fprintf(fid, [...
'<H2>Electrode information</H2>\n' ...
'<Table border="1"><TR>\n' ...
'<TH>Elec #<TH colspan=5>Posn<TH>Nodes<TH>Z<sub>c</sub>\n', ...
'<TR>' ...
'<TH><TH>x<TH>y<TH>z<TH>r<TH>&theta; (&deg;)<TH> <TH>(&Omega;)\n']);

fme= fwd_model.electrode;
n_elec= length(fme);
xyctr= mean(fwd_model.nodes,1);
for i=1:n_elec
    nodes= fme(i).nodes;
    posn= mean(fwd_model.nodes(nodes,:),1);
    posn= [posn, zeros(1,3-length(posn))]; % zero pad
    dxy= posn(1:2) - xyctr(1:2);
    radius= sqrt(norm(dxy));
    angle = 180/pi*atan2(dxy(2),dxy(1));
    nlist= sprintf('%d,', nodes);
    nlist= ['[', nlist(1:end-1),']'];
    fprintf( fid, ['\n<TR>' ...
    '<TD align="right">&nbsp;%d<TD align="right">&nbsp;%1.3f<TD align="right">&nbsp;%1.3f<TD align="right">&nbsp;%1.3f<TD align="right">&nbsp;%1.3f<TD align="right">&nbsp;%1.1f<TD>%s<TD align="right">&nbsp;%1.3f'], ...
    i,posn,radius,angle,nlist ,fme(i).z_contact);
end
fprintf( fid, '\n</Table>\n');

%
% Write out stimulation patterns
%

% Step 1: What are all the meas patterns used?
% It is possible that two different patterns
% use the same electrodes with different amplifications
% we ignore that possibility here
fms= fwd_model.stimulation;
lfms= length(fms);

elec_list = 2.^(0:n_elec-1)';
meas_list = [];
for i=1:lfms
    this_meas_list= (fms(i).meas_pattern~=0)*elec_list;
    meas_list = [meas_list, this_meas_list(:)'];
end
meas_list= unique(meas_list);
lml= length(meas_list);

fprintf(fid, [...
'<H2>Stimulation pattern information</H2>\n' ...
'<Table border="1"><TR>\n' ...
'<TH>Stimulation<TH colspan="%d">Measurement Patterns\n', ...
'<TR><TH>Patterns'],lml);
ne2= floor(n_elec)/2;
for i=1:lml
    d2b= dec2bin( meas_list(i),n_elec);
    for j=length(d2b):-4:1;
        d2b= [d2b(1:j),'<br>',d2b(j+1:end)];
    end
    fprintf(fid,'<TH aligh="center">%s', d2b );
end
fprintf( fid, '\n</Table>\n');


fclose(fid);

if isunix
    system(sprintf('mozilla file:///%s', fname));
else
    system(sprintf('start file:///%s', fname));
end

