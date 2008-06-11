function [eexpi,einsp]= find_frc( imgs, ROI, frate, name, ok, remove_pts)
% FIND_FRC: find candidates for FRC
% points= find_frc( seq, frate)
% ok ==0 (show graph - default)
%      1 (show graph - no click) 
%      2 no graph
%
% Find candidates for FRC from a time seq of EIT data
% frate is framerate (in fps units)

% $Id: find_frc.m,v 1.1 2008-06-11 14:49:06 aadler Exp $

if nargin <5; ok=0; end
if nargin <6; remove_pts=[]; end

seq= ROI*imgs.elem_data;
lseq = length(seq);
% LPF the signal
Fseq= fft(seq);

% Cut off freq
% each point is frate/2/len Hz
% want to cut a 0.25Hz = L *frate/2/len; L=CUTOFF *2*len/frate
L = round( 0.5*2*lseq/frate );
Fseq([1+L+1:end-L])=0; %HPF
Fseq([1,2,end])=0;     %LPF
seq1= ifft(Fseq);

% Test and take real part
if std(imag(seq1))>1e-10; error('FFT code'); end
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

dinout= diff( inout );
fdiff = find( diff(inout) );
fdiff([1,end])= []; % first and last are unreliable
fdiff = setdiff( fdiff, remove_pts);

eexpi= fdiff( (dinout(fdiff)>0) );
einsp= fdiff( (dinout(fdiff)<0) );

window= 3; ww= -window:window;
for i=1:length(eexpi);
  wind= seq( eexpi(i)+ ww );
  ff = find( wind== min(wind) );
  ff= ff(1)-window-1;
  eexpi(i)= eexpi(i) + ff;
end
for i=1:length(einsp);
  wind= seq( einsp(i)+ ww );
  ff = find( wind== max(wind) );
  ff= ff(1)-window -1;
  einsp(i)= einsp(i) + ff;
end

[einsp,eexpi] = remove_some_points( einsp, eexpi, remove_pts );

if ok==2; return; end
  
clf;
axes('position',[0,.00,.3,1.0]);
imgss= imgs; imgss.elem_data= imgs.elem_data(:,[eexpi(1),einsp(1)]);
show_slices(imgss);

axes('position',[.3,.05,.7,.9]);

if ok==1;
   plotpoints( seq, eexpi, einsp, name, ok);
   nname = name(find(name=='/')+1:find(name=='.')-1);
   gg= get(gcf,'paperposition');
   set(gcf,'paperposition',[gg(1:3),gg(3)*.3]);
   print('-dpng','-r50',[nname,'-sig.png']); 
   set(gcf,'paperposition',gg);
   return;
end

while 1;
   xpts=plotpoints( seq, eexpi, einsp, name, ok);
     
   if isempty(xpts); return; end
   [einsp,eexpi] = remove_some_points( einsp, eexpi, xpts);
end

%plot((0:lseq-1)/frate, seq, 'k' );
function x=plotpoints( seq, eexpi, einsp, name, ok);
   plot( seq, 'k' );
   hold on;
   maxs= max(seq)*1.1;
   mins= min(seq)*1.1;
%  plot( [1;1]*eexpi(:)', [mins;maxs]*ones(1,length(eexpi)), 'r')
%  plot( [1;1]*einsp(:)', [mins;maxs]*ones(1,length(einsp)), 'b')
   hh= plot( eexpi, seq(eexpi), 'bo'); set(hh,'LineWidth',4);
   hh= plot( einsp, seq(einsp), 'ro'); set(hh,'LineWidth',4);
   hold off;

if ok==0
   title(['(',name,') CLICK ON LINES TO REMOVE: RETURN TO QUIT']);
   [x,jnk] = ginput;
else 
   title(name);
end

function [einsp,eexpi] = remove_some_points( einsp, eexpi, xpts);
   for xx= xpts(:)'
      d_insp= abs(xx - einsp); 
      d_expi= abs(xx - eexpi); 
      if min(d_insp) < min(d_expi)
         ff= find(d_insp == min(d_insp));
         einsp(ff)= [];
      else
         ff= find(d_expi == min(d_expi));
         eexpi(ff)= [];
      end
   end
