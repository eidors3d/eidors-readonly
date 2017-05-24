function  [c_data, stats]= frame_timing_correct(raw_data, method, fwd_model)
% FRAME_TIMING_CORRECT: Corrects serially collected *.get data so that each
% frame appears to come from an instantaneous point in time.
%
% c_data = frame_timing_correct(filenames, framerates, method)
%   raw_data     = input data
%   method       = 'FFT' or 'linear' (interpolation method)
%   c_data       = output data 208 x N_frames
%   stats.change = change in output
%   fwd_model    = fwd_model (currently ignored)
%
%Corrects serially collected *.get data so that each frame appears 
%to come from an instantaneous point in time.
%outputs *_corrected.get & associated *_corrected.prl files
%intended to be a standalone compilable file.
%
% ISSUES WITH CODE:
%   - FFT method assumes wrap-around. This means the first and last several frames incorrect
%   - Code creates its own data reader - should use eidors_readdata
%   - Code should be able to parse the stim_meas patterns, not assume them
%
%CITATION_REQUEST:
%  AUTHOR: Rebecca Yerworth and Richard Bayford
%  TITLE: The effect of serial data collection on the accuracy of
%  electrical impedance tomography images
%  JOURNAL: Physiological Measurement
%  VOL: 34
%  NUM: 6
%  YEAR: 2013
%  PAGE: 659-669
%  LINK: http://iopscience.iop.org/0967-3334/34/6/659/pdf/0967-3334_34_6_659.pdf
%  DOI: 10.1088/0967-3334/34/6/659
%  PUBMED:  23719130
% 
%%
%try
 %   citeme(serial_correction)
%end

if nargin<=2
   method='fft'; % other choice is linear 
end
if nargin>=3
   eidors_msg('frame_timing_correct: ignoring fwd_model',1);
end

%generate protocol
%imdl.fwd_model.stimulation=mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'},1); 
%n_stim_pattern = size(imdl.fwd_model.stimulation,2);
%n_meas_pattern = size(imdl.fwd_model.stimulation(1).meas_pattern,1);
%prt_len = n_stim_pattern*n_meas_pattern;  %protocol lenght
prt_len=208;
prt_time=0:1/(prt_len):1-1/(prt_len); %subframe timings relative to frame length

    
switch upper(method)
   case 'FFT'
        NFFT=size(raw_data,2);
        c_ave=zeros(prt_len,floor(NFFT/2));
        c_data=zeros(prt_len,NFFT);
        for e=1:prt_len % for each electrode combination
            Y = fft(raw_data(e,1:NFFT)',NFFT)/NFFT;
            a=angle(Y(2:floor(NFFT/2)))-2*pi*(1:NFFT/2-1)'*prt_time(e)/NFFT;
            c_ave(e,1)=Y(1); %Y_phase corrected
            c_ave(e,2:floor(NFFT/2))=2*(cos(a).*abs(Y(2:floor(NFFT/2)))+1i*sin(a).*abs(Y(2:floor(NFFT/2))));
            if isreal(raw_data)
                c_data(e,:)=real(ifft(c_ave(e,:),NFFT))*NFFT;
            else
                c_data(e,:)=(ifft(c_ave(e,:),NFFT))*NFFT;
            end
            c_ave=c_ave(:,1);
        end
   case 'LINEAR'
        %do linear interpolation
        c_data=zeros(prt_len,size(raw_data,2)-1);
        c_data(:,1)=raw_data(:,1); %first time frame not lag corrected.
        for t = 2:size(raw_data,2)
            c_data(:,t)=(raw_data(:,t-1)+(raw_data(:,t)-raw_data(:,t-1)).*(1-prt_time'));
        end
        c_ave=mean(c_data,2); %average frame
   otherwise 
      error('no interpolation method %s available',method)
end
  

stats.change = abs((c_data(1:208,:)-raw_data)./raw_data);
