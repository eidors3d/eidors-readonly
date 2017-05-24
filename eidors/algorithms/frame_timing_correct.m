function  serial_correct()
%Corrects serially collected *.get data so that each frame appears 
%to come from an instantaneous point in time.
%outputs *_corrected.get & associated *_corrected.prl files
%intended to be a standalone compilable file.
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
method='fft'; %to do:give choice of fft or liniear 
%choose file and load date
%Eidors_parsers
global Ehub
%get filenames & store in Ehub
try cd(Ehub.path.data) %set path to last used data file directory (if exists)
end
[filenames,Ehub.path.data, file_type]=uigetfile({...
    %'*.get;*.txt','all EIT data files'; ...
    '*.get','Gotenberg and Visis files, *.get'; ...
    %'*.txt','*.txt files';...
    %'*.*',  'All Files (*.*)' ...
    },'Select filename(s):', ...
    'MultiSelect', 'on');
if filenames == 0
        errordlg('no data selected')
        return
end
cd(Ehub.path.data)
if ischar(filenames)    %if only one file seclected puts as string not cell, so need to convert.
    filenames=cellstr(filenames);
end

%generate protocol
%imdl.fwd_model.stimulation=mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'},1); 
%n_stim_pattern = size(imdl.fwd_model.stimulation,2);
%n_meas_pattern = size(imdl.fwd_model.stimulation(1).meas_pattern,1);
%prt_len = n_stim_pattern*n_meas_pattern;  %protocol lenght
prt_len=208;
prt_time=0:1/(prt_len):1-1/(prt_len); %subframe timings relative to frame length
%scale = 1/13; %untill told otherwise, assume all *.get files collected at 13 frames per second
            %can get this infromation for the *.prl file

%for each file selected....
for i=1:max(size(filenames))
    clear global Ehub;
    global Ehub  %re-declare
    %Ehub.filenames=filenames(i);
%load data
    if file_type == 0
        errordlg('no data selected')
        return
    elseif file_type == 1
        [a,b,ext]=fileparts(filenames{i});%work out filetype
        if strcmp(ext,'.get')
            infile=fopen(filenames{i},'r'); %open file for reading
            if infile==-1;
                errordlg(['Cannot find the file: ',filename,'.get']);
                return
            end;
            fseek(infile,0,'eof');
            eframe=ftell(infile)/256/4;%number of frames in file
            fseek(infile,0,'bof');
            Ehub.raw(:,1:eframe,1)=fread(infile,[208,inf],'208*float32',48*4);%[read 208 boundary voltages, skip 48 pieces of other info] per frame
            %%Loading of addition information.
            fseek(infile,208*4,'bof');
            data2=fread(infile,[48,inf],'48*float32',208*4);
            fclose(infile); %close file
        else
            errordlg('not *.get file');
            return
        end
    end
    
    %correct data
    %[n_dataC n_baseC]=lag_correct(Ehub.raw,Ehub.imdl); %need trimmed version without graphs
    if strcmp(method,'FFT')
        NFFT=size(data,2);  %seems fast enough without truncateing data set to 2^n values
        c_ave=zeros(prt_len,floor(NFFT/2));
        c_data=zeros(prt_len,NFFT);
        for e=1:prt_len % for each electrode combination
            Y = fft(Ehub.raw(e,1:NFFT)',NFFT)/NFFT;
            a=angle(Y(2:floor(NFFT/2)))-2*pi*(1:NFFT/2-1)'*prt_time(e)/NFFT;
            c_ave(e,1)=Y(1); %Y_phase corrected
            c_ave(e,2:floor(NFFT/2))=2*(cos(a).*abs(Y(2:floor(NFFT/2)))+1i*sin(a).*abs(Y(2:floor(NFFT/2))));
            if isreal(Ehub.raw)
                c_data(e,:)=real(ifft(c_ave(e,:),NFFT))*NFFT;
            else
                c_data(e,:)=(ifft(c_ave(e,:),NFFT))*NFFT;
            end
            c_ave=c_ave(:,1);
        end
    else
        %do linear interpolation
        c_data=zeros(prt_len,size(Ehub.raw,2)-1);
        c_data(:,1)=Ehub.raw(:,1); %first time frame not lag corrected.
        for t = 2:size(Ehub.raw,2)
            c_data(:,t)=(Ehub.raw(:,t-1)+(Ehub.raw(:,t)-Ehub.raw(:,t-1)).*(1-prt_time'));
        end
        c_ave=mean(c_data,2); %average frame
    end
    
%output corrected file
   c_data(209:208+48,:)=data2; %attach V, I & auxillary data to cleaned data.  
   outfile=fopen([a b '_corrected' ext],'w'); %open file for writing  %check filename & add '_corrected'
   fseek(outfile,0,'bof');
   fwrite(outfile,c_data,'float32');
   fclose(outfile); %close file
   %copy associated .prl file
   infile=fopen([a b '.prl'],'r'); %open file for reading
   if infile==-1;
       warndlg(['Cannot find the file: ',b,'.prl']);
       fr=13; % default frame rate
   else
       prl=fread(infile);
       fclose(infile);
       outfile=fopen([a b '_corrected.prl'],'w');
       fwrite(outfile,prl);
       fclose(outfile);
       s=importdata([a b '.prl'],'\t');
       line=find(strncmp('Framerate = ',s,12)); %find line in .prl file which gives frame rate.
       fr=str2num(s{line}(13:strfind(s{line},'per')-1)); %framerate
   end
   %statistic/how much did file change?
   temp=abs((c_data(1:208,:)-Ehub.raw)./Ehub.raw*100);% abs(percent diff)
   %[mean(temp(:)) trimmean(temp(:),10) std(temp(:)) prctile(temp(:),[2.5 97.5])]
   h1=msgbox([{strcat('File saved as ''', b, '_corrected.get''')};...
       {strcat('Percent change: ', num2str(mean(temp(:)),3) ,'(', num2str(std(temp(:)),3), '), mean(+/-std)')}],...
       'Sucess');
   end
%h=msgbox('Corrected data saved as *_corrected.get','Finished');
%end 

%function quality
threshold_p=5; %acceptable level of recipocity error (guestimate)
global Ehub
%Ehub.imdl.fwd_model.stimulation=mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current','rotate_meas','do_redundant'},1);%adj16_full_208 equivalent pattern
%    for i = 1:size(Ehub.imdl.fwd_model.stimulation,2) %for each drive combination
%        for j = 1:size(Ehub.imdl.fwd_model.stimulation,2)%for each (other) drive combination
%            [r,c,v]=find(Ehub.imdl.fwd_model.stimulation(j).meas_pattern*Ehub.imdl.fwd_model.stimulation(i).stim_pattern==-2); %which of the jth measurement patternd matches the drive pattern in i
%            [rm,cm,vm]=find(Ehub.imdl.fwd_model.stimulation(i).meas_pattern*Ehub.imdl.fwd_model.stimulation(j).stim_pattern==-2);%and which stimulation patten in i matches this (jth) drive pattern
%            Ehub.imdl.fwd_model.reci(size(Ehub.imdl.fwd_model.stimulation(i).meas_pattern,1)*(i-1)+rm,1)=(size(Ehub.imdl.fwd_model.stimulation(i).meas_pattern,1)*(j-1)+r); %form we we know i,rm has the reci pair j,r
%       end
%    end %matches reci_order in Get_reci 

Ehub.imdl.fwd_model.reci=[39 51 63 75 87 99 111 123 135 147 159 171 183 52 64 76 88 100 112 124 136 148 160 172 184 196 65 77 89 101 113 125 137 149 161 173 185 197 1 78 90 102 114 126 138 150 162 174 186 198 2 14 91 103 115 127 139 151 163 175 187 199 3 15 27 104 116 128 140 152 164 176 188 200 4 16 28 40 117 129 141 153 165 177 189 201 5 17 29 41 53 130   142   154   166   178     190   202     6    18    30    42    54    66   143   155   167   179   191   203     7    19    31    43    55    67    79   156   168   180     192   204     8    20    32    44    56    68    80    92   169   181   193   205     9    21    33    45    57    69    81    93   105   182     194   206    10    22    34    46    58    70    82    94   106   118   195   207    11    23    35    47    59    71    83    95   107   119     131   208    12    24    36    48    60    72    84    96   108   120   132   144    13    25    37    49    61    73    85    97   109   121     133   145   157    26    38    50    62    74    86    98   110   122   134   146   158   170]';

    %reci=(Ehub.raw-Ehub.raw(Ehub.imdl.fwd_model.reci,:,:))./(Ehub.raw+Ehub.raw(Ehub.imdl.fwd_model.reci,:,:))*200;
    reci=(c_data(1:208,:)-c_data(Ehub.imdl.fwd_model.reci,:,:))./(c_data(1:208,:)+c_data(Ehub.imdl.fwd_model.reci,:,:))*200; %use corrected data
Ehub.quality.reci_m=median(reci,2);
Ehub.quality.reci_std=std(reci,0,2);
reci_percent=sum(Ehub.quality.reci_m>threshold_p)/length(Ehub.imdl.fwd_model.reci)*200;%percent of electrode combinations with abs reci>threshold_p%
if max(reci_percent)>threshold_p %threshold may need changing
     i=find(reci_percent>threshold_p);
for j=1:length(i)
    f = questdlg([num2str(round(reci_percent(:,:,i(j)))) '% of electrode combinations have high reciprocity errors. Data may be unreliable.'], ['% electrode combinations with recipocity error >', num2str(threshold_p), '%'],'ignore', 'display', 'display');
if strmatch(f,'display')
    figure
    set(gcf,'position',[400 250 800 450])
    subplot(2,1,1)
    %plot(Ehub.auxiliary.time.data, reci(:,:,i(j)),'.'); xlabel 'time(s)'; ylabel 'recipocity(%)'
    %axis([floor(min(Ehub.auxiliary.time.data(:))) ceil(max(Ehub.auxiliary.time.data(:))) 0 ceil(max(reci(:)))]) %<0 will be a mirror image, so don't display
    plot([1:size(reci,2)]/fr,reci(:,:,i(j))','.'); xlabel 'time(s)'; ylabel 'recipocity(%)';axis tight
    a=axis;
    axis([a(1) a(2) 0 a(4)])
    %axis([floor(min(Ehub.auxiliary.time.data(:))) ceil(max(Ehub.auxiliary.time.data(:))) 0 ceil(max(reci(:)))]) %<0 will be a mirror image, so don't display
    %title(Ehub.filenames{1,1,i(j)})
    title(filenames)
    subplot(2,1,2)
    errorbar(abs(mean(reci(:,:,i(j))')),std(reci(:,:,i(j))'*2),'.')%if I only show the positive errors, I might not see the pattern of where they come from
    axis([1 length(Ehub.imdl.fwd_model.reci) 0 ceil(max(reci(:)))])
    set(gca,'XTick',6:13:208)
    set(gca,'XTickLabel',[' 1-2 ';' 2-3 ';' 3-4 ';' 4-5 ';' 5-6 ';' 6-7 ';' 7-8 ';' 8-9 ';' 9-10';'10-11';'11-12';'12-13';'13-14';'14-15';'15-16';'16-1 '])
    xlabel 'drive electrodes'; ylabel 'recipocity(%)'
    %plot(mean(Ehub.raw')) %should show 16 'u's, if adj 16.
end
end
       
end