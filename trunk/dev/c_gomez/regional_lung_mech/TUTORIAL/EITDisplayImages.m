function EITDisplayImages(eitimages,time,series,clim)
%EITDISPLAYIMAGES   Displays in a new figure the images for the specified
%series and times.
%
% series 'd'   difference
%        'a'   absolute
%        'td'  tidal difference
%        'mtd' average tidal difference
%        'fEIT' functional EIT
%        'C' compliance image
%        'Cmax' maximum compliance map          
%        'Pstar' PEEP at max compliance map
%        'collapse' collapse map
%        'OD' over-distension map
%        'collapseOD' dual lung state map
%
% Signature:
% function EITDisplayImages(eitimages,time,series,clim)
%
% Input:
% eitimages struct      scalar      .structure
%                                   .subject_id
%                                   .subject_type
%                                   .case_date
%                                   .file_name
%                                   .eit_device
%                                   .image_series_length
%                                   .frame_rate
%                                   .difference_image_series
%                                   .absolute_iamge_series
%                                   .image_mask
%
% (time)    double      1xN         vector of image indices
% (series)  char        1xN         [see description above]
% (clim)    double      1           top of color scale
%
% Output:
%
% Copyright C. Gomez-Laberge, November 2010.
% $Id: $

% Set error status to 'no errors present'
s = false;
% Save calling path
callpath = pwd;

calc_colours('defaults');

% Check arguments
switch nargin
    case 1
        % Set time to default value: all images
        time = 1:size(eitimages.difference_image_series,3);
        % Set series to default value: 'difference'
        series = 'd';
        clim=[];
    case 2
        % Set series to default value: 'difference'
        series = 'd';
        clim=[];
    case 3
        clim=[];
    case 4
        % Do nothing
    otherwise
        % Error
        s = true;
        help('EITDisplayImages')
        display('Error EITDisplayImages: invalid arguments');
        error('Error EITDisplayImages: aborting execution');
end

% The figure
hdl = figure;

% Determine imaging times
timelength = length(time);
if timelength > 100
    time = time(1:100);
end
timerange = false;
switch series
    case 'td'
        figname = [eitimages.subject_id ' | '...
            eitimages.case_date ' | '...
            eitimages.file_name ' | TIDAL IMAGES'];
    case 'mtd'
        figname = [eitimages.subject_id ' | '...
            eitimages.case_date ' | '...
            eitimages.file_name ' | AVERAGE TIDAL IMAGE'];
    case 'fEIT'
        figname = [eitimages.subject_id ' | '...
            eitimages.case_date ' | '...
            eitimages.file_name ' | fEIT IMAGE'];
    case 'lungROI'
        figname = [eitimages.subject_id ' | '...
            eitimages.case_date ' | '...
            eitimages.file_name ' | Lung ROI IMAGE'];
    case 'C'
        figname = [eitimages.subject_id ' | '...
            eitimages.case_date ' | '...
            eitimages.file_name ' | LUNG COMPLIANCE'];
    case 'Cmax'
        figname = [eitimages.subject_id ' | '...
            eitimages.case_date ' | MAXIMUM LUNG COMPLIANCE'];
    case 'Pstar'
        figname = [eitimages.subject_id ' | '...
            eitimages.case_date ' | PEEP AT MAX COMPLIANCE'];
    case 'collapse'
        figname = [eitimages.subject_id ' | '...
            eitimages.case_date ' | '...
            eitimages.file_name ' | ALVEOLAR COLLAPSE'];
    case 'OD'
        figname = [eitimages.subject_id ' | '...
            eitimages.case_date ' | '...
            eitimages.file_name ' | ALVEOLAR OVERDISTENSION'];
    case 'collapseOD'
        figname = [eitimages.subject_id ' | '...
            eitimages.case_date ' | '...
            eitimages.file_name ' | ALVEOLAR COLLAPSE & OVERDISTENSION'];        
    case 'test'
        figname = 'TEST MODE';
    otherwise
        if ~isempty(eitimages.frame_rate)
            seconds = (time-1)/eitimages.frame_rate;
            % Round to second decimal
            seconds = round(seconds*100)/100;
            timestart = num2str(seconds(1),'%.2f');
            if timelength>1
                timerange = true;
                timeend = num2str(seconds(end),'%.2f');
            end
            % Figure name
            if timerange
                figname = [eitimages.subject_id ' | '...
                    eitimages.case_date ' | '...
                    eitimages.file_name ' | '...
                    timestart '-' timeend ' s'];
            else
                figname = [eitimages.subject_id ' | '...
                    eitimages.case_date ' | '...
                    eitimages.file_name ' | '...
                    timestart ' s'];
            end
        else
            timestart = num2str(time(1));
            if timelength>1
                timerange = true;
                timeend = num2str(time(end));
            end
            % Figure name
            if timerange
                figname = [eitimages.subject_id ' | '...
                    eitimages.case_date ' | '...
                    eitimages.file_name ' | images '...
                    timestart '-' timeend];
            else
                figname = [eitimages.subject_id ' | '...
                    eitimages.case_date ' | '...
                    eitimages.file_name ' | image '...
                    timestart];
            end
        end
end
set(hdl,'Name',figname);

% Extract frames and reshape into images and display
switch series
    case 'a'
        images = eitimages.absolute_image_series(:,:,time);
        title('Absolute Image Sequence');
    case 'd'
        images = eitimages.difference_image_series(:,:,time);
        unmask = ~eitimages.image_mask;
        if isempty(clim)
            clim = ceil(max(images(:)));
        end
        for i=1:length(time)
            img = images(:,:,i);
            img(unmask) = NaN;
            images(:,:,i) = img;
        end
        calc_colours('backgnd',[.8 .8 .8]);
        calc_colours('ref_level',0);
        calc_colours('cmap_type','blue_yellow');
        calc_colours('clim',clim);
        show_slices(images);
        calc_colours(images,[],clim);
        title('Difference Impedance Images z(t) [% change from reference]');
    case 'td'
        if ~isfield(eitimages,'tidalimages')
            s = true;
            help('EITDisplayImages')
            display('Error EITDisplayImages: tidal images not found');
            error('Error EITDisplayImages: aborting execution');
        end
        images = eitimages.tidalimages(:,:,time);
        unmask = ~eitimages.image_mask;
        if isempty(clim)
            clim = ceil(max(images(:)));
        end
        for i=1:length(time)
            img = images(:,:,i);
            img(unmask) = NaN;
            images(:,:,i) = img;
        end
        calc_colours('ref_level',0);
        calc_colours('backgnd',[.8 .8 .8]);
        calc_colours('cmap_type','blue_yellow');
        calc_colours('clim',clim);
        show_slices(images);
        calc_colours(images,[],clim);
        title('Tidal Images \Deltaz = z_i - z_e [% change from reference]');
    case 'mtd'
        if ~isfield(eitimages,'meantidalimage')
            s = true;
            help('EITDisplayImages')
            display('Error EITDisplayImages: mean tidal image not found');
            error('Error EITDisplayImages: aborting execution');
        end
        img = eitimages.meantidalimage(:,:);
        unmask = ~eitimages.image_mask;
        if isempty(clim)
            clim = ceil(max(img(:)));
        end
        img(unmask) = NaN;
        calc_colours('ref_level',0);
        calc_colours('backgnd',[.8 .8 .8]);
        calc_colours('cmap_type','blue_yellow');
        calc_colours('clim',clim);
        show_slices(img);
        calc_colours(img,[],clim);
        title(['Average Tidal Image \DeltaZ (n=' ...
            num2str(size(eitimages.tidalindices,2))...
            ')  [% change from reference]']);
    case 'fEIT'
        if ~isfield(eitimages,'fEIT')
            s = true;
            help('EITDisplayImages')
            display('Error EITDisplayImages: compliance image not found');
            error('Error EITDisplayImages: aborting execution');
        end
        img = eitimages.fEIT;
        unmask = ~eitimages.image_mask;
        imax = max(img(:));
        clim = imax;
        img(unmask) = NaN;
        calc_colours('ref_level',0);
        calc_colours('backgnd',[.8 .8 .8]);
        calc_colours('image_field',[.6 .6 .6]);        
        calc_colours('image_field_val',-Inf);    %|||cgomez custom field
        calc_colours('image_field_idx',[]);%|||cgomez custom field
        calc_colours('cmap_type','greyscale');
        calc_colours('clim',clim/2);
        show_slices(img-imax/2);
        calc_colours(img,[],-clim); %set to negative clim for greyscale
        title('fEIT Image s_z [ SD(z) ]');
    case 'lungROI'
        if ~isfield(eitimages,'lungROI')
            s = true;
            help('EITDisplayImages')
            display('Error EITDisplayImages: lungROI image not found');
            error('Error EITDisplayImages: aborting execution');
        end
        img = double(eitimages.lungROI);
        unmask = ~eitimages.image_mask;
        fieldmark = 0;
        img(unmask) = NaN;
        image_field = (img==fieldmark);
        imgtemp = img;
        imgtemp(image_field)=[];
        imax = 1;
        if isempty(clim)
            clim  = 1;
        end
        calc_colours('ref_level',0);
        calc_colours('backgnd',[.8 .8 .8]);
        calc_colours('image_field',[.6 .6 .6]);      %new calc_colours field
        calc_colours('image_field_val',fieldmark);   %new calc_colours field
        calc_colours('image_field_idx',image_field); %new calc_colours field
        calc_colours('cmap_type','greyscale');
        calc_colours('clim',clim);
        show_slices(img-imax);
        calc_colours(img,[],-clim); %set to negative clim for greyscale
        title('Lung ROI Image');
    case 'C'
        if ~isfield(eitimages,'complianceimage')
            s = true;
            help('EITDisplayImages')
            display('Error EITDisplayImages: compliance image not found');
            error('Error EITDisplayImages: aborting execution');
        end
        img = eitimages.complianceimage;
        unmask = ~eitimages.image_mask;
        fieldmark = -Inf;
        img(unmask) = NaN;
        image_field = (img==fieldmark);
        imgtemp = img;
        imgtemp(image_field)=[];
        imax = max(imgtemp)-min(imgtemp);
        if isempty(clim)
            OrdOfMag = 10^floor(log10(imax));
            clim  = OrdOfMag * ceil( imax / OrdOfMag );
        end
        calc_colours('ref_level',0);
        calc_colours('backgnd',[.8 .8 .8]);
        calc_colours('image_field',[.6 .6 .6]);      %new calc_colours field
        calc_colours('image_field_val',fieldmark);   %new calc_colours field
        calc_colours('image_field_idx',image_field); %new calc_colours field
        calc_colours('cmap_type','copper');
        calc_colours('clim',clim/2);
        img = img-imax/2-min(imgtemp);
        img(image_field)=fieldmark;
        img(unmask) = NaN;
        show_slices(img);
        calc_colours(img,[],-clim); %set to negative clim for greyscale
        cc = eitimages.complianceimage;
        mm = ~isinf(cc);
        meanC = mean(cc(mm));
        line1 = ['Compliance at PEEP %2.0f cm H_2O\n']; 
        line2 = ['(colour bar units are DeltaZ/cm H_2O)'];
        titleline = [line1 line2];
        title(sprintf(titleline,eitimages.peep),...
            'HorizontalAlignment','center');
        
    case 'Cmax'
        if ~isfield(eitimages,'Cmax')
            s = true;
            help('EITDisplayImages')
            display('Error EITDisplayImages: Cmax image not found');
            error('Error EITDisplayImages: aborting execution');
        end
        img = eitimages.Cmax;
        unmask = ~eitimages.image_mask;
        fieldmark = -Inf;
        img(unmask) = NaN;
        image_field = (img==fieldmark);
        imgtemp = img;
        imgtemp(image_field)=[];
        imax = max(imgtemp)-min(imgtemp);
        if isempty(clim)
            OrdOfMag = 10^floor(log10(imax));
            clim  = OrdOfMag * ceil( imax / OrdOfMag );
        end
        calc_colours('ref_level',0);
        calc_colours('backgnd',[.8 .8 .8]);
        calc_colours('image_field',[.6 .6 .6]);      %new calc_colours field
        calc_colours('image_field_val',fieldmark);   %new calc_colours field
        calc_colours('image_field_idx',image_field); %new calc_colours field
        calc_colours('cmap_type','copper');
        calc_colours('clim',clim/2);
        img = img-imax/2-min(imgtemp);
        img(image_field)=fieldmark;
        img(unmask) = NaN;
        show_slices(img);
        calc_colours(img,[],-clim); %set to negative clim for greyscale
        title('Maximum Compliance During Protocol (% \DeltaZ/cm H_2O)');
    case 'Pstar'
        if ~isfield(eitimages,'Pstar')
            s = true;
            help('EITDisplayImages')
            display('Error EITDisplayImages: Pstar image not found');
            error('Error EITDisplayImages: aborting execution');
        end
        img = eitimages.Pstar;
        uniquevals = unique(img);
        uniquevals(1) = [];
        for i = 1:length(uniquevals)
            valmask = img==uniquevals(i);
            img(valmask)=i;
        end
            
        
        unmask = ~eitimages.image_mask;
        fieldmark = -Inf;
        img(unmask) = NaN;
        image_field = (img==fieldmark);
        imgtemp = img;
        imgtemp(image_field)=[];
        imax = max(imgtemp)-min(imgtemp);
        clim = imax;
        calc_colours('ref_level',0);
        calc_colours('backgnd',[.8 .8 .8]);
        calc_colours('image_field',[.6 .6 .6]);      %new calc_colours field
        calc_colours('image_field_val',fieldmark);   %new calc_colours field
        calc_colours('image_field_idx',image_field); %new calc_colours field
        calc_colours('cmap_type','copper');
        calc_colours('clim',clim/2);
        img = img-imax/2-min(imgtemp);
        img(image_field)=fieldmark;
        img(unmask) = NaN;
        show_slices(img);
        line1 = ['Pressure at Maximum Compliance During Protocol '...
            '(cm H_2O)\n[Dark to light pixels:'];
        line2=[];
        for i = 1:length(uniquevals)-1
            line2 = [line2 ' %d,'];
        end
        titleline = [line1 line2 ' %d PEEP (cm H_2O)]'];
        title(sprintf(titleline, uniquevals));
    case 'collapse'
        if ~isfield(eitimages,'collapse')
            s = true;
            help('EITDisplayImages')
            display('Error EITDisplayImages: collapse map not found');
            error('Error EITDisplayImages: aborting execution');
        end
        img = eitimages.collapse;
        unmask = ~eitimages.image_mask;
        fieldmark = -Inf;
        img(unmask) = NaN;
        image_field = (img==fieldmark);
        imgtemp = img;
        imgtemp(image_field)=[];
        imax = max(imgtemp)-min(imgtemp);
        clim = 1;
        calc_colours('ref_level',0);
        calc_colours('backgnd',[.8 .8 .8]);
        calc_colours('image_field',[.6 .6 .6]);      %new calc_colours field
        calc_colours('image_field_val',fieldmark);   %new calc_colours field
        calc_colours('image_field_idx',image_field); %new calc_colours field
        calc_colours('cmap_type','black_red');
        calc_colours('clim',clim/2);
        img = img-imax/2-min(imgtemp);
        img(image_field)=fieldmark;
        img(unmask) = NaN;
        show_slices(img);
        calc_colours(img,[],-clim); %set to negative clim for greyscale
        title(['Alveolar Collapse: ' num2str(eitimages.collapseSum) ...
            '% [PEEP = ' num2str(eitimages.peep) ' cmH_2O]']);
    case 'OD'
        if ~isfield(eitimages,'overdist')
            s = true;
            help('EITDisplayImages')
            display('Error EITDisplayImages: overdistension map not found');
            error('Error EITDisplayImages: aborting execution');
        end
        img = eitimages.overdist;
        unmask = ~eitimages.image_mask;
        fieldmark = -Inf;
        img(unmask) = NaN;
        image_field = (img==fieldmark);
        imgtemp = img;
        imgtemp(image_field)=[];
        imax = max(imgtemp)-min(imgtemp);
        clim = 1;
        calc_colours('ref_level',0);
        calc_colours('backgnd',[.8 .8 .8]);
        calc_colours('image_field',[.6 .6 .6]);      %new calc_colours field
        calc_colours('image_field_val',fieldmark);   %new calc_colours field
        calc_colours('image_field_idx',image_field); %new calc_colours field
        calc_colours('cmap_type','black_red');
        calc_colours('clim',clim/2);
        img = img-imax/2-min(imgtemp);
        img(image_field)=fieldmark;
        img(unmask) = NaN;
        show_slices(img);
        calc_colours(img,[],-clim); %set to negative clim for greyscale
        title(['Alveolar Overdistension: ' num2str(eitimages.overdistSum) ...
            '% [PEEP = ' num2str(eitimages.peep) ' cmH_2O]']);
    case 'collapseOD'
        if ~isfield(eitimages,'collapse') || ~isfield(eitimages,'overdist')
            s = true;
            help('EITDisplayImages')
            display('Error EITDisplayImages: state maps not found');
            error('Error EITDisplayImages: aborting execution');
        end
        imgC = eitimages.collapse;
        imgOD = eitimages.overdist;
        unmask = ~eitimages.image_mask;
        fieldmark = -Inf;
        imgC(unmask) = NaN;
        image_field = (imgC==fieldmark);
        img = imgOD-imgC;
        clim = 1;
        calc_colours('ref_level',0);
        calc_colours('backgnd',[.8 .8 .8]);
        calc_colours('image_field',[.6 .6 .6]);      %new calc_colours field
        calc_colours('image_field_val',fieldmark);   %new calc_colours field
        calc_colours('image_field_idx',image_field); %new calc_colours field
        calc_colours('cmap_type','blue_black_red');
        calc_colours('clim',clim); 
        img(image_field)=fieldmark;
        img(unmask) = NaN;
        img(1,1) = -1;
        show_slices(img);
        calc_colours(img,[],clim); %set to negative clim for greyscale
        line1=  ['Alveolar Collapse (red): %2.0f%%\n'];
        line2=  ['Alveolar Overdistension (blue): %2.0f%%\n'];
        line3=  ['PEEP: %d cm H2O'];
        titleline = [line1 line2 line3];
        title(sprintf(titleline, eitimages.collapseSum, ...
            eitimages.overdistSum, eitimages.peep),...
            'HorizontalAlignment','center');
    case 'test'
        % Space for testing new images
    otherwise
        % Error
        s = true;
        help('EITDisplayImages')
        display('Error EITDisplayImages: invalid series name');
        error('Error EITDisplayImages: aborting execution');
end

% End of function
end %function
