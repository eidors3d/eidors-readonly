function EITDisplayMovie(eitimages,time,clim,createavi)
%EITDISPLAYMOVIE   Displays in a new figure a ventilation movie for the 
%specified times.
%
% function EITDisplayMovie(eitimages,time,clim,createavi)
%
% Input:
% eitimages struct                   .structure
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
% (time)    double      1xN         vector of image indices\
% (clim)    double      1           top of color scale
% (createavi) bool      1           write movie to avi file in pwd
%
% Output:
%
% Copyright C. Gomez-Laberge, January 2011.
% $Id: $

% Set error status to 'no errors present'
s = false;
% Save calling path
callpath = pwd;

% Check arguments
switch nargin
    case 1
        % Set time to default value: all images
        time = 1:size(eitimages.difference_image_series,3);
        clim = [];
        createavi = false;
    case 2
        clim = [];
        createavi = false;
    case 3
        createavi = false;
    case 4
        % Do nothing
    otherwise
        % Error
        s = true;
        help('EITDisplayMovie')
        display('Error EITDisplayMovie: invalid arguments');
        error('Error EITDisplayMovie: aborting execution');
end

% Determine imaging times
hdl = figure;
timelength = length(time);
timerange = false;
if ~isempty(eitimages.frame_rate)
    frame_rate = eitimages.frame_rate;
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
    frame_rate = 10;
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
set(hdl,'Name',figname);

% Extract frames and reshape into images and display
images = eitimages.difference_image_series(:,:,time);
unmask = ~eitimages.image_mask;
for i=1:length(time)
    image = images(:,:,i);
    image(unmask) = NaN;
    images(:,:,i) = image;
end
if isempty(clim)
    clim = ceil(max(images(:)));
end
calc_colours('backgnd',[.8 .8 .8]);
calc_colours('ref_level',0);
calc_colours('cmap_type','cg_custom');
calc_colours('clim',clim);
for i = 1:length(time)
    show_slices(images(:,:,time(i)));
    F(i) = getframe;
end
% movie(F,1,frame_rate);
close(hdl);
if createavi
    movie2avi(F,[eitimages.subject_id '_mov.avi'],...
        'fps',frame_rate);
end

% End of function
end %function