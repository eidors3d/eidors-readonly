function display_meas(fwd_model,disp_meas,offset,elec_pp);

% DISPLAY_MEAS: display measurements on mesh
% usage: display_meas(fwd_model,disp_meas,offset,elec_pp);
%
% Function that displays the electrodes utilised for each measurement.
% Each electrode is colour-coded depending on its operation
% 
% Current Drive - Red
% Current Sink - Black
% Positive Potential Measurement - Green
% Negative Potential Measurement - Yellow
%
% Pressing any key changes the display to the next measurement
%
% Advanced Options
% "disp_meas" - 'y' or 'n' - selects whether the measuring electrodes are
%                            displayed or not
% "offset" - integer (1 to number of electrodes per plane), default=0
%          - rotates the display by number of electrodes selected
% "elec_pp" - number of electrodes per plane
%
% (C) 2005 by Stephen Murphy. Licensed under GPL version 2.
% $Id: display_meas.m,v 1.1 2005-12-06 14:53:07 aadler Exp $

% old parameters
%function display_meas(vtx,srf,elec,Ib,indH,df,disp_meas,offset,elec_pp);

pp= np_fwd_parameters( fwd_model);
vtx  = pp.vtx;
srf  = pp.srf;
elec = pp.elec;
Ib   = pp.Ib;
df   = pp.df;
indH = pp.indH;

figure
set(gcf,'Name','Wire Mesh Current Injections')

if nargin<2
    disp_meas='n';
end
if nargin<4
    elec_pp= [];
end
    


if disp_meas=='n'

    for loop1=1:size(Ib,2)
    
        trimesh(srf,vtx(:,1),vtx(:,2),vtx(:,3));
        axis image;
        xlabel('x')
        ylabel('y')
        zlabel('z')
        set(gcf,'Colormap',[0 0 0]);
        
        if ~isempty( elec_pp)
            set(gca,'CameraViewAngleMode','Manual');
            camera_angle=360/elec_pp;
            camorbit(camera_angle*offset,0);
        end
            
        hidden off;
        title(['Current Injection ' num2str(loop1)])

        drive_elec=[];
        sink_elec=[];
    
        for loop2=1:size(Ib,1)
            if Ib(loop2,loop1) ~=0 & isreal(sqrt(Ib(loop2,loop1)))
                % Driving electrode
                for loop3=0:(size(elec,2)/3)-1
                
                    if ( elec(loop2,1+(loop3*3)) ~= 0 & elec(loop2,2+(loop3*3)) ~=0 & elec(loop2,3+(loop3*3)) ~=0 )
                        X = [vtx(elec(loop2,1+(loop3*3)),1);vtx(elec(loop2,2+(loop3*3)),1);vtx(elec(loop2,3+(loop3*3)),1)];
                        Y = [vtx(elec(loop2,1+(loop3*3)),2);vtx(elec(loop2,2+(loop3*3)),2);vtx(elec(loop2,3+(loop3*3)),2)];
                        Z = [vtx(elec(loop2,1+(loop3*3)),3);vtx(elec(loop2,2+(loop3*3)),3);vtx(elec(loop2,3+(loop3*3)),3)];
                    end
                    
                    patch(X,Y,Z,'r');
                    drive_elec=[drive_elec;X';Y';Z'];
                end
            elseif Ib(loop2,loop1) ~=0
                %sinking electrode
                for loop3=0:(size(elec,2)/3)-1
            
                    if ( elec(loop2,1+(loop3*3)) ~= 0 & elec(loop2,2+(loop3*3)) ~=0 & elec(loop2,3+(loop3*3)) ~=0 )
                        X = [vtx(elec(loop2,1+(loop3*3)),1);vtx(elec(loop2,2+(loop3*3)),1);vtx(elec(loop2,3+(loop3*3)),1)];
                        Y = [vtx(elec(loop2,1+(loop3*3)),2);vtx(elec(loop2,2+(loop3*3)),2);vtx(elec(loop2,3+(loop3*3)),2)];
                        Z = [vtx(elec(loop2,1+(loop3*3)),3);vtx(elec(loop2,2+(loop3*3)),3);vtx(elec(loop2,3+(loop3*3)),3)];
                    end
                
                    patch(X,Y,Z,'k');
                    sink_elec=[sink_elec;X';Y';Z'];
                end
            end
        end     
    
        pause

    end
    close;
    
elseif disp_meas=='y'
    meas_no=1;
    for loop1=1:size(Ib,2)
    
        while meas_no<=sum(df(1:loop1))
            trimesh(srf,vtx(:,1),vtx(:,2),vtx(:,3));
            axis image;
            xlabel('x')
            ylabel('y')
            zlabel('z')
            set(gcf,'Colormap',[0 0 0]);
            
            if ~isempty(elec_pp)
                set(gca,'CameraViewAngleMode','Manual');
                camera_angle=360/elec_pp;
                camorbit(camera_angle*offset,0);
            end
            
            hidden off;
            title(['Measurement Number ' num2str(meas_no)])

            drive_elec=[];
            sink_elec=[];
        
            for loop2=1:size(Ib,1)
                if Ib(loop2,loop1) ~=0 & isreal(sqrt(Ib(loop2,loop1)))
                    % Driving electrode
                    for loop3=0:(size(elec,2)/3)-1
                    
                        if ( elec(loop2,1+(loop3*3)) ~= 0 & elec(loop2,2+(loop3*3)) ~=0 & elec(loop2,3+(loop3*3)) ~=0 )
                            X = [vtx(elec(loop2,1+(loop3*3)),1);vtx(elec(loop2,2+(loop3*3)),1);vtx(elec(loop2,3+(loop3*3)),1)];
                            Y = [vtx(elec(loop2,1+(loop3*3)),2);vtx(elec(loop2,2+(loop3*3)),2);vtx(elec(loop2,3+(loop3*3)),2)];
                            Z = [vtx(elec(loop2,1+(loop3*3)),3);vtx(elec(loop2,2+(loop3*3)),3);vtx(elec(loop2,3+(loop3*3)),3)];
                        end
                        
                        patch(X,Y,Z,'r');
                        %drive_elec=[drive_elec;X';Y';Z'];
                    end
                elseif Ib(loop2,loop1) ~=0
                    %sinking electrode
                    for loop3=0:(size(elec,2)/3)-1
            
                        if ( elec(loop2,1+(loop3*3)) ~= 0 & elec(loop2,2+(loop3*3)) ~=0 & elec(loop2,3+(loop3*3)) ~=0 )
                            X = [vtx(elec(loop2,1+(loop3*3)),1);vtx(elec(loop2,2+(loop3*3)),1);vtx(elec(loop2,3+(loop3*3)),1)];
                            Y = [vtx(elec(loop2,1+(loop3*3)),2);vtx(elec(loop2,2+(loop3*3)),2);vtx(elec(loop2,3+(loop3*3)),2)];
                            Z = [vtx(elec(loop2,1+(loop3*3)),3);vtx(elec(loop2,2+(loop3*3)),3);vtx(elec(loop2,3+(loop3*3)),3)];
                        end
                
                        patch(X,Y,Z,'k');
                        %sink_elec=[sink_elec;X';Y';Z'];
                    end
                end
            end     
            
            % Patch measuring electrodes
            for loop3=0:(size(elec,2)/3)-1
                if ( elec(indH(meas_no,1),1+(loop3*3)) ~= 0 & elec(indH(meas_no,1),2+(loop3*3)) ~=0 & elec(indH(meas_no,1),3+(loop3*3)) ~=0 )
                    X = [vtx(elec(indH(meas_no,1),1+(loop3*3)),1);vtx(elec(indH(meas_no,1),2+(loop3*3)),1);vtx(elec(indH(meas_no,1),3+(loop3*3)),1)];
                    Y = [vtx(elec(indH(meas_no,1),1+(loop3*3)),2);vtx(elec(indH(meas_no,1),2+(loop3*3)),2);vtx(elec(indH(meas_no,1),3+(loop3*3)),2)];
                    Z = [vtx(elec(indH(meas_no,1),1+(loop3*3)),3);vtx(elec(indH(meas_no,1),2+(loop3*3)),3);vtx(elec(indH(meas_no,1),3+(loop3*3)),3)];
                end
                
                patch(X,Y,Z,'g');
            end
            
            for loop3=0:(size(elec,2)/3)-1
                if ( elec(indH(meas_no,2),1+(loop3*3)) ~= 0 & elec(indH(meas_no,2),2+(loop3*3)) ~=0 & elec(indH(meas_no,2),3+(loop3*3)) ~=0 )
                    X = [vtx(elec(indH(meas_no,2),1+(loop3*3)),1);vtx(elec(indH(meas_no,2),2+(loop3*3)),1);vtx(elec(indH(meas_no,2),3+(loop3*3)),1)];
                    Y = [vtx(elec(indH(meas_no,2),1+(loop3*3)),2);vtx(elec(indH(meas_no,2),2+(loop3*3)),2);vtx(elec(indH(meas_no,2),3+(loop3*3)),2)];
                    Z = [vtx(elec(indH(meas_no,2),1+(loop3*3)),3);vtx(elec(indH(meas_no,2),2+(loop3*3)),3);vtx(elec(indH(meas_no,2),3+(loop3*3)),3)];
                end
                
                patch(X,Y,Z,'y');
            end
            meas_no=meas_no+1;
            pause

        end
    end
    close;
end
    
