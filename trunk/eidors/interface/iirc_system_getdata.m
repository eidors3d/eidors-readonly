function Scan_Data= iirc_system_getdata( Data )
% [sys_config, status]= iirc_system_configure( config_file )
% Read Data from the IIRC EIT system
%
% (C) 2006 tongin oh
% $Id$

Error_Flag=0;

for CurrentFreq = 1 : Data.EIT_Setting_Information.Frequency_Num
    Error_Flag=FreqSet(Data, Error_Flag, CurrentFreq);
    Error_Flag=SetMFCP(Data, Error_Flag, CurrentFreq);
    if Data.EIT_Setting_Information.FreqIndex(CurrentFreq) > 10
        [Scan_Data, ErrorFlag]= Scan_EIT_50(Data, Error_Flag, CurrentFreq);
    else
        [Scan_Data, ErrorFlag]= Scan_EIT_48(Data, Error_Flag, CurrentFreq);
    end
    
end


function Error_Flag=FreqSet(Data, Error_Flag, CurrentFreq)
    Dummy = 51;
    Send_Data = [32 Data.EIT_Setting_Information.FreqInfo(CurrentFreq) 220 Dummy Dummy Dummy Dummy];

    [USB_Status Temp_Data Count] = Multi_USB_Comm(Data.EIT_Setting_Information.Device_Sel, Send_Data);
    if USB_Status == 0
        error('USB Read Error');
    end    
    Receive_Data = Temp_Data(1 : Count);

    Response_Data = [161 Send_Data(2 : 7)];
    if Receive_Data == Response_Data'
    %     msgbox('Frequency Setting OK');
    else    
        error('Frequency Setting Error');
        Error_Flag=1;
    end

function Error_Flag= SetMFCP(Data, Error_Flag, CurrentFreq)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Frequency Setting for Scan(48)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dummy = 51;
    Test_OK = 1;

    Send_Data = [40 0 0 Data.EIT_Setting_Information.FreqInfo(CurrentFreq) Data.EIT_Setting_Information.GICCH(CurrentFreq) Dummy Dummy];

    [USB_Status Temp_Data Count] = Multi_USB_Comm(Data.EIT_Setting_Information.Device_Sel, Send_Data);
    if USB_Status == 0
        error('USB Read Error');
        Test_OK = 0;
    end    
    Receive_Data = Temp_Data(1 : Count);

    Response_Data = [169 Send_Data(2 : 7)];
    if Receive_Data ~= Response_Data'
        error('MFSet Error');
        Test_OK = 0;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CCS & GIC Setting for Scan(48)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Send_Data = [44 0 0 Data.CCS_CAL_Digipot_Value(Data.EIT_Setting_Information.FreqIndex(CurrentFreq), 1 : 4)];

    [USB_Status Temp_Data Count] = Multi_USB_Comm(Data.EIT_Setting_Information.Device_Sel, Send_Data);
    if USB_Status == 0
        error('USB Read Error');
        Test_OK = 0;
    end    
    Receive_Data = Temp_Data(1 : Count);

    Response_Data = [173 Send_Data(2 : 7)];
    if Receive_Data ~= Response_Data'
        error('MFSet Error');
        Test_OK = 0;
    end

    Send_Data = [44 0 1 Data.CCS_CAL_Digipot_Value(Data.EIT_Setting_Information.FreqIndex(CurrentFreq), 6 : 9)];

    [USB_Status Temp_Data Count] = Multi_USB_Comm(Data.EIT_Setting_Information.Device_Sel, Send_Data);
    if USB_Status == 0
        error('USB Read Error');
        Test_OK = 0;
    end    
    Receive_Data = Temp_Data(1 : Count);

    Response_Data = [173 Send_Data(2 : 7)];
    if Receive_Data(1) ~= Response_Data(1)
        msgbox('MFSet Error');
        Test_OK = 0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setting current pattern for Scan(48)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for k = 1 : Data.EIT_Setting_Information.Total_Projection(CurrentFreq)
        Send_Data = [36 0 (Data.EIT_Setting_Information.Total_Projection(CurrentFreq) - 1) (k - 1) ...
             Data.Channel_List(Data.EIT_Setting_Information.Current_Injection_Pattern((CurrentFreq), k, 1)) ...
             Data.Channel_List(Data.EIT_Setting_Information.Current_Injection_Pattern((CurrentFreq), k, 2)) ...
             Data.EIT_Setting_Information.Current_Injection_Pattern((CurrentFreq), k, 3)];

        [USB_Status Temp_Data Count] = Multi_USB_Comm(Data.EIT_Setting_Information.Device_Sel, Send_Data);
        if USB_Status == 0
            msgbox('USB Read Error');
            Test_OK = 0;
        end    
        Receive_Data = Temp_Data(1 : Count);

        Response_Data = [165 Send_Data(2 : 7)];
        if Receive_Data(1) ~= Response_Data(1)
            msgbox('SetCP Error');
            Test_OK = 0;
        end   
    end 


function [Scan_Data, Error_Flag]= Scan_EIT_50(Data, Error_Flag, CurrentFreq)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Complex scan using 16 channel EIT System (50)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Dummy = 51;
    Comm_Status = 1;
    Raw_Data = [];

    Send_Data = [50 1 Dummy Dummy Dummy Dummy Dummy];

    [USB_Status Temp_Data Count] = Multi_USB_Comm(Data.EIT_Setting_Information.Device_Sel, Send_Data);
    if USB_Status == 0
        msgbox('USB Read Error');
    end    
    Receive_Data = Temp_Data(1 : Count);
    if size(Receive_Data) ~= [896 1]
        Comm_Status = 0;
    end   
    Raw_Data = Receive_Data;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data read when 1 scan using 16 channel EIT System (54)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Send_Data = [54 1 Dummy Dummy Dummy Dummy Dummy];
    [USB_Status Temp_Data Count] = Multi_USB_Comm(Data.EIT_Setting_Information.Device_Sel, Send_Data);
    if USB_Status == 0
        msgbox('USB Read Error');
    end    
    Receive_Data = Temp_Data(1 : Count);
    if size(Receive_Data) ~= [896 1]
        Comm_OK = 0;
    end   
    Raw_Data = [Raw_Data; Receive_Data];

    Send_Data = [54 2 Dummy Dummy Dummy Dummy Dummy];
    [USB_Status Temp_Data Count] = Multi_USB_Comm(Data.EIT_Setting_Information.Device_Sel, Send_Data);
    if USB_Status == 0
        msgbox('USB Read Error');
    end    
    Receive_Data = Temp_Data(1 : Count);
    if size(Receive_Data) ~= [896 1]
        Comm_OK = 0;
    end   
    Raw_Data = [Raw_Data; Receive_Data];

    Send_Data = [54 3 Dummy Dummy Dummy Dummy Dummy];
    [USB_Status Temp_Data Count] = Multi_USB_Comm(Data.EIT_Setting_Information.Device_Sel, Send_Data);
    if USB_Status == 0
        msgbox('USB Read Error');
    end    
    Receive_Data = Temp_Data(1 : Count);
    if size(Receive_Data) ~= [896 1]
        Comm_OK = 0;
    end   
    Raw_Data = [Raw_Data; Receive_Data];

    switch (Data.EIT_Setting_Information.FreqIndex(CurrentFreq))
        case 11,
            Comp_Freq_Index = 11;
        case 12,
            Comp_Freq_Index = 14;
        case 13,
            Comp_Freq_Index = 17;
        otherwise
            msgbox('Out of Range !!!');
            Err_Flag = 1;
    end

    if Comm_Status == 1
        for Proj_Num = 1 : Data.EIT_Setting_Information.Total_Projection(CurrentFreq)
            for Data_Step = 1 : Data.EIT_Setting_Information.Total_Channel
                Projection_Data = Raw_Data((14 * Data_Step - 13 : 14 * Data_Step)+(Data.EIT_Setting_Information.Total_Projection(CurrentFreq)*14*(Proj_Num-1)));
                if Projection_Data(1) == ((Data_Step-1)*Channel_Gap + 192)
                    Real_Data_Temp = 256 * Projection_Data(3) + Projection_Data(4);
                    Quad_Data_Temp = 256 * Projection_Data(5) + Projection_Data(6);
                    if Real_Data_Temp >= 32768
                        Real_Data_Freq1 = Real_Data_Temp - 65536;
                    else
                        Real_Data_Freq1 = Real_Data_Temp;
                    end
                    if Quad_Data_Temp >= 32768
                        Quad_Data_Freq1 = Quad_Data_Temp - 65536;
                    else
                        Quad_Data_Freq1 = Quad_Data_Temp;
                    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calibration applied
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %             Theta = 0;
    %             Gain = 1;

                Theta = Data.VMCal(1, Comp_Freq_Index).Phase_Data(1, Data_Step);
                Gain = Data.VMCal(1, Comp_Freq_Index).Gain(1, Data.EIT_Setting_Information.Digipot).Magnitude(Proj_Num, Data_Step);
                Theta_Cal_Real = cos(Theta)*Real_Data_Freq1 + sin(Theta)*Quad_Data_Freq1;
                Theta_Cal_Quad = cos(Theta)*Quad_Data_Freq1 - sin(Theta)*Real_Data_Freq1;
                Amp_Cal_Real = Theta_Cal_Real / Gain;
                Amp_Cal_Quad = Theta_Cal_Quad / Gain;

    %             Scan_Data(Data.EIT_Operating_Information.Cur_Freq).Real_Data(Proj_Num, Data_Step) = Amp_Cal_Real;
    %             Scan_Data(Data.EIT_Operating_Information.Cur_Freq).Quad_Data(Proj_Num, Data_Step) = Amp_Cal_Quad;
    %             Scan_Data(Data.EIT_Operating_Information.Cur_Freq).Magnitude_Data(Proj_Num, Data_Step) = sqrt(Amp_Cal_Real^ 2 + Amp_Cal_Quad^ 2);
    %             Scan_Data(Data.EIT_Operating_Information.Cur_Freq).Theta_Data(Proj_Num, Data_Step) = rad2deg(angle(Amp_Cal_Real + i * Amp_Cal_Quad));

                Scan_Data(Comp_Freq_Index).Meas_Data((Proj_Num-1)*Data.EIT_Setting_Information.Total_Channel+Data_Step, 1) = Amp_Cal_Real;
                Scan_Data(Comp_Freq_Index).Meas_Data((Proj_Num-1)*Data.EIT_Setting_Information.Total_Channel+Data_Step, 2) = Amp_Cal_Quad;
                Scan_Data(Comp_Freq_Index).Meas_Data((Proj_Num-1)*Data.EIT_Setting_Information.Total_Channel+Data_Step, 3) = sqrt(Amp_Cal_Real^ 2 + Amp_Cal_Quad^ 2);
                Scan_Data(Comp_Freq_Index).Meas_Data((Proj_Num-1)*Data.EIT_Setting_Information.Total_Channel+Data_Step, 4) = rad2deg(angle(Amp_Cal_Real + i * Amp_Cal_Quad));

                Real_Data_Temp = 256 * Projection_Data(7) + Projection_Data(8);
                Quad_Data_Temp = 256 * Projection_Data(9) + Projection_Data(10);
                if Real_Data_Temp >= 32768
                    Real_Data_Freq2 = Real_Data_Temp - 65536;
                else
                    Real_Data_Freq2 = Real_Data_Temp;
                end
                if Quad_Data_Temp >= 32768
                    Quad_Data_Freq2 = Quad_Data_Temp - 65536;
                else
                    Quad_Data_Freq2 = Quad_Data_Temp;
                end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calibration applied
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %             Theta = 0;
    %             Gain = 1;

                Theta = Data.VMCal(1, (Comp_Freq_Index)+1).Phase_Data(1, Data_Step);
                Gain = Data.VMCal(1, (Comp_Freq_Index)+1).Gain(1, Data.EIT_Setting_Information.Digipot).Magnitude(Proj_Num, Data_Step);
                Theta_Cal_Real = cos(Theta)*Real_Data_Freq2 + sin(Theta)*Quad_Data_Freq2;
                Theta_Cal_Quad = cos(Theta)*Quad_Data_Freq2 - sin(Theta)*Real_Data_Freq2;
                Amp_Cal_Real = Theta_Cal_Real / Gain;
                Amp_Cal_Quad = Theta_Cal_Quad / Gain;

    %             Scan_Data(Data.EIT_Operating_Information.Cur_Freq+1).Real_Data(Proj_Num, Data_Step) = Amp_Cal_Real;
    %             Scan_Data(Data.EIT_Operating_Information.Cur_Freq+1).Quad_Data(Proj_Num, Data_Step) = Amp_Cal_Quad;
    %             Scan_Data(Data.EIT_Operating_Information.Cur_Freq+1).Magnitude_Data(Proj_Num, Data_Step) = sqrt(Amp_Cal_Real^ 2 + Amp_Cal_Quad^ 2);
    %             Scan_Data(Data.EIT_Operating_Information.Cur_Freq+1).Theta_Data(Proj_Num, Data_Step) = rad2deg(angle(Amp_Cal_Real + i * Amp_Cal_Quad));

                Scan_Data(Comp_Freq_Index+1).Meas_Data((Proj_Num-1)*Data.EIT_Setting_Information.Total_Channel+Data_Step, 1) = Amp_Cal_Real;
                Scan_Data(Comp_Freq_Index+1).Meas_Data((Proj_Num-1)*Data.EIT_Setting_Information.Total_Channel+Data_Step, 2) = Amp_Cal_Quad;
                Scan_Data(Comp_Freq_Index+1).Meas_Data((Proj_Num-1)*Data.EIT_Setting_Information.Total_Channel+Data_Step, 3) = sqrt(Amp_Cal_Real^ 2 + Amp_Cal_Quad^ 2);
                Scan_Data(Comp_Freq_Index+1).Meas_Data((Proj_Num-1)*Data.EIT_Setting_Information.Total_Channel+Data_Step, 4) = rad2deg(angle(Amp_Cal_Real + i * Amp_Cal_Quad));

                Real_Data_Temp = 256 * Projection_Data(11) + Projection_Data(12);
                Quad_Data_Temp = 256 * Projection_Data(13) + Projection_Data(14);
                if Real_Data_Temp >= 32768
                    Real_Data_Freq3 = Real_Data_Temp - 65536;
                else
                    Real_Data_Freq3 = Real_Data_Temp;
                end
                if Quad_Data_Temp >= 32768
                    Quad_Data_Freq3 = Quad_Data_Temp - 65536;
                else
                    Quad_Data_Freq3 = Quad_Data_Temp;
                end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calibration applied
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %             Theta = 0;
    %             Gain = 1;

                Theta = Data.VMCal(1, (Comp_Freq_Index)+2).Phase_Data(1, Data_Step);
                Gain = Data.VMCal(1, (Comp_Freq_Index)+2).Gain(1, Data.EIT_Setting_Information.Digipot).Magnitude(Proj_Num, Data_Step);
                Theta_Cal_Real = cos(Theta)*Real_Data_Freq3 + sin(Theta)*Quad_Data_Freq3;
                Theta_Cal_Quad = cos(Theta)*Quad_Data_Freq3 - sin(Theta)*Real_Data_Freq3;
                Amp_Cal_Real = Theta_Cal_Real / Gain;
                Amp_Cal_Quad = Theta_Cal_Quad / Gain;

                Scan_Data(Comp_Freq_Index+2).Meas_Data((Proj_Num-1)*Data.EIT_Setting_Information.Total_Channel+Data_Step, 1) = Amp_Cal_Real;
                Scan_Data(Comp_Freq_Index+2).Meas_Data((Proj_Num-1)*Data.EIT_Setting_Information.Total_Channel+Data_Step, 2) = Amp_Cal_Quad;
                Scan_Data(Comp_Freq_Index+2).Meas_Data((Proj_Num-1)*Data.EIT_Setting_Information.Total_Channel+Data_Step, 3) = sqrt(Amp_Cal_Real^ 2 + Amp_Cal_Quad^ 2);
                Scan_Data(Comp_Freq_Index+2).Meas_Data((Proj_Num-1)*Data.EIT_Setting_Information.Total_Channel+Data_Step, 4) = rad2deg(angle(Amp_Cal_Real + i * Amp_Cal_Quad));

                else
                    error('Scan Error');
                end
            end
        end
    end


function [Scan_Data, Error_Flag]= Scan_EIT_48(Data, Error_Flag, CurrentFreq)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1 scan using 16 channel EIT System (48)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data.EIT_Operating_Information.Cur_Freq = Data.EIT_Setting_Information.FreqIndex;
    % Data.EIT_Operating_Information.Cur_Gain = Data.EIT_Setting_Information.Digipot;


    Dummy = 51;
    Comm_Status = 1;
    Raw_Data = [];

    Send_Data = [48 1 Dummy Dummy Dummy Dummy Dummy];

    [USB_Status Temp_Data Count] = Multi_USB_Comm(Data.EIT_Setting_Information.Device_Sel, Send_Data);
    if USB_Status == 0
        msgbox('USB Read Error');
    end    
    Receive_Data = Temp_Data(1 : Count);
    if size(Receive_Data) ~= [768 1]
        Comm_Status = 0;
    end   
    Raw_Data = Receive_Data;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data read when 1 scan using 16 channel EIT System (52)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Send_Data = [52 1 Dummy Dummy Dummy Dummy Dummy];
    [USB_Status Temp_Data Count] = Multi_USB_Comm(Data.EIT_Setting_Information.Device_Sel, Send_Data);
    if USB_Status == 0
        msgbox('USB Read Error');
    end    
    Receive_Data = Temp_Data(1 : Count);
    if size(Receive_Data) ~= [768 1]
        Comm_OK = 0;
    end   
    Raw_Data = [Raw_Data; Receive_Data];

    Channel_Gap = 64 / Data.EIT_Setting_Information.Total_Channel;
    if Comm_Status == 1
        for Proj_Num = 1 : Data.EIT_Setting_Information.Total_Projection(CurrentFreq)
            for Data_Step = 1 : Data.EIT_Setting_Information.Total_Channel
                Projection_Data = Raw_Data((6 * Data_Step - 5 : 6 * Data_Step)+(Data.EIT_Setting_Information.Total_Projection(CurrentFreq)*6*(Proj_Num-1)));
                if Projection_Data(1) == ((Data_Step-1)*Channel_Gap + 128)
                    Real_Data_Temp = 256 * Projection_Data(3) + Projection_Data(4);
                    Quad_Data_Temp = 256 * Projection_Data(5) + Projection_Data(6);
                    if Real_Data_Temp >= 32768
                        Real_Data = Real_Data_Temp - 65536;
                    else
                        Real_Data = Real_Data_Temp;
                    end
                    if Quad_Data_Temp >= 32768
                        Quad_Data = Quad_Data_Temp - 65536;
                    else
                        Quad_Data = Quad_Data_Temp;
                    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calibration applied
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %             Theta = 0;
    %             Gain = 1;

                Theta = Data.VMCal(1, Data.EIT_Setting_Information.FreqIndex(CurrentFreq)).Phase_Data(Data_Step);
                Gain = Data.VMCal(1, Data.EIT_Setting_Information.FreqIndex(CurrentFreq)).Gain(1, Data.EIT_Setting_Information.Digipot).Magnitude(Proj_Num, Data_Step);

                Theta_Cal_Real = cos(Theta)*Real_Data + sin(Theta)*Quad_Data;
                Theta_Cal_Quad = cos(Theta)*Quad_Data - sin(Theta)*Real_Data;
                Amp_Cal_Real = Theta_Cal_Real / Gain;
                Amp_Cal_Quad = Theta_Cal_Quad / Gain;

                Scan_Data(Data.EIT_Setting_Information.FreqIndex(CurrentFreq)).Meas_Data((Proj_Num-1)*Data.EIT_Setting_Information.Total_Channel+Data_Step, 1) = Amp_Cal_Real;
                Scan_Data(Data.EIT_Setting_Information.FreqIndex(CurrentFreq)).Meas_Data((Proj_Num-1)*Data.EIT_Setting_Information.Total_Channel+Data_Step, 2) = Amp_Cal_Quad;
                Scan_Data(Data.EIT_Setting_Information.FreqIndex(CurrentFreq)).Meas_Data((Proj_Num-1)*Data.EIT_Setting_Information.Total_Channel+Data_Step, 3) = sqrt(Amp_Cal_Real^ 2 + Amp_Cal_Quad^ 2);
                Scan_Data(Data.EIT_Setting_Information.FreqIndex(CurrentFreq)).Meas_Data((Proj_Num-1)*Data.EIT_Setting_Information.Total_Channel+Data_Step, 4) = rad2deg(angle(Amp_Cal_Real + i * Amp_Cal_Quad));


                else
                    msgbox('Scan Error');
                end
            end
        end
    end
