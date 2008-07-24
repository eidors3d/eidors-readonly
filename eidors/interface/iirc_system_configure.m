function [Data,ErrorFlag]= iirc_system_configure( config_file )
% [sys_config, status]= iirc_system_configure( config_file )
% Configure the IIRC EIT system
%
% (C) 2006 tongin oh
% $Id: iirc_system_configure.m,v 1.1 2006-08-25 08:09:18 tonginoh Exp $

ErrorFlag= 0;
Data= struct;

Data= USBGetNum(Data);
Data= File_Load(config_file, Data);
[Receive_Data, Response_Data] = SysInit(Data);
Data= NVMC( Data );
AGC(Data);
AVG(Data);
% FreqSet
% SetMFCP
% Scan_EIT_48
% Save_ScanData
% 
% Image_EIT(16, 1, Data.EIT_Setting_Information.FreqIndex)

function Data= USBGetNum(Data)

    [USB_Status, Device_Num]=USB_GetNum();
    if USB_Status == 1
        D1 = 'USB Connection : ';
        D2 = int2str(Device_Num);
        MsgTemp = strcat(D1, D2);
        %     msgbox(MsgTemp);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Device select
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if Device_Num == 1
            Data.EIT_Setting_Information.Device_Sel = 0;
        elseif Device_Num == 2
            Data.EIT_Setting_Information.Device_Sel = 1;
        else
            error('Connected too many USB system !!!');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    else
        error('USB No Connection');
    end

function Data= File_Load(File_Name, Data)
    fid = fopen(File_Name, 'r');
% end
    Instruction = fgetl(fid);
    if strcmp(upper(Instruction), upper('File_Start')) == 1
        File_On = 1;
    end    

    while File_On == 1
        Instruction = fgetl(fid);
        if strcmp(upper(Instruction), upper('Total_Ch')) == 1
            Data.EIT_Setting_Information.Total_Channel = fscanf(fid, '%d', 1);
        elseif strcmp(upper(Instruction), upper('Frequency')) == 1
            Data.EIT_Setting_Information.Frequency_Num = fscanf(fid, '%d', 1);
            for j = 1 : Data.EIT_Setting_Information.Frequency_Num
                Temp_Freq_Num = fscanf(fid, '%d', 1);
                Temp_Freq_Unit = fscanf(fid, '%s', 1);
                if strcmp(upper(Temp_Freq_Unit), upper('kHz')) == 1
                    Data.EIT_Setting_Information.Frequency(j) = Temp_Freq_Num * 1000;
                elseif strcmp(upper(Temp_Freq_Unit), upper('Hz')) == 1
                    Data.EIT_Setting_Information.Frequency(j) = Temp_Freq_Num;
                else
                    msgbox('Out of Frequency Range !!!');
                end
                
                switch Data.EIT_Setting_Information.Frequency(j)
                    case 10,
                        Data.EIT_Setting_Information.FreqInfo(j)=0;
                        Data.EIT_Setting_Information.GICCH(j)=0;
                        Data.EIT_Setting_Information.FreqIndex(j)=1;
                    case 50,
                        Data.EIT_Setting_Information.FreqInfo(j)=1;
                        Data.EIT_Setting_Information.GICCH(j)=0;
                        Data.EIT_Setting_Information.FreqIndex(j)=2;
                    case 100,
                        Data.EIT_Setting_Information.FreqInfo(j)=2;
                        Data.EIT_Setting_Information.GICCH(j)=0;
                        Data.EIT_Setting_Information.FreqIndex(j)=3;
                    case 1000,
                        Data.EIT_Setting_Information.FreqInfo(j)=3;
                        Data.EIT_Setting_Information.GICCH(j)=0;
                        Data.EIT_Setting_Information.FreqIndex(j)=4;
                    case 5000,
                        Data.EIT_Setting_Information.FreqInfo(j)=4;
                        Data.EIT_Setting_Information.GICCH(j)=1;
                        Data.EIT_Setting_Information.FreqIndex(j)=5;
                    case 10000,
                        Data.EIT_Setting_Information.FreqInfo(j)=5;
                        Data.EIT_Setting_Information.GICCH(j)=2;
                        Data.EIT_Setting_Information.FreqIndex(j)=6;
                    case 50000,
                        Data.EIT_Setting_Information.FreqInfo(j)=6;
                        Data.EIT_Setting_Information.GICCH(j)=4;
                        Data.EIT_Setting_Information.FreqIndex(j)=7;
                    case 250000,
                        Data.EIT_Setting_Information.FreqInfo(j)=7;
                        Data.EIT_Setting_Information.GICCH(j)=32;
                        Data.EIT_Setting_Information.FreqIndex(j)=8;
                    case 500000,
                        Data.EIT_Setting_Information.FreqInfo(j)=8;
                        Data.EIT_Setting_Information.GICCH(j)=16;
                        Data.EIT_Setting_Information.FreqIndex(j)=9;
                    case 100000,
                        Data.EIT_Setting_Information.FreqInfo(j)=9;
                        Data.EIT_Setting_Information.GICCH(j)=8;
                        Data.EIT_Setting_Information.FreqIndex(j)=10;
                    otherwise,
                        msgbox('Undefined Frequency Range !!!');
                end
                ComplexReadFlag = 1;
                while (ComplexReadFlag == 1)
                    ComplexOnOff = fgetl(fid);
                    if strcmp(upper(ComplexOnOff), upper('Complex_On')) == 1
                        switch Data.EIT_Setting_Information.FreqInfo(j)
                            case 0,
                                Data.EIT_Setting_Information.FreqInfo(j)= 128;
                                Data.EIT_Setting_Information.GICCH(j)=0;
                                Data.EIT_Setting_Information.FreqIndex(j)= 11;
                            case 3,
                                Data.EIT_Setting_Information.FreqInfo(j)= 131;
                                Data.EIT_Setting_Information.GICCH(j)=1;
                                Data.EIT_Setting_Information.FreqIndex(j)= 12;
                            case 6,
                                Data.EIT_Setting_Information.FreqInfo(j)= 134;
                                Data.EIT_Setting_Information.GICCH(j)=32;
                                Data.EIT_Setting_Information.FreqIndex(j)= 13;
                            otherwise
                                msgbox('Undefined Frequency Range !!!');
                        end
                        ComplexReadFlag = 0;
                    elseif strcmp(upper(ComplexOnOff), upper('Complex_Off')) == 1
                        ComplexReadFlag = 0;
                    end
                end
            end
%         elseif strcmp(upper(Instruction), upper('Complex_On')) == 1
%                 Data.EIT_Setting_Information.Complex = 1;
%         elseif strcmp(upper(Instruction), upper('Complex_Off')) == 1
%                 Data.EIT_Setting_Information.Complex = 0;
        elseif strcmp(upper(Instruction), upper('Set_CP')) == 1
            for j = 1 : Data.EIT_Setting_Information.Frequency_Num
                Data.EIT_Setting_Information.Total_Projection(j) = fscanf(fid, '%d', 1);
                for k = 1 : Data.EIT_Setting_Information.Total_Projection(j)
                    for l = 1 : 3
                        Data.EIT_Setting_Information.Current_Injection_Pattern(j, k, l) = fscanf(fid, '%d', 1);
                    end
                end
            end
        elseif strcmp(upper(Instruction), upper('AGC_ON')) == 1
            Data.EIT_Setting_Information.AGC_On_Off = 1;
        elseif strcmp(upper(Instruction), upper('AGC_OFF')) == 1        
            Data.EIT_Setting_Information.AGC_On_Off = 0;
            Data.EIT_Setting_Information.Digipot = fscanf(fid, '%d', 1);
        elseif strcmp(upper(Instruction), upper('AVG_ON')) == 1
            Data.EIT_Setting_Information.AVG_On_Off = 1;
            Data.EIT_Setting_Information.AVG_Num = fscanf(fid, '%d', 1);
        elseif strcmp(upper(Instruction), upper('AVG_OFF')) == 1        
            Data.EIT_Setting_Information.AVG_On_Off = 0;
        elseif strcmp(upper(Instruction), upper('Wait_Time_ON')) == 1
            Data.EIT_Setting_Information.Wait_Time_ON = 1;
            Data.EIT_Setting_Information.Wait_Time1 = fscanf(fid, '%d', 1);
            Data.EIT_Setting_Information.Wait_Time2 = fscanf(fid, '%d', 1);
            Data.EIT_Setting_Information.Wait_Time3 = fscanf(fid, '%d', 1);
        elseif strcmp(upper(Instruction), upper('ADC_Error_Check')) == 1
            Data.EIT_Setting_Information.ADCErrCHK = fscanf(fid, '%d', 1);
        elseif strcmp(upper(Instruction), upper('Calibration_File')) == 1   
            Cal_File = fgetl(fid);
        elseif strcmp(upper(Instruction), upper('Neighboring Injection')) == 1   
            Data.EIT_Setting_Information.Injection_Method = 0;
        elseif strcmp(upper(Instruction), upper('2D Diagonal Injection')) == 1   
            Data.EIT_Setting_Information.Injection_Method = 1;
        elseif strcmp(upper(Instruction), upper('3D Diagonal Injection')) == 1   
            Data.EIT_Setting_Information.Injection_Method = 2;
        elseif strcmp(upper(Instruction), upper('File_End')) == 1
            File_On = 0;  
            fclose(fid);
        end        
    end  

    load VM_Neighboring_Cal.mat
    Data.VMCal = VMCal;

    load CCS_CAL_Digipot_Value.mat
    Data.CCS_CAL_Digipot_Value = CCS_CAL_Digipot_Value;

    load Digipot_Value.mat
    Data.Digipot_Value = Digipot_Value;

function [Receive_Data, Response_Data] = SysInit(Data)
    Dummy=51;
    Send_Data = [0 Data.EIT_Setting_Information.Total_Channel Dummy Dummy Dummy Dummy Dummy];
    [USB_Status Temp_Data Count] = Multi_USB_Comm(Data.EIT_Setting_Information.Device_Sel, Send_Data);
    if USB_Status == 0
        error('USB Read Error');
    end    
    Receive_Data = Temp_Data(1 : Count);
    Response_Data = [129 Send_Data(2 : 7)];
    if Receive_Data == Response_Data'
    %     msgbox('System Initial OK');   
    else
        error('System Initial Error');
        Error_Flag=1;
    end

function Data= NVMC(Data)
    Dummy=51;
    Send_Data = [4 Dummy Dummy Dummy Dummy Dummy Dummy];

    [USB_Status Temp_Data Count] = Multi_USB_Comm(Data.EIT_Setting_Information.Device_Sel, Send_Data);
    if USB_Status == 0
        error('USB Read Error');
    end    
    Receive_Data = Temp_Data(1 : Count);

    if Receive_Data(2) == Data.EIT_Setting_Information.Total_Channel
        Channel_Gap = 64 / Data.EIT_Setting_Information.Total_Channel;
        Data.Channel_List = 0 : Channel_Gap : Channel_Gap * (Data.EIT_Setting_Information.Total_Channel - 1);
        Data.Channel_List_Reverse = zeros(64, 1);
        for i = 1 : Data.EIT_Setting_Information.Total_Channel
            Data.Channel_List_Reverse(Data.Channel_List(i) + 1) = i;
        end
    %     msgbox(['NVMC ' num2str(Receive_Data(2)) ' OK']);
    else
        error(['NVMC Error(' num2str(Receive_Data(2)) ')']);
        Error_Flag=1;
    end

function AGC( Data );
    Dummy = 51;
        
    if Data.EIT_Setting_Information.AGC_On_Off == 1
        Send_Data = [12 1 Dummy Dummy Dummy Dummy Dummy];
    else
        Send_Data(1 : 2) = [12 0];
        Send_Data(3 : 4) = Data.Digipot_Value(Data.EIT_Setting_Information.Digipot, :);
        Send_Data(5 : 7) = [255 Dummy Dummy];
    end    

    [USB_Status Temp_Data Count] = Multi_USB_Comm(Data.EIT_Setting_Information.Device_Sel, Send_Data);
    if USB_Status == 0
        error('USB Read Error');
    end    
    Receive_Data = Temp_Data(1 : Count);

    Response_Data = [141 Send_Data(2 : 7)];
    if Receive_Data == Response_Data'
    %     msgbox('AGC OK');
    else    
        error('AGC Error');
        Error_Flag=1;
    end

function AVG( Data )
    Dummy = 51;
    if Data.EIT_Setting_Information.AVG_On_Off == 1
        Send_Data = [16 (Data.EIT_Setting_Information.AVG_Num - 1) Dummy Dummy Dummy Dummy Dummy];
    else    
        Send_Data = [16 0 Dummy Dummy Dummy Dummy Dummy];
    end

    [USB_Status Temp_Data Count] = Multi_USB_Comm(Data.EIT_Setting_Information.Device_Sel, Send_Data);
    if USB_Status == 0
        msgbox('USB Read Error');
    end    
    Receive_Data = Temp_Data(1 : Count);

    Response_Data = [145 Send_Data(2 : 7)];
    if Receive_Data == Response_Data'
    %     msgbox('AVG OK');
    else    
        msgbox('AVG Error');
        Error_Flag=1;
    end

