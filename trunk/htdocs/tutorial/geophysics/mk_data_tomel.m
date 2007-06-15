function data= mk_data_tomel(data_tomel,name,configuration)
% Convert ERT data from "tomel" to "eidors" format
%
% data_tomel = TOMEL data are given as a matrix with at least 8 colums as,
% for example:
%
%  EZG04 Ring 1 December 2004 ERT data (WEN32_1 set)
%  No    Time      A      B     M      N    mA     Volts
%  Data_Ring1_Dec2004_Wen32_1= [ ...
%   1     223    101    104    102    103   100    0.260163 ...
%
% 160    2285    132    115    105    110   100    0.041110];
% 'Data_Ring1_Dec2004_Wen32_1'
%
% where:
%   col 1 = measurement number
%   col 2 = time elapsed since begining of measurements (in seconds)
%   col 3 = A electrode ( = positive current electrode)
%   col 4 = B electrode ( = negative current electrode)
%   col 5 = M electrode ( = positive potential electrode)
%   col 6 = N electrode ( = negative potential electrode)
%   col 7 = intensity (mA) of the current injected from A to B
%   col 8 = voltage (V) measured V(M) - V(N) for the measurement pattern
%
% name = string containing the name you want to give to the data structure
%
% configuration = string giving a description of the measurement pattern used
%
% EXAMPLE:
%    real_data= mk_data_tomel(data_tomel,'EZG_Ring2_JULY2004','Wenner protocol');
%
% Dominique Gibert, April 2007
%
    data.name= name;
    data.type= 'data';
    data.time= -1;
    data.configuration= configuration;
    data.meas= data_tomel(:,8);
    return;