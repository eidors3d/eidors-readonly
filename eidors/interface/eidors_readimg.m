function img = eidors_readimg( fname, format )
% EIDORS readimg - read reconstructed image files from various EIT equipment
%    manufacturers
%
% Currently the list of supported file formats is:
%    - MCEIT (Goettingen / Viasys) "igt" file format 
%        format = "IGT" or "MCEIT"
%
% Usage
% [vv, auxdata ]= eit_readdata( fname, format )
%     vv      = measurements - data frames in each column
%     auxdata = auxillary data - if provided by system 
%     fname = file name
%
%  if format is unspecified, we attempt to autodetect

% TODO:
% Formats to implement:
%    - Draeger "get" file format
%        format = "GET" or "draeger"
%    - Sheffield MK I "RAW" file format
%        format = "RAW" or "sheffield"
%    - ITS (International Tomography Systems)
%        format = "ITS" or "p2k"
%    - IIRC (Impedance Imaging Research Center, Korea)
%        format = "txt" or "IIRC"