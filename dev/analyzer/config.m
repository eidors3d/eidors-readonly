% Breathing analysis DO show_breaths
%DO manual_breath_select
DO model_breaths
DO show_volume_n_flow
DO stats_volume_n_flow
DO flow_volume_global
DO flow_volume_global_slope
DO TV_image
DO TV_slices
DO flow_volume_image
DO flow_volume_image_change
DO flow_volume_slices
DO flow_volume_components

% Perfusion/Pulsatility analysis 
%DO show_apnoea
%DO show_beats
%DO show_perfusion
%DO show_max_pixel

CONFIG rotate 180
CONFIG min_insp_length 0.2 % seconds
CONFIG min_expi_length 0.5 % seconds
CONFIG FRC_relative_match 0.2 % match start/end FRC
CONFIG FRC_search_window  0.3  % search 100ms for FRC
CONFIG slices 4
CONFIG LP_filter 4.0 % LP filter cut-off frequency [HZ]
CONFIG model_breaths_ncos 4
%CONFIG color_map ocean
%CONFIG color_map gray
 CONFIG color_map viridis
%CONFIG FILE "Example file name.mat" starttime 1 % [s]
%CONFIG FILE "Example file name.mat" endtime 1 % [s]
%CONFIG FILE "Example file name.mat" LPF_fc 1 % LPF cut off [Hz]
%CONFIG FILE "Example file name.mat" SKIPFILE
