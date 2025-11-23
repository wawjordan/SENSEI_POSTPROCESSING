%% Tec2Mat Work 05/04/2024
clc; clear; close all;

folder = 'C:\Users\Will\Documents\MATLAB\VT_Research\Clint_Flat_Plate\working_with_tec2mat\';

%% Example: SENSEI tecplot output (soln.dat, turb.dat, surface_forces.dat)
file_name1 = 'M5TW091_Med-soln.dat';
file_name2 = 'M5TW091_Med-turb.dat';
file_name3 = 'M5TW091_Med-surface_forces.dat';

ZONES_SENSEI1 = tec2mat_MOD(fullfile(folder,file_name1),'debug','safe');
zone_num1  = length(ZONES_SENSEI1);

ZONES_SENSEI2 = tec2mat_MOD(fullfile(folder,file_name2),'debug','safe');
zone_num2  = length(ZONES_SENSEI2);

ZONES_SENSEI3 = tec2mat_MOD(fullfile(folder,file_name3),'debug','safe');
zone_num3  = length(ZONES_SENSEI3);

%%
SENSEI_DATA1 = get_structured_from_sensei_tec2mat(ZONES_SENSEI1,zone_num1);
SENSEI_DATA2 = get_structured_from_sensei_tec2mat(ZONES_SENSEI2,zone_num2);
SENSEI_DATA3 = get_structured_from_sensei_tec2mat(ZONES_SENSEI3,zone_num3);

