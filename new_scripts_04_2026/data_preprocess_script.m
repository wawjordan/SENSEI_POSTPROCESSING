%% Parsing KT-airfoil data Preprocessing (05/14/2026)
clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parent_dir_str = 'SENSEI_POSTPROCESSING';
path_parts = regexp(mfilename('fullpath'), filesep, 'split');
path_idx = find(cellfun(@(s1)strcmp(s1,parent_dir_str),path_parts));
parent_dir = fullfile(path_parts{1:path_idx});
addpath(genpath(parent_dir));
clear parent_dir_str path_idx path_parts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;

inputs = struct();
inputs.epsilon = 0.1;
inputs.kappa   = 0.0;
inputs.tau     = 0.0;
inputs.vinf    = 75.0;
inputs.rhoinf  = 1.0;
inputs.pinf    = 100000.0;
inputs.gamma   = 1.4;

inputs.alpha   = 5; % (degrees)
inputs.nskip   = 4;
inputs.rho_ref = 1.0;
inputs.p_ref   = 100000.0;
inputs.a_ref   = sqrt(inputs.gamma*inputs.p_ref/inputs.rho_ref);

% nondimensionalize inputs
inputs.vinf   = inputs.vinf/inputs.a_ref;
inputs.rhoinf = inputs.rhoinf/inputs.rho_ref;
inputs.pinf   = inputs.pinf/(inputs.rho_ref*inputs.a_ref^2);

airfoil        = kt_airfoil( inputs.epsilon, inputs.kappa, inputs.tau );
airfoil.vinf   = inputs.vinf;
airfoil.rhoinf = inputs.rhoinf;
airfoil.pinf   = inputs.pinf;
airfoil        = airfoil.set_alpha(inputs.alpha);

DATA_DIR='C:\Users\wajordan\Desktop\CASES\';
foldernames1 = {};
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2026-01-15_02.00.10_CURVED_P4_q6_q6_hoex_TT_BC_VARLAY2_ITER2000_IC10'];% 1
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2026-01-15_03.10.03_CURVED_P4_q6_q6_hoex_TT_BC_VARLAY8_ITER2000_IC10'];% 2
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2026-02-13_12.01.06_ORDER_4_bc_IC_10_100_iter_GEO4'];% 3
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_limited_2026-02-13_12.00.09_ORDER_4_bc_IC_10_100_iter_GEO4'];% 4
foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_04_01_2026-04-02_12.20.09'];        % 5
foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_04_01_regress_2026-04-02_12.21.04'];% 6
foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_clustered_04_02_regress_2026-04-02_14.54.19'];% 7
foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_clustered_04_02_2026-04-02_14.54.15'];% 8
foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_clustered_04_02_regress_2026-04-02_19.10.44'];% 9
foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_clustered_04_02_2026-04-02_19.10.39'];% 10
foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_clustered_04_02_regress_2026-04-03_11.01.28'];% 11
foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_clustered_04_02_2026-04-03_10.47.17'];% 12
foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_v2_04_08_regress_2026-04-09_11.38.29']; % 13
foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_v2_04_08_2026-04-09_11.38.25'        ]; % 14
foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_v3_04_09_regress_2026-04-09_18.12.43'        ]; % 15
foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_v3_04_09_2026-04-09_18.12.39'        ]; % 16
foldernames1 = [foldernames1,  'JOUKOWSKI_C_GRID_curved_v3_04_09_regress_2026-04-10_11.40.43'      ]; % 17
foldernames1 = [foldernames1,  'JOUKOWSKI_C_GRID_curved_v3_04_09_2026-04-10_11.40.38'      ]; % 18
foldernames1 = [foldernames1,  'REGRESS_4_JOUKOWSKI_C_GRID_curved_v3_04_09_regress_2026-05-13_12.36.01'      ]; % 19
foldernames1 = [foldernames1,  'ITER_0_JOUKOWSKI_C_GRID_curved_v3_04_09_regress_2026-05-14_20.02.21'      ]; % 20
foldernames1 = cellfun(@(str_b)strcat(DATA_DIR,str_b),foldernames1,UniformOutput=false);

folders = foldernames1([17,18]);

ALL_DATA = struct();
ALL_DATA.inputs = inputs;
ALL_DATA.DATA = struct();
for i = 1:numel(folders)
    folder = folders{i};
    ALL_DATA.DATA(i).folder = folder;
    ALL_DATA.DATA(i).linear = get_airfoil_force_data_from_directory_alt(folder,inputs.alpha,inputs.nskip,airfoil,inputs.rho_ref,inputs.p_ref,inputs.a_ref,false);
    ALL_DATA.DATA(i).curved = get_airfoil_force_data_from_directory_alt(folder,inputs.alpha,inputs.nskip,airfoil,inputs.rho_ref,inputs.p_ref,inputs.a_ref,true);
end
out_dir   = 'C:\Users\wajordan\Desktop\git_MATLAB\SENSEI_POSTPROCESSING\new_scripts_04_2026\data_preprocess';
file_name = 'F_17_18.mat';
save(fullfile(out_dir,file_name),'ALL_DATA');