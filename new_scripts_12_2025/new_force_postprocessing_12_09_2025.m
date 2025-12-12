%% Parsing KT-airfoil data (12/09/2025)
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
inputs.nskip   = 2;
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

% DATA_DIR='C:\Users\wajordan\Desktop\';
% foldernames1 = {             'ALPHA_0_JOUKOWSKI_C_GRID_curved_2025-11-23_13.28.07_K_EXACT'};% 1
% foldernames1 = [foldernames1,'ALPHA_0_JOUKOWSKI_C_GRID_curved_2025-11-23_12.47.49_VAR_REC_8_100'];% 2
% 
% folder = foldernames1{17};
% S = get_soln_data_from_directory(folder);
% G = get_grid_data_from_directory(folder);
% DATA = get_airfoil_force_data_from_directory_alt(folder,inputs.alpha,inputs.nskip,airfoil,inputs.rho_ref,inputs.p_ref,inputs.a_ref,false);
% N1 = [DATA1.F(:).N];
% 
% var = 'CD';
% err_CL_primal_1 = abs([DATA1.H(:).(['primal_',var])]-airfoil.(var));
% err_CL_primal_2 = abs([DATA2.H(:).(['primal_',var])]-airfoil.(var));
% 
% err_CL_ete_1 = abs([DATA1.H(:).(['ete_',var])]-airfoil.(var));
% err_CL_ete_2 = abs([DATA2.H(:).(['ete_',var])]-airfoil.(var));