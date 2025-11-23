%% Force History File Reading (12/17/2024)
clc; clear; close all; fclose('all');

% file_dir1 = 'C:\Users\wajordan\Desktop\Desktop_TEMP\RESULTS\AIRFOIL\CASE_3\JOUKOWSKI_C_GRID_curved_2024-12-16_12.35.28';

file_dir1 = 'C:\Users\wajordan\Desktop\Desktop_TEMP\RESULTS\AIRFOIL\CASE_3_new_investigation\JOUKOWSKI_C_GRID_curved_2025-01-22_18.23.13';
file_dir2 = 'C:\Users\wajordan\Desktop\Desktop_TEMP\RESULTS\AIRFOIL\CASE_3_new_investigation\JOUKOWSKI_C_GRID_curved_2_2025-01-22_18.23.17';

exact_file  = 'kt-001-exact-force_history.dat';
primal_file = 'kt-001-primal-force_history.dat';
ic_file     = 'kt-001-ic-force_history.dat';

% sub_dir = '_kt0129x0033';
sub_dir1 = '_kt1025x0257';
sub_dir2 = '_kt2049x0513';
sub_dir3 = '_kt1025x0257';
sub_dir4 = '_kt2049x0513';
% sub_dir4 = '_kt4097x1025';

alpha = 5;
[CL_ex1,CL_pri1,CL_ic1,~,~,~] = get_data_from_sub_dir(file_dir1,sub_dir1,exact_file,primal_file,ic_file,alpha);
[CL_ex2,CL_pri2,CL_ic2,~,~,~] = get_data_from_sub_dir(file_dir1,sub_dir2,exact_file,primal_file,ic_file,alpha);
[CL_ex3,CL_pri3,CL_ic3,~,~,~] = get_data_from_sub_dir(file_dir2,sub_dir3,exact_file,primal_file,ic_file,alpha);
[CL_ex4,CL_pri4,CL_ic4,~,~,~] = get_data_from_sub_dir(file_dir2,sub_dir4,exact_file,primal_file,ic_file,alpha);

CL_ex = 8*pi*(1.1/(4 + 1/30))*sind(alpha) + 0*CL_ic1;
% CL_ex = 0*CL_ic1;

hold on
plot(CL_ic1,'b-s')
plot(CL_ic2,'b--s')
plot(CL_ic3,'b:o')
plot(CL_ic4,'b-.o')
plot(CL_pri1,'r-s')
plot(CL_pri2,'r--s')
plot(CL_pri3,'r:o')
plot(CL_pri4,'r-.o')
plot(CL_ex1,'k-s')
plot(CL_ex2,'k--s')
plot(CL_ex3,'k:o')
plot(CL_ex4,'k-.o')
plot(CL_ex,'g')

function [CL_ex,CL_pri,CL_ic,CD_ex,CD_pri,CD_ic] = get_data_from_sub_dir(file_dir,sub_dir,exact_file,primal_file,ic_file,alpha)

ic = tec2mat_structured_condensed(fullfile(file_dir,sub_dir,ic_file));
CX_ic = getfield(ic.ZONE,'CX');
CY_ic = getfield(ic.ZONE,'CY');
CL_ic = CY_ic*cosd(alpha) - CX_ic*sind(alpha);
CD_ic = CY_ic*sind(alpha) + CX_ic*cosd(alpha);

ex = tec2mat_structured_condensed(fullfile(file_dir,sub_dir,exact_file));
CX_ex = getfield(ex.ZONE,'CX');
CY_ex = getfield(ex.ZONE,'CY');
CL_ex = CY_ex*cosd(alpha) - CX_ex*sind(alpha) + 0*CL_ic;
CD_ex = CY_ex*sind(alpha) + CX_ex*cosd(alpha) + 0*CD_ic;

pri = tec2mat_structured_condensed(fullfile(file_dir,sub_dir,primal_file));
CX_pri = getfield(pri.ZONE,'CX');
CY_pri = getfield(pri.ZONE,'CY');
CL_pri = CY_pri*cosd(alpha) - CX_pri*sind(alpha) + 0*CL_ic;
CD_pri = CY_pri*sind(alpha) + CX_pri*cosd(alpha) + 0*CD_ic;


end