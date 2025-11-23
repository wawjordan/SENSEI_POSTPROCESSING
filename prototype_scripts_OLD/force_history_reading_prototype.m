%% Force History File Reading (12/14/2024)
clc; clear; close all; fclose('all');

% filedir = 'C:\Users\wajordan\Documents\MATLAB\SENSEI_POSTPROCESSING\prototype scripts';
filedir = 'C:\Users\wajordan\Desktop\Desktop_TEMP\RESULTS\AIRFOIL\CASE_3_new_investigation\JOUKOWSKI_C_GRID_curved_2025-01-22_18.23.13';

% subdir = '_kt0129x0033';
subdir1 = '_kt1025x0257';
subdir2 = '_kt2049x0513';
subdir3 = '_kt1025x0257';
subdir4 = '_kt2049x0513';
% subdir4 = '_kt4097x1025';


% filename1 = 'kt-001-primal-force_history.dat';
% filename2 = 'kt-001-ic-force_history.dat';
% filename3 = 'kt-001-exact-force_history.dat';
% alpha = 5;
% beta  = 0;
% a = 1.1;
% c = 4 + 1/30;
% 
% gamma_nondim = 4*pi*a*sind(alpha+beta);
% 
% CX_exact = -(2/c)*gamma_nondim*sind(alpha);
% CY_exact =  (2/c)*gamma_nondim*cosd(alpha);
% 
% DATA1 = tec2mat_structured_condensed(fullfile(filedir,filename1));
% DATA2 = tec2mat_structured_condensed(fullfile(filedir,filename2));
% DATA3 = tec2mat_structured_condensed(fullfile(filedir,filename3));
% 
% 
% 
% CX = [DATA2(:).ZONE.CX];
% CY = [DATA2(:).ZONE.CY];
% iter = 0:length(CX)-1;
% 
% CX_pri = repmat(DATA1.ZONE.CX,length(CX),1);
% CY_pri = repmat(DATA1.ZONE.CY,length(CY),1);
% 
% CX_ex = repmat(DATA3.ZONE.CX,length(CX),1);
% CY_ex = repmat(DATA3.ZONE.CY,length(CY),1);
% 
% CX_ex1 = CX_ex*0 + CX_exact;
% CY_ex1 = CY_ex*0 + CY_exact;
% 
% CL = CY*cosd(alpha) - CX*sind(alpha);
% CD = CY*sind(alpha) + CX*cosd(alpha);
% 
% CL_pri = CY_pri*cosd(alpha) - CX_pri*sind(alpha);
% CD_pri = CY_pri*sind(alpha) + CX_pri*cosd(alpha);
% 
% CL_ex = CY_ex*cosd(alpha) - CX_ex*sind(alpha);
% CD_ex = CY_ex*sind(alpha) + CX_ex*cosd(alpha);
% 
% CL_ex1 = CY_ex1*cosd(alpha) - CX_ex1*sind(alpha);
% CD_ex1 = CY_ex1*sind(alpha) + CX_ex1*cosd(alpha);
% 
% CL_ex2 = CL_ex*0 + 8*pi*(1.1/(4 + 1/30))*sin(deg2rad(alpha));
% CD_ex2 = CD_ex*0;
% figure;
% hold on;

% plot(iter,CX_pri,'b')
% plot(iter,CX_ex,'k--')
% plot(iter,CX_ex1,'k')
% plot(iter,CX,'r')

% plot(iter,CY_pri,'b')
% plot(iter,CY_ex,'k--')
% plot(iter,CY_ex1,'k')
% plot(iter,CY,'r')

% plot(iter,CL_pri,'b')
% plot(iter,CL_ex,'k--')
% plot(iter,CL_ex1,'k')
% plot(iter,CL_ex2,'g--')
% plot(iter,CL,'r')

% plot(iter,CD_pri,'b')
% plot(iter,CD_ex,'k--')
% plot(iter,CD_ex1,'k')
% plot(iter,CD,'r')


filename1 = 'kt-primal-surface_forces.dat';
filename2 = 'kt-ic-surface_forces.dat';
filename3 = 'kt-exact-surface_forces.dat';

DATA1 = tec2mat_structured_condensed(fullfile(filedir,subdir2,filename1));
DATA2 = tec2mat_structured_condensed(fullfile(filedir,subdir2,filename2));
DATA3 = tec2mat_structured_condensed(fullfile(filedir,subdir2,filename3));

figure;
hold on;
for i = 1:numel(DATA1)
    plot(getfield(DATA1(i).ZONE,'XC'),getfield(DATA1(i).ZONE,'P'),'b','linewidth',2)
end

for i = 1:numel(DATA2)
    plot(getfield(DATA2(i).ZONE,'XC'),getfield(DATA2(i).ZONE,'P'),':')
end
i = 1;
plot(getfield(DATA2(i).ZONE,'XC'),getfield(DATA2(i).ZONE,'P'),'r','linewidth',2)
i = numel(DATA2);
plot(getfield(DATA2(i).ZONE,'XC'),getfield(DATA2(i).ZONE,'P'),'g','linewidth',2)


for i = 1:numel(DATA3)
    plot(getfield(DATA3(i).ZONE,'XC'),getfield(DATA3(i).ZONE,'P'),'k','linewidth',2)
end