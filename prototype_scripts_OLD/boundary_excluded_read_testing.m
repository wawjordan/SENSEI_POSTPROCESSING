%% Prototype reading for files with boundary-excluded error norms (12/05/24)
clc; clear; close all;

file_name = 'C:\Users\wajordan\Desktop\Matlab_Scripts\POST_PROCESSING_NEW\kt-de_corrected_primitive_offset.dat';

[dat,N] = read_err_norm_file_with_ic_and_offsets(file_name);

function [dat,N] = read_err_norm_file_with_ic_and_offsets(file_name)
fid = fopen(file_name,'r');

line1 = fgetl(fid);
line2 = fgetl(fid);
line3 = fgetl(fid);


% Check if lines 1-3 contain a single entry
tmp1 = regexp(line1,'\d*','match');
tmp2 = regexp(line2,'\d*','match');
tmp3 = regexp(line3,'\d*','match');

if (length(tmp1) ~= 1)
    % assume no IC, no offsets
    fclose(fid);
    [dat,N] = read_err_norm_file(file_name);
    return;
elseif (length(tmp2) ~= 1)
    % assume no IC, but with offets
    fclose(fid);
    [dat,N] = read_err_norm_file_with_offsets(file_name);
    return
elseif (length(tmp3) ~= 1)
    % assume IC, but no offets
    fclose(fid);
    [dat,N] = read_err_norm_ic_file(file_name);
    return
end

% assume IC with offsets

% assume that this integer is the number of iterative corrections
n_ic = str2double(tmp1{1});

% assume that this integer is the number of excluded layers
n_layers = str2double(tmp2{1});

% the next line should be the size
line = fgetl(fid);
tmp = sscanf(line,'%d %f');
N = tmp(1);

% this should be the first line with error norms
line = fgetl(fid);
tmp = sscanf(line,'%f');
N_eqns = length(tmp);

% allocate
dat = zeros(N_eqns,3,n_ic+1,n_layers+1);

% rewind
frewind(fid);
% read the first two lines again
fgetl(fid);
fgetl(fid);
% now loop
for j = 0:n_ic
    fgetl(fid); % (iteration) not needed
    for k = 0:n_layers
        fgetl(fid); % (grid size) not needed anymore
        for i = 1:3 % norms (L1, L2, L_inf)
            line = fgetl(fid);
            tmp = sscanf(line,'%f');
            dat(:,i,j+1,k+1) = tmp;
        end
    end
end

fclose(fid);
end


function [dat,N] = read_err_norm_ic_file(file_name)
fid = fopen(file_name,'r');

line1 = fgetl(fid);
line2 = fgetl(fid);


% Check if lines 1-2 contain a single entry
tmp1 = regexp(line1,'\d*','match');
tmp2 = regexp(line2,'\d*','match');

if (length(tmp1) ~= 1)
    % assume no IC, no offsets
    fclose(fid);
    [dat,N] = read_err_norm_file(file_name);
    return;
elseif (length(tmp2) ~= 1)
    % assume no IC, but with offets
    fclose(fid);
    [dat,N] = read_err_norm_file_with_offsets(file_name);
    return
end

% assume IC, no offsets

% assume that this integer is the number of iterative corrections
n_ic = str2double(tmp1{1});

% the next line should be the size
line = fgetl(fid);
tmp = sscanf(line,'%d %f');
N = tmp(1);

% this should be the first line with error norms
line = fgetl(fid);
tmp = sscanf(line,'%f');
N_eqns = length(tmp);

% allocate
dat = zeros(N_eqns,3,n_ic+1);

% rewind
frewind(fid);
% read the first line again
fgetl(fid);
% now loop
for j = 0:n_ic
    fgetl(fid); % (iteration) not needed
    fgetl(fid); % (grid size) not needed anymore
    for i = 1:3 % norms (L1, L2, L_inf)
        line = fgetl(fid);
        tmp = sscanf(line,'%f');
        dat(:,i,j+1) = tmp;
    end
end
fclose(fid);
end



function [dat,N] = read_err_norm_file_with_offsets(file_name)
fid = fopen(file_name,'r');

line = fgetl(fid);

% Check if first line contains single entry
tmp = regexp(line,'\d*','match');

if (length(tmp) ~= 1)
    % try to read using the old method
    fclose(fid);
    [dat,N] = read_err_norm_file(file_name);
    return;
end

% assume that this integer is the number of excluded layers (+1)
n_layers = str2double(tmp{1});

line = fgetl(fid);
tmp = sscanf(line,'%d %f');
N = tmp(1);

line = fgetl(fid);
tmp = sscanf(line,'%f');
N_eqns = length(tmp);

% allocate
dat = zeros(N_eqns,3,n_layers+1);

% rewind
frewind(fid);
% read the first line again
fgetl(fid);
% now loop
for j = 0:n_layers
    fgetl(fid); % not needed anymore
    for i = 1:3 % remaining norms (L2, L_inf)
        line = fgetl(fid);
        tmp = sscanf(line,'%f');
        dat(:,i,j+1) = tmp;
    end
    
end

fclose(fid);
end

function [dat,N] = read_err_norm_file(file_name)
fid = fopen(file_name,'r');

line = fgetl(fid);
tmp = sscanf(line,'%d %f');
N = tmp(1);

line = fgetl(fid);
tmp = sscanf(line,'%f');
N_eqns = length(tmp);
dat = zeros(N_eqns,3);
dat(:,1) = tmp;

for i = 2:3
    line = fgetl(fid);
    tmp = sscanf(line,'%f');
    dat(:,i) = tmp;
end

fclose(fid);
end

% function [dat,N] = read_err_norm_file_with_offsets(file_name)
% fid = fopen(file_name,'r');
% 
% line = fgetl(fid);
% 
% % Check if first line contains single entry
% tmp = regexp(line,'\d*','match');
% 
% if (length(tmp) ~= 1)
%     % try to read using the old method
%     fclose(fid);
%     [dat,N] = read_err_norm_file(file_name);
%     return;
% else
%     % assume that this integer is the number of excluded layers (+1)
%     n_layers = str2double(tmp{1});
%     
%     line = fgetl(fid);
%     tmp = sscanf(line,'%d %f');
%     N = tmp(1);
%     
%     line = fgetl(fid);
%     tmp = sscanf(line,'%f');
%     N_eqns = length(tmp);
%     
%     % allocate
%     dat = zeros(N_eqns,3,n_layers+1);
%     
%     % rewind
%     frewind(fid);
%     % read the first line again
%     fgetl(fid);
%     % now loop
%     for j = 0:n_layers
%         fgetl(fid); % not needed anymore
%         for i = 1:3 % remaining norms (L2, L_inf)
%             line = fgetl(fid);
%             tmp = sscanf(line,'%f');
%             dat(:,i,j+1) = tmp;
%         end
%         
%     end
%     
% end
% 
% fclose(fid);
% end
% 
% function [dat,N] = read_err_norm_file(file_name)
% fid = fopen(file_name,'r');
% 
% line = fgetl(fid);
% tmp = sscanf(line,'%d %f');
% N = tmp(1);
% 
% line = fgetl(fid);
% tmp = sscanf(line,'%f');
% N_eqns = length(tmp);
% dat = zeros(N_eqns,3);
% dat(:,1) = tmp;
% 
% for i = 2:3
%     line = fgetl(fid);
%     tmp = sscanf(line,'%f');
%     dat(:,i) = tmp;
% end
% 
% fclose(fid);
% end