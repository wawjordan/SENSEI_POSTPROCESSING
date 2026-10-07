function [A,lo,hi] = read_array_file( file_name )

fid = fopen(file_name,'r');

line = fgetl(fid);
% match numbers, or groups of numbers separated by colons
vals = regexp(line,'[+-]?[\d:]+','match');

% number of dimensions in the array
ndims = numel(vals);

lo   = ones(1,ndims);
hi   = ones(1,ndims);
dims = ones(1,ndims);
flag = false;
for i = 1:ndims
    tmp1 = str2double(vals{i});
    if isnan(tmp1)
        % assume this is a range
        tmp2 = regexp(vals{i},'[+-]?\d*','match');
        n = numel(tmp2);
        if (n==2)
            tmp3 = cellfun(@str2double,tmp2);
            lo(i) = tmp3(1);
            hi(i) = tmp3(2);
        else
            warning(['unrecognized range format in size descriptor, ',...
                     'outputing as 1D array']);
            flag = true;
        end
    else
        lo(i) = 1;
        hi(i) = str2double(vals{i});
    end
    dims(i) = hi(i) - lo(i) + 1;
end
A = fscanf(fid,'%f');
fclose(fid);
if ( flag )
    lo = 1;
    hi = numel(A);
else
    A = reshape(A,[dims,1]);
end

end

% function A = read_array_file( file_name )
% 
% fid = fopen(file_name,'r');
% 
% 
% line = fgetl(fid);
% dim = regexp(line,'\d*','match');
% dim = cellfun(@str2num,dim);
% A = fscanf(fid,'%f');
% 
% A = reshape(A,[dim,1]);
% fclose(fid);
% 
% end