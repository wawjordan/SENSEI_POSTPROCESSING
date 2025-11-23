function DATA2 = make_named_struct(DATA)
N = numel(DATA);

inputs = cell(1,2*N);

cnt = 0;
for i = 1:N
    cnt = cnt + 1;
    inputs{cnt} = DATA(i).name;
    cnt = cnt + 1;
    inputs{cnt} = DATA(i).data;
end
    
DATA2 = struct(inputs{:});
end