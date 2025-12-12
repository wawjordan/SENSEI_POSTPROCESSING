function write_tecplot_zones_from_DATA(filename,S)
fid = fopen(filename,'w');

write_tecplot_header_from_metadata(fid,S(1).meta_data,true)
write_tecplot_zone_from_DATA(fid,S(1).ZONE,{'XC','YC','ZC','VOL'})
for i = 2:numel(S)
    write_tecplot_header_from_metadata(fid,S(i).meta_data,false)
    write_tecplot_zone_from_DATA(fid,S(i).ZONE,{'XC','YC','ZC','VOL'})
end

fclose(fid);
end


function write_tecplot_zone_from_DATA(fid,ZONE_DATA,exclude_vars)

names = fieldnames(ZONE_DATA);

for v = 1:numel(names)
    if ~any(strcmp(exclude_vars,names{v}))
        fprintf(fid,'%.16E\n', ZONE_DATA.(names{v})(:) );
    end
end
end


function write_tecplot_header_from_metadata(fid,metadata,write_vars)

% VARIABLES = "x", "y", "<greek>r</greek>", "P", "u", "v"
% 
% ZONE
% T="block1"
% ZONETYPE=ORDERED
% I=257
% J=65
% DATAPACKING=BLOCK
% VARLOCATION=([1-2]=NODAL,[3-6]=CELLCENTERED)
% DT=(DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE)
% STRANDID=1
% SOLUTIONTIME=  6.000000000000000E+01
% VARSHARELIST=([1,2]=1)
n_vars = numel(metadata.ZONEVARlist);
if ( write_vars)
    fprintf(fid,['VARIABLES = "%s"',repmat('," %s"',1,n_vars-1),'\n'], metadata.ZONEVARlist{:});
end

if (isfield(metadata,"T"))
    if ~isempty(metadata.T)
        fprintf(fid,'ZONE T="%s"\n',metadata.T);
    else
        fprintf(fid,'ZONE\n');
    end
end

if ( isfield(metadata,"I") )
    fprintf(fid,'I=%d\n',metadata.I);
end
if ( isfield(metadata,"J") )
    fprintf(fid,'J=%d\n',metadata.J);
end
if ( isfield(metadata,"K") )
    fprintf(fid,'K=%d\n',metadata.K);
end

if (isfield(metadata,"ZONETYPE"))
    if ~isempty(metadata.ZONETYPE)
        fprintf(fid,'ZONETYPE=%s\n',metadata.ZONETYPE);
    end
end

if (isfield(metadata,"DATAPACKING"))
    if ~isempty(metadata.DATAPACKING)
        fprintf(fid,'DATAPACKING=%s\n',metadata.DATAPACKING);
    end
end

if (isfield(metadata,"VARLOCATION"))
    if ~isempty(metadata.VARLOCATION)
        fprintf(fid,'VARLOCATION=%s\n',metadata.VARLOCATION);
    end
end

if (isfield(metadata,"DT"))
    if ~isempty(metadata.DT)
        fprintf(fid,'DT=%s\n',metadata.DT);
    end
end

if (isfield(metadata,"STRANDID"))
    if ~isempty(metadata.STRANDID)
        fprintf(fid,'STRANDID=%d\n',metadata.STRANDID);
    end
end

if (isfield(metadata,"SOLUTIONTIME"))
    if ~isempty(metadata.SOLUTIONTIME)
        fprintf(fid,'SOLUTIONTIME=%s\n',metadata.SOLUTIONTIME);
    end
end

if (isfield(metadata,"VARSHARELIST"))
    if ~isempty(metadata.VARSHARELIST)
        fprintf(fid,'VARSHARELIST=(');
        for i = 1:size(metadata.VARSHARELIST,1)-1
            vars  = metadata.VARSHARELIST{i,1};
            zone  = metadata.VARSHARELIST{i,2};
            fprintf(fid,['[%d',repmat(',%d',1,numel(vars)-1),']=%d,'],vars,zone);
        end
        vars  = metadata.VARSHARELIST{end,1};
        zone  = metadata.VARSHARELIST{end,2};
        fprintf(fid,['[%d',repmat(',%d',1,numel(vars)-1),']=%d)\n'],vars,zone);
    end
end

end