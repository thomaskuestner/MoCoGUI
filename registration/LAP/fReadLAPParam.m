function [ SOptions, dDeformRes ] = fReadLAPParam( sParamFile )
% read LAP parameter from txt file (RegGUI)
%
% input:
% sParamFile        file to be read
%
% output:
% SOptions          struct containing algorithm parametrization
% dDeformRes        to be interpolated image resolution (0=no interpolation)
%
% -------------------------------------------------------------------------
% (c) 2015: Thomas Kuestner
% -------------------------------------------------------------------------

try
    cParam = table2cell(readtable(sParamFile,'Delimiter','='));
catch
    fid = fopen(sParamFile,'r');
    tline = fgetl(fid);
    cParam = strsplit(tline,'=');
    while ischar(tline)
        tline = fgetl(fid);
        if(ischar(tline))
            cParam(end+1,:) = strsplit(tline,'=');
        end
    end
    fclose(fid);
    cParam = cParam(2:end,:);
end
for i=1:size(cParam,1)
    cVar = strtok(cParam{i,2}, '%');
    if(strcmp(cVar(1),'''')) % string 
        eval(sprintf('%s = %s;',cParam{i,1},cVar));
    elseif(strcmp(cVar(1),'[')) % empty
        eval(sprintf('%s = [];',cParam{i,1}));
    else % number
        eval(sprintf('%s = %f;',cParam{i,1},str2double(cVar)));
    end
end


end