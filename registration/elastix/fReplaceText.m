function fReplaceText(fFilename, sText, sWhatReplaced)
% replace text (parametrization) for a search string in an external file (RegGUI)
%
% input:
% fFilename         file to be read
% sText             to be replaced string
% sWhatReplaced     search string
%
% -------------------------------------------------------------------------
% (c) 2015: Thomas Kuestner
% -------------------------------------------------------------------------
if(nargin < 3)
    sWhatReplaced = 'InitialTransformParametersFileName';
end
fid = fopen(fFilename,'r'); 
i = 1;                      % row number
tline = fgetl(fid);         % Read line from file, removing newline characters
A{i} = tline;               % A(1) represents first line 
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
    if(ischar(A{i}) && ~isempty(regexp(A{i},[sprintf('^(%s',sWhatReplaced),'\w*'],'once')))
        A{i} = sprintf('(%s "%s")',sWhatReplaced,sText);
    end
end
fclose(fid);

fid = fopen(fFilename, 'w');
for i = 1:numel(A)
    if A{i+1} == -1                 % if next line is empty
        fprintf(fid,'%s', A{i});
        break;
    else
        fprintf(fid,'%s\n', A{i});
    end
end
fclose(fid);
end