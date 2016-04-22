function sFoundtext = fReadText(fFilename, sWhatToLookFor)
% read in string from external file (RegGUI)
%
% input:
% fFilename         file to be read
% sWhatToLookFor    search string
%
% output:
% sFoundtext        found search string with corresponding parametrization
%
% -------------------------------------------------------------------------
% (c) 2015: Thomas Kuestner, Markus Bednarzyk
% -------------------------------------------------------------------------
fid = fopen(fFilename,'r');
i = 1;                      % row number
tline = fgetl(fid);         % Read line from file, removing newline characters
while ischar(tline)    
    if(ischar(tline) && ~isempty(regexp(tline,[sprintf('^(%s',sWhatToLookFor),'\w*'],'once')))
        sFoundtext = tline;
        break;
    end
    i = i+1;
    tline = fgetl(fid);
end
fclose(fid);
end