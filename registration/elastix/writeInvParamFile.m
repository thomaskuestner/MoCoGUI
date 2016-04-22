function writeInvParamFile(paramPath)
% write inverse parameter file for elastix (RegGUI)
%
% input:
% paramPath         parameter file
% sText             to be replaced string
% sWhatReplaced     search string
%
% -------------------------------------------------------------------------
% (c) 2015: Thomas Kuestner, Verena Neumann
% -------------------------------------------------------------------------

%% changes + flags
property_1 = 'HowToCombineTransforms';
value_1 = 'Compose';
flag_1 = 0;                 % shows, if property was found
property_2 = 'Metric';
value_2 = 'DisplacementMagnitudePenalty';
flag_2 = 0;                 % shows, if property was found

%% copy of the original parameterfile
[pathstr,name,ext] = fileparts(paramPath);
invParam = [pathstr,filesep,name,'_inv',ext];
copyfile(paramPath,invParam);

%% initialization
fid = fopen(invParam);
i = 1;                      % row number
tline = fgetl(fid);         % read line from file, removing newline characters
A{i} = tline;

%% read and edit
while ischar(tline)
   i = i+1;
   tline = fgetl(fid);
   A{i} = tline;
   
   % loop for the first change
   if (ischar(A{i}) && ~isempty(regexp(A{i},['(',property_1],'once')) && isempty(regexp(A{i},'//','once')))
      A{i} = sprintf('(%s "%s")',property_1,value_1);
      flag_1 = ~flag_1;
   end
   
    % loop for the second change
   if (ischar(A{i}) && ~isempty(regexp(A{i},['(',property_2],'once')) && isempty(regexp(A{i},'//','once')))
      A{i} = sprintf('(%s "%s")',property_2,value_2);
      flag_2 = ~flag_2;
   end
end

%% control of flags
% When properties not found, they have to be included. Last entry in A hast
% to be -1 to mark the end of the file.

% for property 1
if ~flag_1
    A{length(A)+1} = A{length(A)};
    A{length(A)-1} = sprintf('(%s "%s")',property_1,value_1);
else
end

% for property 2
if ~flag_2
    A{length(A)+1} = A{length(A)};
    A{length(A)-1} = sprintf('(%s "%s")',property_2,value_2);
else
end

%% writing text file
fid = fopen(invParam, 'w');
for i = 1:numel(A)
     if A{i+1} == -1                % if next row the last one
        fprintf(fid,'%s', A{i});    
        break;                      
     else
        fprintf(fid,'%s\r\n', A{i});  % \r important for txt file!
     end
end
fclose(fid);

