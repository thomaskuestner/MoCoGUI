function varargout = LandmarkGUI(varargin)
% GUI to set landmarks in reference and registered images to evaluate the
% registration results
%
% opens from EvalGUI with the current registration results or from the
% command window, then, registration mat-file has to be loaded
%
%
% Before the first start check the path in th LandmarkGUI_OpeningFcn
%
% -------------------------------------------------------------------------
% (c) 2015: Thomas Kuestner, Verena Neumann
% -------------------------------------------------------------------------

% Last Modified by GUIDE v2.5 22-Apr-2016 09:52:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @LandmarkGUI_OpeningFcn, ...
    'gui_OutputFcn',  @LandmarkGUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before LandmarkGUI is made visible.
function LandmarkGUI_OpeningFcn(hObject, ~, handles, varargin)
%% standard paths
currpath = fileparts(mfilename('fullpath'));
if(~exist('GUIPreferences.mat','file')) % first run
    SPaths.sResults = uigetdir(pwd,'Select default result directory');
    SPaths.sData = uigetdir(pwd,'Select default data directory');
    SPaths.sCode = currpath;
    standardVoxelsize = [1,1,1];
    save('GUIPreferences.mat','SPaths','standardVoxelsize');
else   
    load GUIPreferences.mat; 
    try 
        % check if path are reachable
        if(~exist(SPaths.sResults,'dir')), error('path not existing'); end;
        if(~exist(SPaths.sData,'dir')), error('path not existing'); end;
        if(~exist(SPaths.sCode,'dir')), error('path not existing'); end;
    catch
        % loading failed or paths not reachable
        clear 'SPaths' 'standardVoxelsize';
    end
end
if(~exist('SPaths','var'))
    % if no valid GUIPreferences are set use the standard paths:   
    if ~exist([currpath,filesep,'Results'],'dir'); mkdir(currpath,'results'); end
    SPaths.sResults = [currpath,filesep,'results'];    
    SPaths.sData = [currpath,filesep,'example'];
    SPaths.sCode = currpath;
end
handles.SPaths = SPaths;

% set some folders on the path
addpath(genpath([currpath,filesep,'io']));
addpath(genpath([currpath,filesep,'metrics']));
addpath(genpath([currpath,filesep,'registration']));
addpath(genpath([currpath,filesep,'segmentation']));
addpath(genpath([currpath,filesep,'utils']));

%% read inputs
if(nargin > 3)
    for index = 1:2:(nargin-3),
        if nargin-3==index, break, end
        switch lower(varargin{index})
            case 'inarg'
                handles.h = varargin{index+1};
        end
    end
end
% when opening without data, select a data struct from the hard drive
% originating from registration with RegiGUI
if ~isfield(handles, 'h') || isempty(handles.h)
    handles.h = fSelectStruct(handles.SPaths.sResults);
    if ~handles.h
        handles.closeFigure = true;
        % Update handles structure
        guidata(hObject, handles);
        return;
    end
end

%% set icons
try % Try to apply a nice icon
    warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
    jframe = get(hObject, 'javaframe');
    jIcon = javax.swing.ImageIcon([currpath, filesep, 'icons', filesep, 'LandmarkGUI_icon.png']);
    pause(0.001);
    jframe.setFigureIcon(jIcon);
    clear jframe jIcon
catch
    warning('Could not apply a nice icon to the figure :(');
end
dImage = double(imread([currpath,filesep,'icons',filesep,'save.png']))./255;
if size(dImage, 3) == 1, dImage = repmat(dImage, [1 1 3]); end
set(handles.pb_save, 'CData', dImage/max(dImage(:)));

dImage = double(imread([currpath,filesep,'icons',filesep,'open.png']))./255;
if size(dImage, 3) == 1, dImage = repmat(dImage, [1 1 3]); end
set(handles.pb_load, 'CData', dImage/max(dImage(:)));

dImage = double(imread([currpath,filesep,'icons',filesep,'reset.png']))./255;
if size(dImage, 3) == 1, dImage = repmat(dImage, [1 1 3]); end
set(handles.pb_reset_contrast, 'CData', dImage/max(dImage(:)));

dImage = double(imread([currpath,filesep,'icons',filesep,'cogs.png']))./255;
if size(dImage, 3) == 1, dImage = repmat(dImage, [1 1 3]); end
set(handles.pbSettingsMetrics, 'CData', dImage/max(dImage(:)));

x = double(imread([currpath,filesep,'icons',filesep,'papierkorb.jpg']))./255;
I2 = imresize(x, [16 16]);
set(handles.pb_deleteNo, 'CData', I2/max(I2(:)));

%% get voxelsize from preferences, if no SGeo in handles.h!
if ~isfield(handles.h, 'SGeo')|| isempty(handles.h.SGeo.cVoxelsize)
    handles.h.SGeo.cVoxelsize{1} = standardVoxelsize;
end
if ~isfield(handles.h,'nSCT')
    handles.nSCT = 1;           % slice orientation (1 = coronal, 2 = sagittal, 3 = transversal)
else
    handles.nSCT = handles.h.nSCT;
    set(handles.pm_CorSagTrans,'Value',handles.nSCT)
end

if(~isfield(handles.h, 'SDeform'))
    errordlg('No deformation field found. Please provide appropriate file.');
    return;
else
    handles.SDeform = handles.h.SDeform;
end

% Choose default command line output for LandmarkGUI
handles.output = hObject;
handles.slice = floor(size(handles.h.dImg,3)/2);         % slice number
handles.nMovGate = 2;       % compare first Reference Gate with Moving Gate 2
handles.iNGates = size(handles.h.dImg,4);
tmp = repmat({'G'},handles.iNGates,1);
handles.names = strcat(tmp,cellfun(@(x) num2str(x,'%02d'), num2cell(1:handles.iNGates).', 'UniformOutput', false));
for iG=1:handles.iNGates
    handles.fix.(handles.names{iG}) = [];       % arrays for fixed points in different gates
    handles.cPolyCoord.(handles.names{iG}){3,1} = {};    % cell array for Line Polygons
    handles.cPolyROICoord.(handles.names{iG}){3,1} = {};    % cell array for ROI Polygons
end

handles.ROINumber=ones(handles.iNGates,1);       % global number to count lines
handles.lineNumber=ones(handles.iNGates,1);      %                         ROIs
handles.FPointsLines = 1;   % default set-tool is setting points
handles.setActive = 1;      % flag, set function is active for start up
handles.unfixedpoints = 0;  % flag, turns 1 if new landmarks are set
handles.unfixedlines = 0;   % flag, turns 1 if new lines or ROIs are set
handles.axeslink = 1;       % flag for linked axes, important for scrolling
handles.FButtonDown = 0;    % no mouse click at the beginning
handles.exchangePoints = 0;
handles.sImageData='origImage'; % start with original image data (can be changed to interpolated images)
handles.colScale = 4;       % scale default color range
set(handles.markRefAx,'xtick',[],'ytick',[]);
set(handles.markMovAx,'xtick',[],'ytick',[]);
set(handles.axBlackRef, 'xtick', [], 'ytick', []);
set(handles.axBlackMov, 'xtick', [], 'ytick', []);
handles.noScrolling = 0;
handles.cTablePoints = num2cell(repmat('-',1,size(handles.h.dImg,4)));
handles.cTableLines = num2cell(repmat('-',1,size(handles.h.dImg,4)));
handles.cTableROIs = num2cell(repmat('-',1,size(handles.h.dImg,4)));
handles.FCellSelectionActive = 1;
% set(handles.tb_points,'ColumName', handles.names);

% try set(handles.tFilename, 'String', ['Dataset: ',handles.h.sFilename(1:min(end,90))]); catch; end;
methods = {'elastix', 'halar', 'grics', 'LAP', 'demons'};
method = methods{handles.h.nRegMethod};
handles.RegMParam = ['Registration Method: ', method];
if isfield(handles.h, 'sParFile')
    handles.RegMParam = [handles.RegMParam, ', Parameter: ',handles.h.sParFile];
end
set(handles.tRegMParm, 'String', handles.RegMParam)
set(handles.rb_interpImage,'Value',0);
set(handles.rb_origImage,'Value',1);
handles.cMetrics = {'smi', 'snmi', 'efinfo', 'rmi', 'rnmi', 'tmi', 'tnmi', 'ejp', 'gre', 'finfo', 'mssim', 'cc', 'zcc', 'spr', 'ket', 'mse', 'nmse', ...
                    'rmse', 'nrmse', 'psnr', 'tam', 'zca', 'zcr', 'mr', 'ssd', 'msd', 'nssd', 'sad', 'zsad', 'lsad', 'mad', 'shd', 'besov'; ...
                    'Shannon mutual information', 'Shannon normalized mutual information', 'exclusive F-information', 'Renyi mutual information', 'Renyi normalized mutual information', 'Tsallis mutual information', 'Tsallis normalized mutual information', 'energy of joint probability', 'gradient entropy', ...
                    'F-information measures', 'Mean structural similarity', '(Pearson) cross correlation', 'zero-mean normalized (Pearson) cross correlation', 'Spearman rank correlation', ...
                    'Kendall''s tau', 'mean squared error', 'normalized mean squared error', 'root mean squared error', 'normalized root mean squared error', 'peak signal to noise ratio', ...
                    'Tanimoto measure', 'zero crossings (absolute)', 'zero crossings (relative)', 'minimum ratio', 'sum of squared differences', 'median of squared differences', 'normalized sum of squared differences', ...
                    'sum of absolute differences', 'zero-mean sum of absolute differences', 'locally scaled sum of absolute differences', 'median of absolute differences', 'sum of hamming distance', 'besov norm'};
if(~exist('lEvalMetrics','var'))
    handles.lEvalMetrics = logical([0 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0]); % default to be evaluated intensity-based metrics
else
    handles.lEvalMetrics = lEvalMetrics;
end

% normalize images
handles.h.dImg = handles.h.dImg/max(handles.h.dImg(:));
handles.h.dImgReg = handles.h.dImgReg/max(handles.h.dImgReg(:));

handles = plotImages(handles);
axes(handles.axRefGate);
handles.iActive = 0;
% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = LandmarkGUI_OutputFcn(~, ~, handles)

% Get default command line output from handles structure
varargout{1} = handles.output;
try
if(isfield(handles,'closeFigure') && handles.closeFigure)
    LandmarkGUI_CloseRequestFcn(hObject, eventdata, handles)
end
catch
end


function pm_PointsLinesROIs_Callback(hObject, ~, handles)
% select to set points, lines or ROIs
contents = cellstr(get(hObject,'String'));
lead=contents{get(hObject,'value')};

switch lead
    case 'Points'
        handles.FPointsLines = 1;
        set(handles.tb_points,'Data',handles.cTablePoints);
        set(handles.pb_clear_set,'String','clear set points');
        set(handles.pb_clear_all,'String','clear all points');
        
    case 'Lines'
        handles.FPointsLines = 2;
        set(handles.tb_points,'Data',handles.cTableLines);
        set(handles.pb_clear_set,'String','clear set lines');
        set(handles.pb_clear_all,'String','clear all lines');
    case 'ROIs'
        handles.FPointsLines = 3;
        set(handles.tb_points,'Data',handles.cTableROIs);
        set(handles.pb_clear_set,'String','clear set ROIs');
        set(handles.pb_clear_all,'String','clear all ROIs');
    otherwise
        error('Unknown error.')
end
guidata(hObject, handles)


function [exchP, newSlice] = pb_set_lmp_Callback(hObject, eventdata, handles)
% set points, lines and ROIs
handles.noScrolling = 1; guidata(hObject, handles);
if handles.setActive
    
    % set landmark points in the reference image (gate 01)
    if strcmp(handles.sImageData,'interpImage')
        errordlg('Landmark points can only be set in original pixel configuration.')
        exchP = []; newSlice = [];
        return
    end
    
    handles.setActive = 0;
    
    % -------------------------------------------------------------------------
    % Depending on selected type of labeling set points, lines or ROIs
    % ---------------------------
    if handles.FPointsLines == 1
        
        if gca == handles.axRefGate
                        
            % delete old plot of landmarks
            try delete(handles.reflm); delete(handles.reflmLabel); catch; end;
            
            % set new landmarks with fGetPts
            [a, b] = fGetPts(handles.axRefGate);
            
            % get coordinates of pixel
            x = floor(a+0.5); y = floor(b+0.5); % +0.5 necessary for centering of pixel
            
            % delete points outside the image borders
            A = zeros(length(x),2);
            A(:,1) = floor(x);
            A(:,2) = floor(y);
            A(any(A'<0),:) = [];
            A(A(:,1)> size(handles.h.dImg,2),:) = [];
            A(A(:,2)> size(handles.h.dImg,1),:) = [];
            %     if size(A,1) ~= length(x)
            %        errordlg('There were points outside the image. Please check if the labels correspond correctly.')
            %     end
            x = A(:,1);
            y = A(:,2);
            
            % return if no points were selected
            if isempty(x)
                handles.setActive = 1;
                handles.noScrolling = 0;
                try handles.xes = rmfield(handles.xes,'G01'); catch; end;
                guidata(hObject, handles)
                set(handles.tb_points,'Data',handles.cTablePoints);
                exchP=[]; newSlice=[];
                return;
            end
            
            % replace only one point if function is called by the table cell
            % selection callback
            if handles.exchangePoints == 1
                exchP = [x,y];
                newSlice = handles.slice;
                if size(x,1) > 1
                    edlg = errordlg('Please set only one new point.');
                    uiwait(edlg);
                    handles.setActive = 1;
                    guidata(hObject, handles);
                    [exchP, newSlice] = pb_set_lmp_Callback(hObject, eventdata, handles);
                end
                handles.setActive = 1;
                handles.noScrolling = 0;
                handles.FCellSelectionActive = 1;
                handles.exchangePoints = 0;
                guidata(hObject, handles)
                set(handles.tb_points,'Data',handles.cTablePoints);
                return;
            end
            
            % check if number of gates is equal to number of previously set points
            % in other gates
            if isfield(handles,'xes')
                try nPts = length(handles.xes.G02); catch; end;
                try nPts = length(handles.xes.G03); catch; end;
                try nPts = length(handles.xes.G04); catch; end;
                try
                    if nPts ~= length(x)
                        errordlg('Please set the same number of points in every gate.')
                        return
                    end
                catch
                end
            end
            
            
            % pass to global handle
            handles.xes.G01 = x;
            handles.yes.G01 = y;
            handles.zes.G01 = handles.slice.*ones(size(x));
            
            % plot new landmarks
            axes(handles.axRefGate)
            hold on
            handles.reflm = plot(x,y,'yx', 'Linewidth', 1.5,'MarkerSize', 5);
            labels = cellstr( num2str([1:size(x)]'));  % labels corresponding to x,y order
            handles.reflmLabel = text(x(:,1), y(:,1), labels,'Color','y', 'VerticalAlignment','bottom', ...
                'HorizontalAlignment','right');
            hold off
            
            handles.unfixedpoints = 1;
            
        elseif gca == handles.axMovGate
            % set landmark points in current moving image (gate 02-04)
            
            % delete old plot of landmarks
            try delete(handles.movlm.(handles.names{handles.nMovGate}));
                delete(handles.movlmLabel.(handles.names{handles.nMovGate})); catch
            end
            
            axMovGate_ButtonDownFcn(hObject, eventdata, handles)
            % set new landmarks
            [a, b] = fGetPts(handles.axMovGate);            
            
            % get coordinates of pixel
            x = floor(a+0.5); y = floor(b+0.5);
            
            % delete points outside the image borders
            A = zeros(length(x),2);
            A(:,1) = x;
            A(:,2) = y;
            A(any(A'<0),:) = [];
            A(A(:,1)> size(handles.h.dImg,2),:) = [];
            A(A(:,2)> size(handles.h.dImg,1),:) = [];
            %     if size(A,1) ~= length(x)
            %        errordlg('There were points outside the image. Please check if the labels correspond correctly.')
            %     end
            x = A(:,1);
            y = A(:,2);
            
            % return if no points were selected
            if isempty(a)
                handles.setActive = 1;
                handles.noScrolling = 0;
                try handles.xes = rmfield(handles.xes,handles.names{handles.nMovGate}); catch; end;
                guidata(hObject, handles)
                set(handles.tb_points,'Data',handles.cTablePoints);
                exchP=[]; newSlice=[];
                return;
            end
                        
            % replace only one point if function is called by the table cell
            % selection callback
            if handles.exchangePoints == 1
                exchP = [x,y];
                newSlice = handles.slice;
                if size(x,1) > 1
                    edlg = errordlg('Please set only one new point.');
                    uiwait(edlg);
                    handles.setActive = 1;
                    guidata(hObject, handles);
                    [exchP, newSlice] = pb_set_lmp_Callback(hObject, eventdata, handles);
                end
                handles.setActive = 1;
                handles.noScrolling = 0;
                handles.FCellSelectionActive = 1;
                handles.exchangePoints = 0;
                guidata(hObject, handles)
                set(handles.tb_points,'Data',handles.cTablePoints);
                return;
            end
            
            % check if number of gates is equal to number of previously set points
            % in other gates
            if isfield(handles,'xes') && ...
                    any(isfield(handles.xes,handles.names))
%                     (isfield(handles.xes,'G01') ||isfield(handles.xes,'G02') ||isfield(handles.xes,'G03') ||isfield(handles.xes,'G04'))
                for iG=1:length(handles.names)
                    nPts = [];
                    try nPts = length(handles.xes.(handles.names{iG})); catch; end;
%                 try nPts = length(handles.xes.G02); catch; end;
%                 try nPts = length(handles.xes.G03); catch; end;
%                 try nPts = length(handles.xes.G04); catch; end;
                    if(~isempty(nPts) && nPts ~= length(x))
                        errordlg('Please set the same number of points in every gate.')
                        handles.setActive = 1;
                        return
                    end
                end
            end
            
            % save in handles
            handles.xes.(handles.names{handles.nMovGate}) = x;
            handles.yes.(handles.names{handles.nMovGate}) = y;
            handles.zes.(handles.names{handles.nMovGate}) = handles.slice.*ones(size(x));
            
            % plot new landmarks
            axes(handles.axMovGate)
            hold on
            handles.movlm.(handles.names{handles.nMovGate}) = plot(x,y,'yx', 'Linewidth', 1.5,'MarkerSize', 5);
            labels = cellstr( num2str([1:size(x)]'));  % labels corresponding to x,y order
            handles.movlmLabel.(handles.names{handles.nMovGate}) =...
                text(x(:,1), y(:,1), labels,'Color','y', 'VerticalAlignment','bottom', ...
                'HorizontalAlignment','right');
            hold off
            
            handles.unfixedpoints = 1;
        else
            errordlg('Select axes where you want to set points.')
        end
        
        
        % ---------------------------
    elseif handles.FPointsLines == 2; % set lines
        if gca == handles.axRefGate
            if isfield(handles,'hPoly') && isfield(handles.hPoly,'G01')
                delete(handles.hPoly.G01)       % only one line at a time can be set
            end
            % -------------------------------------
            % set landmark lines in reference image
            handles.hPoly.G01 = impoly('Closed',false);
            
            if isempty(handles.hPoly.G01)
                handles.hPoly = rmfield(handles.hPoly,'G01');
                handles.setActive = 1;
                handles.noScrolling = 0;
                guidata(hObject, handles)
                return
            end
            
            setColor(handles.hPoly.G01,'y')
            handles.hPolySlice.G01 = handles.slice;
            handles.xyPos.G01 = getPosition(handles.hPoly.G01);
            handles.unfixedlines = 1;
            
        elseif gca == handles.axMovGate
            if isfield(handles,'hPoly') && isfield(handles.hPoly, handles.names{handles.nMovGate})
                delete(handles.hPoly.(handles.names{handles.nMovGate}))       % only one line at a time can be set
            end
            % -------------------------------------
            % set landmark lines in moving image
            handles.hPoly.(handles.names{handles.nMovGate}) = impoly('Closed',false);
            
            if isempty(handles.hPoly.(handles.names{handles.nMovGate}))
                handles.hPoly = rmfield(handles.hPoly,handles.names{handles.nMovGate});
                handles.setActive = 1;
                handles.noScrolling = 0;
                guidata(hObject, handles)
                return
            end
            
            setColor(handles.hPoly.(handles.names{handles.nMovGate}),'y')
            handles.hPolySlice.(handles.names{handles.nMovGate}) = handles.slice;
            handles.xyPos.(handles.names{handles.nMovGate}) = getPosition(handles.hPoly.(handles.names{handles.nMovGate}));
            
            handles.unfixedlines = 1;
            
        else
            errordlg('Select axes where you want to set points.')
        end
        
        % ---------------------------
    elseif handles.FPointsLines == 3; % set ROIs
        if gca == handles.axRefGate
            if isfield(handles,'hPolyROI.G01')
                delete(handles.hPolyROI.G01)       % only one line at a time can be set
            end
            % -------------------------------------
            % set ROI in reference image, polygon is closed at the end
            handles.hPolyROI.G01 = impoly('Closed',true);
            
            if isempty(handles.hPolyROI.G01)
                handles.hPolyROI = rmfield(handles.hPolyROI,'G01');
                handles.setActive = 1;
                handles.noScrolling = 0;
                guidata(hObject, handles)
                return
            end
            
            setColor(handles.hPolyROI.G01,'y')
            handles.hPolyROISlice.G01 = handles.slice;
            handles.xyPosROI.G01 = getPosition(handles.hPolyROI.G01);
            
            handles.unfixedlines = 1;
            
        elseif gca == handles.axMovGate
            if isfield(handles,'hPolyROI') && isfield(handles.hPolyROI, handles.names{handles.nMovGate})
                delete(handles.hPolyROI.(handles.names{handles.nMovGate}))       % only one line at a time can be set
            end
            % -------------------------------------
            % set landmark lines in moving image
            handles.hPolyROI.(handles.names{handles.nMovGate}) = impoly('Closed',true);
            
            if isempty(handles.hPolyROI.(handles.names{handles.nMovGate}))
                handles.hPolyROI = rmfield(handles.hPolyROI,handles.names{handles.nMovGate});
                handles.setActive = 1;
                handles.noScrolling = 0;
                guidata(hObject, handles)
                return
            end
            
            setColor(handles.hPolyROI.(handles.names{handles.nMovGate}),'y')
            handles.hPolyROISlice.(handles.names{handles.nMovGate}) = handles.slice;
            handles.xyPosROI.(handles.names{handles.nMovGate}) = getPosition(handles.hPolyROI.(handles.names{handles.nMovGate}));
            
            handles.unfixedlines = 1;
            
        else
            errordlg('Select axes where you want to set points.')
        end
        
    else
        error('Unknown error.')
    end
    handles.setActive = 1;
end
handles.noScrolling = 0;
guidata(hObject,handles)


function pb_fix_lmp_Callback(hObject, ~, handles)
% fix the landmark points, lines or ROIs that are set in previous steps

% -------------------------------------------------------------------------
% Depending on selected type of labeling fix points, lines or ROIs
% ---------------------------
if handles.FPointsLines == 1
    
    % check if lmp are set in all 4 gates
    try
        if ~isfield(handles,'xes')||~isfield(handles.xes,'G01')||~isfield(handles.xes,'G02')||~isfield(handles.xes,'G03')||~isfield(handles.xes,'G04')||~isfield(handles,'xes')
            errordlg('Please set corresponding points in all 4 gates.');
            return
        end
    catch
    end
    
    % get the positions of all lmps and save them in fix.G0x by concatenation
    nmb = (1:size(handles.xes.G01))';
    x = handles.xes.G01;
    y = handles.yes.G01;
    z = handles.zes.G01;
    handles.fix.G01 = cat(1,handles.fix.G01,[nmb, x, y, z]);
    for iG=2:handles.iNGates
        x = handles.xes.(handles.names{iG});
        y = handles.yes.(handles.names{iG});
        z = handles.zes.(handles.names{iG});
        handles.fix.(handles.names{iG}) = cat(1,handles.fix.(handles.names{iG}),[nmb, x, y, z]);
    end
%     x = handles.xes.G03;
%     y = handles.yes.G03;
%     z = handles.zes.G03;
%     handles.fix.G03 = cat(1,handles.fix.G03,[nmb, x, y, z]);
%     x = handles.xes.G04;
%     y = handles.yes.G04;
%     z = handles.zes.G04;
%     handles.fix.G04 = cat(1,handles.fix.G04,[nmb, x, y, z]);
    
    handles = rmfield(handles, {'xes','yes','zes'});
    
    % display plotted points in table
    cDataTemp = repmat('x',size(handles.fix.G01,1),4);
    handles.cTablePoints = num2cell(cDataTemp);
    set(handles.tb_points,'Data',handles.cTablePoints);
    
    handles = fPlotFixedMarker(handles);
    
    % remove the landmark points that are set earlier and delete the plot
    try delete(handles.reflm); catch; end;
    try delete(handles.reflmLabel); catch; end;
    try handles=rmfield(handles,'reflm');catch; end
    try handles=rmfield(handles,'reflmLabel');catch; end
    for iI = 1:4
        try delete(handles.movlm.(handles.names{iI})); catch; end;
        try delete(handles.movlmLabel.(handles.names{iI})); catch; end;
        try handles=rmfield(handles,movlm.(handles.names{iI}));catch; end
        try handles=rmfield(handles,movlmLabel.(handles.names{iI}));catch; end
    end
    handles.unfixedpoints = 0; % there are no longer unfixed points
    
    % ---------------------------
elseif handles.FPointsLines == 2
    
    
    %    todo:
    %     % if call from table cell selection put line to correct position in
    %     % cell array
    %     try
    %     if isfield(handles,'lineIndex')
    %         lineNo = handles.lineIndex;
    %         if gca == handles.axRefGate; i=1; else i=handles.nMovGate; end
    %         xyzPoly = [getPosition(handles.hPoly.(handles.names{i})),...
    %             handles.hPolySlice.(handles.names{i}).*ones(size(getPosition(handles.hPoly.(handles.names{i})),1),1)...
    %             lineNo.*ones(size(getPosition(handles.hPoly.(handles.names{i})),1),1)];
    %         xyzPoly = round(xyzPoly);
    %         if handles.nSCT == 2
    %             l = 1; m = 3;
    %             xyzPoly(:,[l,m])= xyzPoly(:,[m,l]);
    %         elseif handles.nSCT == 3
    %             l = 2; m = 3;
    %             xyzPoly(:,[l,m])= xyzPoly(:,[m,l]);
    %         else
    %             % do nothing
    %         end
    %         % find cell array position for the next coordinates
    %         b = handles.cPolyCoord.(handles.names{i})(handles.nSCT,:);
    %         k=1;
    %         % save coordinates
    %         handles.cPolyCoord.(handles.names{i}){handles.nSCT,k} = xyzPoly;
    %         delete(handles.hPoly.(handles.names{i}))
    %         handles = rmfield(handles, 'lineIndex');
    %         guidata(hObject, handles);
    %         return
    %     end
    %     catch; 'here';
    %     end
    
    % save coordinates of the line vertices in a global coordinate system
    % x left-right, y head-foot and z anterior-posterior 
    try
        i = 1;
        xyzPoly = [getPosition(handles.hPoly.(handles.names{i})),...
            handles.hPolySlice.(handles.names{i}).*ones(size(getPosition(handles.hPoly.(handles.names{i})),1),1)...
            handles.lineNumber(i).*ones(size(getPosition(handles.hPoly.(handles.names{i})),1),1)];
%         xyzPoly = round(xyzPoly);
        % increment line number
        handles.lineNumber(i) = handles.lineNumber(i)+1;
        if handles.nSCT == 2
            l = 1; m = 3;
            xyzPoly(:,[l,m])= xyzPoly(:,[m,l]);
        elseif handles.nSCT == 3
            l = 2; m = 3;
            xyzPoly(:,[l,m])= xyzPoly(:,[m,l]);
        else
            % do nothing
        end
        % find cell array position for the next coordinates
        A = cellfun('isempty',handles.cPolyCoord.(handles.names{i}));
        b = find(A(handles.nSCT,:)==0,1,'last');
        if isempty(b); k=1; else k = 1+b; end
        % save coordinates
        handles.cPolyCoord.(handles.names{i}){handles.nSCT,k} = xyzPoly;
        delete(handles.hPoly.(handles.names{i}))
    catch
    end
    try
        i = handles.nMovGate;
        xyzPoly = [getPosition(handles.hPoly.(handles.names{i})),...
            handles.hPolySlice.(handles.names{i}).*ones(size(getPosition(handles.hPoly.(handles.names{i})),1),1),...
            handles.lineNumber(i).*ones(size(getPosition(handles.hPoly.(handles.names{i})),1),1)];
%         xyzPoly = round(xyzPoly);
        % increment line number
        handles.lineNumber(i) = handles.lineNumber(i)+1;
        if handles.nSCT == 2
            l = 1; m = 3;
            xyzPoly(:,[l,m])= xyzPoly(:,[m,l]);
        elseif handles.nSCT == 3
            l = 2; m = 3;
            xyzPoly(:,[l,m])= xyzPoly(:,[m,l]);
        else
            % do nothing
        end
        % find cell array position for the next coordinates
        A = cellfun('isempty',handles.cPolyCoord.(handles.names{i}));
        b = find(A(handles.nSCT,:)==0,1,'last');
        if isempty(b); k=1; else k = 1+b; end
        % save coordinates
        handles.cPolyCoord.(handles.names{i}){handles.nSCT,k} = xyzPoly;
        delete(handles.hPoly.(handles.names{i}))
    catch
    end
    try handles = rmfield(handles,'hPoly'); catch; end
    
    handles = fPlotFixedMarker(handles);
    handles.unfixedlines = 0;
    
    % display plotted lines in table
    for i = 1:size(handles.names,1)
        B = cellfun('isempty',handles.cPolyCoord.(handles.names{i}));
        numLines(i) = sum(~B(:));
    end
    cData = repmat('-',max(numLines),handles.iNGates);
    for i = 1:size(handles.names,1)
        cDataTemp = repmat('L',numLines(i),1);
        cData(1:length(cDataTemp),i) = cDataTemp;
    end
    handles.cTableLines = num2cell(cData);
    set(handles.tb_points,'Data',handles.cTableLines);
    
    
    % ---------------------------
elseif handles.FPointsLines == 3
    % save coordinates of the line vertices in a global coordinate system
    % x left-right, y head-foot and z anterior-posterior 
   try
        i = 1;
        xyzPoly = [getPosition(handles.hPolyROI.(handles.names{i})),...
            handles.hPolyROISlice.(handles.names{i}).*ones(size(getPosition(handles.hPolyROI.(handles.names{i})),1),1)...
            handles.ROINumber(i).*ones(size(getPosition(handles.hPolyROI.(handles.names{i})),1),1)];
        xyzPoly(end+1,:) = xyzPoly(1,:);
        xyzPoly = round(xyzPoly);
        % increment ROI number
        handles.ROINumber(i) = handles.ROINumber(i)+1;
        if handles.nSCT == 2
            l = 1; m = 3;
            xyzPoly(:,[l,m])= xyzPoly(:,[m,l]);
        elseif handles.nSCT == 3
            l = 2; m = 3;
            xyzPoly(:,[l,m])= xyzPoly(:,[m,l]);
        else
            % do nothing
        end
        % find cell array position for the next coordinates
        A = cellfun('isempty',handles.cPolyROICoord.(handles.names{i}));
        b = find(A(handles.nSCT,:)==0,1,'last');
        if isempty(b); k=1; else k = 1+b; end
        % save coordinates
        handles.cPolyROICoord.(handles.names{i}){handles.nSCT,k} = xyzPoly;
        delete(handles.hPolyROI.(handles.names{i}))
    catch
    end
    try
        i = handles.nMovGate;
        xyzPoly = [getPosition(handles.hPolyROI.(handles.names{i})),...
            handles.hPolyROISlice.(handles.names{i}).*ones(size(getPosition(handles.hPolyROI.(handles.names{i})),1),1),...
            handles.ROINumber(i).*ones(size(getPosition(handles.hPolyROI.(handles.names{i})),1),1)];
        xyzPoly(end+1,:) = xyzPoly(1,:);
        xyzPoly = round(xyzPoly);
        % increment ROI number
        handles.ROINumber(i) = handles.ROINumber(i)+1;
        if handles.nSCT == 2
            l = 1; m = 3;
            xyzPoly(:,[l,m])= xyzPoly(:,[m,l]);
        elseif handles.nSCT == 3
            l = 2; m = 3;
            xyzPoly(:,[l,m])= xyzPoly(:,[m,l]);
        else
            % do nothing
        end
        % find cell array position for the next coordinates
        A = cellfun('isempty',handles.cPolyROICoord.(handles.names{i}));
        b = find(A(handles.nSCT,:)==0,1,'last');
        if isempty(b); k=1; else k = 1+b; end
        % save coordinates
        handles.cPolyROICoord.(handles.names{i}){handles.nSCT,k} = xyzPoly;
        delete(handles.hPolyROI.(handles.names{i}))
    catch
    end
    try handles = rmfield(handles,'hPolyROI'); catch; end
    
    handles = fPlotFixedMarker(handles);
    handles.unfixedlines = 0;
    
    % display plotted ROIs in table
    for i = 1:size(handles.names,1)
        B = cellfun('isempty',handles.cPolyROICoord.(handles.names{i}));
        numLines(i) = sum(~B(:));
    end
    cData = repmat('-',max(numLines),handles.iNGates);
    for i = 1:size(handles.names,1)
        cDataTemp = repmat('R',numLines(i),1);
        cData(1:length(cDataTemp),i) = cDataTemp;
    end
    handles.cTableROIs = num2cell(cData);
    set(handles.tb_points,'Data',handles.cTableROIs);
    
else
    error('Unknown error!')
end

guidata(hObject, handles)


function cmp_lmp_Callback(~, ~, handles)
% start comparing set points, lines and ROIs
if handles.nSCT ~= 1
    warning('Computation will be performed in coronal view.');
    set(handles.pm_CorSagTrans,'Value',1); % coronal view
    handles.unfixedlines = 0;
    pm_CorSagTrans_Callback(hObject, [], handles);
end
if size(handles.fix.G01)==[0,0]
    errordlg('There are no fixed points to compare. Please set and fix landmarkpoints first.')
    return;
else
    dVoxSz = handles.h.SGeo.cVoxelsize;
    % get the coordinates of the fixed points that shall be compared
    for iG=1:handles.iNGates
        eval(sprintf('G%02d = handles.fix.G%02d; x%02d = G%02d(:,2); y%02d = G%02d(:,3); z%02d = G%02d(:,4);', iG*ones(1,8)));
    end
    % G02 = handles.fix.G02; x02 = G02(:,2); y02 = G02(:,3); z02 = G02(:,4);
    % G03 = handles.fix.G03; x03 = G03(:,2); y03 = G03(:,3); z03 = G03(:,4);
    % G04 = handles.fix.G04; x04 = G04(:,2); y04 = G04(:,3); z04 = G04(:,4);

    if strcmp(handles.sImageData, 'origImage')
        % transform the coordinates of the moving images according to
        % registration result
        for iI = 2:4
            dBx(:,:,:,iI) = handles.h.SDeform(iI).dBx;
            dBy(:,:,:,iI) = handles.h.SDeform(iI).dBy;
            dBz(:,:,:,iI) = handles.h.SDeform(iI).dBz;
        end
        for iG=2:handles.iNGates
            for i = 1:length(x02)
                eval(sprintf('x%02dn(i) = x%02d(i) + dBx(x%02d(i),y%02d(i),z%02d(i),iG);', iG*ones(1,5)));
                eval(sprintf('y%02dn(i) = y%02d(i) + dBy(x%02d(i),y%02d(i),z%02d(i),iG);', iG*ones(1,5)));
                eval(sprintf('z%02dn(i) = z%02d(i) + dBz(x%02d(i),y%02d(i),z%02d(i),iG);', iG*ones(1,5)));

            end 
        end
    else
        for iG=2:handles.iNGates
            eval(sprintf('x%02dn=x%02d'';y%02dn=y%02d'';z%02dn=z%02d'';', iG*ones(1,6)));
        end

    end
    
    cPoints = cell(2,handles.iNGates-1); % mean+std x images
    for iG=2:handles.iNGates
        % compute the euclidean distance
        eval(sprintf('dist01%02d = sqrt(((x01*dVoxSz{1}(1)-x%02dn''*dVoxSz{iG}(1))).^2+((y01*dVoxSz{1}(2)-y%02dn''*dVoxSz{iG}(2))).^2+((z01*dVoxSz{1}(3)-z%02dn''*dVoxSz{iG}(3))).^2);', iG*ones(1,4)));

        eval(sprintf('cPoints{1,iG-1} = mean(dist01%02d);', iG));
        eval(sprintf('cPoints{2,iG-1} = std(dist01%02d);', iG));
    end
end

% open figure for outputting results
if(isfield(handles, 'cPolyCoord'))
    [xline, yline] = find(~cellfun(@isempty,handles.cPolyCoord.G01));
    iNLinePlots = size(xline,1);
    clear 'xline' 'yline';
else
    iNLinePlots = 0;
end
if(isfield(handles,'cPolyROICoord'))
    [xline, yline] = find(~cellfun(@isempty,handles.cPolyROICoord.G01));
    iNROIPlots = size(xline,1);
    clear 'xline' 'yline';
else
    iNROIPlots = 0;
end
iSimMetrics = get(handles.cbMetrics,'Value');

iNPlots = [handles.iNGates-1, iNLinePlots, iNROIPlots, iSimMetrics];
[hfInsertPoint, hfInsertLineROI, hfPlotLabel, hfPlotMetrics] = fLandmarkCompare(iNPlots, dVoxSz{1});
for iG=2:handles.iNGates
   hfInsertPoint([cPoints{1,iG-1},cPoints{2,iG-1}],iG);
end

if isfield(handles,'cPolyCoord')
    
    if strcmp(handles.sImageData, 'origImage')
        results.cLines = handles.cPolyCoord;
        [diceLine,diceLine0,diceLine1] = fCompareLines(handles.h,results,{hfInsertLineROI, hfPlotLabel}, handles.names, handles.iNGates);
    else
        results.cLines = handles.cPolyCoord;
        [diceLine,diceLine0,diceLine1] = fCompareLinesNoTr(handles.h,results,{hfInsertLineROI, hfPlotLabel}, handles.names, handles.iNGates);
        
    end
    
else
    warning('There are no lines to compare.')
end

if isfield(handles,'cPolyROICoord')
    
    if strcmp(handles.sImageData, 'origImage')
        results.cROIs = handles.cPolyROICoord;
        diceROI = fCompareROIs(handles.h,results,{hfInsertLineROI, hfPlotLabel}, handles.names, handles.iNGates);
        if(iSimMetrics > 0)
            cString = handles.cMetrics(1,handles.lEvalMetrics).';
            metricsROI = fMetricROIs(handles.h,results,hfPlotMetrics, handles.names, handles.iNGates, cString, true);
        end
    else
        results.cROIs = handles.cPolyROICoord;
        diceROI = fCompareROIsNoTr(handles.h,results,{hfInsertLineROI, hfPlotLabel}, handles.names, handles.iNGates);  
        if(iSimMetrics > 0)
            cString = handles.cMetrics(1,handles.lEvalMetrics).';            
            metricsROI = fMetricROIs(handles.h,results,hfPlotMetrics, handles.names, handles.iNGates, cString, false);
        end
    end
else
    warningdlg('There are no ROIs to compare.')
end

try results.EuclidDist = cPoints; catch; end;
try results.diceLine0 = diceLine0; catch; end;
try results.diceLine1 = diceLine1; catch; end;
try results.diceLine = diceLine; catch; end;
try results.diceROI = diceROI; catch; end;
try results.metricsROI = metricsROI; catch; end;


isSave = questdlg('Do you want to save the feature-based results?', 'Save results', 'Yes', 'No','No');
switch isSave
    case 'Yes'
    case 'No'
        return
    case ''
        return
end
methods = {'elastix', 'halar', 'grics', 'lap', 'demons'};
sMethod = methods{handles.h.nRegMethod};
if(~isfield(handles,'SPaths'))
    [sFilename,sPathname] = uiputfile([pwd,filesep,'LandmarkResult.mat'],'Specify file name');
else
    sSavename = [handles.SPaths.sResults,sMethod,filesep,'LandmarkResult.mat'];
    [sFilename,sPathname] = uiputfile(sSavename,'Specify file name');
end
% [sFilename,sPathname] = uiputfile([handles.SPaths.sResults,sMethod,'*.mat'],'Specify file name');
if sFilename == 0; return; end;
save([sPathname, sFilename],'results')
msgbox('Saved results.', 'Save complete.','help')


function pb_nxt_img_Callback(hObject, ~, handles)
% get next image in moving image axes (right image)

% handles.nMovGate = mod(handles.nMovGate+2,3)+2; % is 2,3,4,2,3,4,..
handles.nMovGate = handles.nMovGate + 1;
if(handles.nMovGate > handles.iNGates), handles.nMovGate = 2; end;

% flush image handle in order to force a new drawing (otherwise currently to
% be set but not fixed points will occur in the next image)
delete(handles.hI2); handles.hI2 = [];
handles = showMovImage(handles);
guidata(hObject, handles)


function pb_pre_img_Callback(hObject, ~, handles)
% get previous image in moving image axes (right image)

% handles.nMovGate = mod(handles.nMovGate+3,3)+2; % is 4,3,2,4,3,2,..
handles.nMovGate = handles.nMovGate - 1;
if(handles.nMovGate < 2), handles.nMovGate = handles.iNGates; end;

% flush image handle in order to force a new drawing (otherwise currently to
% be set but not fixed points will occur in the next image)
delete(handles.hI2); handles.hI2 = [];
handles = showMovImage(handles);
guidata(hObject, handles)


function LandmarkGUI_KeyPressFcn(hObject, eventdata, handles)
% collect figure key press events

if eventdata.Key == 's'
    pb_set_lmp_Callback(hObject,eventdata,handles);
elseif eventdata.Key == 'f'
    pb_fix_lmp_Callback(hObject, eventdata, handles);
elseif strcmp(eventdata.Key,'space')
    pb_nxt_img_Callback(hObject, eventdata, handles);
else
%     eventdata.Character;
end


function CBlink_Callback(hObject, ~, handles)
% link reference and moving axes

if get(hObject, 'Value')
    linkaxes([handles.axRefGate, handles.axMovGate])
    handles.axeslink = 1;
else
    linkaxes([handles.axRefGate, handles.axMovGate],'off')
    handles.axeslink = 0;
end
guidata(hObject, handles)


function LandmarkGUI_WindowScrollWheelFcn(hObject, eventdata, handles)
% figure scroll wheel functionality

if handles.noScrolling
    return
end

for i = 1:4
    try handles.xyPos.(handles.names{i}) = getPosition(handles.hPoly.(handles.names{i})); catch; end;
    try handles.xyPosROI.(handles.names{i}) = getPosition(handles.hPolyROI.(handles.names{i})); catch; end;
end

if handles.axeslink
    % delete old landmark point plots
    child_axRefGate = allchild(handles.axRefGate);
    for iI = 1:length(child_axRefGate)
        if child_axRefGate(iI) ~= handles.hI1
            delete(child_axRefGate(iI));
        end
    end
    child_axMovGate = allchild(handles.axMovGate);
    for iI = 1:length(child_axMovGate)
        if child_axMovGate(iI) ~= handles.hI2
            delete(child_axMovGate(iI));
        end
    end
    
    % get image arrays and deformation fields
    if strcmp(handles.sImageData,'origImage')
        dIRef = handles.h.dImg(:,:,:,1);
        dIMove2 = handles.h.dImg(:,:,:,handles.nMovGate);
    elseif strcmp(handles.sImageData,'regImage')
        dIRef = handles.h.dImg(:,:,:,1);
        dIMove2 = handles.h.dImgReg(:,:,:,handles.nMovGate);
    elseif strcmp(handles.sImageData,'interpImage')
        dIRef = handles.dImgInterp(:,:,:,1);
        dIMove2 = handles.dImgInterp(:,:,:,handles.nMovGate);
    else
        error('Unknown Error!')
    end
    
    % get scroll events
    if eventdata.VerticalScrollCount < 0
        handles.slice = max([1 handles.slice - 1]);
    elseif eventdata.VerticalScrollCount > 0
        handles.slice = min([size(dIRef, 3) handles.slice + 1]);
    else
        %nothing
    end
    
    % set new image data in axes
    set(handles.hI1, 'CData', dIRef(:,:,handles.slice));
    set(handles.hI2, 'CData', dIMove2(:,:,handles.slice));
    
    set(handles.tSliceRef, 'String',[num2str(handles.slice),'/',num2str(size(dIRef,3))]);
    set(handles.tSliceMov, 'String',[num2str(handles.slice),'/',num2str(size(dIRef,3))]);
    
    % check if there are any fixed marker in this slice and plot them
    handles = fPlotFixedMarker(handles);
    if handles.unfixedpoints
        try
            if handles.slice == handles.zes.G01(1,1) %unfixed points are all in the same slice
                axes(handles.axRefGate)
                x = handles.xes.G01;
                y = handles.yes.G01;
                hold on
                handles.reflm = plot(x,y,'yx', 'Linewidth', 1.5,'MarkerSize', 5);
                labels = cellstr( num2str([1:size(x)]'));  %' # labels correspond to their order
                handles.reflmLabel = text(x(:,1), y(:,1), labels,'Color','y', 'VerticalAlignment','bottom', ...
                    'HorizontalAlignment','right');
                hold off
            end
        catch
        end
        try
            if handles.slice == handles.zes.(handles.names{handles.nMovGate})(1,1) %unfixed points are all in the same slice
                axes(handles.axMovGate)
                x = handles.xes.(handles.names{handles.nMovGate});
                y = handles.yes.(handles.names{handles.nMovGate});
                labels = cellstr( num2str([1:size(x)]'));
                hold on
                handles.movlm.(handles.names{handles.nMovGate})...
                    = plot(x,y,'yx', 'Linewidth', 1.5,'MarkerSize', 5);
                handles.movlmLabel.(handles.names{handles.nMovGate})...
                    = text(x(:,1), y(:,1), labels,'Color','y', 'VerticalAlignment','bottom', ...
                    'HorizontalAlignment','right');
                hold off
            end
        catch
        end
    else
        try delete(handles.reflm); catch; end;
        try delete(handles.reflmLabel); catch; end;
        try delete(handles.movlm); catch; end;
        try delete(handles.movlmLabel); catch; end;
    end
    
    % plot lines and ROIs in reference image
    try
        if isfield(handles,'hPoly') && isfield(handles.hPoly,'G01') && handles.hPolySlice.G01 == handles.slice
            handles.hPoly.G01 = impoly(handles.axRefGate, handles.xyPos.G01,'Closed', false);
            setColor(handles.hPoly.G01,'y')
        end
    catch
    end
    try
        if isfield(handles,'hPolyROI') && isfield(handles.hPolyROI,'G01') && handles.hPolyROISlice.G01 == handles.slice
            handles.hPolyROI.G01 = impoly(handles.axRefGate, handles.xyPosROI.G01,'Closed', true);
            setColor(handles.hPolyROI.G01,'y')
        end
    catch
    end
    % plot lines and ROIs in moving image
    iI = handles.nMovGate;
    try
        if isfield(handles,'hPoly') && isfield(handles.hPoly,handles.names{iI}) && handles.hPolySlice.(handles.names{iI}) == handles.slice
            handles.hPoly.(handles.names{iI}) = impoly(handles.axMovGate, handles.xyPos.(handles.names{iI}),'Closed', false);
            setColor(handles.hPoly.(handles.names{iI}),'y')
        end
    catch
    end
    try
        if isfield(handles,'hPolyROI') && isfield(handles.hPolyROI,handles.names{iI}) && handles.hPolyROISlice.(handles.names{iI}) == handles.slice
            handles.hPolyROI.(handles.names{iI}) = impoly(handles.axMovGate, handles.xyPosROI.(handles.names{iI}),'Closed', true);
            setColor(handles.hPolyROI.(handles.names{iI}),'y')
        end
    catch
    end
        
else
    if gca == handles.axRefGate
        %         disp(num2str(gca))
        % delete old landmark point plots
        child_axRefGate = allchild(handles.axRefGate);
        for iI = 1:length(child_axRefGate)
            if child_axRefGate(iI) ~= handles.hI1
                delete(child_axRefGate(iI));
            end
        end
        
        % get image array
        dIRef = handles.h.dImg(:,:,:,1);
        set(handles.tSliceRef, 'String',[num2str(handles.slice),'/',num2str(size(dIRef,3))]);
        
        % get scroll events
        if eventdata.VerticalScrollCount < 0
            handles.slice = max([1 handles.slice - 1]);
        else
            handles.slice = min([size(dIRef, 3) handles.slice + 1]);
        end
        
        % set new image data in axes
        set(handles.hI1, 'CData', dIRef(:,:,handles.slice));
        
        % check if there are any fixed marker in this slice, and plot them
        handles = fPlotFixedMarker(handles);
        
        % if there are any unfixed points, plot them in the right slice
        if handles.unfixedpoints
            try
                if handles.slice == handles.zes.G01(1,1)
                    axes(handles.axRefGate)
                    x = handles.xes.G01;
                    y = handles.yes.G01;
                    hold on
                    handles.reflm = plot(x,y,'yx', 'Linewidth', 1.5,'MarkerSize', 5);
                    labels = cellstr( num2str([1:size(x)]'));  %' # labels correspond to their order
                    handles.reflmLabel = text(x(:,1), y(:,1), labels,'Color','y', 'VerticalAlignment','bottom', ...
                        'HorizontalAlignment','right');
                    hold off
                end
            catch
            end
        else
            try delete(handles.reflm); catch; end;
            try delete(handles.reflmLabel); catch; end;
        end
        
        % plot lines and ROIs in reference image
        try
            if isfield(handles,'hPoly') && isfield(handles.hPoly,'G01') && handles.hPolySlice.G01 == handles.slice
                handles.hPoly.G01 = impoly(handles.axRefGate, handles.xyPos.G01,'Closed', false);
                setColor(handles.hPoly.G01,'y')
            end
        catch
        end
        try
            if isfield(handles,'hPolyROI') && isfield(handles.hPolyROI,'G01') && handles.hPolyROISlice.G01 == handles.slice
                handles.hPolyROI.G01 = impoly(handles.axRefGate, handles.xyPosROI.G01,'Closed', true);
                setColor(handles.hPolyROI.G01,'y')
            end
        catch
        end
        
    else % do the same stuff only in moving image
        %         disp(num2str(gca))
        child_axMovGate = allchild(handles.axMovGate);
        for iI = 1:length(child_axMovGate)
            if child_axMovGate(iI) ~= handles.hI2
                delete(child_axMovGate(iI));
            end
        end
        
        % get image arrays and deformation fields
        dIMove2 = handles.h.dImg(:,:,:,handles.nMovGate);
        set(handles.tSliceMov, 'String',[num2str(handles.slice),'/',num2str(size(dIMove2,3))]);
        
        % get scroll events
        if eventdata.VerticalScrollCount < 0
            handles.slice = max([1 handles.slice - 1]);
        else
            handles.slice = min([size(dIMove2, 3) handles.slice + 1]);
        end
        set(handles.hI2, 'CData', dIMove2(:,:,handles.slice));
        handles = fPlotFixedMarker(handles);
        
        if handles.unfixedpoints
            try
                if handles.slice == handles.zes.(handles.names{handles.nMovGate})(1,1)
                    axes(handles.axMovGate)
                    x = handles.xes.(handles.names{handles.nMovGate});
                    y = handles.yes.(handles.names{handles.nMovGate});
                    labels = cellstr( num2str([1:size(x)]'));
                    hold on
                    handles.movlm.(handles.names{handles.nMovGate})...
                        = plot(x,y,'yx', 'Linewidth', 1.5,'MarkerSize', 5);
                    handles.movlmLabel.(handles.names{handles.nMovGate})...
                        = text(x(:,1), y(:,1), labels,'Color','y', 'VerticalAlignment','bottom', ...
                        'HorizontalAlignment','right');
                    hold off
                end
            catch
            end%display('Error: failed to plot marker in moving image while scrolling');
        else
            try delete(handles.reflm); catch; end;
            try delete(handles.reflmLabel); catch; end;
        end
        
        % plot lines and ROIs in moving image
        iI = handles.nMovGate;
        try
            if isfield(handles,'hPoly') && isfield(handles.hPoly,handles.names{iI}) && handles.hPolySlice.(handles.names{iI}) == handles.slice
                handles.hPoly.(handles.names{iI}) = impoly(handles.axMovGate, handles.xyPos.(handles.names{iI}),'Closed', false);
                setColor(handles.hPoly.(handles.names{iI}),'y')
            end
        catch
        end
        try
            if isfield(handles,'hPolyROI') && isfield(handles.hPolyROI,handles.names{iI}) && handles.hPolyROISlice.(handles.names{iI}) == handles.slice
                handles.hPolyROI.(handles.names{iI}) = impoly(handles.axMovGate, handles.xyPosROI.(handles.names{iI}),'Closed', true);
                setColor(handles.hPolyROI.(handles.names{iI}),'y')
            end
        catch
        end
        
    end
    
end

if(handles.iActive == 0)
    axes(handles.axRefGate);
    set(handles.markRefAx,'Color','white')
    set(handles.markMovAx,'Color','black')
elseif(handles.iActive == 1)
    axes(handles.axMovGate);
    set(handles.markRefAx,'Color','black')
    set(handles.markMovAx,'Color','white')
end

guidata(hObject, handles)


function LandmarkGUI_WindowButtonDownFcn(hObject, ~, handles)
% figure left mouse button click
if gca == handles.axRefGate
    set(handles.markRefAx,'Color','white')
    set(handles.markMovAx,'Color','black')
    handles.iActive = 0;
elseif gca == handles.axMovGate
    set(handles.markRefAx,'Color','black')
    set(handles.markMovAx,'Color','white')
    handles.iActive = 1;
else
    
end
% Save starting parameters
handles.dPosStart = get(gca, 'CurrentPoint');

handles.FButtonDown = 1;
guidata(hObject, handles)


function LandmarkGUI_WindowButtonUpFcn(hObject, ~, handles)
% figure mouse button release
handles.FButtonDown = 0;
handles.colMin = handles.colRange(1);
handles.colMax = handles.colRange(2);
guidata(hObject, handles)


function LandmarkGUI_WindowButtonMotionFcn(hObject, ~, handles)
% figure mouse movement
if isfield(handles,'FButtonDown')
    if handles.FButtonDown
        
        iD = get(gca, 'CurrentPoint') - handles.dPosStart; % Mouse distance travelled since button down
        
        switch get(hObject, 'SelectionType')
            
            case 'normal' % left mouse button
                % move image
                % to be implemented
                
            case 'alt' % right mouse button
                % zoom image
                % to be implemented
                
            case 'extend' % middle mouse button or 'shift' key
                % contrast and brightness
                handles.colWidth  = handles.colMax-handles.colMin;
                handles.colWidth  = handles.colMax.*exp(-iD(1,2)*0.02);
                handles.colCenter = (handles.colMax+handles.colMin)/2;
                handles.colCenter = handles.colCenter.*exp(iD(1,1)*0.02);
                handles.colRange  = [handles.colCenter-handles.colWidth/2, handles.colCenter+handles.colWidth/2];
                
                if handles.axeslink
                    caxis(handles.axMovGate,handles.colRange);
                    caxis(handles.axRefGate,handles.colRange);
                    guidata(hObject, handles);
                else
                    try caxis(gca,handles.colRange); guidata(hObject, handles); catch; end
                end
        end
    else
        return
    end
end
guidata(hObject, handles)


function LandmarkGUI_CloseRequestFcn(hObject, eventdata, handles)
% close LandmarkGUI
selection = questdlg({'Close this figure?', '(Points, Lines and ROIs that were not saved will get lost.)'},...
    'Close Figure - Save results first!',...
    'Yes (I saved the results.)','No','Yes (I saved the results.)');
switch selection,
    case 'Yes (I saved the results.)'
        % save GUIPreference
        currpath = fileparts(mfilename('fullpath'));
        if(isfield(handles,'SPaths'))
            SPaths = handles.SPaths;
            standardVoxelsize = [1 1 1];
            if(exist([currpath, filesep, 'GUIPreferences.mat'],'file'))
                save([currpath, filesep, 'GUIPreferences.mat'],'SPaths','standardVoxelsize', 'lEvalMetrics', '-append');
            else
                save([currpath, filesep, 'GUIPreferences.mat'],'SPaths','standardVoxelsize', 'lEvalMetrics');
            end
        end
        delete(gcf)
    case 'No'
        return
end


% --- Nothing shall happen, if axes background is clicked
% the following axes are only for  the white frame
function markRefAx_ButtonDownFcn(hObject, eventdata, handles)
function markMovAx_ButtonDownFcn(hObject, eventdata, handles)
function axRefGate_ButtonDownFcn(hObject, eventdata, handles)
function axMovGate_ButtonDownFcn(hObject, eventdata, handles)

% --- Executes when LandmarkGUI is resized.
function LandmarkGUI_ResizeFcn(hObject, eventdata, handles)
function tDeleteNo_Callback(hObject, eventdata, handles)
function tDeleteNo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pb_reset_contrast_Callback(hObject, ~, handles)
% reset all contrasts
dIRef = handles.h.dImg(:,:,:,1);
handles.colMax = max(dIRef(:))/handles.colScale;
handles.colMin = 0;
handles.colRange  = [handles.colMin, handles.colMax];
if handles.axeslink
    caxis(handles.axMovGate,handles.colRange);
    caxis(handles.axRefGate,handles.colRange);
else
    if gca == handles.axMovGate
        caxis(gca,handles.colRange);
    else
        caxis(gca,handles.colRange);
    end
end
guidata(hObject, handles)


function pb_clear_all_Callback(hObject, eventdata, handles)
% clear all fixed landmarkpoints

if handles.FPointsLines == 1
    isSure = questdlg('Are you sure you want to clear all fixed landmark points?', 'clear all points', 'Yes, delete all', 'Cancel','Cancel');
    switch isSure
        case 'Yes, delete all'
        case 'Cancel'
            return
        case ''
            return
    end
    
    % remove the x and y values
    try handles=rmfield(handles,'xes'); catch; end;
    try handles=rmfield(handles,'yes'); catch; end;
    try handles=rmfield(handles,'zes');
    catch
    end
    % delete the plot of the set lmps in reference and moving image
    try delete(handles.reflm);catch; end
    try delete(handles.reflmLabel);catch; end
    for iI = 2:4
        try delete(handles.movlm.(handles.names{iI}));catch; end
        try delete(handles.movlmLabel.(handles.names{iI}));catch; end
    end
    % delete the plot of the fixed lmps in reference and moving image
    for iI = 1:4
        try delete(handles.fixlm.(handles.names{iI}));catch; end
        try delete(handles.fixlmLabel.(handles.names{iI}));catch; end
        handles.fix.(handles.names{iI}) = [];
    end
    % remove the fixed values and create the empty fields
    try handles=rmfield(handles,'fix');catch;end
%     handles.fix.G01 = [];
%     handles.fix.G02 = [];
%     handles.fix.G03 = [];
%     handles.fix.G04 = [];
    
%     set(handles.textG02,'String','0');
%     set(handles.textG03,'String','0');
%     set(handles.textG04,'String','0');
    
    % update table
    cData = repmat('-',1,handles.iNGates);
    handles.cTablePoints = num2cell(cData);
    set(handles.tb_points,'Data',handles.cTablePoints);
    
elseif handles.FPointsLines == 2
    isSure = questdlg('Are you sure you want to clear all lines?', 'clear all lines and ROIs', 'Yes, delete all', 'Cancel','Cancel');
    switch isSure
        case 'Yes, delete all'
        case 'Cancel'
            return
        case ''
            return
    end
    % reset cell array for Line and ROI Polygons
    handles.cPolyCoord = rmfield(handles.cPolyCoord,handles.names);
    for iI=1:handles.iNGates
        eval(sprintf('handles.cPolyCoord.G%02d{3,1} = {};', iG));
    end
%     handles.cPolyCoord.G01{3,1} = {};
%     handles.cPolyCoord.G02{3,1} = {};
%     handles.cPolyCoord.G03{3,1} = {};
%     handles.cPolyCoord.G04{3,1} = {};
    
    handles = fPlotFixedMarker(handles);
    
    % update table
    cData = repmat('-',1,handles.iNGates);
    handles.cTableLines = num2cell(cData);
    set(handles.tb_points,'Data',handles.cTableLines);
    
elseif handles.FPointsLines == 3
    isSure = questdlg('Are you sure you want to clear all ROIs?', 'clear all lines and ROIs', 'Yes, delete all', 'Cancel','Cancel');
    switch isSure
        case 'Yes, delete all'
        case 'Cancel'
            return
        case ''
            return
    end
    % reset cell array for Line and ROI Polygons
    handles.cPolyROICoord = rmfield(handles.cPolyROICoord,handles.names);
    for iI=1:handles.iNGates
        eval(sprintf('handles.cPolyROICoord.G%02d{3,1} = {};',iG));
    end
%     handles.cPolyROICoord.G01{3,1} = {};
%     handles.cPolyROICoord.G02{3,1} = {};
%     handles.cPolyROICoord.G03{3,1} = {};
%     handles.cPolyROICoord.G04{3,1} = {};
    
    handles = fPlotFixedMarker(handles);
    
    % update table
    cData = repmat('-',1,handles.iNGates);
    handles.cTableROIs = num2cell(cData);
    set(handles.tb_points,'Data',handles.cTableROIs);
else
    error('Unknown error!')
end
eventdata.VerticalScrollCount = 0;
figure1_WindowScrollWheelFcn(hObject, eventdata, handles)

guidata(hObject, handles)


function pb_deleteNo_Callback(hObject, eventdata, handles)
% delete row number X of the uitable

% get row numbers
sDelNo = get(handles.tDeleteNo,'String');
if handles.FPointsLines == 1
    
    isSure = questdlg(['Do you want to delete landmark points in row '...
        sDelNo,'?'], 'clear corresponding points',...
        'Yes, delete!', 'Cancel','Cancel');
    switch isSure
        case 'Yes, delete!'
        case 'Cancel'
            return
        case ''
            return
    end
    
    dDelNo = sort(str2num(sDelNo),'descend');
    for iG = 1:handles.iNGates
        for iI = 1:length(dDelNo)
            try eval(sprintf('handles.fix.G%02d(dDelNo(iI),:)=[];',iG)); catch; end
%             try handles.fix.G02(dDelNo(iI),:)=[]; catch; end
%             try handles.fix.G03(dDelNo(iI),:)=[]; catch; end
%             try handles.fix.G04(dDelNo(iI),:)=[]; catch; end
        end
    end
    % display plotted points in table
    cData = repmat('x',size(handles.fix.G01,1),handles.iNGates);
    handles.cTablePoints = num2cell(cData);
    set(handles.tb_points,'Data',handles.cTablePoints);
    
    % plot the fixed landmark point in red in the reference image and moving
    % image
    % delete old landmark point plots
    child_axRefGate = allchild(handles.axRefGate);
    for iI = 1:length(child_axRefGate)
        if child_axRefGate(iI) ~= handles.hI1
            delete(child_axRefGate(iI));
        end
    end
    child_axMovGate = allchild(handles.axMovGate);
    for iI = 1:length(child_axMovGate)
        if child_axMovGate(iI) ~= handles.hI2
            delete(child_axMovGate(iI));
        end
    end
    
elseif handles.FPointsLines == 2
    
    isSure = questdlg(['Do you want to delete the lines with number '...
        sDelNo,'?'], 'clear corresponding points',...
        'Yes, delete!', 'Cancel','Cancel');
    switch isSure
        case 'Yes, delete!'
        case 'Cancel'
            return
        case ''
            return
    end
    
    dDelNo = sort(str2num(sDelNo),'descend');
    for iI = 1:length(dDelNo)
        for i = 1:size(handles.names,1);
            cPolyCoord = handles.cPolyCoord.(handles.names{i})(handles.nSCT,:);
            for j = 1:size(cPolyCoord,2)
                try
                    if dDelNo(iI) == cPolyCoord{j}(1,handles.iNGates)
                        cPolyCoord{j} = [];
                    end
                catch
                end
            end
            handles.cPolyCoord.(handles.names{i})(handles.nSCT,:) = cPolyCoord;
        end
    end
    % display plotted lines in table
    for i = 1:size(handles.names,1)
        B = cellfun('isempty',handles.cPolyCoord.(handles.names{i}));
        numLines(i) = sum(~B(:));
    end
    cData = repmat('-',max(numLines),handles.iNGates);
    for i = 1:size(handles.names,1)
        cDataTemp = repmat('L',numLines(i),1);
        cData(1:length(cDataTemp),i) = cDataTemp;
    end
    handles.cTableLines = num2cell(cData);
    set(handles.tb_points,'Data',handles.cTableLines);
    
elseif handles.FPointsLines == 3
    
    isSure = questdlg(['Do you want to delete the ROIs with number '...
        sDelNo,'?'], 'clear corresponding points',...
        'Yes, delete!', 'Cancel','Cancel');
    switch isSure
        case 'Yes, delete!'
        case 'Cancel'
            return
        case ''
            return
    end
    
    dDelNo = sort(str2num(sDelNo),'descend');
    for iI = 1:length(dDelNo)
        for i = 1:size(handles.names,1);
            cPolyCoord = handles.cPolyROICoord.(handles.names{i})(handles.nSCT,:);
            for j = 1:size(cPolyCoord,2)
                try
                    if dDelNo(iI) == cPolyCoord{j}(1,handles.iNGates)
                        cPolyCoord{j} = [];
                    end
                catch
                end
            end
            handles.cPolyROICoord.(handles.names{i})(handles.nSCT,:) = cPolyCoord;
        end
    end
    % display plotted lines in table
    for i = 1:size(handles.names,1)
        B = cellfun('isempty',handles.cPolyROICoord.(handles.names{i}));
        numLines(i) = sum(~B(:));
    end
    cData = repmat('-',max(numLines),handles.iNGates);
    for i = 1:size(handles.names,1)
        cDataTemp = repmat('R',numLines(i),1);
        cData(1:length(cDataTemp),i) = cDataTemp;
    end
    handles.cTableROIs = num2cell(cData);
    set(handles.tb_points,'Data',handles.cTableROIs);
end
eventdata.VerticalScrollCount = 0;
figure1_WindowScrollWheelFcn(hObject, eventdata, handles)
handles = fPlotFixedMarker(handles);
guidata(hObject, handles)


function pb_clear_set_Callback(hObject, ~, handles)
% clear currently set points
if handles.FPointsLines == 1
    
    isSure = questdlg('Are you sure you want to delete the set (yellow) points in all gates?', 'clear corresponding points',...
        'Yes, delete!', 'Cancel','Cancel');
    switch isSure
        case 'Yes, delete!'
        case 'Cancel'
            return
        case ''
            return
    end
    
    try handles = rmfield(handles,{'xes','yes','zes'}); catch; end;
    try delete(handles.reflm); catch; end;
    try delete(handles.reflmLabel); catch; end;
    for iG=1:handles.iNGates
        try delete(handles.movlm.(handles.names{iG})); catch; end;
        try delete(handles.movlmLabel.(handles.names{iG})); catch; end;
    end
%     try delete(handles.movlm.G03); catch; end;
%     try delete(handles.movlm.G04); catch; end;
%     try delete(handles.movlmLabel.G02); catch; end;
%     try delete(handles.movlmLabel.G03); catch; end;
%     try delete(handles.movlmLabel.G04); catch; end;
    handles.unfixedpoints = 0;
    
elseif handles.FPointsLines == 2
    
    isSure = questdlg('Are you sure you want to delete the set (yellow) lines?', 'clear lines and ROIs',...
        'Yes, delete!', 'Cancel','Cancel');
    switch isSure
        case 'Yes, delete!'
        case 'Cancel'
            return
        case ''
            return
    end
    
    for iG=1:handles.iNGates
        try delete(handles.hPoly.(handles.names{iG})); catch; end;
    end
%     try delete(handles.hPoly.G01); catch; end;
%     try delete(handles.hPoly.G02); catch; end;
%     try delete(handles.hPoly.G03); catch; end;
%     try delete(handles.hPoly.G04); catch; end;
    try handles = rmfield(handles, 'hPoly'); catch; end;
    
elseif handles.FPointsLines == 3
    
    isSure = questdlg('Are you sure you want to delete the set (yellow) ROIs?', 'clear lines and ROIs',...
        'Yes, delete!', 'Cancel','Cancel');
    switch isSure
        case 'Yes, delete!'
        case 'Cancel'
            return
        case ''
            return
    end
    
    for iG=1:handles.iNGates
        try delete(handles.hPolyROI.(handles.names{iG})); catch; end;
    end
%     try delete(handles.hPolyROI.G01); catch; end;
%     try delete(handles.hPolyROI.G02); catch; end;
%     try delete(handles.hPolyROI.G03); catch; end;
%     try delete(handles.hPolyROI.G04); catch; end;
    try handles = rmfield(handles, 'hPolyROI'); catch; end;
    
else
    error('Unknown error!')
end
guidata(hObject, handles)


function tb_points_CellSelectionCallback(hObject, eventdata, handles)
% change set landmarks
if handles.FCellSelectionActive
    handles.FCellSelectionActive = 0;
    
    if strcmp(handles.sImageData,'interpImage')
        errordlg('Landmark points can only be set in original pixel configuration.')
        return
    end
    
    % coordinates of selected cell
    coords = eventdata.Indices;
    
    if handles.FPointsLines == 1
        
        if ~isempty(coords)
            % set corresponding cell to '-'
            handles.cTablePoints{coords(1),coords(2)} = '-';
            set(handles.tb_points,'Data',handles.cTablePoints);
            
            % go to gate in which the point shall be deleted
            if coords(2) == 1
                % mark the point
                axes(handles.axRefGate)
                
                set(handles.markRefAx,'Color','white')
                set(handles.markMovAx,'Color','black')
                
                x = handles.fix.G01(coords(1),2);
                y = handles.fix.G01(coords(1),3);
                
                handles.fixlmLabel.G01...
                    = text(x, y, num2str(coords(1)),'Color','c','VerticalAlignment','bottom', ...
                    'HorizontalAlignment','right');
                handles.exchangePoints = 1;
                % activate set function
                [exchP, newSlice] = pb_set_lmp_Callback(hObject,eventdata,handles);
                if ~isempty(exchP)
                    % fix point directly
                    handles.fix.G01(coords(1),2) = exchP(1);
                    handles.fix.G01(coords(1),3) = exchP(2);
                    handles.fix.G01(coords(1),4) = newSlice;
                    guidata(hObject, handles)
                else
                    handles.exchangePoints = 0;
                    handles.fixlmLabel.G01...
                        = text(x, y, num2str(coords(1)),'Color','r','VerticalAlignment','bottom', ...
                        'HorizontalAlignment','right');
                    handles.cTablePoints{coords(1),coords(2)} = 'x';
                    set(handles.tb_points,'Data',handles.cTablePoints);
                    handles.FCellSelectionActive = 1;
                    guidata(hObject, handles)
                    return;
                end
            else
                % mark the point
                if handles.nMovGate ~= coords(2);
                    handles.nMovGate = coords(2);
                    handles = showMovImage(handles);
                end
                axes(handles.axMovGate)
                
                set(handles.markRefAx,'Color','black')
                set(handles.markMovAx,'Color','white')
                
                x = handles.fix.(handles.names{handles.nMovGate})(coords(1),2);
                y = handles.fix.(handles.names{handles.nMovGate})(coords(1),3);
                
                handles.fixlmLabel.(handles.names{handles.nMovGate})...
                    = text(x, y, num2str(coords(1)),'Color','c','VerticalAlignment','bottom', ...
                    'HorizontalAlignment','right');
                handles.exchangePoints = 1;
                % activate set function
                [exchP, newSlice] = pb_set_lmp_Callback(hObject,eventdata,handles);
                if ~isempty(exchP)
                    % fix point directly
                    handles.fix.(handles.names{handles.nMovGate})(coords(1),2) = exchP(1);
                    handles.fix.(handles.names{handles.nMovGate})(coords(1),3) = exchP(2);
                    handles.fix.(handles.names{handles.nMovGate})(coords(1),4) = newSlice;
                else
                    handles.exchangePoints = 0;
                    handles.fixlmLabel.(handles.names{handles.nMovGate})...
                        = text(x, y, num2str(coords(1)),'Color','r','VerticalAlignment','bottom', ...
                        'HorizontalAlignment','right');
                    handles.cTablePoints{coords(1),coords(2)} = 'x';
                    set(handles.tb_points,'Data',handles.cTablePoints);
                    handles.FCellSelectionActive = 1;
                    return;
                end
            end
            
            % plot the fixed landmark point in red in the reference image and moving
            % image
            % delete old landmark point plots
            child_axRefGate = allchild(handles.axRefGate);
            for iI = 1:length(child_axRefGate)
                if child_axRefGate(iI) ~= handles.hI1
                    delete(child_axRefGate(iI));
                end
            end
            child_axMovGate = allchild(handles.axMovGate);
            for iI = 1:length(child_axMovGate)
                if child_axMovGate(iI) ~= handles.hI2
                    delete(child_axMovGate(iI));
                end
            end
            handles = fPlotFixedMarker(handles);
            
            handles.cTablePoints{coords(1),coords(2)} = 'x';
            set(handles.tb_points,'Data',handles.cTablePoints);
            
        end
        handles.exchangePoints = 0;
        
    elseif handles.FPointsLines == 2
        
        %     % go to corresponding gate
        %     if coords(2) ~= 1 && handles.nMovGate ~= coords(2)
        %             handles.nMovGate = coords(2);
        %             handles = showMovImage(handles);
        %     end
        %
        %     % highlight line
        %     if coords(2) == 1
        %         set(handles.fixLineRefR,'Color','c','MarkerFaceColor','c')
        %     else
        %         set(handles.fixLineMovR,'Color','c','MarkerFaceColor','c')
        %     end
        %
        %     isSure = questdlg('Do you want to exchange the selected line?',...
        %         'Set new line', 'Yes, set new line. (Not yet working)', 'Cancel','Cancel');
        %     switch isSure
        %         case 'Yes, set new line. (Not yet working)'
        %             set(handles.fixLineRefR,'Color','r','MarkerFaceColor','r')
        %             set(handles.fixLineRefR,'Color','r','MarkerFaceColor','r')
        %             return
        %         case 'Cancel'
        %             set(handles.fixLineRefR,'Color','r','MarkerFaceColor','r')
        %             set(handles.fixLineRefR,'Color','r','MarkerFaceColor','r')
        %             return
        %         case ''
        %             set(handles.fixLineRefR,'Color','r','MarkerFaceColor','r')
        %             set(handles.fixLineRefR,'Color','r','MarkerFaceColor','r')
        %             return
        %     end
        % to be implemented: set new line and fix it to the right place in the
        % cell array (see pb_fix_lmp_Callback)
        %     handles.lineIndex = coords(1);
        %     pb_set_lmp_Callback(hObject,[],handles);
        
    elseif handles.FPointsLines == 3
        
    end
    handles.FCellSelectionActive = 1;
end
guidata(hObject, handles)


function pb_save_Callback(~, ~, handles)
% save set landmarks
if handles.nSCT ~= 1
    errordlg('Please save the landmark points from the coronal view!')
    return;
end

% if isempty(handles.fix.G01)
%     errordlg('Set and fix landmarkpoints!')
% else
    results.fixedLMP = handles.fix;
    results.cLines = handles.cPolyCoord;
    results.cROIs = handles.cPolyROICoord;
    results.SGeo.dImgSize = size(handles.h.dImg(:,:,:,1));
    results.SGeo.dVoxelsize = handles.h.SGeo.cVoxelsize;
    methods = {'elastix', 'halar', 'grics', 'lap', 'demons'};
    sMethod = methods{handles.h.nRegMethod};
    
    [sFilename,sPathname] = uiputfile([handles.SPaths.sResults,filesep,sMethod,filesep,'*.mat'],'Specify file name');
    if sFilename == 0; return; end;
    save([sPathname, sFilename],'results')
    msgbox('Saved results.', 'Save complete.','help')
% end


function pb_load_Callback(hObject, ~, handles)
% load preset landmarks
isSure = questdlg({'Loading landmarks will delete current fixed points, lines and ROIs.',' '...
    'Do you want to continue?'},...
    'Load landmarks.', 'Yes, clear fixed lines.', 'Cancel','Cancel');
switch isSure
    case 'Yes, clear fixed lines.'
    case 'Cancel'
        return
    case ''
        return
end

if handles.nSCT ~= 1
    errordlg('Please load the landmark points in the coronal view!')
    return;
end

[sFilename, sPathname] = uigetfile(handles.SPaths.sResults);
if sFilename == 0; return; end
load([sPathname, sFilename]);
if ~exist('results','var') || ~isfield(results,'fixedLMP')
    errordlg('The file does not contain the proper data structure!')
    return;
end
if all(results.SGeo.dImgSize == size(handles.h.dImg(:,:,:,1)))
    % if all dimensions agree take the results as they are
    try handles.fix = results.fixedLMP; catch; end;
    try handles.cPolyCoord = results.cLines; catch; end;
    try handles.cPolyROICoord = results.cROIs; catch; end;
else
    % transform landmark points to new grid
    % assumption: physical space of image and lm-point field is the same
    if ~isfield(results, 'SGeo')|| isempty(results.SGeo.dVoxelsize)
        voxSizePoints = standardVoxelsize;
    else
        voxSizePoints = results.SGeo.dVoxelsize;
    end
    voxSizeImage = handles.h.SGeo.dVoxelsize;
    for iI = 1:size(handles.names,1)
        xCoords = round((results.fixedLMP.(handles.names{iI})(:,2)-0.5)*voxSizePoints{iI}(1)/voxSizeImage(1)+0.5);
        yCoords = round((results.fixedLMP.(handles.names{iI})(:,3)-0.5)*voxSizePoints{iI}(2)/voxSizeImage(2)+0.5);
        zCoords = round((results.fixedLMP.(handles.names{iI})(:,4)-0.5)*voxSizePoints{iI}(3)/voxSizeImage(3)+0.5);
        fix = [results.fixedLMP.(handles.names{iI})(:,1), xCoords, yCoords, zCoords];
        handles.fix.(handles.names{iI}) = fix;
    end
    
    
    % transform also ROIs and lines
    B = ~(cellfun('isempty',results.cROIs.(handles.names{iI})));
    [a,b]=find(B);
    for iI = 1:size(handles.names,1)
        cROIs = results.cROIs.(handles.names{iI});
        for rR = 1:sum(B(:))
            dROI = cROIs{a(rR),b(rR)};
            xCoords = round((dROI(:,1)-0.5)*voxSizePoints{iI}(1)/voxSizeImage(1)+0.5);
            yCoords = round((dROI(:,2)-0.5)*voxSizePoints{iI}(2)/voxSizeImage(2)+0.5);
            zCoords = round((dROI(:,3)-0.5)*voxSizePoints{iI}(3)/voxSizeImage(3)+0.5);
            cROIs{a(rR),b(rR)} = [xCoords,yCoords,zCoords, cROIs{a(rR),b(rR)}(:,4)];
        end
        handles.cPolyROICoord.(handles.names{iI}) = cROIs;
    end
    
    B = ~(cellfun('isempty',results.cLines.(handles.names{iI})));
    [a,b]=find(B);
    for iI = 1:size(handles.names,1)
        cLines = results.cLines.(handles.names{iI});
        for lL = 1:sum(B(:))
            dLine = cLines{a(lL),b(lL)};
            xCoords = round((dLine(:,1)-0.5)*voxSizePoints{iI}(1)/voxSizeImage(1)+0.5);
            yCoords = round((dLine(:,2)-0.5)*voxSizePoints{iI}(2)/voxSizeImage(2)+0.5);
            zCoords = round((dLine(:,3)-0.5)*voxSizePoints{iI}(3)/voxSizeImage(3)+0.5);
            cLines{a(lL),b(lL)} = [xCoords,yCoords,zCoords, cLines{a(lL),b(lL)}(:,4)];
        end
        handles.cPolyCoord.(handles.names{iI}) = cLines;
    end
    
    
    
end
% display plotted lines in table
for i = 1:size(handles.names,1)
    B = cellfun('isempty',handles.cPolyCoord.(handles.names{i}));
    numLines(i) = sum(~B(:));
end
cData = repmat('-',max(numLines),handles.iNGates);
for i = 1:size(handles.names,1)
    cDataTemp = repmat('L',numLines(i),1);
    cData(1:length(cDataTemp),i) = cDataTemp;
end
handles.cTableLines = num2cell(cData);
% display plotted ROIs in table
for i = 1:size(handles.names,1)
    B = cellfun('isempty',handles.cPolyROICoord.(handles.names{i}));
    numLines(i) = sum(~B(:));
end
cData = repmat('-',max(numLines),handles.iNGates);
for i = 1:size(handles.names,1)
    cDataTemp = repmat('R',numLines(i),1);
    cData(1:length(cDataTemp),i) = cDataTemp;
end
handles.cTableROIs = num2cell(cData);
    
cData = repmat('x',size(handles.fix.G01,1),handles.iNGates);
handles.cTablePoints = num2cell(cData);
set(handles.tb_points,'Data',handles.cTablePoints);
guidata(hObject, handles)


function viewPanel_SelectionChangeFcn(hObject, eventdata, handles)
% change viewed images
dImg = handles.h.dImg(:,:,:,1);
dVoxSz = handles.h.SGeo.cVoxelsize{1}; % just scale to reference

if eventdata.NewValue == handles.rb_origImage
    handles.sImageData = 'origImage';
    handles = plotImages(handles);
elseif eventdata.NewValue == handles.rb_regImage
    handles.sImageData = 'regImage';
    handles = plotImages(handles);
elseif eventdata.NewValue == handles.rb_interpImage
    handles.sImageData = 'interpImage';
    if handles.nSCT == 2
        X1 = 1:size(dImg,1); X2 = (1:size(dImg,2))'; X3 = 1:size(dImg,3);
        Y1 = 1:size(dImg,1); Y2 = (1:dVoxSz(2)/dVoxSz(3):size(dImg,2))'; Y3 = 1:size(dImg,3);
        try
            for iI = 1:size(handles.names,1)
                xCoords = round((handles.fix.(handles.names{iI})(:,2)-0.5)*dVoxSz(3)/dVoxSz(2)+0.5);
                yCoords = round((handles.fix.(handles.names{iI})(:,3)-0.5)*1+0.5);
                zCoords = round((handles.fix.(handles.names{iI})(:,4)-0.5)*1+0.5);
                fix = [handles.fix.(handles.names{iI})(:,1), xCoords, yCoords, zCoords];
                handles.fixInterp.(handles.names{iI}) = fix;
            end
        catch; disp('There are no fixed points')
        end
        
        % to do:
        % if there are any lines or ROIs:
        try
            
        catch
        end
        
        
    elseif handles.nSCT == 3
        X1 = (1:size(dImg,1))'; X2 = 1:size(dImg,2); X3 = 1:size(dImg,3);
        Y1 = (1:dVoxSz(2)/dVoxSz(3):size(dImg,1))'; Y2 = 1:size(dImg,2); Y3 = 1:size(dImg,3);
        try
            for iI = 1:size(handles.names,1)
                xCoords = round((handles.fix.(handles.names{iI})(:,2)-0.5)*1+0.5);
                yCoords = round((handles.fix.(handles.names{iI})(:,3)-0.5)*dVoxSz(3)/dVoxSz(2)+0.5);
                zCoords = round((handles.fix.(handles.names{iI})(:,4)-0.5)*1+0.5);
                fix = [handles.fix.(handles.names{iI})(:,1), xCoords, yCoords, zCoords];
                handles.fixInterp.(handles.names{iI}) = fix;
            end
        catch; disp('There are no fixed points')
        end
        
        % to do:
        % if there are any lines or ROIs:
        try
            
        catch
        end
        
    else
        set(hObject, 'Value', 0);
        set(eventdata.OldValue,'Value',1)
        return;
    end
    h = waitbar(0,'Please wait...'); st=1;steps = 4;
    try handles = rmfield(handles,'dImgInterp');catch; end
    
    for iI = 1:size(handles.h.dImg,4)
        handles.dImgInterp(:,:,:,iI) = ...
            interpn(X1,X2,X3,handles.h.dImg(:,:,:,iI),Y1,Y2,Y3); waitbar(st/steps);st=st+1;
    end
    try close(h); catch; end;
    
    handles.colRange = [0 handles.colMax];
    slice = round(handles.slice);
    dImg1 = handles.dImgInterp(:,:,:,1);
    dIMove2 = handles.dImgInterp(:,:,:,handles.nMovGate);
    
    % plot reference image
    axes(handles.axRefGate)
    cla
    handles.hI1 = imshow(dImg1(:,:,slice),'DisplayRange',handles.colRange,'Parent', handles.axRefGate);
    
    % plot moving image
    axes(handles.axMovGate)
    cla
    handles.hI2 = imshow(dIMove2(:,:,slice),'DisplayRange',handles.colRange,'Parent', handles.axMovGate);
    
    % remove the landmark points that are set (not fixed) and delete the plot
    try
        handles=rmfield(handles,'xes');
        handles=rmfield(handles,'yes');
    catch
    end
    try delete(handles.reflm); catch; end;
    try delete(handles.reflmLabel); catch; end;
    try handles=rmfield(handles,'reflm');catch; end
    try handles=rmfield(handles,'reflmLabel');catch; end
    for iI = 1:handles.iNGates
        try handles=rmfield(handles,movlm.(handles.names{iI}));catch; end
        try handles=rmfield(handles,movlmLabel.(handles.names{iI}));catch; end
    end
end
handles = fPlotFixedMarker(handles);
guidata(hObject, handles)


function pm_CorSagTrans_Callback(hObject, ~, handles)
% change image orientation (cor, sag, tra)
if handles.unfixedlines == 1
    isSure = questdlg('Changing orientation deletes unfixed lines and ROIs. Do you want to continue?',...
        'Change orientation', 'Yes, clear set lines.', 'Cancel','Cancel');
    switch isSure
        case 'Yes, clear set lines.'
        case 'Cancel'
            return
        case ''
            return
    end
end
handles.unfixedlines = 0;
contents = cellstr(get(handles.pm_CorSagTrans,'String'));
lead=contents{get(handles.pm_CorSagTrans,'value')};
% delete old landmark point plots
try delete(handles.reflm);catch; end
try delete(handles.reflmLabel);catch; end
try delete(handles.movlm.(handles.names{handles.nMovGate}));catch; end
try delete(handles.movlmLabel.(handles.names{handles.nMovGate}));catch; end
try delete(handles.fixlm.G01);catch; end
try delete(handles.fixlmLabel.G01);catch; end
try delete(handles.fixlm.(handles.names{handles.nMovGate}));catch; end
try delete(handles.fixlmLabel.(handles.names{handles.nMovGate}));catch; end

% get image array and fixed points
dImg = handles.h.dImg;
dImgReg = handles.h.dImgReg;
if ~isempty(handles.fix.G01)
    for iG=1:handles.iNGates
        eval(sprintf('fixedG%02d = handles.fix.(handles.names{iG})(:,2:4);', iG));
    end
%     fixedG01 = handles.fix.G01(:,2:4);
%     fixedG02 = handles.fix.G02(:,2:4);
%     fixedG03 = handles.fix.G03(:,2:4);
%     fixedG04 = handles.fix.G04(:,2:4);
end
% permutation vector for orientation changes
CorSag = [1,3,2,4];
CorTra = [3,2,1,4];
SagCor = [1,3,2,4];
SagTra = [2,3,1,4];
TraCor = [3,2,1,4];
TraSag = [3,1,2,4];

switch lead                                             % for all combinations of orientation changes:
    case 'coronal'
        nSCT=1;                                         % 1. set: - flag for slice orientation
        if handles.nSCT == 2 % sag -> cor
            dImgRot = permute(dImg,SagCor);             % 2. a) permute image data according to permutation vector
            dImgRegRot = permute(dImgReg,SagCor);
            i = 1; j = 3;
            for iG=1:handles.iNGates
                try
                    eval(sprintf('fixedG%02d(:,[i,j]) = fixedG%02d(:,[j,i]);', iG, iG)); %    b) swap coordinates for the fixed landmark points
                    handles.SDeform(iG).dBx = permute(handles.SDeform(iG).dBx,SagCor);
                    handles.SDeform(iG).dBy = permute(handles.SDeform(iG).dBy,SagCor);
                    handles.SDeform(iG).dBz = permute(handles.SDeform(iG).dBz,SagCor);
%                     fixedG01(:,[i,j])= fixedG01(:,[j,i]);       
%                     fixedG02(:,[i,j])= fixedG02(:,[j,i]);
%                     fixedG03(:,[i,j])= fixedG03(:,[j,i]);
%                     fixedG04(:,[i,j])= fixedG04(:,[j,i]);
                catch
                end
            end
        elseif handles.nSCT == 3 % tra -> cor
            dImgRot = permute(dImg,TraCor);
            dImgRegRot = permute(dImgReg,TraCor);
            i = 2; j = 3;
            for iG=1:handles.iNGates
                try
                    eval(sprintf('fixedG%02d(:,[i,j]) = fixedG%02d(:,[j,i]);', iG, iG)); %    b) swap coordinates for the fixed landmark points
                    handles.SDeform(iG).dBx = permute(handles.SDeform(iG).dBx,TraCor);
                    handles.SDeform(iG).dBy = permute(handles.SDeform(iG).dBy,TraCor);
                    handles.SDeform(iG).dBz = permute(handles.SDeform(iG).dBz,TraCor);
%                     fixedG01(:,[i,j])= fixedG01(:,[j,i]);
%                     fixedG02(:,[i,j])= fixedG02(:,[j,i]);
%                     fixedG03(:,[i,j])= fixedG03(:,[j,i]);
%                     fixedG04(:,[i,j])= fixedG04(:,[j,i]);
                catch
                end
            end
        else return                                     % 3. do nothing, if slice orientation did not change
        end
    case 'sagittal'
        nSCT=2;
        if handles.nSCT == 1 % cor -> sag
            dImgRot = permute(dImg,CorSag);
            dImgRegRot = permute(dImgReg,CorSag);
            i = 1; j = 3;
            for iG=1:handles.iNGates
                try
                    eval(sprintf('fixedG%02d(:,[i,j]) = fixedG%02d(:,[j,i]);', iG, iG)); %    b) swap coordinates for the fixed landmark points
                    handles.SDeform(iG).dBx = permute(handles.SDeform(iG).dBx,CorSag);
                    handles.SDeform(iG).dBy = permute(handles.SDeform(iG).dBy,CorSag);
                    handles.SDeform(iG).dBz = permute(handles.SDeform(iG).dBz,CorSag);
%                     fixedG01(:,[i,j])= fixedG01(:,[j,i]);
%                     fixedG02(:,[i,j])= fixedG02(:,[j,i]);
%                     fixedG03(:,[i,j])= fixedG03(:,[j,i]);
%                     fixedG04(:,[i,j])= fixedG04(:,[j,i]);
                catch
                end
            end
        elseif handles.nSCT == 3 % tra -> sag
                dImgRot = permute(dImg,TraSag);
                dImgRegRot = permute(dImgReg,TraSag);
                for iG=1:handles.iNGates
                    try
                        eval(sprintf('fixedG%02d = circshift(fixedG%02d,[0,-1]);', iG, iG)); %    b) swap coordinates for the fixed landmark points
                        handles.SDeform(iG).dBx = permute(handles.SDeform(iG).dBx,TraSag);
                        handles.SDeform(iG).dBy = permute(handles.SDeform(iG).dBy,TraSag);
                        handles.SDeform(iG).dBz = permute(handles.SDeform(iG).dBz,TraSag);

%                         fixedG01 = circshift(fixedG01,[0,-1]);
%                         fixedG02 = circshift(fixedG02,[0,-1]);
%                         fixedG03 = circshift(fixedG03,[0,-1]);
%                         fixedG04 = circshift(fixedG04,[0,-1]);
                    catch
                    end
                end
        else return
        end
    case 'transverse'
        nSCT=3;
        if handles.nSCT == 1 % cor -> tra
            dImgRot = permute(dImg,CorTra);
            dImgRegRot = permute(dImgReg,CorTra);
            i = 2; j = 3;
            for iG=1:handles.iNGates
                try
                    eval(sprintf('fixedG%02d(:,[i,j]) = fixedG%02d(:,[j,i]);', iG, iG)); %    b) swap coordinates for the fixed landmark points
                    handles.SDeform(iG).dBx = permute(handles.SDeform(iG).dBx,CorTra);
                    handles.SDeform(iG).dBy = permute(handles.SDeform(iG).dBy,CorTra);
                    handles.SDeform(iG).dBz = permute(handles.SDeform(iG).dBz,CorTra);

%                     fixedG01(:,[i,j])= fixedG01(:,[j,i]);
%                     fixedG02(:,[i,j])= fixedG02(:,[j,i]);
%                     fixedG03(:,[i,j])= fixedG03(:,[j,i]);
%                     fixedG04(:,[i,j])= fixedG04(:,[j,i]);
                catch
                end
            end
        elseif handles.nSCT == 2 % sag -> tra
            dImgRot = permute(dImg,SagTra);
            dImgRegRot = permute(dImgReg,SagTra);
            for iG=1:handles.iNGates
                try
                    eval(sprintf('fixedG%02d = circshift(fixedG%02d,[0,1]);', iG, iG)); %    b) swap coordinates for the fixed landmark points
                    handles.SDeform(iG).dBx = permute(handles.SDeform(iG).dBx,SagTra);
                    handles.SDeform(iG).dBy = permute(handles.SDeform(iG).dBy,SagTra);
                    handles.SDeform(iG).dBz = permute(handles.SDeform(iG).dBz,SagTra);

%                     fixedG01 = circshift(fixedG01,[0,1]);
%                     fixedG02 = circshift(fixedG02,[0,1]);
%                     fixedG03 = circshift(fixedG03,[0,1]);
%                     fixedG04 = circshift(fixedG04,[0,1]);
                catch
                end
            end
        else return
        end
    otherwise
        error('unexpected error');
end
% set new values in handle structure.
handles.nSCT = nSCT;
handles.lead = lead;
handles.h.dImg = dImgRot;
handles.h.dImgReg = dImgRegRot;
handles.slice = floor(size(handles.h.dImg,3)/2);
try
    for iG=1:handles.iNGates
        eval(sprintf('handles.fix.(handles.names{iG})(:,2:4) = fixedG%02d;', iG));
    end
%     handles.fix.G01(:,2:4) = fixedG01;
%     handles.fix.G02(:,2:4) = fixedG02;
%     handles.fix.G03(:,2:4) = fixedG03;
%     handles.fix.G04(:,2:4) = fixedG04;
catch
end
set(handles.rb_interpImage,'Value',0);
set(handles.rb_origImage,'Value',1);
handles.sImageData = 'origImage';
handles = plotImages(handles);
handles = fPlotFixedMarker(handles);
guidata(hObject, handles)


function SOut = fSelectStruct(sResultpath)
% select to be loaded registration result
[sFilename, sPathname] = uigetfile(sResultpath);
if sFilename == 0
    SOut = 0;
    return
end
h = waitbar(0,'Load image data ...');
load([sPathname, sFilename]); waitbar(1)
try close(h); catch; end;
if exist('SRegiResult','var')
    SOut = SRegiResult;
else
    SOut=0;
    return
end


function handles = plotImages(handles)
% plot reference and moving images

% get images
dIRef = handles.h.dImg(:,:,:,1);
if strcmp(handles.sImageData, 'origImage')
    dIMove2 = handles.h.dImg(:,:,:,handles.nMovGate);
elseif strcmp(handles.sImageData, 'regImage')
    dIMove2 = handles.h.dImgReg(:,:,:,handles.nMovGate);
end
dAspectRatio = handles.h.SGeo.cVoxelsize{1};

% determine range of colours
handles.colMax = max(dIRef(:))/handles.colScale;
handles.colMin = 0;
handles.colCenter = handles.colMax/2;
handles.colWidth = handles.colMax - handles.colMin;
handles.colRange = [0 handles.colMax];

% get slice of interest
slice = handles.slice;%round(size(dIMove2,3)/2);

% plot reference image
axes(handles.axRefGate)
if(~isfield(handles,'hI1') || isempty(handles.hI1))
    cla
    handles.hI1 = imshow(dIRef(:,:,slice),'DisplayRange',handles.colRange,'Parent', handles.axRefGate);
else
    set(handles.hI1, 'CData', dIRef(:,:,slice));
    set(handles.axRefGate, 'CLim', handles.colRange);
end
if strcmp(handles.sImageData,'origImage')||strcmp(handles.sImageData,'regImage')
    if handles.nSCT == 2
        daspect([dAspectRatio(2) dAspectRatio(3) 1]);
    elseif handles.nSCT == 3
        daspect([dAspectRatio(3) dAspectRatio(1) 1]);
    else
        daspect([dAspectRatio(1) dAspectRatio(2) 1]);
    end
end
xlim(handles.axRefGate,[0.5 size(dIRef,2)+0.5]);
ylim(handles.axRefGate,[0.5 size(dIRef,1)+0.5]);
colormap(handles.axRefGate,'gray');
freezeColors;

% plot moving image
dAspectRatio = handles.h.SGeo.cVoxelsize{handles.nMovGate};
axes(handles.axMovGate)
if(~isfield(handles,'hI2') || isempty(handles.hI2))
    cla
    handles.hI2 = imshow(dIMove2(:,:,slice),'DisplayRange',handles.colRange,'Parent', handles.axMovGate);
else
    set(handles.hI2, 'CData', dIMove2(:,:,slice));
    set(handles.axMovGate, 'CLim', handles.colRange);
end
if strcmp(handles.sImageData,'origImage')||strcmp(handles.sImageData,'regImage')
    if handles.nSCT == 2
        daspect([dAspectRatio(2) dAspectRatio(3) 1]);
    elseif handles.nSCT == 3
        daspect([dAspectRatio(3) dAspectRatio(1) 1]);
    else
        daspect([dAspectRatio(1) dAspectRatio(2) 1]);
    end
end
xlim(handles.axMovGate,[0.5 size(dIMove2,2)+0.5]);
ylim(handles.axMovGate,[0.5 size(dIMove2,1)+0.5]);
linkaxes([handles.axRefGate, handles.axMovGate])
set(handles.CBlink, 'Value',1);
handles.axeslink = 1;


function handles = showMovImage(handles)
% show moving images
dAspectRatio = handles.h.SGeo.cVoxelsize{handles.nMovGate};
if strcmp(handles.sImageData,'origImage')
    dIMove2 = handles.h.dImg(:,:,:,handles.nMovGate);
elseif strcmp(handles.sImageData,'regImage')
    dIMove2 = handles.h.dImgReg(:,:,:,handles.nMovGate);
elseif strcmp(handles.sImageData,'origImage')
    dIMove2 = handles.dImgInterp(:,:,:,handles.nMovGate);
else
    error('Unknown Error!')
end
% plot moving image
axes(handles.axMovGate)
if(~isfield(handles,'hI2') || isempty(handles.hI2))
    cla
    handles.hI2 = imshow(dIMove2(:,:,handles.slice),handles.colRange);
else
    set(handles.hI2, 'CData', dIMove2(:,:,handles.slice));
    set(handles.axMovGate, 'CLim', handles.colRange);
end
if strcmp(handles.sImageData,'origImage')||strcmp(handles.sImageData,'regImage')
    if handles.nSCT == 2
        daspect([dAspectRatio(2) dAspectRatio(3) 1]);
    elseif handles.nSCT == 3
        daspect([dAspectRatio(3) dAspectRatio(1) 1]);
    else
        daspect([dAspectRatio(1) dAspectRatio(2) 1]);
    end
end
xlim(handles.axMovGate,[0.5 size(dIMove2,2)+0.5]);
ylim(handles.axMovGate,[0.5 size(dIMove2,1)+0.5]);
set(handles.tGateMov, 'String',['Image 0',num2str(handles.nMovGate)]);
% if there are some set landmark points, plot them
try
    if handles.slice == handles.zes.(handles.names{handles.nMovGate})(1,1)
        axes(handles.axMovGate)
        x = handles.xes.(handles.names{handles.nMovGate});
        y = handles.yes.(handles.names{handles.nMovGate});
        labels = cellstr( num2str([1:size(x)]'));  %' # labels correspond to their order
        hold on
        handles.movlm.(handles.names{handles.nMovGate})...
            = plot(x,y,'yx', 'Linewidth', 1.5,'MarkerSize', 5);
        handles.movlmLabel.(handles.names{handles.nMovGate})...
            = text(x(:,1), y(:,1), labels,'Color','y', 'VerticalAlignment','bottom', ...
            'HorizontalAlignment','right');
        hold off
    end
catch %display('Error: failed to plot marker in moving image while scrolling');
end

handles = fPlotFixedMarker(handles);


function handles = fPlotFixedMarker(handles)
% plot the fixed landmark point in red in the reference image and moving
% image
for iI = 1:handles.iNGates
    try delete(handles.fixlm.(handles.names{iI})); delete(handles.fixlmLabel.(handles.names{iI})); catch; end;
end

if gca == handles.axRefGate || handles.axeslink == 1
    try delete(handles.fixLineRefW,handles.fixLineRefR); catch; end;
    
    if strcmp(handles.sImageData,'origImage')||strcmp(handles.sImageData,'regImage')
        % if there are some fixed landmark points in the current slice, plot them
        try R = (handles.fix.G01(:,4) == handles.slice); catch; end;
        try
            axes(handles.axRefGate)
            x = handles.fix.G01(:,2).*R;
            y = handles.fix.G01(:,3).*R;
            labels = cellstr(num2str([1:size(x)]'));
            labels(x==0)=[];x(x==0)=[];y(y==0)=[];
            hold on
            handles.fixlm.G01...
                = plot(x,y,'rx', 'Linewidth', 1.5,'MarkerSize', 5);
            handles.fixlmLabel.G01...
                = text(x(:,1), y(:,1), labels,'Color','r', 'VerticalAlignment','bottom', ...
                'HorizontalAlignment','right');
            hold off
        catch
        end
    elseif strcmp(handles.sImageData,'interpImage')
        % if there are some fixed landmark points, plot them
        try R = (handles.fixInterp.G01(:,4) == handles.slice); catch; end;
        try
            axes(handles.axRefGate)
            x = handles.fixInterp.G01(:,2).*R;
            y = handles.fixInterp.G01(:,3).*R;
            labels = cellstr(num2str([1:size(x)]'));
            labels(x==0)=[];x(x==0)=[];y(y==0)=[];
            hold on
            handles.fixlm.G01...
                = plot(x,y,'rx', 'Linewidth', 1.5,'MarkerSize', 5);
            handles.fixlmLabel.G01...
                = text(x(:,1), y(:,1), labels,'Color','r', 'VerticalAlignment','bottom', ...
                'HorizontalAlignment','right');
            hold off
        catch
        end
    end
    for i = 1:size(handles.cPolyCoord.G01,2)
        % for sagittal and transverse: permute coordinates back!
        xyzPoly = handles.cPolyCoord.G01{handles.nSCT,i};
        if handles.nSCT == 2
            try l = 1; m = 3; xyzPoly(:,[l,m])= xyzPoly(:,[m,l]); catch; end;
        elseif handles.nSCT == 3
            try l = 2; m = 3; xyzPoly(:,[l,m])= xyzPoly(:,[m,l]); catch; end;
        end
        try
            if xyzPoly(1,3) == handles.slice
                x = xyzPoly(:,1);
                y = xyzPoly(:,2);
                lineNo = num2str(xyzPoly(1,4));
                hold on
                handles.fixLineRefW = plot(x,y,'w','LineWidth',3);
                handles.fixLineRefR = plot(x,y,'or-','LineWidth',1.5,...
                    'MarkerFaceColor','r','MarkerEdgeColor','w');
                handles.fixLineLabel = text(x(1)+2,y(1)+2,lineNo,'Color','r');
                hold off
            end
        catch
        end
    end
    for i = 1:size(handles.cPolyROICoord.G01,2)
        % for sagittal and transverse: permute coordinates back!
        xyzPoly = handles.cPolyROICoord.G01{handles.nSCT,i};
        if handles.nSCT == 2
            try l = 1; m = 3; xyzPoly(:,[l,m])= xyzPoly(:,[m,l]); catch; end;
        elseif handles.nSCT == 3
            try l = 2; m = 3; xyzPoly(:,[l,m])= xyzPoly(:,[m,l]); catch; end;
        end
        try
            if xyzPoly(1,3) == handles.slice
                x = xyzPoly(:,1);
                y = xyzPoly(:,2);
                lineNo = num2str(xyzPoly(1,4));
                hold on
                handles.fixROIRefW = plot(x,y,'w','LineWidth',3);
                handles.fixROIRefR = plot(x,y,'or-','LineWidth',1.5,...
                    'MarkerFaceColor','r','MarkerEdgeColor','w');
                handles.fixROILabel = text(x(1)+2,y(1)+2,lineNo,'Color','r');
                hold off
            end
        catch
        end
    end
end
if gca == handles.axMovGate || handles.axeslink == 1
    try delete(handles.fixLineMovW,handles.fixLineMovR); catch; end;
    % plot also marker in the moving image
    if strcmp(handles.sImageData,'origImage')||strcmp(handles.sImageData,'regImage')
        % if there are some fixed landmark points, plot them
        try R = (handles.fix.(handles.names{handles.nMovGate})(:,4) == handles.slice); catch; end;
        try
            axes(handles.axMovGate)
            x = handles.fix.(handles.names{handles.nMovGate})(:,2).*R;
            y = handles.fix.(handles.names{handles.nMovGate})(:,3).*R;
            labels = cellstr(num2str([1:size(x)]'));
            labels(x==0)=[];x(x==0)=[];y(y==0)=[];%
            hold on
            handles.fixlm.(handles.names{handles.nMovGate})...
                = plot(x,y,'rx', 'Linewidth', 1.5,'MarkerSize', 5);
            handles.fixlmLabel.(handles.names{handles.nMovGate})...
                = text(x(:,1), y(:,1), labels,'Color','r', 'VerticalAlignment','bottom', ...
                'HorizontalAlignment','right');
            hold off
        catch
        end
    elseif strcmp(handles.sImageData,'interpImage')
        % if there are some fixed landmark points, plot them
        try R = (handles.fixInterp.(handles.names{handles.nMovGate})(:,4) == handles.slice); catch; end;
        try
            axes(handles.axMovGate)
            x = handles.fixInterp.(handles.names{handles.nMovGate})(:,2).*R;
            y = handles.fixInterp.(handles.names{handles.nMovGate})(:,3).*R;
            labels = cellstr(num2str([1:size(x)]'));
            labels(x==0)=[];x(x==0)=[];y(y==0)=[];
            hold on
            handles.fixlm.(handles.names{handles.nMovGate})...
                = plot(x,y,'rx', 'Linewidth', 1.5,'MarkerSize', 5);
            handles.fixlmLabel.(handles.names{handles.nMovGate})...
                = text(x(:,1), y(:,1), labels,'Color','r', 'VerticalAlignment','bottom', ...
                'HorizontalAlignment','right');
            hold off
        catch
        end
    end
    iI = handles.nMovGate;
    for i = 1:size(handles.cPolyCoord.(handles.names{iI}),2)
        % for sagittal and transverse: permute coordinates back!
        xyzPoly = handles.cPolyCoord.(handles.names{iI}){handles.nSCT,i};
        if handles.nSCT == 2
            try l = 1; m = 3; xyzPoly(:,[l,m])= xyzPoly(:,[m,l]); catch; end;
        elseif handles.nSCT == 3
            try l = 2; m = 3; xyzPoly(:,[l,m])= xyzPoly(:,[m,l]); catch; end;
        end
        try
            if xyzPoly(1,3) == handles.slice
                x = xyzPoly(:,1);
                y = xyzPoly(:,2);
                lineNo = num2str(xyzPoly(1,4));
                hold on
                handles.fixLineMovW = plot(x,y,'w','LineWidth',3);
                handles.fixLineMovR = plot(x,y,'or-','LineWidth',1.5,...
                    'MarkerFaceColor','r','MarkerEdgeColor','w');
                handles.fixLineLabel = text(x(1)+2,y(1)+2,lineNo,'Color','r');
                hold off
            end
        catch
        end
    end
    iI = handles.nMovGate;
    for i = 1:size(handles.cPolyROICoord.(handles.names{iI}),2)
        % for sagittal and transverse: permute coordinates back!
        xyzPoly = handles.cPolyROICoord.(handles.names{iI}){handles.nSCT,i};
        if handles.nSCT == 2
            try l = 1; m = 3; xyzPoly(:,[l,m])= xyzPoly(:,[m,l]); catch; end;
        elseif handles.nSCT == 3
            try l = 2; m = 3; xyzPoly(:,[l,m])= xyzPoly(:,[m,l]); catch; end;
        end
        try
            if xyzPoly(1,3) == handles.slice
                x = xyzPoly(:,1);
                y = xyzPoly(:,2);
                lineNo = num2str(xyzPoly(1,4));
                hold on
                handles.fixROIMovW = plot(x,y,'w','LineWidth',3);
                handles.fixROIMovR = plot(x,y,'or-','LineWidth',1.5,...
                    'MarkerFaceColor','r','MarkerEdgeColor','w');
                handles.fixROILabel = text(x(1)+2,y(1)+2,lineNo,'Color','r');
                hold off
            end
        catch
        end
    end
    
%     set(handles.markRefAx,'Color','black')
%     set(handles.markMovAx,'Color','white')
end


% change settings for similarity metrics
function pbSettingsMetrics_Callback(hObject, eventdata, handles)
% open small window to decide similarity metrics
cShow = cell(size(handles.cMetrics,2),1);
for iI=1:length(cShow)
    cShow{iI} = [handles.cMetrics{2,iI}, ' (', upper(handles.cMetrics{1,iI}),')'];
end
iInd = fSmallPopup(cShow,'Choose intensity-based similarity metrics', length(cShow), find(handles.lEvalMetrics));
if(isempty(iInd))
    return;
end
handles.lEvalMetrics = false(1,size(handles.cMetrics,2));
handles.lEvalMetrics(iInd) = true;
guidata(hObject, handles);
