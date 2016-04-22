function varargout = EvalGUI(varargin)
% GUI to evaluate the registration results
%
% opens from RegiGUI with the current registration results or from the
% command window, then, registration mat-file has to be loaded
%
% For detailed documentation see "Documentation of GUIs for MR image
% registration" (pdf).
%
% Before the first start check the path in th EvalGUI_OpeningFcnResult
%
% -------------------------------------------------------------------------
% (c) 2015: Thomas Kuestner, Verena Neumann
% -------------------------------------------------------------------------

% Last Modified by GUIDE v2.5 06-Dec-2015 14:26:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @EvalGUI_OpeningFcn, ...
    'gui_OutputFcn',  @EvalGUI_OutputFcn, ...
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

% -------------------------------------------------------------------------
% basic GUI functions
% -------------------------------------------------------------------------

% --- Executes just before EvalGUI is made visible.
function EvalGUI_OpeningFcn(hObject, eventdata, handles, varargin)
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
            case 'egmode'
                handles.EGmode = varargin{index+1};
            case 'path'
                h = fwaitbar(0,'Loading image data ...');
                [sPathname, sFilename, sExt] = fileparts(varargin{index+1});
                load([sPathname, filesep, sFilename, sExt]); fwaitbar(1,h);
                try close(h); catch; end;
                handles.h = SRegiResult;
                if isempty(sPathname)
                    handles.closeFigure = true;
                    % Update handles structure
                    guidata(hObject, handles);
                    return;
                end
                handles.h.sFilename = sFilename;
                handles.SPaths.sData = sPathname;
        end
    end
end
% when opening without data, select a data struct from the hard drive
% originating from registration with RegiGUI
if ~isfield(handles, 'h') || isempty(handles.h)
    [handles.h, sPathname, sFilename] = fSelectStruct(handles);
    if isempty(sPathname)
        handles.closeFigure = true;
        % Update handles structure
        guidata(hObject, handles);
        error('Missing data for viewing!');
    end
    handles.h.sFilename = sFilename;
    handles.SPaths.sData = sPathname;
end

%% set icons
try % Try to apply a nice icon
    warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
    jframe = get(hObject, 'javaframe');
    jIcon = javax.swing.ImageIcon([currpath, filesep, 'icons', filesep, 'EvalGUI_icon.png']);
    pause(0.001);
    jframe.setFigureIcon(jIcon);
    clear jframe jIcon
catch
    warning('Could not apply a nice icon to the figure :(');
end
dImage = double(imread(['icons',filesep,'layers.png']))./255;
if size(dImage, 3) == 1, dImage = repmat(dImage, [1 1 3]); end
set(handles.pb_DispGates, 'CData', dImage/max(dImage(:)));

dImage = double(imread(['icons',filesep,'reset.png']))./255;
if size(dImage, 3) == 1, dImage = repmat(dImage, [1 1 3]); end
set(handles.pb_resetContrast, 'CData', dImage/max(dImage(:)));

dImage = double(imread(['icons',filesep,'save.png']))./255;
if size(dImage, 3) == 1, dImage = repmat(dImage, [1 1 3]); end
set(handles.pb_save, 'CData', dImage/max(dImage(:)));

dImage = double(imread(['icons',filesep,'doc_export.png']))./255;
if size(dImage, 3) == 1, dImage = repmat(dImage, [1 1 3]); end
set(handles.pb_SaveAs, 'CData', dImage/max(dImage(:)));

%% check inputs
if ~isfield(handles.h, 'SGeo') || isempty(handles.h.SGeo.cVoxelsize)
    handles.h.SGeo.cVoxelsize{1} = standardVoxelsize;
    handles.h.SGeo.cOrientation{1} = 'coronal';
end
if(~isfield(handles.h,'nSCT'))
    if(~isfield(handles.h.SGeo, 'cOrientation'))
        handles.nSCT = 1;           % slice orientation (default: coronal)
    else
        cont = {'coronal', 'sagittal', 'transverse'};
        handles.nSCT = find(strcmp(cont,handles.h.SGeo.cOrientation{1}));
    end
else
    handles.nSCT = handles.h.nSCT;
end
set(handles.pop_SagCorTrans,'Value',handles.nSCT)

if(~isfield(handles.h,'sRegMethod'))
    handles.h.sRegMethod = 'N/A';
end
if(~isfield(handles.h, 'sParFile'))
    handles.h.sParFile = 'N/A';
end
if(~isfield(handles.h, 'sFilename'))
    handles.h.sFilename = 'default';
end
if(~isfield(handles.h,'sShownames'))
    handles.h.sShownames = cell(1,size(handles.h.dImg,4));
end
if(~isfield(handles.h,'nRegMethod'))
    handles.h.nRegMethod = 1;
end

handles.iActualGates = size(handles.h.dImg,4);
if(handles.iActualGates < 4) % if input has less than four gates
    % pad with zeros
    handles.h.dImg = cat(4,handles.h.dImg,zeros(size(handles.h.dImg,1),size(handles.h.dImg,2),size(handles.h.dImg,3),4-size(handles.h.dImg,4)));
    handles.h.dImgReg = cat(4,handles.h.dImg,zeros(size(handles.h.dImg,1),size(handles.h.dImg,2),size(handles.h.dImg,3),4-size(handles.h.dImg,4)));
    for iI=size(handles.h.dImg,4)+1:4
        handles.h.SDeform(iI).dBx = zeros(size(handles.h.SDeform(2).dBx));
        handles.h.SDeform(iI).dBy = zeros(size(handles.h.SDeform(2).dBy));
        handles.h.SDeform(iI).dBz = zeros(size(handles.h.SDeform(2).dBz));
        handles.h.SDeform(iI).dFx = zeros(size(handles.h.SDeform(2).dFx));
        handles.h.SDeform(iI).dFy = zeros(size(handles.h.SDeform(2).dFy));
        handles.h.SDeform(iI).dFz = zeros(size(handles.h.SDeform(2).dFz));
    end
end


%% Choose default command line output for EvalGUI
handles.output = hObject;
handles.slice = floor(size(handles.h.dImg,3)/2);         % slice number
handles.quiverScale = 1;
handles.quiverFactor = 8;
handles.dStep = [0.25, 1]; % quiverScale, quiverFactor
handles.level_set = 2;
handles.level_set_incr = 0.25;
handles.tcmap=[ 0 0 0; 7 0 17; 14 0 33; 21 0 50; 29 0 67; 36 0 84; 43 0 100; 50 0 117;
57 0 134; 64 0 150; 72 0 167; 80 3 164; 89 7 156; 97 11 149; 106 15 142; 114 19 134;
123 23 127; 131 27 119; 140 31 112; 149 35 105; 157 39 97; 166 43 90; 174 47 82;
183 51 75; 192 55 68; 200 59 60; 209 63 53; 217 67 45; 226 71 38; 234 75 31;
243 79 23; 252 83 16; 255 88 12; 255 95 12; 255 102 11; 255 109 11; 255 116 10;
255 123 10; 255 130 9; 255 137 9; 255 144 8; 255 151 8; 255 158 7; 255 165 7;
255 172 6; 255 179 6; 255 186 5; 255 193 4; 255 200 4; 255 207 3; 255 214 3; 255 221 2;
255 228 2; 255 235 1; 255 242 1; 255 249 0; 255 252 22; 255 252 55; 255 253 88;
255 253 122; 255 254 155; 255 254 188; 255 255 222; 255 255 255]/255;
handles.dView = [-30,60];
handles.dCamlight = [-45,45];

handles.jacDiv = 0;         % 'gradient' method to measure change of defField
handles.sImageData='origImage'; % start with original image data (can be changed to registered images)
handles.FButtonDown = 0;    % start with no pressed button
handles.newGates = 0;       % the gates have not been changed yet
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
if(~isfield(handles.h.SDeform(2),'dBx'))
    set(handles.rbBackward, 'Visible', 'off');
end
handles.iShown = 1:4;
handles.methods = {'elastix', 'halar', 'grics', 'LAP', 'demons'};

% normalize image data to [0,1]
for k = 1:size(handles.h.dImg,4)
    dImgData = handles.h.dImg(:,:,:,k);
    handles.h.dImg(:,:,:,k)    = (dImgData-min(dImgData(:)))./(max(dImgData(:))-min(dImgData(:)));
    dImgData = handles.h.dImgReg(:,:,:,k);
    dImgData(dImgData<0) = 0;
    handles.h.dImgReg(:,:,:,k) = (dImgData-min(dImgData(:)))./(max(dImgData(:))-min(dImgData(:)));
end

if isfield(handles, 'EGmode') && strcmp(handles.EGmode,'auto')
    fEvalAuto(hObject,handles);
    handles.closeFigure = true;
    % Update handles structure
    guidata(hObject, handles);
    
else    
    set(handles.cb_showQuiver,'Value',1)
%     set(handles.textFilename, 'String', ['Dataset: ',handles.h.sFilename(1:min(end,94))])
%     methods = {'elastix', 'halar', 'grics', 'LAP', 'demons'};
%     method = methods{handles.h.nRegMethod};
    handles.RegMParam = ['Registration Method: ', handles.h.sRegMethod];
    handles.RegMParam = [handles.RegMParam, ', Parameter: ',handles.h.sParFile];
    set(handles.tRegMParm, 'String', handles.RegMParam)
    sStandardLabel = {'Image 01 (reference/fixed)'; 'Image 02'; 'Image 03'; 'Image 04'};
    for iI=1:4
        if(~isempty(handles.h.sShownames{iI}))
            eval(sprintf('set(handles.txGate0%d,''String'',{sStandardLabel{iI}, handles.h.sShownames{iI}});',iI));
        else % if iNGates < 4
            eval(sprintf('set(handles.txGate0%d,''String'',sStandardLabel(iI));',iI));
        end
    end    
    
    % compute det(Jac), divergence and arrow length of deformation field
    handles = fComputeDiv(handles);
    
    % Update handles structure
    guidata(hObject, handles);
    
    % plot all gates
    fPlotImageAndDefField(hObject, handles);
    drawnow update;
end


function EvalGUI_ResizeFcn(hObject, eventdata, handles)


function varargout = EvalGUI_OutputFcn(hObject, eventdata, handles)
try varargout{1}= handles.output; catch; end;
if(isfield(handles,'closeFigure') && handles.closeFigure)
    EvalGUI_CloseRequestFcn(hObject, eventdata, handles)
end


function EvalGUI_CloseRequestFcn(hObject, eventdata, handles)
% save GUIPreference
try
currpath = fileparts(mfilename('fullpath'));
if(isfield(handles,'SPaths'))
    SPaths = handles.SPaths;
    standardVoxelsize = [1 1 1];
    lEvalMetrics = handles.lEvalMetrics;
    save([currpath, filesep, 'GUIPreferences.mat'],'SPaths','standardVoxelsize', 'lEvalMetrics', '-append');
end
delete(hObject);
catch
end


%% ------------------------------------------------------------------------
% callback functions
% -------------------------------------------------------------------------
function pop_SagCorTrans_Callback(hObject, eventdata, handles)
% change view between sagital, coronal and transverse slices

% get selected slice orientation
contents = cellstr(get(handles.pop_SagCorTrans,'String'));
lead=contents{get(handles.pop_SagCorTrans,'value')};

% get image array and deformation fields
dImg = handles.h.dImg;
dImgReg = handles.h.dImgReg;
dFx2 = handles.h.SDeform(2).dFx;
dFy2 = handles.h.SDeform(2).dFy;
dFz2 = handles.h.SDeform(2).dFz;
dFx3 = handles.h.SDeform(3).dFx;
dFy3 = handles.h.SDeform(3).dFy;
dFz3 = handles.h.SDeform(3).dFz;
dFx4 = handles.h.SDeform(4).dFx;
dFy4 = handles.h.SDeform(4).dFy;
dFz4 = handles.h.SDeform(4).dFz;
if(strcmp(get(handles.rbForward,'Visible'),'on'))
    dBx2 = handles.h.SDeform(2).dBx;
    dBy2 = handles.h.SDeform(2).dBy;
    dBz2 = handles.h.SDeform(2).dBz;
    dBx3 = handles.h.SDeform(3).dBx;
    dBy3 = handles.h.SDeform(3).dBy;
    dBz3 = handles.h.SDeform(3).dBz;
    dBx4 = handles.h.SDeform(4).dBx;
    dBy4 = handles.h.SDeform(4).dBy;
    dBz4 = handles.h.SDeform(4).dBz;
    lBackward = true;
else
    lBackward = false;
end

% permutation vector for orientation changes
CorSag = [1,3,2,4];
CorTra = [3,2,1,4];
SagCor = [1,3,2,4];
SagTra = [2,3,1,4];
TraCor = [3,2,1,4];
TraSag = [3,1,2,4];

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% very sensitive code: be attentive when changing it!

% especially Fx/Fy/Fz permutation is tricky!

switch lead                                             % for all combinations of orientation changes:
    case 'coronal'
        nSCT=1;                                         % 1. set: - flag for slice orientation
        handles.quiverFactor = 8;                       %         - standard quiver densitiy
        handles.quiverScale = 1;                        %         - standard quiver length
%         set(handles.slider_quiverScale,'value', handles.quiverScale);   %         - slider positions
%         set(handles.slider_quiverFac,'value', handles.quiverFactor);
        if handles.nSCT == 2 % sag -> cor
            dImgRot = permute(dImg,SagCor);             % 2. a) permute all image data according to permutation vector
            dImgRegRot = permute(dImgReg,SagCor);
            handles.h.SDeform(2).dFx = permute(dFz2,SagCor);             %    b) permute all deformation fields (x, y and z coordinates)
            handles.h.SDeform(2).dFz = permute(dFx2,SagCor);             %    c) set permuted deformFields into old ones
            handles.h.SDeform(2).dFy = permute(dFy2,SagCor);             %      ! swap names (e.g. y <-> z) according to permutation
            handles.h.SDeform(3).dFx = permute(dFz3,SagCor);
            handles.h.SDeform(3).dFz = permute(dFx3,SagCor);
            handles.h.SDeform(3).dFy = permute(dFy3,SagCor);
            handles.h.SDeform(4).dFx = permute(dFz4,SagCor);
            handles.h.SDeform(4).dFz = permute(dFx4,SagCor);
            handles.h.SDeform(4).dFy = permute(dFy4,SagCor);
            if(lBackward)
                handles.h.SDeform(2).dBx = permute(dBz2,SagCor);             %    b) permute all deformation fields (x, y and z coordinates)
                handles.h.SDeform(2).dBz = permute(dBx2,SagCor);             %    c) set permuted deformFields into old ones
                handles.h.SDeform(2).dBy = permute(dBy2,SagCor);             %      ! swap names (e.g. y <-> z) according to permutation
                handles.h.SDeform(3).dBx = permute(dBz3,SagCor);
                handles.h.SDeform(3).dBz = permute(dBx3,SagCor);
                handles.h.SDeform(3).dBy = permute(dBy3,SagCor);
                handles.h.SDeform(4).dBx = permute(dBz4,SagCor);
                handles.h.SDeform(4).dBz = permute(dBx4,SagCor);
                handles.h.SDeform(4).dBy = permute(dBy4,SagCor);
            end
            if isfield(handles,'divF2')
                handles.divF2 = permute(handles.divF2,SagCor);
                handles.divF3 = permute(handles.divF3,SagCor);
                handles.divF4 = permute(handles.divF4,SagCor);
            end
            if isfield(handles,'detJ2')
                handles.detJ2 = permute(handles.detJ2,SagCor);
                handles.detJ3 = permute(handles.detJ3,SagCor);
                handles.detJ4 = permute(handles.detJ4,SagCor);
            end
            if isfield(handles,'arrL2')
                handles.arrL2 = permute(handles.arrL2,SagCor);
                handles.arrL3 = permute(handles.arrL3,SagCor);
                handles.arrL4 = permute(handles.arrL4,SagCor);
            end
            if isfield(handles,'gradL2')
                handles.gradL2 = permute(handles.gradL2,SagCor);
                handles.gradL3 = permute(handles.gradL3,SagCor);
                handles.gradL4 = permute(handles.gradL4,SagCor);
            end
        elseif handles.nSCT == 3 % tra -> cor
            dImgRot = permute(dImg,TraCor);
            dImgRegRot = permute(dImgReg,TraCor);
            handles.h.SDeform(2).dFz = permute(dFy2,TraCor);
            handles.h.SDeform(2).dFy = permute(dFz2,TraCor);
            handles.h.SDeform(2).dFx = permute(dFx2,TraCor);
            handles.h.SDeform(3).dFz = permute(dFy3,TraCor);
            handles.h.SDeform(3).dFy = permute(dFz3,TraCor);
            handles.h.SDeform(3).dFx = permute(dFx3,TraCor);
            handles.h.SDeform(4).dFz = permute(dFy4,TraCor);
            handles.h.SDeform(4).dFy = permute(dFz4,TraCor);
            handles.h.SDeform(4).dFx = permute(dFx4,TraCor);
            if(lBackward)
                handles.h.SDeform(2).dBz = permute(dBy2,TraCor);
                handles.h.SDeform(2).dBy = permute(dBz2,TraCor);
                handles.h.SDeform(2).dBx = permute(dBx2,TraCor);
                handles.h.SDeform(3).dBz = permute(dBy3,TraCor);
                handles.h.SDeform(3).dBy = permute(dBz3,TraCor);
                handles.h.SDeform(3).dBx = permute(dBx3,TraCor);
                handles.h.SDeform(4).dBz = permute(dBy4,TraCor);
                handles.h.SDeform(4).dBy = permute(dBz4,TraCor);
                handles.h.SDeform(4).dBx = permute(dBx4,TraCor);
            end
            if isfield(handles,'divF2')
                handles.divF2 = permute(handles.divF2,TraCor);
                handles.divF3 = permute(handles.divF3,TraCor);
                handles.divF4 = permute(handles.divF4,TraCor);
            end
            if isfield(handles,'detJ2')
                handles.detJ2 = permute(handles.detJ2,TraCor);
                handles.detJ3 = permute(handles.detJ3,TraCor);
                handles.detJ4 = permute(handles.detJ4,TraCor);
            end
            if isfield(handles,'arrL2')
                handles.arrL2 = permute(handles.arrL2,TraCor);
                handles.arrL3 = permute(handles.arrL3,TraCor);
                handles.arrL4 = permute(handles.arrL4,TraCor);
            end
            if isfield(handles,'gradL2')
                handles.gradL2 = permute(handles.gradL2,TraCor);
                handles.gradL3 = permute(handles.gradL3,TraCor);
                handles.gradL4 = permute(handles.gradL4,TraCor);
            end
        else return                                     % 3. do nothing, if slice orientation did not change
        end
    case 'sagittal'
        nSCT=2;
        handles.quiverFactor = 6;
        handles.quiverScale = 1;
%         set(handles.slider_quiverScale,'value', 4);
%         set(handles.slider_quiverFac,'value', 5);
        if handles.nSCT == 1 % cor -> sag
            dImgRot = permute(dImg,CorSag);
            dImgRegRot = permute(dImgReg,CorSag);
            handles.h.SDeform(2).dFx = permute(dFz2,CorSag);
            handles.h.SDeform(2).dFz = permute(dFx2,CorSag);
            handles.h.SDeform(2).dFy = permute(dFy2,CorSag);
            handles.h.SDeform(3).dFx = permute(dFz3,CorSag);
            handles.h.SDeform(3).dFz = permute(dFx3,CorSag);
            handles.h.SDeform(3).dFy = permute(dFy3,CorSag);
            handles.h.SDeform(4).dFx = permute(dFz4,CorSag);
            handles.h.SDeform(4).dFz = permute(dFx4,CorSag);
            handles.h.SDeform(4).dFy = permute(dFy4,CorSag);
            if(lBackward)
                handles.h.SDeform(2).dBx = permute(dBz2,CorSag);
                handles.h.SDeform(2).dBz = permute(dBx2,CorSag);
                handles.h.SDeform(2).dBy = permute(dBy2,CorSag);
                handles.h.SDeform(3).dBx = permute(dBz3,CorSag);
                handles.h.SDeform(3).dBz = permute(dBx3,CorSag);
                handles.h.SDeform(3).dBy = permute(dBy3,CorSag);
                handles.h.SDeform(4).dBx = permute(dBz4,CorSag);
                handles.h.SDeform(4).dBz = permute(dBx4,CorSag);
                handles.h.SDeform(4).dBy = permute(dBy4,CorSag);
            end
            if isfield(handles,'divF2')
                handles.divF2 = permute(handles.divF2,CorSag);
                handles.divF3 = permute(handles.divF3,CorSag);
                handles.divF4 = permute(handles.divF4,CorSag);
            end
            if isfield(handles,'detJ2')
                handles.detJ2 = permute(handles.detJ2,CorSag);
                handles.detJ3 = permute(handles.detJ3,CorSag);
                handles.detJ4 = permute(handles.detJ4,CorSag);
            end
            if isfield(handles,'arrL2')
                handles.arrL2 = permute(handles.arrL2,CorSag);
                handles.arrL3 = permute(handles.arrL3,CorSag);
                handles.arrL4 = permute(handles.arrL4,CorSag);
            end
            if isfield(handles,'gradL2')
                handles.gradL2 = permute(handles.gradL2,CorSag);
                handles.gradL3 = permute(handles.gradL3,CorSag);
                handles.gradL4 = permute(handles.gradL4,CorSag);
            end
        elseif handles.nSCT == 3 % tra -> sag
            dImgRot = permute(dImg,TraSag);
            dImgRegRot = permute(dImgReg,TraSag);
            handles.h.SDeform(2).dFy = permute(dFz2,TraSag);
            handles.h.SDeform(2).dFz = permute(dFx2,TraSag);
            handles.h.SDeform(2).dFx = permute(dFy2,TraSag);
            handles.h.SDeform(3).dFy = permute(dFz3,TraSag);
            handles.h.SDeform(3).dFz = permute(dFx3,TraSag);
            handles.h.SDeform(3).dFx = permute(dFy3,TraSag);
            handles.h.SDeform(4).dFy = permute(dFz4,TraSag);
            handles.h.SDeform(4).dFz = permute(dFx4,TraSag);
            handles.h.SDeform(4).dFx = permute(dFy4,TraSag);
            if(lBackward)
                handles.h.SDeform(2).dBy = permute(dBz2,TraSag);
                handles.h.SDeform(2).dBz = permute(dBx2,TraSag);
                handles.h.SDeform(2).dBx = permute(dBy2,TraSag);
                handles.h.SDeform(3).dBy = permute(dBz3,TraSag);
                handles.h.SDeform(3).dBz = permute(dBx3,TraSag);
                handles.h.SDeform(3).dBx = permute(dBy3,TraSag);
                handles.h.SDeform(4).dBy = permute(dBz4,TraSag);
                handles.h.SDeform(4).dBz = permute(dBx4,TraSag);
                handles.h.SDeform(4).dBx = permute(dBy4,TraSag);
            end
            if isfield(handles,'divF2')
                handles.divF2 = permute(handles.divF2,TraSag);
                handles.divF3 = permute(handles.divF3,TraSag);
                handles.divF4 = permute(handles.divF4,TraSag);
            end
            if isfield(handles,'detJ2')
                handles.detJ2 = permute(handles.detJ2,TraSag);
                handles.detJ3 = permute(handles.detJ3,TraSag);
                handles.detJ4 = permute(handles.detJ4,TraSag);
            end
            if isfield(handles,'arrL2')
                handles.arrL2 = permute(handles.arrL2,TraSag);
                handles.arrL3 = permute(handles.arrL3,TraSag);
                handles.arrL4 = permute(handles.arrL4,TraSag);
            end
            if isfield(handles,'gradL2')
                handles.gradL2 = permute(handles.gradL2,TraSag);
                handles.gradL3 = permute(handles.gradL3,TraSag);
                handles.gradL4 = permute(handles.gradL4,TraSag);
            end
        else return
        end
    case 'transverse'
        nSCT=3;
        handles.quiverFactor = 5;
        handles.quiverScale = 1;
%         set(handles.slider_quiverScale,'value', 4);
%         set(handles.slider_quiverFac,'value', 6);
        if handles.nSCT == 1 % cor -> tra
            dImgRot = permute(dImg,CorTra);
            dImgRegRot = permute(dImgReg,CorTra);
            handles.h.SDeform(2).dFz = permute(dFy2,CorTra);
            handles.h.SDeform(2).dFy = permute(dFz2,CorTra);
            handles.h.SDeform(2).dFx = permute(dFx2,CorTra);
            handles.h.SDeform(3).dFz = permute(dFy3,CorTra);
            handles.h.SDeform(3).dFy = permute(dFz3,CorTra);
            handles.h.SDeform(3).dFx = permute(dFx3,CorTra);
            handles.h.SDeform(4).dFz = permute(dFy4,CorTra);
            handles.h.SDeform(4).dFy = permute(dFz4,CorTra);
            handles.h.SDeform(4).dFx = permute(dFx4,CorTra);
            if(lBackward)
                handles.h.SDeform(2).dBz = permute(dBy2,CorTra);
                handles.h.SDeform(2).dBy = permute(dBz2,CorTra);
                handles.h.SDeform(2).dBx = permute(dBx2,CorTra);
                handles.h.SDeform(3).dBz = permute(dBy3,CorTra);
                handles.h.SDeform(3).dBy = permute(dBz3,CorTra);
                handles.h.SDeform(3).dBx = permute(dBx3,CorTra);
                handles.h.SDeform(4).dBz = permute(dBy4,CorTra);
                handles.h.SDeform(4).dBy = permute(dBz4,CorTra);
                handles.h.SDeform(4).dBx = permute(dBx4,CorTra);
            end
            if isfield(handles,'divF2')
                handles.divF2 = permute(handles.divF2,CorTra);
                handles.divF3 = permute(handles.divF3,CorTra);
                handles.divF4 = permute(handles.divF4,CorTra);
            end
            if isfield(handles,'detJ2')
                handles.detJ2 = permute(handles.detJ2,CorTra);
                handles.detJ3 = permute(handles.detJ3,CorTra);
                handles.detJ4 = permute(handles.detJ4,CorTra);
            end
            if isfield(handles,'arrL2')
                handles.arrL2 = permute(handles.arrL2,CorTra);
                handles.arrL3 = permute(handles.arrL3,CorTra);
                handles.arrL4 = permute(handles.arrL4,CorTra);
            end
            if isfield(handles,'gradL2')
                handles.gradL2 = permute(handles.gradL2,CorTra);
                handles.gradL3 = permute(handles.gradL3,CorTra);
                handles.gradL4 = permute(handles.gradL4,CorTra);
            end
        elseif handles.nSCT == 2 % sag -> tra
            dImgRot = permute(dImg,SagTra);
            dImgRegRot = permute(dImgReg,SagTra);
            handles.h.SDeform(2).dFz = permute(dFy2,SagTra);
            handles.h.SDeform(2).dFx = permute(dFz2,SagTra);
            handles.h.SDeform(2).dFy = permute(dFx2,SagTra);
            handles.h.SDeform(3).dFz = permute(dFy3,SagTra);
            handles.h.SDeform(3).dFx = permute(dFz3,SagTra);
            handles.h.SDeform(3).dFy = permute(dFx3,SagTra);
            handles.h.SDeform(4).dFz = permute(dFy4,SagTra);
            handles.h.SDeform(4).dFx = permute(dFz4,SagTra);
            handles.h.SDeform(4).dFy = permute(dFx4,SagTra);
            if(lBackward)
                handles.h.SDeform(2).dBz = permute(dBy2,SagTra);
                handles.h.SDeform(2).dBx = permute(dBz2,SagTra);
                handles.h.SDeform(2).dBy = permute(dBx2,SagTra);
                handles.h.SDeform(3).dBz = permute(dBy3,SagTra);
                handles.h.SDeform(3).dBx = permute(dBz3,SagTra);
                handles.h.SDeform(3).dBy = permute(dBx3,SagTra);
                handles.h.SDeform(4).dBz = permute(dBy4,SagTra);
                handles.h.SDeform(4).dBx = permute(dBz4,SagTra);
                handles.h.SDeform(4).dBy = permute(dBx4,SagTra); 
            end
            if isfield(handles,'divF2')
                handles.divF2 = permute(handles.divF2,SagTra);
                handles.divF3 = permute(handles.divF3,SagTra);
                handles.divF4 = permute(handles.divF4,SagTra);
            end
            if isfield(handles,'detJ2')
                handles.detJ2 = permute(handles.detJ2,SagTra);
                handles.detJ3 = permute(handles.detJ3,SagTra);
                handles.detJ4 = permute(handles.detJ4,SagTra);
            end
            if isfield(handles,'arrL2')
                handles.arrL2 = permute(handles.arrL2,SagTra);
                handles.arrL3 = permute(handles.arrL3,SagTra);
                handles.arrL4 = permute(handles.arrL4,SagTra);
            end
            if isfield(handles,'gradL2')
                handles.gradL2 = permute(handles.gradL2,SagTra);
                handles.gradL3 = permute(handles.gradL3,SagTra);
                handles.gradL4 = permute(handles.gradL4,SagTra);
            end
        else return
        end
    otherwise
        error('unexpected error');
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
set(handles.slider_quiverScale,'value', handles.quiverScale);   %         - slider positions
set(handles.tQuiverSize,'String',num2str(handles.quiverScale));
dVal = get(handles.slider_quiverFac, 'Max') - handles.quiverFactor + 1;
set(handles.slider_quiverFac,'value', dVal);
set(handles.txtDens,'String',num2str(1/handles.quiverFactor,'%.3f'));
% set new values in handle structure.
handles.h.dImg = dImgRot;
handles.h.dImgReg = dImgRegRot;
handles.nSCT = nSCT;
handles.slice= floor(size(handles.h.dImg,3)/2);
if(isfield(handles,'hI1'))
    if(ishandle(handles.hI1))
        delete(handles.hI1);
    end
    handles.hI1 = [];
end
if(isfield(handles,'hI'))
    for n=2:4
        if(ishandle(handles.hI(n)))
            delete(handles.hI(n));
        end
        handles.hI(n) = -1;
    end
end
if(isfield(handles,'hQ'))
    for n=2:4
        if(ishandle(handles.hQ(n)))
            delete(handles.hQ(n));
        end
        handles.hQ(n) = -1;
    end
end

% plot new slice orientation
fPlotImageAndDefField(hObject, handles)


function pb_SegLung_Callback(hObject, eventdata, handles)
% segment lung
if handles.nSCT == 2 || handles.nSCT == 3
    errordlg('Please segment the lungs in the coronal view.')
    return
end

% which type of image shall be segmented?
if(isfield(eventdata,'autoMode'))
    qSegImg = 'transformed images';
else
    qSegImg = questdlg('Which images shall be segmented?', ' ', 'original images', 'transformed images', 'transformed images');
end

switch qSegImg
    case 'original images'
        dImg2Seg = handles.h.dImg(:,:,:,1:handles.iActualGates);
        isOrigImg = 1;
    case 'transformed images'
        dImg2Seg = handles.h.dImgReg(:,:,:,1:handles.iActualGates);
        dImg2Seg(:,:,:,1) = handles.h.dImg(:,:,:,1);
        isOrigImg = 0;
    case ''
        return
end

% segment lungs in the 4 gates and display overlap measures
h = fwaitbar(0,'Segmenting lungs. Please wait...'); st = 0; steps = 2*size(dImg2Seg,4)-1;
% check if there are already masks and if the gates have changed since the last computation
if (~isfield(handles,'lung3d') || handles.newGates == 1 || ~strcmp(qSegImg, handles.qSegImg))
    for iI=1:size(dImg2Seg,4)
        handles.lung3d{iI} = fLungSegmKohlm(dImg2Seg,iI,1); st= st+1;fwaitbar(st/steps,h);
    end
else
    st= size(dImg2Seg,4);fwaitbar(st/steps,h);
end

if isOrigImg % transform images according to deformation field
    for iI = 2:size(dImg2Seg,4)
        uB = {handles.h.SDeform(iI).dBy,handles.h.SDeform(iI).dBx,handles.h.SDeform(iI).dBz};
        handles.lung3d{iI} = imshift_3D(handles.lung3d{iI},uB,'bilinear');
    end
end
handles.newGate = 0;  % flag reset to 0
hData.dImg = dImg2Seg;
hData.lung3d = handles.lung3d;
hData.nGate = 1;
CutOffValue = CutLungGUI('inarg', hData);
if CutOffValue == 1; try close(h); catch; end; return; end;
handles.lung3dcut = handles.lung3d;
for iI=1:size(dImg2Seg,4)
    handles.lung3dcut{iI}(1:CutOffValue,:,:)= 0;
end

dice = zeros(size(dImg2Seg,4)-1,1);
jac = zeros(size(dImg2Seg,4)-1,1);
for k = 1:size(dImg2Seg,4)-1
    [dice(k), jac(k)] = fEvalOverlap(handles.lung3dcut{1},handles.lung3dcut{k+1}); st= st+1;fwaitbar(st/steps,h);
end
try close(h); catch; end;

handles.qSegImg = qSegImg;
if(isOrigImg) % starting from original image
    if(isfield(handles,'lung'))
        handles.lung.isPerformed = handles.lung.isPerformed | [true, false]; % original images, transformed images
    else
        handles.lung.isPerformed = [true, false];
    end
    handles.lung.diceOri = dice;
    handles.lung.jacOri = jac;
else
    if(isfield(handles,'lung'))
        handles.lung.isPerformed = handles.lung.isPerformed | [false, true]; % original images, transformed images
    else
        handles.lung.isPerformed = [false, true];
    end
    handles.lung.diceReg = dice;
    handles.lung.jacReg = jac;
end

iInd = get(handles.pm_Gate,'Value');
lMap = repmat(handles.iShown,2,1);
lMap = lMap(:);
if(mod(iInd,2) == 0) % even -> transformed
    if(handles.lung.isPerformed(2))
        set(handles.textDice,'String',num2str(handles.lung.diceReg(lMap(iInd)-1),'%6.4g'));
        set(handles.textJac,'String',num2str(handles.lung.jacReg(lMap(iInd)-1),'%6.4g'));
    else
        set(handles.textDice,'String', '');
        set(handles.textJac,'String', '');
    end
else % odd -> original
    if(handles.lung.isPerformed(1))
        set(handles.textDice,'String',num2str(handles.lung.diceOri(lMap(iInd)-1),'%6.4g'));
        set(handles.textJac,'String',num2str(handles.lung.jacOri(lMap(iInd)-1),'%6.4g'));
    else
        set(handles.textDice,'String', '');
        set(handles.textJac,'String', '');
    end 
end
guidata(hObject, handles);


function pb_lmPoints_Callback(hObject, eventdata, handles)
% open landmarkpoint GUI
h.dImg = handles.h.dImg(:,:,:,1:handles.iActualGates);
h.dImgReg = handles.h.dImgReg(:,:,:,1:handles.iActualGates);
h.SDeform = handles.h.SDeform(1:handles.iActualGates);
h.nRegMethod = handles.h.nRegMethod;
h.sFilename = handles.h.sFilename;
h.sParFile = handles.h.sParFile;
h.nSCT = handles.nSCT;
try h.SPaths = handles.h.SPaths; catch; end;
try h.SPaths = handles.SPaths; catch; end;
h.SGeo = handles.h.SGeo;
LandmarkGUI('inarg', h);


function slider_quiverFac_Callback(hObject, eventdata, handles)
% Slider for quiver Factor manipulation (density)

% get slider position and compute corresponding standard quiver Factor
x = get(hObject,'Max') - get(hObject,'Value') + 1;
% x = get(hObject,'Value');
k = get(hObject,'Min'):handles.dStep(2):get(hObject,'Max');
y = interp1(k,k,x,'nearest');
% dStep = get(hObject,'SliderStep');
% y = get(hObject,'Min'):dStep(1):get(hObject,'Max');
% y = 20:-1:1;
% y = [12,10,8,7,6,5,4,3,2,1];

% get deformation field arrays
if(get(handles.rbForward,'Value') > 0)
    dFx2 = handles.h.SDeform(2).dFx;
    dFy2 = handles.h.SDeform(2).dFy;
    dFx3 = handles.h.SDeform(3).dFx;
    dFy3 = handles.h.SDeform(3).dFy;
    dFx4 = handles.h.SDeform(4).dFx;
    dFy4 = handles.h.SDeform(4).dFy;
else % backward
    dFx2 = handles.h.SDeform(2).dBx;
    dFy2 = handles.h.SDeform(2).dBy;
    dFx3 = handles.h.SDeform(3).dBx;
    dFy3 = handles.h.SDeform(3).dBy;
    dFx4 = handles.h.SDeform(4).dBx;
    dFy4 = handles.h.SDeform(4).dBy;
end

% for different slice orientation adjust quiver Factor
% handles.quiverFactor = y(round(x));
% handles.quiverFactor = round(1/x);
handles.quiverFactor = y;
set(handles.txtDens,'String',num2str(1/handles.quiverFactor,'%.3f'));

% get new positions for quiver arrows
iX = 1:handles.quiverFactor:size(dFx2, 2);
iY = 1:handles.quiverFactor:size(dFy2, 1);

% draw new quiver arrows
% if(isfield(handles,'hQ') && handles.hQ(2) > 0 && ishandle(handles.hQ(2)))
%     delete(handles.hQ(2));
% end
try
    delete(handles.hQ(2));
catch
end
axes(handles.Gate02)
hold on
handles.hQ(2)=quiver(iX, iY, handles.quiverScale*dFx2(1:handles.quiverFactor:end, 1:handles.quiverFactor:end, handles.slice),...
    handles.quiverScale*dFy2(1:handles.quiverFactor:end, 1:handles.quiverFactor:end, handles.slice), 0);
set(handles.hQ(2), 'Linewidth', 1.5, 'Color', 'y');
hold off

% if(isfield(handles,'hQ') && handles.hQ(3) > 0 && ishandle(handles.hQ(3)))
%     delete(handles.hQ(3));
% end
try
    delete(handles.hQ(3));
catch
end
axes(handles.Gate03)
hold on
handles.hQ(3)=quiver(iX, iY, handles.quiverScale*dFx3(1:handles.quiverFactor:end, 1:handles.quiverFactor:end, handles.slice),...
    handles.quiverScale*dFy3(1:handles.quiverFactor:end, 1:handles.quiverFactor:end, handles.slice), 0);
set(handles.hQ(3), 'Linewidth', 1.5, 'Color', 'y');
hold off

% if(isfield(handles,'hQ') && handles.hQ(4) > 0 && ishandle(handles.hQ(4)))
%     delete(handles.hQ(4));
% end
try
    delete(handles.hQ(4));
catch
end
axes(handles.Gate04)
hold on
handles.hQ(4)=quiver(iX, iY, handles.quiverScale*dFx4(1:handles.quiverFactor:end, 1:handles.quiverFactor:end, handles.slice),...
    handles.quiverScale*dFy4(1:handles.quiverFactor:end, 1:handles.quiverFactor:end, handles.slice), 0);
set(handles.hQ(4), 'Linewidth', 1.5, 'Color', 'y');
hold off

guidata(hObject, handles)


function slider_quiverScale_Callback(hObject, eventdata, handles)
% Slider for quiver size

% get slider position and set standard scale value
x = get(hObject,'Value');
% dStep = get(hObject,'SliderStep');
k = get(hObject,'Min'):handles.dStep(1):get(hObject,'Max');
% k = [0.25 0.5 0.75 1 1.5 2 2.5 3 4 5 6];
% y = k(round(x));
y = interp1(k,k,x,'nearest');
set(handles.tQuiverSize, 'String', num2str(y))
handles.quiverScale = y;

% get deformation field arrays
if(get(handles.rbForward,'Value') > 0)
    dFx2 = handles.h.SDeform(2).dFx;
    dFy2 = handles.h.SDeform(2).dFy;
    dFx3 = handles.h.SDeform(3).dFx;
    dFy3 = handles.h.SDeform(3).dFy;
    dFx4 = handles.h.SDeform(4).dFx;
    dFy4 = handles.h.SDeform(4).dFy;
else % backward
    dFx2 = handles.h.SDeform(2).dBx;
    dFy2 = handles.h.SDeform(2).dBy;
    dFx3 = handles.h.SDeform(3).dBx;
    dFy3 = handles.h.SDeform(3).dBy;
    dFx4 = handles.h.SDeform(4).dBx;
    dFy4 = handles.h.SDeform(4).dBy;
end

% get new positions for quiver arrows
iX = 1:handles.quiverFactor:size(dFx2, 2);
iY = 1:handles.quiverFactor:size(dFy2, 1);

% draw new quiver arrows
% if(isfield(handles,'hQ') && handles.hQ(2) > 0 && ishandle(handles.hQ(2)))
%     delete(handles.hQ(2));
% end
try
    delete(handles.hQ(2));
catch
end
axes(handles.Gate02)
hold on
handles.hQ(2)=quiver(iX, iY, handles.quiverScale*dFx2(1:handles.quiverFactor:end, 1:handles.quiverFactor:end, handles.slice),...
    handles.quiverScale*dFy2(1:handles.quiverFactor:end, 1:handles.quiverFactor:end, handles.slice), 0);
set(handles.hQ(2), 'Linewidth', 1.5, 'Color', 'y');
hold off

% if(isfield(handles,'hQ') && handles.hQ(3) > 0 && ishandle(handles.hQ(3)))
%     delete(handles.hQ(3));
% end
try
    delete(handles.hQ(3));
catch
end
axes(handles.Gate03)
hold on
handles.hQ(3)=quiver(iX, iY, handles.quiverScale*dFx3(1:handles.quiverFactor:end, 1:handles.quiverFactor:end, handles.slice),...
    handles.quiverScale*dFy3(1:handles.quiverFactor:end, 1:handles.quiverFactor:end, handles.slice), 0);
set(handles.hQ(3), 'Linewidth', 1.5, 'Color', 'y');
hold off

% if(isfield(handles,'hQ') && handles.hQ(4) > 0 && ishandle(handles.hQ(4)))
%     delete(handles.hQ(4));
% end
try
    delete(handles.hQ(4));
catch
end
axes(handles.Gate04)
hold on
handles.hQ(4)=quiver(iX, iY, handles.quiverScale*dFx4(1:handles.quiverFactor:end, 1:handles.quiverFactor:end, handles.slice),...
    handles.quiverScale*dFy4(1:handles.quiverFactor:end, 1:handles.quiverFactor:end, handles.slice), 0);
set(handles.hQ(4), 'Linewidth', 1.5, 'Color', 'y');
hold off
guidata(hObject, handles)


function cb_showQuiver_Callback(hObject, eventdata, handles)
% checkbox for quivers
if(get(handles.rbForward,'Value') > 0)
    dFx2 = handles.h.SDeform(2).dFx;
    dFy2 = handles.h.SDeform(2).dFy;
    dFx3 = handles.h.SDeform(3).dFx;
    dFy3 = handles.h.SDeform(3).dFy;
    dFx4 = handles.h.SDeform(4).dFx;
    dFy4 = handles.h.SDeform(4).dFy;
else % backward
    dFx2 = handles.h.SDeform(2).dBx;
    dFy2 = handles.h.SDeform(2).dBy;
    dFx3 = handles.h.SDeform(3).dBx;
    dFy3 = handles.h.SDeform(3).dBy;
    dFx4 = handles.h.SDeform(4).dBx;
    dFy4 = handles.h.SDeform(4).dBy;
end
iX = 1:handles.quiverFactor:size(dFx2, 2);
iY = 1:handles.quiverFactor:size(dFy2, 1);

if get(hObject, 'Value')
    axes(handles.Gate02)
    hold on
    handles.hQ(2)=quiver(iX, iY, handles.quiverScale*dFx2(1:handles.quiverFactor:end, 1:handles.quiverFactor:end, handles.slice),...
        dFy2(1:handles.quiverFactor:end, 1:handles.quiverFactor:end, handles.slice), 0);
    set(handles.hQ(2), 'Linewidth', 1.5, 'Color', 'y');
    hold off
    
    axes(handles.Gate03)
    hold on
    handles.hQ(3)=quiver(iX, iY, handles.quiverScale*dFx3(1:handles.quiverFactor:end, 1:handles.quiverFactor:end, handles.slice),...
        dFy3(1:handles.quiverFactor:end, 1:handles.quiverFactor:end, handles.slice), 0);
    set(handles.hQ(3), 'Linewidth', 1.5, 'Color', 'y');
    hold off
    
    axes(handles.Gate04)
    hold on
    handles.hQ(4)=quiver(iX, iY, handles.quiverScale*dFx4(1:handles.quiverFactor:end, 1:handles.quiverFactor:end, handles.slice),...
        dFy4(1:handles.quiverFactor:end, 1:handles.quiverFactor:end, handles.slice), 0);
    set(handles.hQ(4), 'Linewidth', 1.5, 'Color', 'y');
    hold off
else
    delete(handles.hQ(2))
    delete(handles.hQ(3))
    delete(handles.hQ(4))
end
guidata(hObject, handles)


function pb_DispGates_Callback(hObject, eventdata, handles)
% change shown images
handles.newGates = 1; % the gates have changed
% what are the gates?
nGates = fLargePopup(handles.h.sShownames);
if isempty(nGates) % cancel was pressed
    return;
end
if length(nGates) ~= 4
    errordlg('Please choose 4 gates!')
    return;
end
if max(nGates) > size(handles.h.dImg,4)
    errordlg(sprintf('There are only %g gates in the current data set.',...
        size(handles.h.dImg,4)));
    return;
end
if ~isempty(find(nGates(2:end)==1,1))
    errordlg(sprintf('You cannot take the first gate at position %g. There exists no deformation field.',...
        find(nGates(2:end),1)+1));
    return;
end
cont = {'coronal', 'sagittal', 'transverse'};
nOldSCT = handles.nSCT;
if handles.nSCT ~= find(strcmp(cont,handles.h.SGeo.cOrientation{1}));
    % turn images into loaded in orientation    
    set(handles.pop_SagCorTrans,'Value', find(strcmp(cont,handles.h.SGeo.cOrientation{1})));
    
    % change view between sagital, coronal and transverse slices
    pop_SagCorTrans_Callback(hObject, eventdata, handles);

    lTurn = true;
else
    lTurn = false;
end

% store original data for the possibility to change the gates that are displayed
handles.h.dImgSave = handles.h.dImg;
handles.h.dImgRegSave = handles.h.dImgReg;
handles.h.SDeformSave = handles.h.SDeform;

handles.iShown = nGates;
% get the data from the global handle and store it to the used dImg,
% dImgReg and SDeform variables
for k = 1:4
    handles.h.dImg(:,:,:,k) = handles.h.dImgSave(:,:,:,nGates(k));
    handles.h.dImgReg(:,:,:,k) = handles.h.dImgRegSave(:,:,:,nGates(k));
    handles.h.SDeform(k) = handles.h.SDeformSave(nGates(k));
end
if(nGates(1) == 1)
    set(handles.txGate01,'String',{['Image ',num2str(nGates(1),'%02u'), ' (reference/fixed)'],handles.h.sShownames{1}});
    cStringDD = cell(6,1); l=1;
else
    set(handles.txGate01,'String',{['Image ',num2str(nGates(1),'%02u')],handles.h.sShownames{nGates(1)}});
    cStringDD = cell(8,1); l=1;
    cStringDD{l} = [num2str(nGates(1),'%02u'), ' (original)']; l=l+1;
    cStringDD{l} = [num2str(nGates(1),'%02u'), ' (transformed)']; l=l+1;
end
set(handles.txGate02,'String',{['Image ',num2str(nGates(2),'%02u')],handles.h.sShownames{nGates(2)}});
cStringDD{l} = [num2str(nGates(2),'%02u'), ' (original)']; l=l+1;
cStringDD{l} = [num2str(nGates(2),'%02u'), ' (transformed)']; l=l+1;
set(handles.txGate03,'String',{['Image ',num2str(nGates(3),'%02u')],handles.h.sShownames{nGates(3)}});
cStringDD{l} = [num2str(nGates(3),'%02u'), ' (original)']; l=l+1;
cStringDD{l} = [num2str(nGates(3),'%02u'), ' (transformed)']; l=l+1;
set(handles.txGate04,'String',{['Image ',num2str(nGates(4),'%02u')],handles.h.sShownames{nGates(4)}});
cStringDD{l} = [num2str(nGates(4),'%02u'), ' (original)']; l=l+1;
cStringDD{l} = [num2str(nGates(4),'%02u'), ' (transformed)']; % l=l+1;

set(handles.pm_Gate,'String', cStringDD);

% compute det(Jac), divergence and arrow length of deformation field
handles = fComputeDiv(handles);

guidata(hObject, handles)
% plot all gates
fPlotImageAndDefField(hObject, handles);
% plot metrics
if(isfield(handles,'cSimiReg'))
    cString = handles.cMetrics(1,handles.lEvalMetrics);
    iInd = get(handles.pm_Gate,'Value');
    lMap = repmat(handles.iShown,2,1);
    lMap = lMap(:);
    lTmp = mat2cell(repmat(': ',length(cString),1),ones(length(cString),1),2);
    if(mod(iInd,2) == 0) % even -> transformed
        cString = strcat(cString,lTmp,cellfun(@(x) sprintf('%.4f',x), handles.cSimiReg{lMap(iInd)-1}));
    else % odd -> original
        cString = strcat(cString,lTmp,cellfun(@(x) sprintf('%.4f',x), handles.cSimiOri{lMap(iInd)-1}));
    end
    set(handles.lb_metrics, 'String', cString);
end
if(isfield(handles,'lung'))
    iInd = get(handles.pm_Gate,'Value');
    lMap = repmat(handles.iShown,2,1);
    lMap = lMap(:);
    if(mod(iInd,2) == 0) % even -> transformed
        if(handles.lung.isPerformed(2))
            set(handles.textDice,'String',num2str(handles.lung.diceReg(lMap(iInd)-1),'%6.4g'));
            set(handles.textJac,'String',num2str(handles.lung.jacReg(lMap(iInd)-1),'%6.4g'));
        else
            set(handles.textDice,'String', '');
            set(handles.textJac,'String', '');
        end
    else % odd -> original
        if(handles.lung.isPerformed(1))
            set(handles.textDice,'String',num2str(handles.lung.diceOri(lMap(iInd)-1),'%6.4g'));
            set(handles.textJac,'String',num2str(handles.lung.jacOri(lMap(iInd)-1),'%6.4g'));
        else
            set(handles.textDice,'String', '');
            set(handles.textJac,'String', '');
        end 
    end
end

% turn back to desired view
if(lTurn)
    set(handles.pop_SagCorTrans,'Value', nOldSCT);
    
    % change view between sagital, coronal and transverse slices
    pop_SagCorTrans_Callback(hObject, eventdata, handles);
end
guidata(hObject, handles);


function pb_resetContrast_Callback(hObject, eventdata, handles)
% reset all contrasts and views
dIMove2 = handles.h.dImg(:,:,:,2);
colRange01 = [0, max(dIMove2(:))/2];
if strcmp(handles.sImageData, 'origImage');
%     handles.colRange = [0, max(dIMove2(:))/2];
    handles.colRange = [0, max(handles.h.dImg(:))];
elseif strcmp(handles.sImageData, 'regImage');
%     dIMoveReg2 = handles.h.dImgReg(:,:,:,2);
%     handles.colRange = [0, max(dIMoveReg2(:))/2];
    handles.colRange = [0, max(handles.h.dImgReg(:))];
elseif strcmp(handles.sImageData,'diffOri')
    handles.colRange = [0, max(dIMove2(:))/2];
elseif strcmp(handles.sImageData,'diffReg')
    dIMoveReg2 = handles.h.dImgReg(:,:,:,2);
    handles.colRange = [0, max(dIMoveReg2(:))/2];
elseif strcmp(handles.sImageData,'detJac');
    handles.colRange = [min(handles.detJ2(:))/1.5, max(handles.detJ2(:))/1.5];
elseif strcmp(handles.sImageData,'divergence');
    handles.colRange = [min(handles.divF2(:))/1.5, max(handles.divF2(:))/1.5];
elseif strcmp(handles.sImageData,'arrLength');
    handles.colRange = [min(handles.arrL2(:))/1.5, max(handles.arrL2(:))/1.5];
elseif strcmp(handles.sImageData,'isosurface')
    handles.colRange = [];
    view(handles.dView);
    camlight(handles.dCamlight(1),handles.dCamlight(2));
else
    errordlg('Unknown error in reset contrast function.')
end

caxis(handles.Gate01,colRange01);
if(~isempty(handles.colRange))
    caxis(handles.Gate02,handles.colRange);
    caxis(handles.Gate03,handles.colRange);
    caxis(handles.Gate04,handles.colRange);
else % isosurface
    colormap(handles.Gate02, handles.tcmap);
    colormap(handles.Gate03, handles.tcmap);
    colormap(handles.Gate04, handles.tcmap);
end
event.VerticalScrollCount = 0;
EvalGUI_WindowScrollWheelFcn(hObject, event, handles); % faked workaround so that reference image is also directly scaled
guidata(hObject, handles);


function uipanel8_SelectionChangeFcn(hObject, eventdata, handles)
% behaviour for change of view radio button group
fPlotImageAndDefField(hObject, handles);


function pm_Gate_Callback(hObject, eventdata, handles)
% dropdown menu of shown image metrics

% get current to be viewed gate
if(isfield(handles,'cSimiReg'))
    cString = handles.cMetrics(1,handles.lEvalMetrics).';
    iInd = get(handles.pm_Gate,'Value');
    lMap = repmat(handles.iShown(2:end),2,1);
    lMap = lMap(:);
    lTmp = mat2cell(repmat(': ',length(cString),1),ones(length(cString),1),2);
    if(mod(iInd,2) == 0) % even -> transformed
        cString = strcat(cString,lTmp,cellfun(@(x) sprintf('%.4f',x), handles.cSimiReg{lMap(iInd)-1}, 'UniformOutput', false).');
    else % odd -> original
        cString = strcat(cString,lTmp,cellfun(@(x) sprintf('%.4f',x), handles.cSimiOri{lMap(iInd)-1}, 'UniformOutput', false).');
    end
    set(handles.lb_metrics, 'String', cString);
end
if(isfield(handles,'lung'))
    iInd = get(handles.pm_Gate,'Value');
    lMap = repmat(handles.iShown(2:end),2,1);
    lMap = lMap(:);
    if(mod(iInd,2) == 0) % even -> transformed
        if(handles.lung.isPerformed(2))
            set(handles.textDice,'String',num2str(handles.lung.diceReg(lMap(iInd)-1),'%6.4g'));
            set(handles.textJac,'String',num2str(handles.lung.jacReg(lMap(iInd)-1),'%6.4g'));
        else
            set(handles.textDice,'String', '');
            set(handles.textJac,'String', '');
        end
    else % odd -> original
        if(handles.lung.isPerformed(1))
            set(handles.textDice,'String',num2str(handles.lung.diceOri(lMap(iInd)-1),'%6.4g'));
            set(handles.textJac,'String',num2str(handles.lung.jacOri(lMap(iInd)-1),'%6.4g'));
        else
            set(handles.textDice,'String', '');
            set(handles.textJac,'String', '');
        end 
    end
end


function pb_Intensity_Callback(hObject, eventdata, handles)
% evaluate intensity-based metrics

cString = handles.cMetrics(1,handles.lEvalMetrics).';
lTmp = num2cell(repmat('''',length(cString),1));
lComma = num2cell(repmat(',',length(cString),1));
sString = cell2mat(strcat(lTmp,cString,lTmp,lComma).');
sString = sString(1:end-1);

cValues = cell(1,size(handles.h.dImg,4));
h = fwaitbar(0, 'Evaluating similarity metrics'); steps = handles.iActualGates-1; st = 0;
for k=1:handles.iActualGates-1
    eval(sprintf('cValuesOri{k} = similarity_measure(handles.h.dImg(:,:,:,1), handles.h.dImg(:,:,:,k+1), %s);',sString));
    eval(sprintf('cValuesReg{k} = similarity_measure(handles.h.dImg(:,:,:,1), handles.h.dImgReg(:,:,:,k+1), %s);',sString));
    st = st+1; fwaitbar(st/steps,h);
end
try close(h); catch; end;

% get current to be viewed gate
iInd = get(handles.pm_Gate,'Value');
lMap = repmat(handles.iShown(2:end),2,1);
lMap = lMap(:);
lTmp = mat2cell(repmat(': ',length(cString),1),ones(length(cString),1),2);
if(mod(iInd,2) == 0) % even -> transformed
    cString = strcat(cString,lTmp,cellfun(@(x) sprintf('%.4f',x), cValuesReg{lMap(iInd)-1}, 'UniformOutput', false).');
else % odd -> original
    cString = strcat(cString,lTmp,cellfun(@(x) sprintf('%.4f',x), cValuesOri{lMap(iInd)-1}, 'UniformOutput', false).');
end
handles.cSimiOri = cValuesOri;
handles.cSimiReg = cValuesReg;
set(handles.lb_metrics, 'String', cString);
guidata(hObject, handles);


function pb_Intensity_ButtonDownFcn(hObject, eventdata, handles)
% changes/selects intensity-based metrics

% RIGHTCLICK
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


function pb_isoLevelVal_Callback(hObject, eventdata, handles)
% change isolevel of isosurface plot
handles.level_set = str2double(get(handles.edLevelSet,'String'));
handles = fPlotIsosurface(handles, handles.Gate02, handles.level_set, 2, false);
set(handles.text21, 'String', sprintf('|df| = %.2f', handles.level_set ));
handles = fPlotIsosurface(handles, handles.Gate03, handles.level_set, 3, false);
set(handles.text22, 'String', sprintf('|df| = %.2f', handles.level_set ));
handles = fPlotIsosurface(handles, handles.Gate04, handles.level_set, 4, false);
set(handles.text23, 'String', sprintf('|df| = %.2f', handles.level_set ));
linkprop([handles.Gate02, handles.Gate03, handles.Gate04],{'CameraPosition','CameraUpVector'});
guidata(hObject,handles);


function cb_cut_Callback(hObject, eventdata, handles)
% cut away upper part of isosurface to see inside
handles = fPlotIsosurface(handles, handles.Gate02, handles.level_set, 2, false);
handles = fPlotIsosurface(handles, handles.Gate03, handles.level_set, 3, false);
handles = fPlotIsosurface(handles, handles.Gate04, handles.level_set, 4, false);
linkprop([handles.Gate02, handles.Gate03, handles.Gate04],{'CameraPosition','CameraUpVector'});
guidata(hObject,handles);


function edLevelSet_Callback(hObject, eventdata, handles)
% changes isolevel (after pressing enter)
pb_isoLevelVal_Callback(hObject, eventdata, handles);


function pop_Grad_Callback(hObject, eventdata, handles)
% change gradient direction
iInd = get(hObject,'Value');
dVoxSize = handles.h.SGeo.cVoxelsize(handles.iShown);
if(get(handles.rbForward,'Value') > 0)
    dFx2 = handles.h.SDeform(2).dFx;
    dFy2 = handles.h.SDeform(2).dFy;
    dFz2 = handles.h.SDeform(2).dFz;
    dFx3 = handles.h.SDeform(3).dFx;
    dFy3 = handles.h.SDeform(3).dFy;
    dFz3 = handles.h.SDeform(3).dFz;
    dFx4 = handles.h.SDeform(4).dFx;
    dFy4 = handles.h.SDeform(4).dFy;
    dFz4 = handles.h.SDeform(4).dFz;
else % backward
    dFx2 = handles.h.SDeform(2).dBx;
    dFy2 = handles.h.SDeform(2).dBy;
    dFz2 = handles.h.SDeform(2).dBz;
    dFx3 = handles.h.SDeform(3).dBx;
    dFy3 = handles.h.SDeform(3).dBy;
    dFz3 = handles.h.SDeform(3).dBz;
    dFx4 = handles.h.SDeform(4).dBx;
    dFy4 = handles.h.SDeform(4).dBy;
    dFz4 = handles.h.SDeform(4).dBz;
end
[dGradFx_x, dGradFx_y, dGradFx_z] = gradient(dFx2,dVoxSize{2}(1),dVoxSize{2}(2),dVoxSize{2}(3));
[dGradFy_x, dGradFy_y, dGradFy_z] = gradient(dFy2,dVoxSize{2}(1),dVoxSize{2}(2),dVoxSize{2}(3));
[dGradFz_x, dGradFz_y, dGradFz_z] = gradient(dFz2,dVoxSize{2}(1),dVoxSize{2}(2),dVoxSize{2}(3));
switch iInd
    case 1 % x
        handles.gradL2 = sqrt(dGradFx_x.^2 + dGradFy_x.^2 + dGradFz_x.^2);
    case 2 % y 
        handles.gradL2 = sqrt(dGradFx_y.^2 + dGradFy_y.^2 + dGradFz_y.^2);
    case 3 % z
        handles.gradL2 = sqrt(dGradFx_z.^2 + dGradFy_z.^2 + dGradFz_z.^2);
end

[dGradFx_x, dGradFx_y, dGradFx_z] = gradient(dFx3,dVoxSize{3}(1),dVoxSize{3}(2),dVoxSize{3}(3));
[dGradFy_x, dGradFy_y, dGradFy_z] = gradient(dFy3,dVoxSize{3}(1),dVoxSize{3}(2),dVoxSize{3}(3));
[dGradFz_x, dGradFz_y, dGradFz_z] = gradient(dFz3,dVoxSize{3}(1),dVoxSize{3}(2),dVoxSize{3}(3));
switch iInd
    case 1 % x
        handles.gradL3 = sqrt(dGradFx_x.^2 + dGradFy_x.^2 + dGradFz_x.^2);
    case 2 % y 
        handles.gradL3 = sqrt(dGradFx_y.^2 + dGradFy_y.^2 + dGradFz_y.^2);
    case 3 % z
        handles.gradL3 = sqrt(dGradFx_z.^2 + dGradFy_z.^2 + dGradFz_z.^2);
end

[dGradFx_x, dGradFx_y, dGradFx_z] = gradient(dFx4,dVoxSize{4}(1),dVoxSize{4}(2),dVoxSize{4}(3));
[dGradFy_x, dGradFy_y, dGradFy_z] = gradient(dFy4,dVoxSize{4}(1),dVoxSize{4}(2),dVoxSize{4}(3));
[dGradFz_x, dGradFz_y, dGradFz_z] = gradient(dFz4,dVoxSize{4}(1),dVoxSize{4}(2),dVoxSize{4}(3));
switch iInd
    case 1 % x
        handles.gradL4 = sqrt(dGradFx_x.^2 + dGradFy_x.^2 + dGradFz_x.^2);
    case 2 % y 
        handles.gradL4 = sqrt(dGradFx_y.^2 + dGradFy_y.^2 + dGradFz_y.^2);
    case 3 % z
        handles.gradL4 = sqrt(dGradFx_z.^2 + dGradFy_z.^2 + dGradFz_z.^2);
end

try axes(handles.Gate02); colorbar('off'); catch; end
try axes(handles.Gate03); colorbar('off'); catch; end
try axes(handles.Gate04); colorbar('off'); catch; end
handles.colRange = [min(handles.gradL2(:)),max(handles.gradL2(:))];
handles = fPlotDiv(handles, handles.Gate02, handles.gradL2, dFx2, dFy2, 2);
handles = fPlotDiv(handles, handles.Gate03, handles.gradL3, dFx3, dFy3, 3);
handles = fPlotDiv(handles, handles.Gate04, handles.gradL4, dFx4, dFy4, 4);
guidata(hObject,handles);


function pb_save_Callback(hObject, eventdata, handles)
% save results as mat file

sMethod = handles.methods{handles.h.nRegMethod};
if(~isfield(handles.h,'SPaths'))
    [sFilename,sPathname] = uiputfile([pwd,filesep,handles.h.sFilename(1:end-4),'_EvalResult.mat'],'Specify file name');
else
    sSavename = [handles.h.SPaths.sResults,sMethod,filesep,handles.h.sFilename(1:end-4),'_EvalResult.mat'];
    [sFilename,sPathname] = uiputfile(sSavename,'Specify file name');
end

sSavename = [sPathname,filesep,sFilename];
if ~isfield(handles, 'EGmode') || ~strcmp(handles.EGmode,'auto')
    h = fwaitbar(0,'Saving Results...'); st = 0; steps = 19;
    try results.dDivF2 = handles.divF2; catch; end; st=st+1; fwaitbar(st/steps,h);
    try results.dDivF3 = handles.divF3; catch; end; st=st+1; fwaitbar(st/steps,h);
    try results.dDivF4 = handles.divF4; catch; end; st=st+1; fwaitbar(st/steps,h);
    try results.dDetJ2 = handles.detJ2; catch; end; st=st+1; fwaitbar(st/steps,h);
    try results.dDetJ3 = handles.detJ3; catch; end; st=st+1; fwaitbar(st/steps,h);
    try results.dDetJ4 = handles.detJ4; catch; end; st=st+1; fwaitbar(st/steps,h);
    try results.darrL2 = handles.arrL2; catch; end; st=st+1; fwaitbar(st/steps,h);
    try results.darrL3 = handles.arrL3; catch; end; st=st+1; fwaitbar(st/steps,h);
    try results.darrL4 = handles.arrL4; catch; end; st=st+1; fwaitbar(st/steps,h);
    try results.dgradL2 = handles.gradL2; catch; end; st=st+1; fwaitbar(st/steps,h);
    try results.dgradL3 = handles.gradL3; catch; end; st=st+1; fwaitbar(st/steps,h);
    try results.dgradL4 = handles.gradL4; catch; end; st=st+1; fwaitbar(st/steps,h);
    try results.lLung3d = handles.lung3d; catch; end; st=st+1; fwaitbar(st/steps,h);
    try results.dDiceOri = handles.lung.diceOri; catch; end; st=st+1; fwaitbar(st/steps,h);
    try results.dJaccOri = handles.lung.jacOri; catch; end; st=st+1; fwaitbar(st/steps,h);
    try results.dDiceReg = handles.lung.diceReg; catch; end; st=st+1; fwaitbar(st/steps,h);
    try results.dJaccReg = handles.lung.jacReg; catch; end; st=st+1; fwaitbar(st/steps,h);
    try results.SSimiOri = handles.cSimiOri; catch; end; st=st+1; fwaitbar(st/steps,h);
    try results.SSimiReg = handles.cSimiReg; catch; end; st=st+1; fwaitbar(st/steps,h);
    save(sSavename,'results'); try close(h); catch; end;
    msgbox('Saved results.', 'Save complete.','help')
else
    sSavename = [sSavename(1:end-4),'_SEvalResult_auto.mat'];
    try results.dDice = handles.lung.dice; catch; end;
    try results.dJacc = handles.lung.jac; catch; end;
    try results.dDiceOri = handles.lung.diceOri; catch; end; 
    try results.dJaccOri = handles.lung.jacOri; catch; end; 
    try results.dDiceReg = handles.lung.diceReg; catch; end;
    try results.dJaccReg = handles.lung.jacReg; catch; end; 
    try results.SSimiOri = handles.cSimiOri; catch; end;
    try results.SSimiReg = handles.cSimiReg; catch; end;
    save(sSavename,'results'); try close(h); catch; end;
end


function pb_SaveAs_Callback(hObject, eventdata, handles)
% export results to xls/csv file
fExport(handles,'normal');


function pb_SaveAs_ButtonDownFcn(hObject, eventdata, handles)
% export results to base workspace
fExport(handles,'right');


function fExport(handles,sType)
% common result export function

sMethod = handles.methods{handles.h.nRegMethod};
if(strcmp(sType,'normal'))
    if(~isfield(handles,'SPaths'))
        [sFilename,sPathname] = uiputfile(pwd,'Specify file name');
    else
        [sFilename,sPathname] = uiputfile({'*.xls;*.xlsx','Microsoft Excel (*.xls, *.xlsx)'; '*.csv', 'Comma-separated file (*.csv)' }, 'Specify file name', [handles.SPaths.sResults,filesep,sMethod]);
    end

    if sFilename == 0; return; end
end

h = fwaitbar(0,'Saving Results...'); st = 0; steps = 19;
try results.dDivF2 = handles.divF2; catch; end; st=st+1; fwaitbar(st/steps,h);
try results.dDivF3 = handles.divF3; catch; end; st=st+1; fwaitbar(st/steps,h);
try results.dDivF4 = handles.divF4; catch; end; st=st+1; fwaitbar(st/steps,h);
try results.dDetJ2 = handles.detJ2; catch; end; st=st+1; fwaitbar(st/steps,h);
try results.dDetJ3 = handles.detJ3; catch; end; st=st+1; fwaitbar(st/steps,h);
try results.dDetJ4 = handles.detJ4; catch; end; st=st+1; fwaitbar(st/steps,h);
try results.darrL2 = handles.arrL2; catch; end; st=st+1; fwaitbar(st/steps,h);
try results.darrL3 = handles.arrL3; catch; end; st=st+1; fwaitbar(st/steps,h);
try results.darrL4 = handles.arrL4; catch; end; st=st+1; fwaitbar(st/steps,h);
try results.dgradL2 = handles.gradL2; catch; end; st=st+1; fwaitbar(st/steps,h);
try results.dgradL3 = handles.gradL3; catch; end; st=st+1; fwaitbar(st/steps,h);
try results.dgradL4 = handles.gradL4; catch; end; st=st+1; fwaitbar(st/steps,h);
try results.lLung3d = handles.lung3d; catch; end; st=st+1; fwaitbar(st/steps,h);
try results.dDiceOri = handles.lung.diceOri; catch; end; st=st+1; fwaitbar(st/steps,h);
try results.dJaccOri = handles.lung.jacOri; catch; end; st=st+1; fwaitbar(st/steps,h);
try results.dDiceReg = handles.lung.diceReg; catch; end; st=st+1; fwaitbar(st/steps,h);
try results.dJaccReg = handles.lung.jacReg; catch; end; st=st+1; fwaitbar(st/steps,h);
try results.SSimiOri = handles.cSimiOri; catch; end; st=st+1; fwaitbar(st/steps,h);
try results.SSimiReg = handles.cSimiReg; catch; end; st=st+1; fwaitbar(st/steps,h);
if(strcmp(sType,'normal'))
    if(isfield(results,'SSimiOri') || isfield(results,'dDiceOri') || isfield(results,'dDiceReg')) % TODO: width and length!!!! zuerst in lines: original dann reg, in spalten: sim dann lung
        if(isfield(results,'SSimiOri') && ~isfield(results,'SSimiReg'))
            if(~isfield(results,'dDiceOri') && ~isfield(results,'dDiceReg'))
                iLength = length(results.SSimiOri);
                iWidth = length(results.SSimiOri{1});
                cLung = {};
                lWrite = [true, false, false, false]; % simOri, simReg, lungOri, lungReg
            elseif(~isfield(results,'dDiceOri') && isfield(results,'dDiceReg'))
                iLength = length(results.SSimiOri);
                iWidth = length(results.SSimiOri{1}) + 2;
                cLung = {'Dice','Jaccard'};
                lWrite = [true, false, false, true]; % simOri, simReg, lungOri, lungReg
            elseif(isfield(results,'dDiceOri') && ~isfield(results,'dDiceReg'))
                iLength = length(results.SSimiOri);
                iWidth = length(results.SSimiOri{1}) + 2;
                cLung = {'Dice','Jaccard'};
                lWrite = [true, false, true, false]; % simOri, simReg, lungOri, lungReg
            else
                iLength = 2*length(results.SSimiOri);
                iWidth = length(results.SSimiOri{1}) + 2;
                cLung = {'Dice','Jaccard'};
                lWrite = [true, false, true, true]; % simOri, simReg, lungOri, lungReg
            end
            cHeader = cat(2,handles.cMetrics(1,handles.lEvalMetrics),cLung);
            
        elseif(~isfield(results,'SSimiOri') && isfield(results,'SSimiReg'))
            if(~isfield(results,'dDiceOri') && ~isfield(results,'dDiceReg'))                
                iLength = length(results.SSimiReg);
                iWidth = length(results.SSimiReg{1});
                cLung = {};
                lWrite = [false, true, false, false]; % simOri, simReg, lungOri, lungReg
            elseif(~isfield(results,'dDiceOri') && isfield(results,'dDiceReg'))
                iLength = length(results.SSimiReg);
                iWidth = length(results.SSimiReg{1}) + 2;
                cLung = {'Dice','Jaccard'};
                lWrite = [false, true, false, true]; % simOri, simReg, lungOri, lungReg
            elseif(isfield(results,'dDiceOri') && ~isfield(results,'dDiceReg'))
                iLength = length(results.SSimiReg);
                iWidth = length(results.SSimiReg{1}) + 2;
                cLung = {'Dice','Jaccard'};
                lWrite = [false, true, true, false]; % simOri, simReg, lungOri, lungReg
            else
                iLength = 2*length(results.SSimiReg);
                iWidth = length(results.SSimiReg{1}) + 2;
                cLung = {'Dice','Jaccard'};
                lWrite = [false, true, true, true]; % simOri, simReg, lungOri, lungReg
            end  
            cHeader = cat(2,handles.cMetrics(1,handles.lEvalMetrics),cLung);
            
        elseif(isfield(results,'SSimiOri') && isfield(results,'SSimiReg'))
            iLength = 2*length(results.SSimiReg);
            if(~isfield(results,'dDiceOri') && ~isfield(results,'dDiceReg'))                
                iWidth = length(results.SSimiReg{1});
                cLung = {};
                lWrite = [true, true, false, false]; % simOri, simReg, lungOri, lungReg
            elseif(~isfield(results,'dDiceOri') && isfield(results,'dDiceReg'))
                iWidth = length(results.SSimiOri{1}) + 2;
                cLung = {'Dice','Jaccard'};
                lWrite = [true, true, false, true]; % simOri, simReg, lungOri, lungReg
            elseif(isfield(results,'dDiceOri') && ~isfield(results,'dDiceReg'))
                iWidth = length(results.SSimiOri{1}) + 2;
                cLung = {'Dice','Jaccard'};
                lWrite = [true, true, true, false]; % simOri, simReg, lungOri, lungReg
            else
                iWidth = length(results.SSimiOri{1}) + 2;
                cLung = {'Dice','Jaccard'};
                lWrite = [true, true, true, true]; % simOri, simReg, lungOri, lungReg
            end  
            cHeader = cat(2,handles.cMetrics(1,handles.lEvalMetrics),cLung);
        else
            if(~isfield(results,'dDiceOri') && ~isfield(results,'dDiceReg'))
                iLength = 0; % not possible
                iWidth = 0;
                cHeader = {};
                lWrite = [false, false, false, false]; % simOri, simReg, lungOri, lungReg
            elseif(~isfield(results,'dDiceOri') && isfield(results,'dDiceReg'))
                iLength = length(results.dDiceReg);
                iWidth = 2;
                cHeader = {'Dice','Jaccard'};
                lWrite = [false, false, false, true]; % simOri, simReg, lungOri, lungReg
            elseif(isfield(results,'dDiceOri') && ~isfield(results,'dDiceReg'))
                iLength = length(results.dDiceOri);
                iWidth = 2;
                cHeader = {'Dice','Jaccard'};
                lWrite = [false, false, true, false]; % simOri, simReg, lungOri, lungReg
            else
                iLength = length(results.dDiceOri);
                iWidth = 2;
                cHeader = {'Dice','Jaccard'};
                lWrite = [false, false, true, true]; % simOri, simReg, lungOri, lungReg
            end
        end
        cMatrix = cell(iLength+2,1+iWidth);
        cMatrix{1,1} = [sprintf('filename: %s', handles.h.sFilename), ' registration method: ',sMethod, ' parameter file: ', handles.h.sParFile];
        cMatrix(2,2:end) = cHeader;
        for iI=1:iLength
            if((all(lWrite(1:2)) || all(lWrite(3:4))) && iI > iLength/2) % just for Ori+Reg iLength is twice as long
                iIdx = iI-iLength/2;
                cMatrix{2+iI,1} = sprintf('Image %02d (transformed)', iIdx+1);
                if(lWrite(2) && ~lWrite(4))
                    cMatrix(2+iI,2:2+length(results.SSimiReg{iIdx})-1) = results.SSimiReg{iIdx};
                elseif(~lWrite(2) && lWrite(4))
                    cMatrix(2+iI,2:end) = {results.dDiceReg(iIdx),results.dJaccReg(iIdx)};
                elseif(lWrite(2) && lWrite(4))
                    cMatrix(2+iI,2:end) = cat(2,results.SSimiReg{iIdx},{results.dDiceReg(iIdx),results.dJaccReg(iIdx)});
                end
            else
                cMatrix{2+iI,1} = sprintf('Image %02d (original)', iI+1);
                if(lWrite(1) && ~lWrite(3))
                    cMatrix(2+iI,2:2+length(results.SSimiOri{iI})-1) = results.SSimiOri{iI};
                elseif(~lWrite(1) && lWrite(3))
                    cMatrix(2+iI,2:end) = {results.dDiceOri(iI),results.dJaccOri(iI)};
                elseif(lWrite(1) && lWrite(3))
                    cMatrix(2+iI,2:end) = cat(2,results.SSimiOri{iI},{results.dDiceOri(iI),results.dJaccOri(iI)});
                end
            end
        end
        [~,~,sExt] = fileparts(sFilename);
        switch sExt
            case {'.xls', '.xlsx'}
                xlswrite([sPathname,sFilename],cMatrix);
            case {'.csv'}
                fid = fopen([sPathname,sFilename],'w');
                fprintf(fid, '%s\n', cMatrix{1,1});
                formatSpec = repmat('%s,', 1, size(cMatrix,2));
                formatSpec = [formatSpec(1:end-1), '\n'];
                fprintf(fid, formatSpec, cMatrix{2,:});
                for iI=3:size(cMatrix,1)
                    formatSpec = ['%s,', repmat('%f,', 1, size(cMatrix,2)-1)];
                    formatSpec = [formatSpec(1:end-1), '\n'];
                    fprintf(fid, formatSpec, cMatrix{iI,:});
                end
                fclose(fid);
        end
    end
else
    assignin('base','sEvalGUI',results);
end
try close(h); catch; end;
msgbox('Results exported.', 'Export complete.','help')


function viewPanel_SelectionChangeFcn(hObject, eventdata, handles)
% change between original images, registered images and deformation field
% changes (det(Jac), divergence, arrow length, isosurface)

set(handles.uitoggletool4,'Enable','off');
set(handles.cb_showQuiver,'Enable','on');
set(handles.cb_cut,'Visible','off');
set(handles.pb_isoLevelVal,'Visible','off');
set(handles.edLevelSet,'Visible','off');
set(handles.pop_Grad,'Visible','off');
if eventdata.NewValue ~= handles.rb_isosurface % going back from isosurface
    if(isfield(handles,'hIso'))
        for iI=2:4
            if(ishandle(handles.hIso(1,iI))); delete(handles.hIso(1,iI)); end;
            if(ishandle(handles.hIso(2,iI))); delete(handles.hIso(2,iI)); end;
        end
    end
end

if eventdata.NewValue == handles.rb_origImage
    try axes(handles.Gate02); colorbar('off'); catch; end
    try axes(handles.Gate03); colorbar('off'); catch; end
    try axes(handles.Gate04); colorbar('off'); catch; end
    dIMove2 = handles.h.dImg(:,:,:,2);
    dIMove3 = handles.h.dImg(:,:,:,3);
    dIMove4 = handles.h.dImg(:,:,:,4);
%     handles.colRange = [nanmin([min(dIMove2(:)),min(dIMove3(:)),min(dIMove4(:))]) nanmax([max(dIMove2(:))/3,max(dIMove3(:))/3,max(dIMove4(:))/3])];
    handles.colRange = [nanmin(handles.h.dImg(:)), nanmax(handles.h.dImg(:))];
    axes(handles.Gate01); unfreezeColors; caxis(handles.colRange); colormap(handles.Gate01,'gray'); freezeColors;
    handles = fPlotSelectionChange(handles, handles.Gate02, dIMove2(:,:,handles.slice), handles.colRange, 2);
    handles = fPlotSelectionChange(handles, handles.Gate03, dIMove3(:,:,handles.slice), handles.colRange, 3);
    handles = fPlotSelectionChange(handles, handles.Gate04, dIMove4(:,:,handles.slice), handles.colRange, 4);
    handles.sImageData = 'origImage';
elseif eventdata.NewValue == handles.rb_regImage
    try axes(handles.Gate02); colorbar('off'); catch; end
    try axes(handles.Gate03); colorbar('off'); catch; end
    try axes(handles.Gate04); colorbar('off'); catch; end
    dIMoveReg2 = handles.h.dImgReg(:,:,:,2);
    dIMoveReg3 = handles.h.dImgReg(:,:,:,3);
    dIMoveReg4 = handles.h.dImgReg(:,:,:,4);
%     handles.colRange = [nanmin([min(dIMoveReg2(:)),min(dIMoveReg3(:)),min(dIMoveReg4(:))]) nanmax([max(dIMoveReg2(:))/3,max(dIMoveReg3(:))/3,max(dIMoveReg4(:))/3])];
    handles.colRange = [nanmin(handles.h.dImgReg(:)), nanmax(handles.h.dImgReg(:))];
    axes(handles.Gate01); unfreezeColors; caxis(handles.colRange); colormap(handles.Gate01,'gray'); freezeColors;
    handles = fPlotSelectionChange(handles, handles.Gate02, dIMoveReg2(:,:,handles.slice), handles.colRange, 2);
    handles = fPlotSelectionChange(handles, handles.Gate03, dIMoveReg3(:,:,handles.slice), handles.colRange, 3);
    handles = fPlotSelectionChange(handles, handles.Gate04, dIMoveReg4(:,:,handles.slice), handles.colRange, 4);
    handles.sImageData = 'regImage';
    
elseif eventdata.NewValue == handles.rb_DiffOri
    try axes(handles.Gate02); colorbar('off'); catch; end
    try axes(handles.Gate03); colorbar('off'); catch; end
    try axes(handles.Gate04); colorbar('off'); catch; end
    dIRef = handles.h.dImg(:,:,:,1);
    dIMove2 = handles.h.dImg(:,:,:,2);
    dIMove3 = handles.h.dImg(:,:,:,3);
    dIMove4 = handles.h.dImg(:,:,:,4);
    handles = fPlotDiff(handles, handles.Gate02, dIRef, dIMove2, [], [], 2);
    handles = fPlotDiff(handles, handles.Gate03, dIRef, dIMove3, [], [], 3);
    handles = fPlotDiff(handles, handles.Gate04, dIRef, dIMove4, [], [], 4);
    handles.sImageData = 'diffOri';
elseif eventdata.NewValue == handles.rb_DiffReg
    try axes(handles.Gate02); colorbar('off'); catch; end
    try axes(handles.Gate03); colorbar('off'); catch; end
    try axes(handles.Gate04); colorbar('off'); catch; end
    dIRef = handles.h.dImg(:,:,:,1);
    dIMoveReg2 = handles.h.dImgReg(:,:,:,2);
    dIMoveReg3 = handles.h.dImgReg(:,:,:,3);
    dIMoveReg4 = handles.h.dImgReg(:,:,:,4);
    handles = fPlotDiff(handles, handles.Gate02, dIRef, dIMoveReg2, [], [], 2);
    handles = fPlotDiff(handles, handles.Gate03, dIRef, dIMoveReg3, [], [], 3);
    handles = fPlotDiff(handles, handles.Gate04, dIRef, dIMoveReg4, [], [], 4);
    handles.sImageData = 'diffReg';    
elseif eventdata.NewValue == handles.rb_detJac
    handles = fSetCDataDiv(handles,handles.detJ2,handles.detJ3,handles.detJ4);
    handles.sImageData = 'detJac';
    
elseif eventdata.NewValue == handles.rb_divergence
    handles = fSetCDataDiv(handles,handles.divF2,handles.divF3,handles.divF4);
    handles.sImageData = 'divergence';
    
elseif   eventdata.NewValue == handles.rb_arrlength
    handles = fSetCDataDiv(handles,handles.arrL2,handles.arrL3,handles.arrL4);
    handles.sImageData = 'arrLength';
    
elseif   eventdata.NewValue == handles.rb_nablaLength
    set(handles.pop_Grad,'Visible','on');
    handles = fSetCDataDiv(handles,handles.gradL2,handles.gradL3,handles.gradL4);
    handles.sImageData = 'nablaLength';
    
elseif   eventdata.NewValue == handles.rb_isosurface
    set(handles.uitoggletool4,'Enable','on');
    set(handles.cb_showQuiver,'Value',0);
    set(handles.cb_showQuiver,'Enable','off');
    set(handles.cb_cut,'Visible','on');
    set(handles.pb_isoLevelVal,'Visible','on');
    set(handles.edLevelSet,'Visible','on');
    handles = fPlotIsosurface(handles, handles.Gate02, handles.level_set, 2, true);
    handles = fPlotIsosurface(handles, handles.Gate03, handles.level_set, 3, true);
    handles = fPlotIsosurface(handles, handles.Gate04, handles.level_set, 4, true);
    handles.sImageData = 'isosurface';
else
    errordlg('Unknown error')
    return;
end

if(strcmp(handles.sImageData,'isosurface'))
    linkprop([handles.Gate02, handles.Gate03, handles.Gate04],{'CameraPosition','CameraUpVector'});
    set(handles.text20, 'String',[num2str(handles.slice),'/',num2str(size(handles.h.dImg,3))]);
    set(handles.text21, 'String', sprintf('|df| = %.2f', handles.level_set ));
    set(handles.text22, 'String', sprintf('|df| = %.2f', handles.level_set ));
    set(handles.text23, 'String', sprintf('|df| = %.2f', handles.level_set ));
else
    linkaxes([handles.Gate01, handles.Gate02, handles.Gate03, handles.Gate04]);
    set(handles.text20, 'String',[num2str(handles.slice),'/',num2str(size(handles.h.dImg,3))]);
    set(handles.text21, 'String',[num2str(handles.slice),'/',num2str(size(handles.h.dImg,3))]);
    set(handles.text22, 'String',[num2str(handles.slice),'/',num2str(size(handles.h.dImg,3))]);
    set(handles.text23, 'String',[num2str(handles.slice),'/',num2str(size(handles.h.dImg,3))]);
end
guidata(hObject, handles)

%% ------------------------------------------------------------------------
% end callback functions
% -------------------------------------------------------------------------


%% ------------------------------------------------------------------------
% Window button functions
% -------------------------------------------------------------------------
function EvalGUI_WindowScrollWheelFcn(hObject, eventdata, handles)
% executes on scroll wheel click while the figure is in focus

% get image arrays and deformation fields
dIRef = handles.h.dImg(:,:,:,1);
dIMove2 = handles.h.dImg(:,:,:,2);
dIMove3 = handles.h.dImg(:,:,:,3);
dIMove4 = handles.h.dImg(:,:,:,4);
if(get(handles.rbForward,'Value') > 0)
    dFx2 = handles.h.SDeform(2).dFx;
    dFy2 = handles.h.SDeform(2).dFy;
    dFx3 = handles.h.SDeform(3).dFx;
    dFy3 = handles.h.SDeform(3).dFy;
    dFx4 = handles.h.SDeform(4).dFx;
    dFy4 = handles.h.SDeform(4).dFy;
else % backward
    dFx2 = handles.h.SDeform(2).dBx;
    dFy2 = handles.h.SDeform(2).dBy;
    dFx3 = handles.h.SDeform(3).dBx;
    dFy3 = handles.h.SDeform(3).dBy;
    dFx4 = handles.h.SDeform(4).dBx;
    dFy4 = handles.h.SDeform(4).dBy;
end

% get scroll events
lMoving = true;
if eventdata.VerticalScrollCount < 0
    handles.slice = max([1 handles.slice - 1]);
elseif(eventdata.VerticalScrollCount > 0)
    handles.slice = min([size(dIRef, 3) handles.slice + 1]);
else 
    lMoving = false;
end

% reference image
% set new image data in axes => always keep reference image shown
set(handles.hI1, 'CData', dIRef(:,:,handles.slice));
axes(handles.Gate01); colormap(handles.Gate01,'gray'); freezeColors;

% moving images
if(lMoving)
    if strcmp(handles.sImageData,'origImage')
        set(handles.hI(2), 'CData', dIMove2(:,:,handles.slice)); colormap(handles.Gate02,'gray'); %freezeColors;
        set(handles.hI(3), 'CData', dIMove3(:,:,handles.slice)); colormap(handles.Gate03,'gray'); %freezeColors;
        set(handles.hI(4), 'CData', dIMove4(:,:,handles.slice)); colormap(handles.Gate04,'gray'); %freezeColors;
    elseif strcmp(handles.sImageData,'regImage')
        dIMoveReg2 = handles.h.dImgReg(:,:,:,2); 
        dIMoveReg3 = handles.h.dImgReg(:,:,:,3); 
        dIMoveReg4 = handles.h.dImgReg(:,:,:,4); 
        set(handles.hI(2), 'CData', dIMoveReg2(:,:,handles.slice)); colormap(handles.Gate02,'gray'); %freezeColors;
        set(handles.hI(3), 'CData', dIMoveReg3(:,:,handles.slice)); colormap(handles.Gate03,'gray'); %freezeColors;
        set(handles.hI(4), 'CData', dIMoveReg4(:,:,handles.slice)); colormap(handles.Gate04,'gray'); %freezeColors;
    elseif strcmp(handles.sImageData,'diffOri')
        dIZero = zeros(size(dIRef,1),size(dIRef,2));
        set(handles.hI(2), 'CData', cat(3, dIRef(:,:,handles.slice), dIMove2(:,:,handles.slice), dIZero))
        set(handles.hI(3), 'CData', cat(3, dIRef(:,:,handles.slice), dIMove3(:,:,handles.slice), dIZero))
        set(handles.hI(4), 'CData', cat(3, dIRef(:,:,handles.slice), dIMove4(:,:,handles.slice), dIZero))
    elseif strcmp(handles.sImageData,'diffReg')
        dIMoveReg2 = handles.h.dImgReg(:,:,:,2);
        dIMoveReg3 = handles.h.dImgReg(:,:,:,3);
        dIMoveReg4 = handles.h.dImgReg(:,:,:,4);
        dIZero = zeros(size(dIRef,1),size(dIRef,2));
        set(handles.hI(2), 'CData', cat(3, dIRef(:,:,handles.slice), dIMoveReg2(:,:,handles.slice), dIZero))
        set(handles.hI(3), 'CData', cat(3, dIRef(:,:,handles.slice), dIMoveReg3(:,:,handles.slice), dIZero))
        set(handles.hI(4), 'CData', cat(3, dIRef(:,:,handles.slice), dIMoveReg4(:,:,handles.slice), dIZero))
    elseif strcmp(handles.sImageData,'detJac')
        set(handles.hI(2) ,'CData', handles.detJ2(:,:,handles.slice)); colormap(handles.Gate02,'jet'); %freezeColors;
        set(handles.hI(3) ,'CData', handles.detJ3(:,:,handles.slice)); colormap(handles.Gate03,'jet'); %freezeColors;
        set(handles.hI(4) ,'CData', handles.detJ4(:,:,handles.slice)); colormap(handles.Gate04,'jet'); %freezeColors;
    elseif strcmp(handles.sImageData,'divergence')
        set(handles.hI(2) ,'CData', handles.divF2(:,:,handles.slice)); colormap(handles.Gate02,'jet'); %freezeColors;
        set(handles.hI(3) ,'CData', handles.divF3(:,:,handles.slice)); colormap(handles.Gate03,'jet'); %freezeColors;
        set(handles.hI(4) ,'CData', handles.divF4(:,:,handles.slice)); colormap(handles.Gate04,'jet'); %freezeColors;
    elseif strcmp(handles.sImageData,'arrLength')
        set(handles.hI(2) ,'CData', handles.arrL2(:,:,handles.slice)); colormap(handles.Gate02,'jet'); %freezeColors;
        set(handles.hI(3) ,'CData', handles.arrL3(:,:,handles.slice)); colormap(handles.Gate03,'jet'); %freezeColors;
        set(handles.hI(4) ,'CData', handles.arrL4(:,:,handles.slice)); colormap(handles.Gate04,'jet'); %freezeColors;        
    elseif strcmp(handles.sImageData,'nablaLength')
        set(handles.hI(2) ,'CData', handles.gradL2(:,:,handles.slice)); colormap(handles.Gate02,'jet'); %freezeColors;
        set(handles.hI(3) ,'CData', handles.gradL3(:,:,handles.slice)); colormap(handles.Gate03,'jet'); %freezeColors;
        set(handles.hI(4) ,'CData', handles.gradL4(:,:,handles.slice)); colormap(handles.Gate04,'jet'); %freezeColors;
    elseif strcmp(handles.sImageData,'isosurface')
        % cut contour accordingly
        if(get(handles.cb_cut,'Value') > 0)
            handles = fPlotIsosurface(handles, handles.Gate02, handles.level_set, 2, false);
            handles = fPlotIsosurface(handles, handles.Gate03, handles.level_set, 3, false);
            handles = fPlotIsosurface(handles, handles.Gate04, handles.level_set, 4, false);
        else
            colormap(handles.Gate02, handles.tcmap);
            colormap(handles.Gate03, handles.tcmap);
            colormap(handles.Gate04, handles.tcmap);
        end
    else
        errordlg('Unknown error in scroll function')
    end

    % set new quiver data in axes
    try
        set(handles.hQ(2), 'UData', dFx2(1:handles.quiverFactor:end, 1:handles.quiverFactor:end, handles.slice), ...
            'VData', dFy2(1:handles.quiverFactor:end, 1:handles.quiverFactor:end, handles.slice));
        set(handles.hQ(3), 'UData', dFx3(1:handles.quiverFactor:end, 1:handles.quiverFactor:end, handles.slice), ...
            'VData', dFy3(1:handles.quiverFactor:end, 1:handles.quiverFactor:end, handles.slice));
        set(handles.hQ(4), 'UData', dFx4(1:handles.quiverFactor:end, 1:handles.quiverFactor:end, handles.slice), ...
            'VData', dFy4(1:handles.quiverFactor:end, 1:handles.quiverFactor:end, handles.slice));
    catch
    end

end

% display current slice number
if(strcmp(handles.sImageData,'isosurface'))
    set(handles.text20, 'String',[num2str(handles.slice),'/',num2str(size(dIRef,3))]);
    set(handles.text21, 'String', sprintf('|df| = %.2f', handles.level_set ));
    set(handles.text22, 'String', sprintf('|df| = %.2f', handles.level_set ));
    set(handles.text23, 'String', sprintf('|df| = %.2f', handles.level_set ));
else
    set(handles.text20, 'String',[num2str(handles.slice),'/',num2str(size(dIRef,3))]);
    set(handles.text21, 'String',[num2str(handles.slice),'/',num2str(size(dIRef,3))]);
    set(handles.text22, 'String',[num2str(handles.slice),'/',num2str(size(dIRef,3))]);
    set(handles.text23, 'String',[num2str(handles.slice),'/',num2str(size(dIRef,3))]);
end
guidata(hObject, handles)


function EvalGUI_WindowButtonDownFcn(hObject, eventdata, handles)
% Save starting parameters for brightness/contrast scaling
handles.dPosStart = get(gca, 'CurrentPoint');

handles.FButtonDown = 1;
handles.colMin = handles.colRange(1);
handles.colMax = handles.colRange(2);
guidata(hObject, handles)


function EvalGUI_WindowButtonMotionFcn(hObject, eventdata, handles)
% apply brightness/contrast scaling
try
    if handles.FButtonDown
        switch get(hObject, 'SelectionType')
            case 'extend'
                iD = get(gca, 'CurrentPoint') - handles.dPosStart; % Mouse distance travelled since button down
                
                % contrast and brightness
                handles.colWidth  = handles.colMax-handles.colMin;
                handles.colWidth  = handles.colMax.*exp(-iD(1,2)*0.02);
                handles.colCenter = (handles.colMax+handles.colMin)/2;
                handles.colCenter = handles.colCenter.*exp(iD(1,1)*0.02);
                handles.colRange  = [handles.colCenter-handles.colWidth/2, handles.colCenter+handles.colWidth/2];    
                caxis(handles.Gate01,handles.colRange);
                if(~strcmp(handles.sImageData, 'origImage') || ~strcmp(handles.sImageData, 'regImage'))
                    if strcmp(handles.sImageData, 'isosurface')
                        handles.colRange = [];
                    end
                    if(~isempty(handles.colRange))
                        caxis(handles.Gate02,handles.colRange);
                        caxis(handles.Gate03,handles.colRange);
                        caxis(handles.Gate04,handles.colRange);
                    else % isosurface 
                        colormap(handles.Gate02, handles.tcmap);
                        colormap(handles.Gate03, handles.tcmap);
                        colormap(handles.Gate04, handles.tcmap);
                    end
                end
                event.VerticalScrollCount = 0;
                EvalGUI_WindowScrollWheelFcn(hObject, event, handles); % faked workaround so that reference image is also directly scaled
                guidata(hObject, handles);
        end
    else
        return
    end
catch
end
guidata(hObject, handles)


function EvalGUI_WindowButtonUpFcn(hObject, eventdata, handles)
% get button release event
handles.FButtonDown = 0;
guidata(hObject, handles)


%% ------------------------------------------------------------------------
% end Window button functions
% -------------------------------------------------------------------------


%% ------------------------------------------------------------------------
% non-Callback functions
% -------------------------------------------------------------------------

function fPlotImageAndDefField(hObject, handles)
% plot images and deformation fields

% get images and deformation fields
dIRef = handles.h.dImg(:,:,:,1);
dIMove2 = handles.h.dImg(:,:,:,2);
dIMove3 = handles.h.dImg(:,:,:,3);
dIMove4 = handles.h.dImg(:,:,:,4);
dIMoveReg2 = handles.h.dImgReg(:,:,:,2);
dIMoveReg3 = handles.h.dImgReg(:,:,:,3);
dIMoveReg4 = handles.h.dImgReg(:,:,:,4);
if(get(handles.rbForward,'Value') > 0)
    dFx2 = handles.h.SDeform(2).dFx;
    dFy2 = handles.h.SDeform(2).dFy;
    dFx3 = handles.h.SDeform(3).dFx;
    dFy3 = handles.h.SDeform(3).dFy;
    dFx4 = handles.h.SDeform(4).dFx;
    dFy4 = handles.h.SDeform(4).dFy;
else % backward
    dFx2 = handles.h.SDeform(2).dBx;
    dFy2 = handles.h.SDeform(2).dBy;
    dFx3 = handles.h.SDeform(3).dBx;
    dFy3 = handles.h.SDeform(3).dBy;
    dFx4 = handles.h.SDeform(4).dBx;
    dFy4 = handles.h.SDeform(4).dBy;
end

% determine range of colours
% handles.colMax = nanmax([max(dIMove2(:))/5,max(dIMove3(:))/5,max(dIMove4(:))/5]);
handles.colMax = nanmax([max(dIRef(:)),max(dIMove2(:)),max(dIMove3(:)),max(dIMove4(:))]);
handles.colMin = 0;
handles.colRange = [handles.colMin, handles.colMax];
dAspectRatio = handles.h.SGeo.cVoxelsize{1};
% plot reference image
axes(handles.Gate01)
if(~isfield(handles,'hI1') || isempty(handles.hI1))
    cla
    handles.hI1 = imshow(dIRef(:,:,handles.slice),handles.colRange);
else
    set(handles.hI1, 'CData', dIRef(:,:,handles.slice));
    set(handles.Gate01, 'CLim', handles.colRange);
end
colormap(handles.Gate01,'gray');
if handles.nSCT == 2 % sag
    daspect([dAspectRatio(2) dAspectRatio(3) 1]);
elseif handles.nSCT == 3 % tra
    daspect([dAspectRatio(3) dAspectRatio(1) 1]);
else % cor
    daspect([dAspectRatio(1) dAspectRatio(2) 1]);
end
xlim(handles.Gate01,[0.5 size(dIRef,2)+0.5]);
ylim(handles.Gate01,[0.5 size(dIRef,1)+0.5]);
freezeColors;

% plot moving images and corresponding quiver arrows
if strcmp(handles.sImageData, 'origImage')
    handles.colMax = max(dIMove2(:))/3;
    handles = fPlotImg(handles, handles.Gate02, dIMove2, dFx2, dFy2, 2);
    handles = fPlotImg(handles, handles.Gate03, dIMove3, dFx3, dFy3, 3);
    handles = fPlotImg(handles, handles.Gate04, dIMove4, dFx4, dFy4, 4);
elseif strcmp(handles.sImageData, 'regImage')
    handles.colMax = max(dIMoveReg2(:))/3;
    handles = fPlotImg(handles, handles.Gate02, dIMoveReg2, dFx2, dFy2, 2);
    handles = fPlotImg(handles, handles.Gate03, dIMoveReg3, dFx3, dFy3, 3);
    handles = fPlotImg(handles, handles.Gate04, dIMoveReg4, dFx4, dFy4, 4);
elseif strcmp(handles.sImageData, 'diffOri')
    handles.colRange = [nanmin([dIRef(:); dIMove2(:)]), nanmax([dIRef(:); dIMove2(:)])];
    handles = fPlotDiff(handles, handles.Gate02, dIRef, dIMove2, dFx2, dFy2, 2);
    handles = fPlotDiff(handles, handles.Gate03, dIRef, dIMove3, dFx3, dFy3, 3);
    handles = fPlotDiff(handles, handles.Gate04, dIRef, dIMove4, dFx4, dFy4, 4);
elseif strcmp(handles.sImageData, 'diffReg')
    handles.colRange = [nanmin([dIRef(:); dIMoveReg2(:)]), nanmax([dIRef(:); dIMoveReg2(:)])];
    handles = fPlotDiff(handles, handles.Gate02, dIRef, dIMoveReg2, dFx2, dFy2, 2);
    handles = fPlotDiff(handles, handles.Gate03, dIRef, dIMoveReg3, dFx3, dFy3, 3);
    handles = fPlotDiff(handles, handles.Gate04, dIRef, dIMoveReg4, dFx4, dFy4, 4);
elseif strcmp(handles.sImageData, 'detJac')
    % clear old color bars
    try axes(handles.Gate02); colorbar('off'); catch; end
    try axes(handles.Gate03); colorbar('off'); catch; end
    try axes(handles.Gate04); colorbar('off'); catch; end
    handles.colRange = [min(handles.detJ2(:)),max(handles.detJ2(:))];
    handles = fPlotDiv(handles, handles.Gate02, handles.detJ2, dFx2, dFy2, 2);
    handles = fPlotDiv(handles, handles.Gate03, handles.detJ3, dFx3, dFy3, 3);
    handles = fPlotDiv(handles, handles.Gate04, handles.detJ4, dFx4, dFy4, 4);
elseif strcmp(handles.sImageData, 'divergence')
    try axes(handles.Gate02); colorbar('off'); catch; end
    try axes(handles.Gate03); colorbar('off'); catch; end
    try axes(handles.Gate04); colorbar('off'); catch; end
    handles.colRange = [min(handles.divF2(:)),max(handles.divF2(:))];
    handles = fPlotDiv(handles, handles.Gate02, handles.divF2, dFx2, dFy2, 2);
    handles = fPlotDiv(handles, handles.Gate03, handles.divF3, dFx3, dFy3, 3);
    handles = fPlotDiv(handles, handles.Gate04, handles.divF4, dFx4, dFy4, 4);
elseif   strcmp(handles.sImageData, 'arrLength')
    try axes(handles.Gate02); colorbar('off'); catch; end
    try axes(handles.Gate03); colorbar('off'); catch; end
    try axes(handles.Gate04); colorbar('off'); catch; end
    handles.colRange = [min(handles.arrL2(:)),max(handles.arrL2(:))];
    handles = fPlotDiv(handles, handles.Gate02, handles.arrL2, dFx2, dFy2, 2);
    handles = fPlotDiv(handles, handles.Gate03, handles.arrL3, dFx3, dFy3, 3);
    handles = fPlotDiv(handles, handles.Gate04, handles.arrL4, dFx4, dFy4, 4);
elseif   strcmp(handles.sImageData, 'nablaLength')
    try axes(handles.Gate02); colorbar('off'); catch; end
    try axes(handles.Gate03); colorbar('off'); catch; end
    try axes(handles.Gate04); colorbar('off'); catch; end
    handles.colRange = [min(handles.gradL2(:)),max(handles.gradL2(:))];
    handles = fPlotDiv(handles, handles.Gate02, handles.gradL2, dFx2, dFy2, 2);
    handles = fPlotDiv(handles, handles.Gate03, handles.gradL3, dFx3, dFy3, 3);
    handles = fPlotDiv(handles, handles.Gate04, handles.gradL4, dFx4, dFy4, 4);
    
elseif   strcmp(handles.sImageData, 'isosurface')    
    set(handles.uitoggletool4,'Enable','on');
    set(handles.cb_showQuiver,'Value',0);
    set(handles.cb_showQuiver,'Enable','off');
    set(handles.cb_cut,'Visible','on');
    set(handles.pb_isoLevelVal,'Visible','on');
    set(handles.edLevelSet,'Visible','on');
    handles = fPlotIsosurface(handles, handles.Gate02, handles.level_set, 2, true);
    handles = fPlotIsosurface(handles, handles.Gate03, handles.level_set, 3, true);
    handles = fPlotIsosurface(handles, handles.Gate04, handles.level_set, 4, true);    
else
    errordlg('Unknown error in plotImageAndDefField function')
end
linkaxes([handles.Gate01, handles.Gate02, handles.Gate03, handles.Gate04])
handles.delImg = 0;
% display current slice number
if(strcmp(handles.sImageData, 'isosurface'))
    set(handles.text20, 'String',[num2str(handles.slice),'/',num2str(size(dIRef,3))]);
    set(handles.text21, 'String', sprintf('|df| = %.2f', handles.level_set ));
    set(handles.text22, 'String', sprintf('|df| = %.2f', handles.level_set ));
    set(handles.text23, 'String', sprintf('|df| = %.2f', handles.level_set ));
else
    set(handles.text20, 'String',[num2str(handles.slice),'/',num2str(size(dIRef,3))]);
    set(handles.text21, 'String',[num2str(handles.slice),'/',num2str(size(dIRef,3))]);
    set(handles.text22, 'String',[num2str(handles.slice),'/',num2str(size(dIRef,3))]);
    set(handles.text23, 'String',[num2str(handles.slice),'/',num2str(size(dIRef,3))]);
end
guidata(hObject, handles)


function handles = fComputeDiv(handles)
% precompute deformation field values (divergence, det(jac), arrowLength)
if(get(handles.rbForward,'Value') > 0)
    dFx2 = handles.h.SDeform(2).dFx;
    dFy2 = handles.h.SDeform(2).dFy;
    dFz2 = handles.h.SDeform(2).dFz;
    dFx3 = handles.h.SDeform(3).dFx;
    dFy3 = handles.h.SDeform(3).dFy;
    dFz3 = handles.h.SDeform(3).dFz;
    dFx4 = handles.h.SDeform(4).dFx;
    dFy4 = handles.h.SDeform(4).dFy;
    dFz4 = handles.h.SDeform(4).dFz;
else % backward
    dFx2 = handles.h.SDeform(2).dBx;
    dFy2 = handles.h.SDeform(2).dBy;
    dFz2 = handles.h.SDeform(2).dBz;
    dFx3 = handles.h.SDeform(3).dBx;
    dFy3 = handles.h.SDeform(3).dBy;
    dFz3 = handles.h.SDeform(3).dBz;
    dFx4 = handles.h.SDeform(4).dBx;
    dFy4 = handles.h.SDeform(4).dBy;
    dFz4 = handles.h.SDeform(4).dBz;
end
dVoxSize = handles.h.SGeo.cVoxelsize(handles.iShown);

% determinant of 'jacobian' (numeric derivative: diff)
h = fwaitbar(0,'Preparing deformation fields. Please wait...'); st = 0; steps = 6;
[Fx2x, Fx2y, Fx2z] = gradient(dFx2,dVoxSize{2}(1),dVoxSize{2}(2),dVoxSize{2}(3));
[Fy2x, Fy2y, Fy2z] = gradient(dFy2,dVoxSize{2}(1),dVoxSize{2}(2),dVoxSize{2}(3));
[Fz2x, Fz2y, Fz2z] = gradient(dFz2,dVoxSize{2}(1),dVoxSize{2}(2),dVoxSize{2}(3));
handles.detJ2 = (Fx2x+1).*(Fy2y+1).*(Fz2z+1) + Fx2y.*Fy2z.*Fz2x + Fx2z.*Fz2y.*Fy2x ...
    - Fx2y.*Fy2x.*(Fz2z+1) - Fx2z.*Fz2x.*(Fy2y+1) - Fy2z.*Fz2y.*(Fx2x+1);
st= st+1;fwaitbar(st/steps,h);

[Fx3x, Fx3y, Fx3z] = gradient(dFx3,dVoxSize{3}(1),dVoxSize{3}(2),dVoxSize{3}(3));
[Fy3x, Fy3y, Fy3z] = gradient(dFy3,dVoxSize{3}(1),dVoxSize{3}(2),dVoxSize{3}(3));
[Fz3x, Fz3y, Fz3z] = gradient(dFz3,dVoxSize{3}(1),dVoxSize{3}(2),dVoxSize{3}(3));
handles.detJ3 = (Fx3x+1).*(Fy3y+1).*(Fz3z+1) + Fx3y.*Fy3z.*Fz3x + Fx3z.*Fz3y.*Fy3x ...
    - Fx3y.*Fy3x.*(Fz3z+1) - Fx3z.*Fz3x.*(Fy3y+1) - Fy3z.*Fz3y.*(Fx3x+1);
 
st= st+1;fwaitbar(st/steps,h);

[Fx4x, Fx4y, Fx4z] = gradient(dFx4,dVoxSize{4}(1),dVoxSize{4}(2),dVoxSize{4}(3));
[Fy4x, Fy4y, Fy4z] = gradient(dFy4,dVoxSize{4}(1),dVoxSize{4}(2),dVoxSize{4}(3));
[Fz4x, Fz4y, Fz4z] = gradient(dFz4,dVoxSize{4}(1),dVoxSize{4}(2),dVoxSize{4}(3));
handles.detJ4 = (Fx4x+1).*(Fy4y+1).*(Fz4z+1) + Fx4y.*Fy4z.*Fz4x + Fx4z.*Fz4y.*Fy4x ...
    - Fx4y.*Fy4x.*(Fz4z+1) - Fx4z.*Fz4x.*(Fy4y+1) - Fy4z.*Fz4y.*(Fx4x+1);

st= st+1;fwaitbar(st/steps,h);

handles.colRange = [min(handles.detJ2(:)/1), max(handles.detJ2(:)/1)];

% divergence
x = 1:dVoxSize{2}(1):dVoxSize{2}(1)*size(dFx2,1);
y = 1:dVoxSize{2}(2):dVoxSize{2}(2)*size(dFx2,2);
z = 1:dVoxSize{2}(3):dVoxSize{2}(3)*size(dFx2,3);
handles.divF2 = divergence(x,y,z, dFx2, dFy2, dFz2)*(1); st= st+1;fwaitbar(st/steps,h);
x = 1:dVoxSize{3}(1):dVoxSize{3}(1)*size(dFx3,1);
y = 1:dVoxSize{3}(2):dVoxSize{3}(2)*size(dFx3,2);
z = 1:dVoxSize{3}(3):dVoxSize{3}(3)*size(dFx3,3);
handles.divF3 = divergence(x,y,z, dFx3, dFy3, dFz3)*(1); st= st+1;fwaitbar(st/steps,h);
x = 1:dVoxSize{4}(1):dVoxSize{4}(1)*size(dFx4,1);
y = 1:dVoxSize{4}(2):dVoxSize{4}(2)*size(dFx4,2);
z = 1:dVoxSize{4}(3):dVoxSize{4}(3)*size(dFx4,3);
handles.divF4 = divergence(x,y,z, dFx4, dFy4, dFz4)*(1); st= st+1;fwaitbar(st/steps,h);
handles.colRange = [min(handles.divF2(:)/1.5),max(handles.divF2(:)/1.5)];  try close(h); catch; end;

% arrow length
F2 = dFx2.^2 + dFy2.^2 + dFz2.^2;
handles.arrL2 = F2.^(1/2);
F3 = dFx3.^2 + dFy3.^2 + dFz3.^2;
handles.arrL3 = F3.^(1/2);
F4 = dFx4.^2 + dFy4.^2 + dFz4.^2;
handles.arrL4 = F4.^(1/2);
handles.colRange = [min(handles.arrL2(:)/1.5),max(handles.arrL2(:)/1.5)];

% gradient of deformation field
iInd = get(handles.pop_Grad,'Value');
switch iInd
    case 1 % x
        handles.gradL2 = sqrt(Fx2x.^2 + Fy2x.^2 + Fz2x.^2);
    case 2 % y 
        handles.gradL2 = sqrt(Fx2y.^2 + Fy2y.^2 + Fz2y.^2);
    case 3 % z
        handles.gradL2 = sqrt(Fx2z.^2 + Fy2z.^2 + Fz2z.^2);
end

switch iInd
    case 1 % x
        handles.gradL3 = sqrt(Fx3x.^2 + Fy3x.^2 + Fz3x.^2);
    case 2 % y 
        handles.gradL3 = sqrt(Fx3y.^2 + Fy3y.^2 + Fz3y.^2);
    case 3 % z
        handles.gradL3 = sqrt(Fx3z.^2 + Fy3z.^2 + Fz3z.^2);
end

switch iInd
    case 1 % x
        handles.gradL4 = sqrt(Fx4x.^2 + Fy4x.^2 + Fz4x.^2);
    case 2 % y 
        handles.gradL4 = sqrt(Fx4y.^2 + Fy4y.^2 + Fz4y.^2);
    case 3 % z
        handles.gradL4 = sqrt(Fx4z.^2 + Fy4z.^2 + Fz4z.^2);
end
% handles.colRange = [min(handles.gradL2(:)/1.5),max(handles.gradL2(:)/1.5)];

% isosurface
handles.max_amp = ceil(max([max(dFx2(:)), max(dFx3(:)), max(dFx4(:)), max(dFy2(:)), max(dFy3(:)), max(dFy4(:)), max(dFz2(:)), max(dFz3(:)), max(dFz4(:))]));


function handles = fPlotDiv(handles, GateAx, dData, dFx, dFy, n)
% plot divergence, det(jac) and arrowLength

% example: plotDiv(handles, handles.Gate02, handles.detJ2, 2)
% get position, size and density of quiver arrows
if(~isempty(dFx))
    quiverFactor = handles.quiverFactor;
    quiverScale = handles.quiverScale;
    iX = 1:quiverFactor:size(dFx, 2);
    iY = 1:quiverFactor:size(dFy, 1);
end

dAspectRatio = handles.h.SGeo.cVoxelsize{handles.iShown(n)};
axes(GateAx)
if(~isfield(handles,'hI') || length(handles.hI) < n || handles.hI(n) < 0)
    if(isfield(handles,'hI') && ishandle(handles.hI(n)))
        delete(handles.hI(n));
    end        
    handles.hI(n) = imshow(dData(:,:,handles.slice),handles.colRange);
    xLimits = [0.5 size(dData,2)+0.5];
    yLimits = [0.5 size(dData,1)+0.5];
else
    set(handles.hI(n), 'CData', dData(:,:,handles.slice));
    set(GateAx, 'CLim', handles.colRange);
    xLimits = get(GateAx,'XLim');
    yLimits = get(GateAx,'YLim');
end
divPos = get(GateAx,'position');
colormap(GateAx,jet);
colorbar('eastoutside','Color',[1 1 1]);
set(GateAx,'position',divPos);
if handles.nSCT == 2
    daspect([dAspectRatio(2) dAspectRatio(3) 1]);
elseif handles.nSCT == 3
    daspect([dAspectRatio(3) dAspectRatio(1) 1]);
else
    daspect([dAspectRatio(1) dAspectRatio(2) 1]);
end
xlim(GateAx,xLimits);
ylim(GateAx,yLimits);
if(get(handles.cb_showQuiver,'Value') && ~isempty(dFx))
    hold on
    if(~isfield(handles, 'hQ') || length(handles.hQ) < n || handles.hQ(n) < 0)
        if(isfield(handles,'hQ') && length(handles.hQ) >= n && ishandle(handles.hQ(n)))
            delete(handles.hQ(n));
        end 
        handles.hQ(n) = quiver(iX, iY, quiverScale*dFx(1:quiverFactor:end, 1:quiverFactor:end, handles.slice),...
            dFy(1:quiverFactor:end, 1:quiverFactor:end, handles.slice), 0);
    else
        set(handles.hQ(n), 'XData', iX, 'YData', iY, 'UData', quiverScale*dFx(1:quiverFactor:end, 1:quiverFactor:end, handles.slice), 'VData', dFy(1:quiverFactor:end, 1:quiverFactor:end, handles.slice));
    end        
    set(handles.hQ(n), 'Linewidth', 1.5, 'Color', 'y');
    hold off
end
freezeColors;


function handles = fPlotDiff(handles, GateAx, dData1, dData2, dFx, dFy, n)
% plot difference images

% get position, size and density of quiver arrows
if(~isempty(dFx))
    quiverFactor = handles.quiverFactor;
    quiverScale = handles.quiverScale;
    iX = 1:quiverFactor:size(dFx, 2);
    iY = 1:quiverFactor:size(dFy, 1);
end

dAspectRatio = handles.h.SGeo.cVoxelsize{handles.iShown(n)};
dIZero = zeros(size(dData1,1),size(dData1,2));
axes(GateAx)
if(~isfield(handles,'hI') || length(handles.hI) < n || handles.hI(n) < 0)
    if(isfield(handles,'hI') && ishandle(handles.hI(n)))
        delete(handles.hI(n));
    end 
    handles.hI(n) = imshow(cat(3,dData1(:,:,handles.slice), dData2(:,:,handles.slice), dIZero), handles.colRange);
    xLimits = [0.5 size(dData1,2)+0.5];
    yLimits = [0.5 size(dData1,1)+0.5];
    xlim(GateAx,xLimits);
    ylim(GateAx,yLimits);
else
    set(handles.hI(n), 'CData', cat(3,dData1(:,:,handles.slice), dData2(:,:,handles.slice), dIZero));
    set(GateAx, 'CLim', handles.colRange);
%     xLimits = get(GateAx,'XLim');
%     yLimits = get(GateAx,'YLim');
end
% divPos = get(GateAx,'position'); 
% set(GateAx,'position',divPos);
if handles.nSCT == 2
    daspect([dAspectRatio(2) dAspectRatio(3) 1]);
elseif handles.nSCT == 3
    daspect([dAspectRatio(3) dAspectRatio(1) 1]);
else
    daspect([dAspectRatio(1) dAspectRatio(2) 1]);
end
if(get(handles.cb_showQuiver,'Value') && ~isempty(dFx))
    hold on
    if(~isfield(handles, 'hQ') || length(handles.hQ) < n || handles.hQ(n) < 0)
        if(isfield(handles,'hQ') && length(handles.hQ) >= n && ishandle(handles.hQ(n)))
            delete(handles.hQ(n));
        end 
        handles.hQ(n) = quiver(iX, iY, quiverScale*dFx(1:quiverFactor:end, 1:quiverFactor:end, handles.slice),...
            dFy(1:quiverFactor:end, 1:quiverFactor:end, handles.slice), 0);
    else
        set(handles.hQ(n), 'XData', iX, 'YData', iY, 'UData', quiverScale*dFx(1:quiverFactor:end, 1:quiverFactor:end, handles.slice), 'VData', dFy(1:quiverFactor:end, 1:quiverFactor:end, handles.slice));
    end        
    set(handles.hQ(n), 'Linewidth', 1.5, 'Color', 'y');
    hold off
end
freezeColors;


function handles = fPlotImg(handles, GateAx, dData, dFx, dFy, n)
% plot images and deformation field
% example: handles = plotImg(handles, handles.Gate02, dIMove2, dFx2, dFy2, 2)

% get position, size and density of quiver arrows
quiverFactor = handles.quiverFactor;
quiverScale = handles.quiverScale;
iX = 1:quiverFactor:size(dFx, 2);
iY = 1:quiverFactor:size(dFy, 1);

dAspectRatio = handles.h.SGeo.cVoxelsize{handles.iShown(n)};
axes(GateAx)
if(~isfield(handles,'hI') || length(handles.hI) < n || handles.hI(n) < 0)
    if(isfield(handles,'hI') && length(handles.hI) >= n && ishandle(handles.hI(n)))
        delete(handles.hI(n));
    end 
    handles.hI(n) = imshow(dData(:,:,handles.slice),handles.colRange);
    xLimits = [0.5 size(dData,2)+0.5];
    yLimits = [0.5 size(dData,1)+0.5];
    xlim(GateAx,xLimits);
    ylim(GateAx,yLimits);
else
    set(handles.hI(n), 'CData', dData(:,:,handles.slice));
    set(GateAx, 'CLim', handles.colRange);
%     xLimits = get(GateAx,'XLim');
%     yLimits = get(GateAx,'YLim');
end
if handles.nSCT == 2
    daspect([dAspectRatio(2) dAspectRatio(3) 1]);
elseif handles.nSCT == 3
    daspect([dAspectRatio(3) dAspectRatio(1) 1]);
else
    daspect([dAspectRatio(1) dAspectRatio(2) 1]);
end

if get(handles.cb_showQuiver,'Value')
    hold on
    if(~isfield(handles, 'hQ') || length(handles.hQ) < n || handles.hQ(n) < 0)
        if(isfield(handles,'hQ') && length(handles.hQ) >= n && ishandle(handles.hQ(n)))
            delete(handles.hQ(n));
        end 
        handles.hQ(n) = quiver(iX, iY, quiverScale*dFx(1:quiverFactor:end, 1:quiverFactor:end, handles.slice),...
            dFy(1:quiverFactor:end, 1:quiverFactor:end, handles.slice), 0);
    else
        set(handles.hQ(n), 'XData', iX, 'YData', iY, 'UData', quiverScale*dFx(1:quiverFactor:end, 1:quiverFactor:end, handles.slice), 'VData', dFy(1:quiverFactor:end, 1:quiverFactor:end, handles.slice));
    end        
    set(handles.hQ(n), 'Linewidth', 1.5, 'Color', 'y');
    hold off
end
freezeColors;


function handles = fPlotSelectionChange(handles, haxes, dData, colRange, n)
% plot images for selection change
axes(haxes);
if(~isfield(handles,'hI') || length(handles.hI) < n || handles.hI(n) < 0)
    if(isfield(handles,'hIso'))
        if(ishandle(handles.hIso(1,n))); delete(handles.hIso(1,n)); end
        if(ishandle(handles.hIso(2,n))); delete(handles.hIso(2,n)); end
    end
    handles.hI(n) = imshow(dData, colRange); colormap(haxes,'gray'); freezeColors;
else
    set(handles.hI(n), 'CData', dData); caxis(colRange); colormap(haxes,'gray'); freezeColors;
end


function handles = fSetCDataDiv(handles,dData2, dData3, dData4)
% update divergence, det(jac) and arrow length plots

% example: handles = fSetCDataDiv(handles,handles.detJ2,handles.detJ3,handles.detJ4)

try axes(handles.Gate02); colorbar('off'); catch; end
try axes(handles.Gate03); colorbar('off'); catch; end
try axes(handles.Gate04); colorbar('off'); catch; end
% draw det(Jac) / divergence / arrow length
handles.colRange = [nanmin([min(dData2(:)/1.0),min(dData3(:)/1.0),min(dData4(:)/1.0)]),nanmax([max(dData2(:)/1.0),max(dData3(:)/1.0),max(dData4(:)/1.0)])];
axes(handles.Gate02); n = 2;
if(~isfield(handles,'hI') || length(handles.hI) < n || handles.hI(n) < 0)
    if(isfield(handles,'hI') && ishandle(handles.hI(n)))
        delete(handles.hI(n));
    end 
    handles.hI(n) = imshow(dData2(:,:,handles.slice), handles.colRange); colormap(handles.Gate02,'jet');
else
    set(handles.hI(n), 'CData', dData2(:,:,handles.slice)); caxis(handles.colRange); colormap(handles.Gate02,'jet');
end
divPos = get(handles.Gate02,'position'); colorbar('eastoutside'); set(handles.Gate02,'position',divPos); freezeColors;

axes(handles.Gate03); n = 3; 
if(~isfield(handles,'hI') || length(handles.hI) < n || handles.hI(n) < 0)
    if(isfield(handles,'hI') && ishandle(handles.hI(n)))
        delete(handles.hI(n));
    end 
    handles.hI(n) = imshow(dData3(:,:,handles.slice), handles.colRange); colormap(handles.Gate03,'jet');
else
    set(handles.hI(n), 'CData', dData3(:,:,handles.slice)); caxis(handles.colRange); colormap(handles.Gate03,'jet');
end
divPos = get(handles.Gate03,'position'); colorbar('eastoutside'); set(handles.Gate03,'position',divPos); freezeColors;

axes(handles.Gate04); n = 4;
if(~isfield(handles,'hI') || length(handles.hI) < n || handles.hI(n) < 0)
    if(isfield(handles,'hI') && ishandle(handles.hI(n)))
        delete(handles.hI(n));
    end 
    handles.hI(n) = imshow(dData4(:,:,handles.slice), handles.colRange); colormap(handles.Gate04,'jet');
else
    set(handles.hI(n), 'CData', dData4(:,:,handles.slice)); caxis(handles.colRange); colormap(handles.Gate04,'jet');
end
divPos = get(handles.Gate04,'position'); colorbar('eastoutside'); set(handles.Gate04,'position',divPos); freezeColors;


function handles = fPlotIsosurface(handles, haxes, level_set, n, lNew)
% plot isosurface
if(n == 2)
    handles.hWaitbar = fwaitbar(0, 'Calculating isosurfaces');
end
steps = 9; st = (n-2)* 3 + 1;

% Swap ux and uy due to our notation:
if( get(handles.rbForward,'Value') > 0) % forward
    u = {handles.h.SDeform(n).dFy, handles.h.SDeform(n).dFx, handles.h.SDeform(n).dFz};
else % backward
    u = {handles.h.SDeform(n).dBy, handles.h.SDeform(n).dBx, handles.h.SDeform(n).dBz};
end

if(get(handles.cb_cut,'Value') > 0) % cut accordingly
	% shown orientation -> already included in SDeform.dXX
    u = {u{1}(:,:,handles.slice:end), u{2}(:,:,handles.slice:end), u{3}(:,:,handles.slice:end)};
end
[M,N,P] = size(u{1});
[x,y,z] = meshgrid(1:M,1:N,1:P);

amp = sqrt(u{1}.^2 + u{2}.^2 + u{3}.^2);

% max_amp = ceil(max(amp(:)));
max_amp = handles.max_amp;
C_face = round(63*level_set/max_amp)+1;

axes(haxes);
% delete prior object handles
if(isfield(handles,'hQ') && handles.hQ(n) > 0 && ishandle(handles.hQ(n)))
    delete(handles.hQ(n));
end
handles.hQ(n) = -1; % invalidate handle
if(isfield(handles,'hI') && handles.hI(n) > 0 && ishandle(handles.hI(n)))
    delete(handles.hI(n));
end
handles.hI(n) = -1; % invalidate handle
if(isfield(handles,'hIso'))
    if(size(handles.hIso,2) >= n)
        if(handles.hIso(1,n) > 0 && ishandle(handles.hIso(1,n)))
            delete(handles.hIso(1,n));
        end
        if(handles.hIso(2,n) > 0 && ishandle(handles.hIso(2,n)))
            delete(handles.hIso(2,n));
        end
    end
end
handles.hIso(1,n) = patch(isosurface(x,y,z,amp,level_set));
fwaitbar(st/steps,handles.hWaitbar); st = st+1;
set(handles.hIso(1,n),'FaceColor',handles.tcmap(C_face,:),'EdgeColor','none', 'FaceAlpha', 1);
isonormals(x,y,z,amp,handles.hIso(1,n));
fwaitbar(st/steps,handles.hWaitbar); st = st+1;
set(haxes, 'YDir', 'reverse','ZDir', 'reverse');
daspect(haxes,[1 1 1]);
handles.hIso(2,n) = patch(isocaps(x,y,z,amp,level_set),'FaceColor','interp','EdgeColor','none');
fwaitbar(st/steps,handles.hWaitbar); % st = st+1;
colormap(haxes,handles.tcmap); freezeColors();
if(lNew)
    divPos = get(haxes,'position'); colorbar('eastoutside'); set(haxes,'position',divPos);
end
caxis([0,max_amp]);

xlabel(haxes, 'x', 'Color', 'w');
ylabel(haxes, 'y', 'Color', 'w');
zlabel(haxes, 'z', 'Color', 'w');

% axis tight
view(handles.dView);
camlight(handles.dCamlight(1),handles.dCamlight(2));
lighting gouraud;

if(n == 4) % delete waitbar
    if(isfield(handles, 'hWaitbar'))
        close(handles.hWaitbar);
        handles = rmfield(handles,'hWaitbar');
    end
end


function fEvalAuto(hObject, handles)
% automatic evaluation mode
h = fwaitbar(0,'Evaluating registration results in automatic mode'); steps = 3; st=1;

pb_Intensity_Callback(hObject, [], handles);
fwaitbar(st/steps, h); st = st + 1;
eventdata.autoMode = true;
pb_SegLung_Callback(hObject, eventdata, handles)
fwaitbar(st/steps, h); st = st + 1;
pb_save_Callback(hObject, 0, handles);
fwaitbar(st/steps, h);
try close(h); catch; end;
return


function [SOut, sPathname, sFilename] = fSelectStruct(handles)
% load data from mat file
[sFilename, sPathname] = uigetfile(handles.SPaths.sResults, 'Select a result mat-file from RegGUI.');
if sFilename == 0
    SOut=[];sPathname=[];sFilename=[];
    return
end
h = fwaitbar(0,'Loading image data ...');
load([sPathname, sFilename]); fwaitbar(1,h);
try close(h); catch; end;
try 
    SOut = SRegiResult;
catch
    errordlg('This is not a result mat-file from RegGUI');
    return;
end
