function varargout = RegGUI(varargin)
% GUI framework to perform registration of gated MR images with different
% registration algorithms.
% Batch Processing.
% input data to be loaded in the GUI: mat-files containing dImg variable
% with image data (4D) and
% optional SGeo with dVoxelsize (3x1) and sOrientation (e.g. 'coronal')
%
% For detailed documentation see "Documentation of GUIs for MR image
% registration" (pdf).
%
% Before the first start check the paths in th RegGUI_OpeningFcn
%
% -------------------------------------------------------------------------
% (c) 2015: Thomas Kuestner, Verena Neumann
% -------------------------------------------------------------------------

% Last Modified by GUIDE v2.5 22-Apr-2016 09:53:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @RegGUI_OpeningFcn, ...
    'gui_OutputFcn',  @RegGUI_OutputFcn, ...
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


% --- Executes just before RegGUI is made visible.
function RegGUI_OpeningFcn(hObject, eventdata, handles, varargin)
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

if(~exist('standardVoxelsize','var'))
    handles.dVoxelsize = [1,1,1];
else
    handles.dVoxelsize = standardVoxelsize;
    set(handles.voxel_x,'String',sprintf('%.2f',handles.dVoxelsize(1)));
    set(handles.voxel_y,'String',sprintf('%.2f',handles.dVoxelsize(2)));
    set(handles.voxel_z,'String',sprintf('%.2f',handles.dVoxelsize(3)));
end

handles.SPaths = SPaths;
handles.currpath = currpath;
handles.sDemoDataPath = [currpath,filesep,'example'];
set(handles.txtOutput,'String',handles.SPaths.sResults);

% set some folders on the path
addpath(genpath([currpath,filesep,'io']));
addpath(genpath([currpath,filesep,'metrics']));
addpath(genpath([currpath,filesep,'registration']));
addpath(genpath([currpath,filesep,'segmentation']));
addpath(genpath([currpath,filesep,'utils']));


%% icons
try % Try to apply a nice icon
    warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
    jframe = get(hObject, 'javaframe');
    jIcon = javax.swing.ImageIcon([currpath, filesep, 'icons', filesep, 'REGGUI_icon.png']);
    pause(0.001);
    jframe.setFigureIcon(jIcon);
    clear jframe jIcon
catch
    warning('Could not apply a nice icon to the figure :(');
end
% rightarrow
set(handles.axArrowRight,'xtick',[],'ytick',[]);
set(handles.axArrowRight2,'xtick',[],'ytick',[]);
dIcon = 240/255*zeros(100,100);
dIcon(45:55,1:85) = 1;
dIcon(50,100) = 1;
for i=1:15
   dIcon(50-i:50+i,100-i) = 1; 
end
axes(handles.axArrowRight);
imshow(dIcon);
axes(handles.axArrowRight2);
imshow(dIcon);

% ok button
dImage = double(imread(['icons',filesep,'checkmark.png']))./255;
if size(dImage, 3) == 1, dImage = repmat(dImage, [1 1 3]); end
set(handles.buttonOK, 'CData', dImage/max(dImage(:)));
set(handles.buttonOK, 'Visible','off');

% cancel button
dImage = double(imread(['icons',filesep,'cancel.png']))./255;
if size(dImage, 3) == 1, dImage = repmat(dImage, [1 1 3]); end
set(handles.buttonCancel, 'CData', dImage/max(dImage(:)));
set(handles.buttonCancel, 'Visible','off');

% images button 4D
dImage = double(imread(['icons',filesep,'folder_plus.png']))./255;
if size(dImage, 3) == 1, dImage = repmat(dImage, [1 1 3]); end
set(handles.button4D, 'CData', dImage/max(dImage(:)));

% images button workspace
dImage = double(imread(['icons',filesep,'doc_plus.png']))./255;
if size(dImage, 3) == 1, dImage = repmat(dImage, [1 1 3]); end
set(handles.buttonWorkspace, 'CData', dImage/max(dImage(:)));

% load registration button
dImage = double(imread(['icons',filesep,'open.png']))./255;
if size(dImage, 3) == 1, dImage = repmat(dImage, [1 1 3]); end
set(handles.pb_LoadReg, 'CData', dImage/max(dImage(:)));

% clear images button
dImage = double(imread(['icons',filesep,'round_delete.png']))./255;
if size(dImage, 3) == 1, dImage = repmat(dImage, [1 1 3]); end
set(handles.clearFigures, 'CData', dImage/max(dImage(:)));

% edit parameter file button
dImage = double(imread(['icons',filesep,'edit.png']))./255;
if size(dImage, 3) == 1, dImage = repmat(dImage, [1 1 3]); end
set(handles.pb_editParam, 'CData', dImage/max(dImage(:)));

%% prestore axes
handles.quiverScale = 1;                % quiver arrows
handles.quiverFac = 10;
dFx = zeros(256,256,1);
iX = 1:handles.quiverFac:size(dFx, 2);
iY = 1:handles.quiverFac:size(dFx, 1);

axes(handles.Gate01);
cla;
handles.hI1 = imshow(zeros(256,256,1),[0 1]);

axes(handles.Gate02);
cla;
handles.hI2 = imshow(zeros(256,256,1),[0 1]);
hold on;
handles.hQ2=quiver(iX, iY, dFx(1:handles.quiverFac:end, 1:handles.quiverFac:end, 1),...
    dFx(1:handles.quiverFac:end, 1:handles.quiverFac:end, 1), handles.quiverScale);
set(handles.hQ2, 'Linewidth', 1.5, 'Color', 'y', 'Visible', 'off');

axes(handles.Gate03);
cla;
handles.hI3 = imshow(zeros(256,256,1),[0 1]);
hold on;
handles.hQ3=quiver(iX, iY, dFx(1:handles.quiverFac:end, 1:handles.quiverFac:end, 1),...
    dFx(1:handles.quiverFac:end, 1:handles.quiverFac:end, 1), handles.quiverScale);
set(handles.hQ3, 'Linewidth', 1.5, 'Color', 'y', 'Visible', 'off');

axes(handles.Gate04);
cla;
handles.hI4 = imshow(zeros(256,256,1),[0 1]);
hold on;
handles.hQ4=quiver(iX, iY, dFx(1:handles.quiverFac:end, 1:handles.quiverFac:end, 1),...
    dFx(1:handles.quiverFac:end, 1:handles.quiverFac:end, 1), handles.quiverScale);
set(handles.hQ4, 'Linewidth', 1.5, 'Color', 'y', 'Visible', 'off');

%% check for existing registration algorithms
sPaths = dir([currpath,filesep,'registration']);
lPath = cell2mat({sPaths(:).isdir});
cPaths = {sPaths(lPath).name};
cPaths = cPaths(~strcmp(cPaths,'parameterFiles') & ~strcmp(cPaths,'.') & ~strcmp(cPaths,'..'));
set(handles.selectRegiMethod,'String',cPaths);
handles.dRegMapping = {'elastix', 'halar', 'GRICS', 'LAP', 'demons'; ...
                        1       , 2      , 3      , 4    , 5      };

%% set some default values
handles.slice = ones(2,4); % first line: shown slice, second line: max slice
handles.lOpen = true(1,4);
handles.iIndGlobal = zeros(1,4); % global index into handles.loadedImg
handles.iMinShow = 15;

handles.output = hObject;
handles.nRegMethod = handles.dRegMapping{2,strcmp(handles.dRegMapping(1,:),cPaths{get(handles.selectRegiMethod,'Value')})};
handles.sRegMethod = handles.dRegMapping{1,strcmp(handles.dRegMapping(1,:),cPaths{get(handles.selectRegiMethod,'Value')})};
handles.sParFile = 'NA';
handles.FAutoEval = 0;
handles.doubleClick = false(2,1); % listImages, listRegs
handles.allReg = []; % containing all registration information
handles.loadedImg = []; % containing all info about loaded images

%% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = RegGUI_OutputFcn(hObject, eventdata, handles)
%% customize edit textboxes
try
    jScrollPane = findjobj(handles.listImages);
%     jViewPort = jScrollPane.getViewport;
%     jEditbox = jViewPort.getComponent(0);
%     jEditbox.setWrapping(false);    
    set(jScrollPane,'HorizontalScrollBarPolicy',30); % jScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED;

    jScrollPane = findjobj(handles.listRegs);
%     jViewPort = jScrollPane.getViewport;
%     jEditbox = jViewPort.getComponent(0);
%     jEditbox.setWrapping(false);   
    set(jScrollPane,'HorizontalScrollBarPolicy',30); % jScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED;

catch 
end

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes when user attempts to close RegGUI.
function RegGUI_CloseRequestFcn(hObject, eventdata, handles)
% save GUIPreference
try
currpath = fileparts(mfilename('fullpath'));
if(isfield(handles,'SPaths'))
    SPaths = handles.SPaths;
    standardVoxelsize = [1 1 1];
    save([currpath, filesep, 'GUIPreferences.mat'],'SPaths','standardVoxelsize','-append');
end
delete(hObject);
catch 
end


% --- Executes on scroll wheel click while the figure is in focus.
function RegGUI_WindowScrollWheelFcn(hObject, eventdata, handles)

for iJ=1:4
    if(handles.iIndGlobal(iJ) == 0 || ~(isempty(handles.loadedImg) || isempty(handles.allReg))) % TODO: erlaube, dass man auch durch die regErgebnisse scrollen kann! => negativ für regi und positiv für loadedImg
        continue;
    end
    
    lPlotDF = false;
    if(handles.iIndGlobal(iJ) > 0) % from loadedImg
        dImg = scaleImg(handles.loadedImg{handles.iIndGlobal(iJ)}.dImg,[0 1]);            
    else % from listReg
        dImg = scaleImg(handles.allReg{-handles.iIndGlobal(iJ),1}{iJ}.dImg,[0 1]);
        if(isfield(handles.allReg{-handles.iIndGlobal(iJ),3},'SDeform') && ~isempty(handles.allReg{-handles.iIndGlobal(iJ),3}.SDeform(iJ).dFx))
            dFx = handles.allReg{-handles.iIndGlobal(iJ),3}.SDeform(iJ).dFx;
            dFy = handles.allReg{-handles.iIndGlobal(iJ),3}.SDeform(iJ).dFy;
            iX = 1:handles.quiverFac:size(dFx, 2);
            iY = 1:handles.quiverFac:size(dFy, 1);
            lPlotDF = true;
        end  
    end
    % set new slice number
    if eventdata.VerticalScrollCount < 0
        handles.slice(1,iJ) = max([1 handles.slice(1,iJ)-1]);
    else
        handles.slice(1,iJ) = min([size(dImg, 3) handles.slice(1,iJ)+1]);
    end
    eval(sprintf('set(handles.hI%d, ''CData'', dImg(:,:,handles.slice(1,iJ)));', iJ));
    eval(sprintf('set(handles.sliceNo%d, ''String'', ''%02d/%02d'');', iJ, handles.slice(1,iJ), handles.slice(2,iJ)));
    if(lPlotDF)
        eval(sprintf('set(handles.hQ%d, ''XData'', iX, ''YData'', iY, ''UData'', dFx(1:handles.quiverFac:end, 1:handles.quiverFac:end, handles.slice(1,iJ)), ''VData'', dFy(1:handles.quiverFac:end, 1:handles.quiverFac:end, handles.slice(1,iJ)), ''Color'', ''y'', ''Visible'', ''on'');',iJ));
    end    
%     eval(sprintf('set(handles.hI%d, ''CData'', dImg(:,:,handles.slice(1,iJ)));', iJ));
end
guidata(hObject, handles)

%% ------------------------------------------------------------------------
% GUI functionalities (callbacks)
% -------------------------------------------------------------------------



function fLoadImg(hObject,handles,sSelectionType)
% load images from files and folders (parent)
if strcmp(sSelectionType, 'normal')
   [sFilename, sPath, sExt] = uigetfile( ...
            {'*.*', 'All Files'; ...
            '*.dcm; *.DCM; *.ima; *.IMA; *.mat; *.MAT; *.jpg; *.jpeg; *.JPG; *.JPEG; *.tif; *.tiff; *.TIF; *.TIFF; *.gif; *.GIF; *.bmp; *.BMP; *.png; *.PNG; *.nii; *.NII; *.gipl; *.GIPL; *.mhd; *.MHD', 'All images'; ...
            '*.mat; *.MAT', 'Matlab File (*.mat)'; ...
            '*.jpg; *.jpeg; *.JPG; *.JPEG', 'JPEG-Image (*.jpg)'; ...
            '*.tif; *.tiff; *.TIF; *.TIFF;', 'TIFF-Image (*.tif)'; ...
            '*.gif; *.GIF', 'Gif-Image (*.gif)'; ...
            '*.bmp; *.BMP', 'Bitmaps (*.bmp)'; ...
            '*.png; *.PNG', 'Portable Network Graphics (*.png)'; ...
            '*.dcm; *.DCM; *.ima; *.IMA', 'DICOM Files (*.dcm, *.ima)'; ...
            '*.nii; *.NII', 'NifTy Files (*.nii)'; ...
            '*.gipl; *.GIPL', 'Guys Image Processing Lab Files (*.gipl)'; ...
            '*.mhd; *.MHD', 'MetaImages (*.mhd)'}, ...
            'OpenLocation'  , handles.SPaths.sData, ...
            'Multiselect'   , 'on');
        if isnumeric(sPath), return, end;   % Dialog aborted 
        if(length(unique(sExt)) > 1)
            errordlg('no mixed file extension loading allowed!');
            return;
        end 
    csFilenames{1} = sFilename;
else
    sPath = uigetdir(handles.SPaths.sData);
    if isnumeric(sPath), return, end;

    sPath = [sPath, filesep];
    SFiles = dir(sPath);
    SFiles = SFiles(~[SFiles.isdir]);
    csFilenames = cell(length(SFiles), 1);
    for i = 1:length(SFiles), csFilenames{i} = SFiles(i).name; end
end

sExtFirst = 0;
for iI=1:length(csFilenames)      
    [~,~,sExtCurr] = fileparts(csFilenames{iI});   
    
    switch sExtCurr
        case {'.mat', '.MAT'}
            sExtCurr = 1;
            lCont = true;
            
        case {'.dcm', '.DCM', '.ima', '.IMA'}
            sExtCurr = 2;
            lCont = false;
            
        case {'.mhd', '.MHD'}
            sExtCurr = 3;
            lCont = true;
            
        case {'.gipl', '.GIPL'}
            sExtCurr = 4;
            lCont = true;
            
        case {'.nii', '.NII'}
            sExtCurr = 5;
            lCont = false;
            
        case {'.jpg', '.jpeg', '.tif', '.tiff', '.gif', '.GIF', '.bmp', '.BMP', '.png', '.PNG'}
            sExtCurr = 6;
            lCont = true;
    end
    
    if(iI == 1), sExtFirst = sExtCurr; end    
    if(sExtCurr ~= sExtFirst)
        errordlg('no mixed file extension loading allowed!');
        return;
    end
    
    imgInfo =  fGetImg(csFilenames{iI},sPath,handles);
    handles.loadedImg = cat(1, handles.loadedImg, imgInfo);
    if(~lCont) % for DICOM: just necessary to load first file -> remaining ones belong to same 3D volume
        break;
    end
end

% get current content
cList = cellfun(@(x) x.sShow, handles.loadedImg,'UniformOutput',false);
set(handles.listImages,'String',cList);

guidata(hObject, handles)


function cimgInfo = fGetImg(csFilename,sPath,handles)
% load images from files and folders (child) => perform the actual loading

% no mixed file extension loading allowed -> sufficient to check 1st file
[~,sFile,sExt] = fileparts(csFilename);
    
switch sExt
    case '.workspace'
        imgInfo.sType = 0;
        imgInfo.sData = ['workspace',filesep,sFile];
        imgInfo.dVoxelsize = [str2num(get(handles.voxel_x,'String')), str2num(get(handles.voxel_y,'String')), str2num(get(handles.voxel_z,'String'))];
        contents = cellstr(get(handles.pm_corsagtrans,'String'));
        imgInfo.orientation = contents{get(handles.pm_corsagtrans,'value')};
        
        evalin('base',sprintf('if exist(''%s'',''var''); assignin(''caller'',''dImg'',%s); end;',sFile, sFile));
        imgInfo.size = size(dImg);
        if(length(imgInfo.size) > 3)
            dSize4D = imgInfo.size(4);
        else
            dSize4D = 1;
        end
        sShowname = sFile;
        
    case {'.mat', '.MAT'}  
        imgInfo.sType = 1;
        imgInfo.sData = [sPath,csFilename];
        warning('off','MATLAB:load:variableNotFound');
        try
            load([sPath,csFilename],'SGeo');
        catch
        end
        if exist('SGeo','var')
            imgInfo.dVoxelsize = SGeo.dVoxelsize;
            imgInfo.orientation = SGeo.sOrientation;
        else
            imgInfo.dVoxelsize = [str2num(get(handles.voxel_x,'String')), str2num(get(handles.voxel_y,'String')), str2num(get(handles.voxel_z,'String'))];
            contents = cellstr(get(handles.pm_corsagtrans,'String'));
            imgInfo.orientation = contents{get(handles.pm_corsagtrans,'value')};
        end
        
        tmp = whos('-file', [sPath, filesep, csFilename]);
%             idx = arrayfun(@(x) strcmp(x.name,'dImg'), tmp); % image must be named: dImg if a 4D array shall be loaded in      
        for iJ = 1:length(tmp) % find 1st double variable
            if(strcmp(tmp(iJ).class, 'double') || strcmp(tmp(iJ).class, 'single'))
                imgInfo.idxMat = iJ;
                break;
            end
        end
                       
        if(length(tmp(imgInfo.idxMat).size) > 3) % 4D array found
            dSize4D = tmp(imgInfo.idxMat).size(4);
        else
            dSize4D = 1;
        end
        imgInfo.size = tmp(imgInfo.idxMat).size;
        sShowname = csFilename;
        % load image
        hWait = fwaitbar(0.5, 'Loading MAT file');
        sLoad = load(imgInfo.sData);
        cNames = {tmp(:).name};
        eval(sprintf('dImg = double(sLoad.%s);', cNames{iJ}));
        close(hWait);  
        
    case {'.dcm', '.DCM', '.ima', '.IMA'} % DICOM
        imgInfo.sType = 2;
        imgInfo.sData = sPath;
        sTag  = dicominfo([sPath,filesep,sFile,sExt]);
        imgInfo.size = [sTag.Rows, sTag.Columns, 0];
        imgInfo.dVoxelsize = [sTag.PixelSpacing; sTag.SliceThickness];
        dOrient = reshape(sTag.ImageOrientationPatient, [3, 2])';        
        dOrientIndicator = sum(abs(dOrient));
        [~, iInd] = min(dOrientIndicator);
        switch iInd
            case 1
                imgInfo.orientation = 'sagittal';
            case 2
                imgInfo.orientation = 'coronal';
            case 3
                imgInfo.orientation = 'transverse';
        end
        if(strcmp(sPath(end),filesep))
            [~,sShowname] = fileparts(sPath(1:end-1));
        else
            [~,sShowname] = fileparts(sPath);
        end
        dSize4D = 1;
        % load complete image folder directly
        try 
            [dImg, sInfo] = fReadDICOM(imgInfo.sData);
        catch
            errordl('No mixed DICOM image loading supported! Please place images in a separate folders!');
        end
        imgInfo.size(3) = size(dImg,3);
        % check if regrouping is necessary
        dSlicePos = zeros(length(sInfo),1);
        for iS = 1:length(sInfo)
            dSlicePos(iS) = sInfo(iS).STag.SliceLocation;
        end
        if(length(unique(dSlicePos)) ~= length(dSlicePos))
            dUnSlicePos = unique(dSlicePos);
            dImgNew = zeros(size(dImg,1),size(dImg,2),length(dUnSlicePos),size(dImg,3)/length(dUnSlicePos));
            for iS = 1:length(dUnSlicePos)
                dImgNew(:,:,iS,:) = permute(dImg(:,:,dSlicePos == dUnSlicePos(iS)),[1 2 4 3]);
            end
            dImg = dImgNew; clear 'dImgNew';
            dSize4D = size(dImg,4); 
            imgInfo.size(3) = size(dImg,3);
            imgInfo.size(4) = dSize4D;
        end
         
    case {'.mhd', '.MHD'} % MHD/RAW
        imgInfo.sType = 3;
        imgInfo.sData = [sPath,csFilename];
        currPath = pwd;
        cd(sPath);
        [dImg, info]=read_mhd(csFilename);
        cd(currPath);
        imgInfo.dVoxelsize = info.spacing;
        imgInfo.size = info.size;
        [yaw, pitch, roll] = dcm2angle(info.orientation);
        if(yaw == 0 && pitch == 0 && roll == 0)
            imgInfo.orientation = 'transverse';
        elseif(yaw == 0 && abs(pitch) == pi/2 && roll == 0)
            imgInfo.orientation = 'sagittal';
        elseif(yaw == 0 && pitch == 0 && abs(roll) == pi/2)
            imgInfo.orientation = 'coronal';
        else
            contents = cellstr(get(handles.pm_corsagtrans,'String'));
            imgInfo.orientation = contents{get(handles.pm_corsagtrans,'value')};
        end
        dSize4D = 1;
        
    case {'.gipl', '.GIPL'}
        imgInfo.sType = 4;
        imgInfo.sData = [sPath,csFilename];
        [imgInfo.dImg, imgInfo.dVoxelsize, imgInfo.size] = giplread(imgInfo.sData);
        contents = cellstr(get(handles.pm_corsagtrans,'String'));
        imgInfo.orientation = contents{get(handles.pm_corsagtrans,'value')};
        dSize4D = 1;
        
    case {'.nii', '.NII'}
        imgInfo.sType = 5;
        imgInfo.sData = sPath;
        info = nii_read_header(sPath,csFilename);
        imgInfo.dImg = nii_read_volume(info);
        
        % TODO: derive information from header
        imgInfo.dVoxelsize = [1 1 1];
        contents = cellstr(get(handles.pm_corsagtrans,'String'));
        imgInfo.orientation = contents{get(handles.pm_corsagtrans,'value')};
        dSize4D = 1;
        
    case {'.jpg', '.jpeg', '.tif', '.tiff', '.gif', '.GIF', '.bmp', '.BMP', '.png', '.PNG'}
        imgInfo.sType = 6;
        imgInfo.sData = [sPath,csFilename];
        imgInfo.dImg = imread([sPath,csFilename]);
        imgInfo.dVoxelsize = [1 1 1];
        contents = cellstr(get(handles.pm_corsagtrans,'String'));
        imgInfo.orientation = contents{get(handles.pm_corsagtrans,'value')};
        dSize4D = 1;
        
    otherwise
        errordlg('Sorry file extension not supported yet :-(');
end

cimgInfo = cell(dSize4D,1);
for iI = 1:dSize4D % split up 4D data and always just store (max.) 3D data
    imgInfo.i4D = iI;
    imgInfo.dImg = dImg(:,:,:,iI);
    imgInfo.sShowname = sShowname;
    if(length(sShowname) < handles.iMinShow )
        imgInfo.sShow = sprintf('%s - %d x %d x %d x %d (%.2f x %.2f x %.2f) - %s', sShowname(1:min([handles.iMinShow ,length(sShowname)])), imgInfo.size(1), imgInfo.size(2), imgInfo.size(3), iI, imgInfo.dVoxelsize(1), imgInfo.dVoxelsize(2), imgInfo.dVoxelsize(3), imgInfo.orientation);
    else
        imgInfo.sShow = sprintf('%s... - %d x %d x %d x %d (%.2f x %.2f x %.2f) - %s', sShowname(1:min([handles.iMinShow ,length(sShowname)])), imgInfo.size(1), imgInfo.size(2), imgInfo.size(3), iI, imgInfo.dVoxelsize(1), imgInfo.dVoxelsize(2), imgInfo.dVoxelsize(3), imgInfo.orientation);
    end
    cimgInfo{iI} = imgInfo;
end
set(handles.voxel_x,'String',imgInfo.dVoxelsize(1));
set(handles.voxel_y,'String',imgInfo.dVoxelsize(2));
set(handles.voxel_z,'String',imgInfo.dVoxelsize(3));
sCST = {'sagittal' 'coronal' 'transverse'};
idx = find(strcmp(sCST, imgInfo.orientation));
set(handles.pm_corsagtrans,'Value',idx);


function button4D_Callback(hObject, eventdata, handles)
% load images from file
fLoadImg(hObject,handles,'normal');


function button4D_ButtonDownFcn(hObject, eventdata, handles)
% load images from folder
fLoadImg(hObject,handles,'rightclick');


function buttonWorkspace_Callback(hObject, eventdata, handles)
% load images from workspace
cVars = evalin('base','whos');
lInd = {cVars(:).class};
lSize = {cVars(:).size};
lInd = (strcmp(lInd,'double') | strcmp(lInd,'single') | strcmp(lInd,'int8') | strcmp(lInd,'int16') | strcmp(lInd,'int32') | strcmp(lInd,'int64') | ...
    strcmp(lInd,'uint8') | strcmp(lInd,'uint16') | strcmp(lInd,'uint32') | strcmp(lInd,'uint64')) & cellfun(@(x) length(x) > 1, lSize);
if(nnz(lInd) == 0)
    return;
end
cVars = {cVars(lInd).name};
iIndLoad = fSmallPopup(cVars,'Load image(s) from workspace',length(cVars));
if(isempty(iIndLoad))
    return;
end
for iI = 1:length(iIndLoad)    
    imgInfo =  fGetImg([cVars{iIndLoad(iI)},'.workspace'],'',handles);
    handles.loadedImg = cat(1, handles.loadedImg, imgInfo);
end

% get current content
cList = cellfun(@(x) x.sShow, handles.loadedImg,'UniformOutput',false);
set(handles.listImages,'String',cList);
guidata(hObject, handles);


function loadParamFile_Callback(hObject, eventdata, handles)
% load parameter file
[Filename, Pathname] = uigetfile(('*.txt'), 'Select parameter file', [handles.currpath,filesep,'registration',filesep,'parameterFiles']);
if Filename==0, return, end
handles.sParFile=[Pathname,Filename];
set(handles.txtParafile,'String',Filename);
guidata(hObject,handles)


function standardParameterFile_Callback(hObject, eventdata, handles)
% set default parameter files

state = get(hObject,'Value'); % is standard parameter box checked?
if state == 1
    if handles.nRegMethod == 1 %elastix
        handles.sParFile = [handles.SPaths.sCode,filesep,'registration',filesep,'parameterFiles',filesep,'elastix_DEFAULT.txt'];
        set(handles.txtParafile,'String','elastix_DEFAULT.txt')
        
    elseif handles.nRegMethod == 2 %halar
        handles.sParFile = [handles.SPaths.sCode,filesep,'registration',filesep,'parameterFiles',filesep,'lreg_DEFAULT.txt'];
        set(handles.txtParafile,'String','lreg_DEFAULT.txt') 
        
    elseif handles.nRegMethod == 4 %LAP
        handles.sParFile = [handles.SPaths.sCode,filesep,'registration',filesep,'parameterFiles',filesep,'LAP_DEFAULT.txt'];
        set(handles.txtParafile,'String','LAP_DEFAULT.txt')
        
    elseif handles.nRegMethod == 5 %demons
        handles.sParFile = [handles.SPaths.sCode,filesep,'registration',filesep,'parameterFiles',filesep,'demons_DEFAULT.txt'];
        set(handles.txtParafile,'String','demons_DEFAULT.txt')
    else
        set(handles.txtParafile,'String',' ')
        errordlg('No Parameter File available or necessary for chosen registration method')
        set(hObject,'Value',0)
    end
else
    handles.sParFile = 'NA';
    set(handles.txtParafile,'String',' ')
end
guidata(hObject,handles)


function selectRegiMethod_Callback(hObject, eventdata, handles)
% select registration method
% define different Registration Methods and set flag
contents = cellstr(get(hObject,'String'));
sRegMethod=contents{get(hObject,'value')};

nRegMethod = handles.dRegMapping{2,strcmp(handles.dRegMapping(1,:),sRegMethod)};
% switch sRegMethod
%     case 'Select Registration Method'
%         sRegMethod = 'NA';
%         nRegMethod = 0;
%     case 'elastix'
%         nRegMethod=1;
%     case 'halar'
%         nRegMethod=2;
%     case 'GRICS'
%         nRegMethod=3;
%     case 'LAP'
%         nRegMethod=4;
%     case 'demons'
%         nRegMethod=5;
%     otherwise
%         error('unexpected error');
% end
set(handles.standardParameterFile,'Value',0)
set(handles.txtParafile,'String',[])
% Save the handles structure.
handles.nRegMethod=nRegMethod;
handles.sRegMethod=sRegMethod;
guidata(hObject,handles)


function pushbutton22_Callback(hObject, eventdata, handles)
% set result directory
sPath = uigetdir(handles.SPaths.sResults);
if isnumeric(sPath), return, end;
handles.SPaths.sResults = sPath;
set(handles.txtOutput,'String',sPath);
guidata(hObject,handles);


function runRegistration_Callback(hObject, eventdata, handles)
% perform registration

% get selected registrations
iInd = get(handles.listRegs,'Value');
if(isempty(iInd) || isempty(handles.allReg))
    return;
end

for iI=1:length(iInd)
    if(handles.allReg{iInd(iI),2}.lDone)
        continue;
    end
    cSavename = unique(cellfun(@(x) x.sShowname, handles.allReg{iInd(iI),1},'UniformOutput',false));
    if(length(cSavename) > 1)
        sSavename = '';
        for iSave=1:length(cSavename)
            sSavename = [sSavename, cSavename{iSave}, '_'];
        end
        sSavename = sSavename(1:end-1);
    else
        sSavename = cSavename{1};
    end
    % check if already a result file is existing => yes: print a warning
    % and ask how to proceed
    if(exist([handles.SPaths.sResults,filesep,handles.allReg{iInd(iI),2}.sRegMethod,filesep,sSavename,'_registration.mat'],'file'))
        sAnswer = questdlg('Registration result does already exist and will be overwritten!', 'File already existing', 'Continue', 'Abort', 'Continue');
        if(strcmp(sAnswer,'Abort'))
            continue;
        end
    end
    if(~exist([handles.SPaths.sResults,filesep,handles.allReg{iInd(iI),2}.sRegMethod],'dir'))
        mkdir([handles.SPaths.sResults,filesep,handles.allReg{iInd(iI),2}.sRegMethod]);
    end
    
    % check if chosen registration is ready to run (all parameters set)
    if(strcmp(handles.allReg{iInd(iI),2}.sRegMethod, 'NA') || strcmp(handles.allReg{iInd(iI),2}.sParafile,'NA'))
        errordlg(sprintf('Not all necessary parameters are set for registration %02d',iInd(iI)), 'Parameters missing');
        continue;
    end
    
    % prepare images
    iShift = 2-length(handles.allReg{iInd(iI),1}{1}.size);
    dFix = handles.allReg{iInd(iI),1}{1}.dImg;
    dMove = cell2mat(shiftdim(cellfun(@(x) x.dImg, handles.allReg{iInd(iI),1}(2:end),'UniformOutput',false),iShift));
    SGeo.cVoxelsize = cellfun(@(x) x.dVoxelsize, handles.allReg{iInd(iI),1}, 'UniformOutput', false);
    SGeo.sOrientation = cellfun(@(x) x.orientation, handles.allReg{iInd(iI),1}, 'UniformOutput', false);
    
    sParafile = handles.allReg{iInd(iI),2}.sParafile;
    iDim = handles.allReg{iInd(iI),2}.iDim; % 1=2D, 2=3D
    
    % prepare output struct
%     sDim = {'2D','3D'};    
%     SRegiResult.sFilename = [handles.allReg{iInd(iI),1}{1}.sShowname];
%     SRegiResult.sParFile = handles.allReg{iInd(iI),2}.sParafile;
%     SRegiResult.SGeo.dVoxelsize = SGeo.cVoxelsize;
%     SRegiResult.sRegMethod = [handles.allReg{iInd(iI),2}.sRegMethod,' (',sDim(iDim),')'];
%     SRegiResult.sShownames = strcat(cellfun(@(x) x.sShowname, handles.allReg{iInd(iI),1}(:), 'UniformOutput', false), cellfun(@(x) sprintf(' (%02d)',x.i4D), handles.allReg{iInd(iI),1}(:), 'UniformOutput', false));
    
    % perform actual registration of the images
    switch handles.allReg{iInd(iI),2}.sRegMethod
        case 'elastix'
            ela = tic;
            [SDeform, dImgReg] = fRegElastix(dFix, dMove, sParafile, iDim, SGeo, [handles.SPaths.sResults,filesep,'elastix'], handles.SPaths.sCode);
            dComptime = toc(ela); 
            SRegiResult.nRegMethod = 1; % elastix  
    
        case 'halar'
            hal = tic;
            [SDeform, dImgReg] = fRegHalar(dFix, dMove, sParafile, iDim, SGeo, [handles.SPaths.sResults,filesep,'halar'], handles.SPaths.sCode);  %handles.sData, handles.sParFile);
            dComptime = toc(hal);
            SRegiResult.nRegMethod = 2; % halar          
    
        case 'GRICS'
            close(h);
            errordlg('GRICS is not implemented yet!')
            return;
    
        case 'LAP'
            lap3d = tic;
            [dFix, dMove, SDeform, dImgReg, cVoxelInterp] = fRegLAP(dFix, dMove, sParafile, iDim, SGeo);
            dComptime = toc(lap3d);           
            
            SRegiResult.nRegMethod = 4; % LAP
            SRegiResult.SGeo.cVoxelsize = cVoxelInterp;
    
        case 'demons'
            dem = tic;
            [dFix, dMove, SDeform, dImgReg, cVoxelInterp] = fRegDemons(dFix, dMove, sParafile, iDim, SGeo);
            dComptime = toc(dem);

            SRegiResult.nRegMethod = 5; % Demons
            SRegiResult.SGeo.cVoxelsize = cVoxelInterp;
            
        otherwise
            errordlg('Choose a registration method first!')
            return
    end
    
    disp(['Time taken for ',  handles.allReg{iInd(iI),2}.sRegMethod,' registration: ', num2str(dComptime)]);
    SRegiResult = fGenerateRegStruct(SRegiResult,dImgReg, SDeform, handles.allReg(iInd(iI),:), dComptime, dFix, dMove, handles.dRegMapping);
    save([handles.SPaths.sResults,filesep,handles.allReg{iInd(iI),2}.sRegMethod,filesep,sSavename,'_registration.mat'],...
        'SRegiResult', '-v7.3');
            
    % save performed registration
    handles.allReg{iInd(iI),2}.lDone = true;
    handles.allReg{iInd(iI),3} = SRegiResult;
end

% change appearance in listbox
bColors = {'FFFFFF','0d9f1e'};
pre1 = '<HTML><BODY bgcolor="';
pre2 = '"><FONT color="';
post = '</FONT></BODY></HTML>'; 
listboxStr = cell(size(handles.allReg,1),1);
for i = 1:size(handles.allReg,1)
    str = [pre1 bColors{handles.allReg{i,2}.lDone+1} pre2 '000000">' handles.allReg{i,2}.sShow post];
    listboxStr{i} = str;
end
set( handles.listRegs, 'String', listboxStr );

% automatic evaluation if flag is set
if handles.FAutoEval == 1
    evalRegiResults_Callback(hObject, [], handles);
end
guidata(hObject,handles)


function SRegiResult = fGenerateRegStruct(SRegiResult,dImgReg, SDeform, cReg, dComptime, dFix, dMove, dRegMapping)
% generate a registration result struct
if(nargin < 6)
    iShift = 1-length(cReg{1,1}{1}.size);
    dFix = cReg{1,1}{1}.dImg;
    dMove = cell2mat(shiftdim(cellfun(@(x) x.dImg, cReg{1,1}(2:end),'UniformOutput',false),iShift));
    dRegMapping = {'elastix', 'halar', 'GRICS', 'LAP', 'demons'; ...
                        1       , 2      , 3      , 4    , 5      };
end
if(nargin < 4), dComptime = 0; end;

SGeo.cVoxelsize = cellfun(@(x) x.dVoxelsize, cReg{1,1}, 'UniformOutput', false);
SGeo.sOrientation = cellfun(@(x) x.orientation, cReg{1,1}, 'UniformOutput', false);

iDim = cReg{1,2}.iDim; % 1=2D, 2=3D
sDim = {'2D','3D'};
% prepare output struct
SRegiResult.dComptime = dComptime;
SRegiResult.dImg = cat(4,dFix,dMove);
SRegiResult.dImgReg = dImgReg;
SRegiResult.SDeform = SDeform;
SRegiResult.sFilename = [cReg{1,1}{1}.sShowname];
SRegiResult.sParFile = cReg{1,2}.sParafile;
if(isfield(SRegiResult,'SGeo'))
    if(~isfield(SRegiResult.SGeo,'cVoxelsize'))
        SRegiResult.SGeo.cVoxelsize = SGeo.cVoxelsize;
    end
else
    SRegiResult.SGeo.cVoxelsize = SGeo.cVoxelsize;
end
SRegiResult.SGeo.cOrientation = SGeo.sOrientation;
SRegiResult.iDim = iDim;
SRegiResult.sRegMethod = [cReg{1,2}.sRegMethod,' (',sDim{iDim},') -> ',cReg{1,1}{1}.sShowname(1:min([length(cReg{1,1}{1}.sShowname), 20])),' (',num2str(cReg{1,1}{1}.i4D),')'];
SRegiResult.sShownames = strcat(cellfun(@(x) x.sShowname, cReg{1,1}(:), 'UniformOutput', false).', cellfun(@(x) sprintf(' (%02d)',x.i4D), cReg{1,1}(:), 'UniformOutput', false).');

SRegiResult.nRegMethod = dRegMapping{2,strcmp(dRegMapping(1,:),cReg{1,2}.sRegMethod)};


function loadDemoData_Callback(hObject, eventdata, handles)
% load Demo data
h = fwaitbar(0,'Loading demo data. Please wait...'); st = 0; steps = 1;
% load([handles.sDemoDataPath,filesep,'Example',filesep,'data_la_Ph4_Tol100_t000_Ext00_EspOff_closest_recon_SRegiResult.mat']);
if(exist('K:\CS_Retro\ISBI16\RESULTS\elastix\la\real_optimized.mat','file'))
    sDemoPath = 'K:\CS_Retro\ISBI16\RESULTS\elastix\la\real_optimized.mat';
else
    [sFile,sPath] = uigetfile();
    if(isempty(sFile))
        return;
    end
    sDemoPath = [sPath,filesep,sFile];
end
load(sDemoPath);

st= st+1;fwaitbar(st/steps,h);try close(h); catch; end;

cRegImg = cell(1,size(SRegiResult.dImg,4));
for iI=1:size(SRegiResult.dImg,4)
    cRegImg{iI}.dImg = SRegiResult.dImg(:,:,:,iI);
    
    cRegImg{iI}.i4D = iI;
    cRegImg{iI}.size = [size(SRegiResult.dImg,1),size(SRegiResult.dImg,2),size(SRegiResult.dImg,3)];
%     if(iI > 1)
%         cRegImg{iI}.SDeform = SRegiResult.SDeform(iI);
%     end
    cRegImg{iI}.sShowname = 'la';
    cRegImg{iI}.sShow = sprintf('la 256x256x72x%d (1.95x1.95x5)',iI);
    cRegImg{iI}.sType = 1; % mat file
    cRegImg{iI}.dVoxelsize = [1.95, 1.95, 5];
    cRegImg{iI}.orientation = 'coronal';
    cRegImg{iI}.sData = sDemoPath;
end

handles.sParFile = SRegiResult.sParFile;
handles.sPathname = 'DemoData';
handles.nRegMethod = 1;
handles.sRegMethod = 'elastix';
set(handles.selectRegiMethod,'Value',2);
set(handles.popupDim,'Value',2);

set(handles.voxel_x,'String',1.95);
handles.SGeo.dVoxelsize(1) = 1.95;
set(handles.voxel_y,'String',1.95);
handles.SGeo.dVoxelsize(2) = 1.95;
set(handles.voxel_z,'String',5);
handles.SGeo.dVoxelsize(3) = 5;
set(handles.pm_corsagtrans,'Value',1);

set(handles.txtParafile,'String',handles.sParFile);
%set(handles.text3,'String',handles.sFilename); update listbox!
handles = fInsertRegistration(handles, cRegImg, true);
SRegRes = fGenerateRegStruct([],SRegiResult.dImgReg, SRegiResult.SDeform, handles.allReg(end,:), handles.dRegMapping); % needed for old example data
handles.allReg{end,3} = SRegRes;
% else for new example: handles.allReg{end,3} = SRegiResult;
handles = plotImageAndDefField(handles, size(handles.allReg,1));
guidata(hObject,handles)


function pb_LoadReg_Callback(hObject, eventdata, handles)
% load already processed registration
[sFilename, sPath, sExt] = uigetfile('*.mat', 'Select registration result', handles.SPaths.sResults);
if isnumeric(sPath), return, end;   % Dialog aborted 

h = fwaitbar(0, 'Loading registration'); st = 0;
load([sPath,sFilename]);

if(~exist('SRegiResult','var'))
    try close(h); catch; end;
    errordlg('Not a valid registration result file');
    return;
end
cRegImg = cell(1,size(SRegiResult.dImg,4));
steps = size(SRegiResult.dImg,4) + 1;
st= st+1;fwaitbar(st/steps,h);
try
    for iI=1:size(SRegiResult.dImg,4)
        cRegImg{iI}.dImg = SRegiResult.dImg(:,:,:,iI);

        cRegImg{iI}.i4D = iI;
        cRegImg{iI}.size = [size(SRegiResult.dImg,1),size(SRegiResult.dImg,2),size(SRegiResult.dImg,3)];
        cRegImg{iI}.dVoxelsize = SRegiResult.SGeo.cVoxelsize{iI};
        cRegImg{iI}.orientation = SRegiResult.SGeo.cOrientation{iI};
        cRegImg{iI}.sType = 1; % mat file => no backtracking to actual image files
        cRegImg{iI}.sData = [sPath,sFilename];

        cRegImg{iI}.sShowname = SRegiResult.sShownames{iI};
        if(length(cRegImg{iI}.sShowname) < handles.iMinShow)
            cRegImg{iI}.sShow = sprintf('%s - %d x %d x %d x %d (%.2f x %.2f x %.2f) - %s', cRegImg{iI}.sShowname(1:min([handles.iMinShow,length(cRegImg{iI}.sShowname)])), cRegImg{iI}.size(1), cRegImg{iI}.size(2), cRegImg{iI}.size(3), iI, cRegImg{iI}.dVoxelsize(1), cRegImg{iI}.dVoxelsize(2), cRegImg{iI}.dVoxelsize(3), cRegImg{iI}.orientation);
        else
            cRegImg{iI}.sShow = sprintf('%s... - %d x %d x %d x %d (%.2f x %.2f x %.2f) - %s', cRegImg{iI}.sShowname(1:min([handles.iMinShow,length(cRegImg{iI}.sShowname)])), cRegImg{iI}.size(1), cRegImg{iI}.size(2), cRegImg{iI}.size(3), iI, cRegImg{iI}.dVoxelsize(1), cRegImg{iI}.dVoxelsize(2), cRegImg{iI}.dVoxelsize(3), cRegImg{iI}.orientation);
        end
        st= st+1;fwaitbar(st/steps,h);
    end
    handles.sParFile = SRegiResult.sParFile;
    handles.nRegMethod = SRegiResult.nRegMethod;
    contents = cellstr(get(handles.selectRegiMethod,'String'));
    contents = contents(2:end);
    handles.sRegMethod = contents{handles.nRegMethod};
    set(handles.selectRegiMethod,'Value',handles.nRegMethod+1);
    set(handles.popupDim,'Value',SRegiResult.iDim);
    set(handles.txtParafile,'String',handles.sParFile);
catch
    try close(h); catch; end;
    errordlg('Not a valid registration result file');
    return;
end
handles = fInsertRegistration(handles, cRegImg, true);
handles.allReg{end,3} = SRegiResult;
try close(h); catch; end;
handles = plotImageAndDefField(handles, size(handles.allReg,1));

guidata(hObject,handles)


function clearFigures_Callback(hObject, eventdata, handles)
% clear all axes
handles = fClearAllAxes(handles);
guidata(hObject, handles)


function evalRegiResults_Callback(hObject, eventdata, handles)
% evaluate registration result in EvalGUI
iInd = get(handles.listRegs,'Value');
if(isempty(handles.allReg)); return; end;
if(~isempty(iInd))
    for iI=1:length(iInd)
        if(~handles.allReg{iInd(iI),2}.lDone)
            continue;
        end
        
        cSavename = unique(cellfun(@(x) x.sShowname, handles.allReg{iInd(iI),1},'UniformOutput',false));
        if(length(cSavename) > 1)
            sSavename = '';
            for iSave=1:length(cSavename)
                sSavename = [sSavename, cSavename{iSave}, '_'];
            end
            sSavename = sSavename(1:end-1);
        else
            sSavename = cSavename{1};
        end
        if(length(handles.allReg(iInd(iI),:)) < 3 || isempty(handles.allReg{iInd(iI),3})) % not stored in struct -> load from file
            if(~exist([handles.SPaths.sResults,filesep,handles.allReg{iInd(iI),2}.sRegMethod,filesep,sSavename,'_registration.mat'],'file'))
                errordlg('Registration result file not existing');
                return;
            end
            load([handles.SPaths.sResults,filesep,handles.allReg{iInd(iI),2}.sRegMethod,filesep,sSavename,'_registration.mat']);
            if(~exist('SRegiResult','var'))
                errordl('Not a valid registration result file');
                return;
            end
            h = SRegiResult;
        else
            h = handles.allReg{iInd(iI),3};
        end
       
        if(handles.allReg{iInd(iI),1}{1}.sType == 1) % mat
            h.sPathname = fileparts(handles.allReg{iInd(iI),1}{1}.sData);
        elseif(handles.allReg{iInd(iI),1}{1}.sType == 2) % DICOM
            h.sPathname = handles.allReg{iInd(iI),1}{1}.sData;
        end        
        h.SPaths = handles.SPaths;
        cSCT = cellstr(get(handles.pm_corsagtrans,'String'));
        iMap = [2 1 3];
        h.nSCT = iMap(strcmp(cSCT,handles.allReg{iInd(iI),1}{1}.orientation));

        if ~handles.FAutoEval
            EvalGUI('inarg', h);
        else
            EvalGUI('inarg', h, 'EGmode', 'auto')
        end  
    end
end


function cb_autoEval_Callback(hObject, eventdata, handles)
% check for automatic evaluation
if get(hObject, 'Value')
    handles.FAutoEval = 1;
else
    handles.FAutoEval = 0;
end
guidata(hObject, handles)


function pb_landmarkGUI_Callback(hObject, eventdata, handles)
% evaluate registration result in LandmarkGUI
iInd = get(handles.listRegs,'Value');
if(isempty(handles.allReg)); return; end;
if(~isempty(iInd))
    for iI=1:length(iInd)
        if(~handles.allReg{iInd(iI),2}.lDone)
            continue;
        end
        
        cSavename = unique(cellfun(@(x) x.sShowname, handles.allReg{iInd(iI),1},'UniformOutput',false));
        if(length(cSavename) > 1)
            sSavename = '';
            for iSave=1:length(cSavename)
                sSavename = [sSavename, cSavename{iSave}, '_'];
            end
            sSavename = sSavename(1:end-1);
        else
            sSavename = cSavename{1};
        end
        if(length(handles.allReg(iInd(iI),:)) < 3 || isempty(handles.allReg{iInd(iI),3})) % not stored in struct -> load from file
            if(~exist([handles.SPaths.sResults,filesep,handles.allReg{iInd(iI),2}.sRegMethod,filesep,sSavename,'_registration.mat'],'file'))
                errordlg('Registration result file not existing');
                return;
            end
            load([handles.SPaths.sResults,filesep,handles.allReg{iInd(iI),2}.sRegMethod,filesep,sSavename,'_registration.mat']);
            if(~exist('SRegiResult','var'))
                errordl('Not a valid registration result file');
                return;
            end
            h = SRegiResult;
        else
            h = handles.allReg{iInd(iI),3};
        end
        
        if(handles.allReg{iInd(iI),1}{1}.sType == 1) % mat
            h.sPathname = fileparts(handles.allReg{iInd(iI),1}{1}.sData);
        elseif(handles.allReg{iInd(iI),1}{1}.sType == 2) % DICOM
            h.sPathname = handles.allReg{iInd(iI),1}{1}.sData;
        end        
        h.SPaths = handles.SPaths;
        cSCT = cellstr(get(handles.pm_corsagtrans,'String'));
        iMap = [2 1 3];
        h.nSCT = iMap(strcmp(cSCT,handles.allReg{iInd(iI),1}{1}.orientation));

        LandmarkGUI('inarg', h);        
    end
end


function listRegs_Callback(hObject, eventdata, handles)
% show or change registration information
iInd = get(hObject, 'Value');
if(~isempty(handles.allReg))
    regInfo = handles.allReg{iInd,2};
    sRegiMethod = get(handles.selectRegiMethod,'String');
    idx = find(strcmp(sRegiMethod, regInfo.sRegMethod));
    if(isempty(idx))
        set(handles.selectRegiMethod,'Value',1);
    else
        set(handles.selectRegiMethod,'Value',idx);
    end
    
    set(handles.popupDim,'Value',regInfo.iDim);
    [~,sTmpParafile] = fileparts(regInfo.sParafile);
    set(handles.txtParafile,'String',sTmpParafile);
    if(handles.allReg{iInd,2}.lDone) % show images and df
        handles = plotImageAndDefField(handles, iInd);
    end
end % show information nevertheless
if strcmp(get(handles.RegGUI,'SelectionType'),'open') % double click
    handles.doubleClick(2) = true;
    set(handles.buttonOK,'Visible','on');
    set(handles.buttonCancel,'Visible','on');
    bColors = {'FFFFFF','0d9f1e'};
    pre1 = '<HTML><BODY bgcolor="';
    pre2 = '"><FONT color="';
    post = '</FONT></BODY></HTML>'; 
    listboxStr = cell(size(handles.allReg,1),1);
    for i = 1:size(handles.allReg,1)
        if(i == iInd)
            str = [pre1 bColors{handles.allReg{i,2}.lDone+1} pre2 '000000">' handles.allReg{i,2}.sShow post];
        else
            str = [pre1 bColors{handles.allReg{i,2}.lDone+1} pre2 'c3c3c3">' handles.allReg{i,2}.sShow post];
        end
        listboxStr{i} = str;
    end
    set( handles.listRegs, 'String', listboxStr );
else
    set(handles.buttonOK,'Visible','off');
    set(handles.buttonCancel,'Visible','off');
end
guidata(hObject,handles);


function listImages_Callback(hObject, eventdata, handles)
% show or change loaded images
iInd = get(hObject, 'Value');
if(~isempty(handles.loadedImg))
    imgInfo = handles.loadedImg{iInd};
    set(handles.voxel_x,'String',imgInfo.dVoxelsize(1));
    set(handles.voxel_y,'String',imgInfo.dVoxelsize(2));
    set(handles.voxel_z,'String',imgInfo.dVoxelsize(3));
    sCST = {'sagittal' 'coronal' 'transverse'};
    idx = find(strcmp(sCST, imgInfo.orientation));
    set(handles.pm_corsagtrans,'Value',idx);
end % show information nevertheless
if strcmp(get(handles.RegGUI,'SelectionType'),'open') % double click
    handles.doubleClick(1) = true;
    set(handles.buttonOK,'Visible','on');
    set(handles.buttonCancel,'Visible','on');
    pre = '<HTML><FONT color="';
    post = '</FONT></HTML>'; 
    listboxStr = cell(length(handles.loadedImg),1);
    for i = 1:length(handles.loadedImg)
        if(i == iInd)
            str = [pre '000000">' handles.loadedImg{i}.sShow post];
        else
            str = [pre 'c3c3c3">' handles.loadedImg{i}.sShow post];
        end
        listboxStr{i} = str;
    end
    set( handles.listImages, 'String', listboxStr );
else % single click
    if(~isempty(handles.loadedImg)) 
        for iJ=1:length(iInd) % show images on axis
            if(all(~handles.lOpen))
                break;
            end
            if(~any(ismember(handles.iIndGlobal,iInd(iJ))))
                iOpen = find(handles.lOpen,1,'first');
                handles.lOpen(iOpen) = false;
                dImg = scaleImg(handles.loadedImg{iInd(iJ)}.dImg,[0 1]);
                handles.slice(:,iOpen) = [round(size(dImg,3)/2); size(dImg,3)];
                eval(sprintf('set(handles.hI%d, ''CData'', dImg(:,:,round(size(dImg,3)/2)));',iOpen));
                eval(sprintf('set(handles.sliceNo%d, ''String'', ''%02d/%02d'');', iOpen, handles.slice(1,iOpen), handles.slice(2,iOpen)));
                eval(sprintf('set(handles.showName%d, ''String'', ''%s - %02d'');', iOpen, handles.loadedImg{iInd(iJ)}.sShowname, handles.loadedImg{iInd(iJ)}.i4D));
                handles.iIndGlobal(iOpen) = iInd(iJ);
            end
        end
    end
    set(handles.buttonOK,'Visible','off');
    set(handles.buttonCancel,'Visible','off');
end
guidata(hObject,handles);


function listRegs_KeyPressFcn(hObject, eventdata, handles)
% delete registration from list
if(strcmp(eventdata.Key,'delete'))
    iInds = get(handles.listRegs,'Value');
    if(~isempty(handles.allReg))
        lInd = false(size(handles.allReg,1),1);
        lInd(iInds) = true;
        handles.allReg = handles.allReg(~lInd);
        if(any(ismember(handles.iIndGlobal, -iInds)))
            handles = fClearAllAxes(handles);
        end
        iNextVal = min([max(iInds(:)) + 1, max([length(handles.allReg), 1])]);
        if(~isempty(handles.allReg)) % single element deleted
            cList = cellfun(@(x) x.sShow, handles.allReg(:,2),'UniformOutput',false);
        else
            cList = 'insert registrations';
        end
        set(handles.listRegs,'String',cList, 'Value', iNextVal);
    end
end
guidata(hObject,handles);


function listImages_KeyPressFcn(hObject, eventdata, handles)
% delete images from list
if(strcmp(eventdata.Key,'delete'))
    iInds = get(handles.listImages,'Value');
    if(~isempty(handles.loadedImg))
        lInd = false(length(handles.loadedImg),1);
        lInd(iInds) = true;
        handles.loadedImg = handles.loadedImg(~lInd);
        if(any(ismember(handles.iIndGlobal, iInds)))
            handles = fClearAllAxes(handles);
        end
        iNextVal = min([max(iInds(:)) + 1, max([length(handles.loadedImg), 1])]);
        if(~isempty(handles.loadedImg))
            cList = cellfun(@(x) x.sShow, handles.loadedImg,'UniformOutput',false);
        else
            cList = {'insert images'};
        end
        set(handles.listImages,'String',cList, 'Value', iNextVal);
    end
end
guidata(hObject,handles);


function buttonAddReg_Callback(hObject, ~, handles)
% add registration to list

% get selected images
iInds = get(handles.listImages,'Value');
if(~isempty(handles.loadedImg))
    if(length(iInds) < 2)
        errordlg('Please select at least two images!', 'Insufficient amount of images for registration');
        return;
    end
    lInd = false(length(handles.loadedImg),1);
    lInd(iInds) = true;
    % open a popup window to select reference image
    dSelImg = handles.loadedImg(lInd);
    cShow = cell(length(dSelImg),1);
    for iI=1:length(dSelImg)
        cShow{iI} = dSelImg{iI}.sShow;
    end
    iIndRefImg = fSmallPopup(cShow,'Choose reference (fixed) image',1);
    if(isempty(iIndRefImg))
        return
    end
    lInd = true(length(dSelImg),1);
    lInd(iIndRefImg) = false;
    cRegImg = {dSelImg{iIndRefImg},dSelImg{lInd}};
    handles = fInsertRegistration(handles, cRegImg, false);
end

% write into listRegs
guidata(hObject,handles);

function handles = fInsertRegistration(handles, cRegImg, lDone)
% insert a new registration with all necessary information

% check for same size of all images
if(size(unique(cell2mat(cellfun(@(x) x.size, cRegImg, 'UniformOutput', false).'),'rows'),1) > 1)
    errordlg('Images have different size. Please select other images or interpolate to same size first!', 'Different image size');
    return;
end
handles.allReg{end+1,1} = cRegImg;
sReg.sRegMethod = handles.sRegMethod;
%     sReg.sParafile = get(handles.txtParafile,'String');
sReg.sParafile = handles.sParFile;
[~,sTmpParafile] = fileparts(sReg.sParafile);
if(length(handles.allReg{end,1}{1}.size) == 2) % just 2D image
    sReg.iDim = 1;
else
    sReg.iDim = get(handles.popupDim,'value'); % 1=2D, 2=3D
end
sTmp = {'2D','3D'};
sShownames = cellfun(@(x) x.sShowname, handles.allReg{end,1},'UniformOutput',false);
if(length(unique(sShownames)) < length(sShownames)) % same datasets are included
    if(length(handles.allReg{end,1}) > 4)
        sReg.sShow = sprintf('reg %02d: %s %s (%s) %s (%02d) <-- %s (%02d), %s (%02d), %s (%02d), ...',size(handles.allReg,1), sReg.sRegMethod, sTmp{sReg.iDim}, sTmpParafile, handles.allReg{end,1}{1}.sShowname(1:min([length(handles.allReg{end,1}{1}.sShowname) handles.iMinShow ])), handles.allReg{end,1}{1}.i4D, handles.allReg{end,1}{2}.sShowname(1:min([length(handles.allReg{end,1}{2}.sShowname) handles.iMinShow ])), handles.allReg{end,1}{2}.i4D, handles.allReg{end,1}{3}.sShowname(1:min([length(handles.allReg{end,1}{3}.sShowname) handles.iMinShow ])), handles.allReg{end,1}{3}.i4D, handles.allReg{end,1}{4}.sShowname(1:min([length(handles.allReg{end,1}{4}.sShowname) handles.iMinShow ])), handles.allReg{end,1}{4}.i4D);
    elseif(length(handles.allReg{end,1}) == 4)
        sReg.sShow = sprintf('reg %02d: %s %s (%s) %s (%02d) <-- %s (%02d), %s (%02d), %s (%02d)',size(handles.allReg,1), sReg.sRegMethod, sTmp{sReg.iDim}, sTmpParafile, handles.allReg{end,1}{1}.sShowname(1:min([length(handles.allReg{end,1}{1}.sShowname) handles.iMinShow ])), handles.allReg{end,1}{1}.i4D, handles.allReg{end,1}{2}.sShowname(1:min([length(handles.allReg{end,1}{2}.sShowname) handles.iMinShow ])), handles.allReg{end,1}{2}.i4D, handles.allReg{end,1}{3}.sShowname(1:min([length(handles.allReg{end,1}{3}.sShowname) handles.iMinShow ])), handles.allReg{end,1}{3}.i4D, handles.allReg{end,1}{4}.sShowname(1:min([length(handles.allReg{end,1}{4}.sShowname) handles.iMinShow ])), handles.allReg{end,1}{4}.i4D);
    elseif(length(handles.allReg{end,1}) == 3)
        sReg.sShow = sprintf('reg %02d: %s %s (%s) %s (%02d) <-- %s (%02d), %s (%02d)',size(handles.allReg,1), sReg.sRegMethod, sTmp{sReg.iDim}, sTmpParafile, handles.allReg{end,1}{1}.sShowname(1:min([length(handles.allReg{end,1}{1}.sShowname) handles.iMinShow ])), handles.allReg{end,1}{1}.i4D, handles.allReg{end,1}{2}.sShowname(1:min([length(handles.allReg{end,1}{2}.sShowname) handles.iMinShow ])), handles.allReg{end,1}{2}.i4D, handles.allReg{end,1}{3}.sShowname(1:min([length(handles.allReg{end,1}{3}.sShowname) handles.iMinShow ])), handles.allReg{end,1}{3}.i4D);
    else
        sReg.sShow = sprintf('reg %02d: %s %s (%s) %s (%02d) <-- %s (%02d)',size(handles.allReg,1), sReg.sRegMethod, sTmp{sReg.iDim}, sTmpParafile, handles.allReg{end,1}{1}.sShowname(1:min([length(handles.allReg{end,1}{1}.sShowname) handles.iMinShow ])), handles.allReg{end,1}{1}.i4D, handles.allReg{end,1}{2}.sShowname(1:min([length(handles.allReg{end,1}{2}.sShowname) handles.iMinShow ])), handles.allReg{end,1}{2}.i4D);
    end
else       
    if(length(handles.allReg{end,1}) > 4)
        sReg.sShow = sprintf('reg %02d: %s %s (%s) %s <-- %s, %s, %s, ...',size(handles.allReg,1), sReg.sRegMethod, sTmp{sReg.iDim}, sTmpParafile, handles.allReg{end,1}{1}.sShowname(1:min([length(handles.allReg{end,1}{1}.sShowname) handles.iMinShow ])), handles.allReg{end,1}{2}.sShowname(1:min([length(handles.allReg{end,1}{2}.sShowname) handles.iMinShow ])), handles.allReg{end,1}{3}.sShowname(1:min([length(handles.allReg{end,1}{3}.sShowname) handles.iMinShow ])), handles.allReg{end,1}{4}.sShowname(1:min([length(handles.allReg{end,1}{4}.sShowname) handles.iMinShow ])));
    elseif(length(handles.allReg{end,1}) == 4)
        sReg.sShow = sprintf('reg %02d: %s %s (%s) %s <-- %s, %s, %s',size(handles.allReg,1), sReg.sRegMethod, sTmp{sReg.iDim}, sTmpParafile, handles.allReg{end,1}{1}.sShowname(1:min([length(handles.allReg{end,1}{1}.sShowname) handles.iMinShow ])), handles.allReg{end,1}{2}.sShowname(1:min([length(handles.allReg{end,1}{2}.sShowname) handles.iMinShow ])), handles.allReg{end,1}{3}.sShowname(1:min([length(handles.allReg{end,1}{3}.sShowname) handles.iMinShow ])), handles.allReg{end,1}{4}.sShowname(1:min([length(handles.allReg{end,1}{4}.sShowname) handles.iMinShow ])));
    elseif(length(handles.allReg{end,1}) == 3)
        sReg.sShow = sprintf('reg %02d: %s %s (%s) %s <-- %s, %s',size(handles.allReg,1), sReg.sRegMethod, sTmp{sReg.iDim}, sTmpParafile, handles.allReg{end,1}{1}.sShowname(1:min([length(handles.allReg{end,1}{1}.sShowname) handles.iMinShow ])), handles.allReg{end,1}{2}.sShowname(1:min([length(handles.allReg{end,1}{2}.sShowname) handles.iMinShow])), handles.allReg{end,1}{3}.sShowname(1:min([length(handles.allReg{end,1}{3}.sShowname) handles.iMinShow])));
    else
        sReg.sShow = sprintf('reg %02d: %s %s (%s) %s <-- %s',size(handles.allReg,1), sReg.sRegMethod, sTmp{sReg.iDim}, sTmpParafile, handles.allReg{end,1}{1}.sShowname(1:min([length(handles.allReg{end,1}{1}.sShowname) handles.iMinShow])), handles.allReg{end,1}{2}.sShowname(1:min([length(handles.allReg{end,1}{2}.sShowname) handles.iMinShow])));
    end
end
sReg.lDone = lDone;
handles.allReg{end,2} = sReg;
%     cList = cellfun(@(x) x.sShow, handles.allReg(:,2),'UniformOutput',false);
%     set(handles.listRegs,'String',cList);   
bColors = {'FFFFFF','0d9f1e'};
pre1 = '<HTML><BODY bgcolor="';
pre2 = '"><FONT color="';
post = '</FONT></BODY></HTML>'; 
listboxStr = cell(size(handles.allReg,1),1);
for i = 1:size(handles.allReg,1)
    str = [pre1 bColors{handles.allReg{i,2}.lDone+1} pre2 '000000">' handles.allReg{i,2}.sShow post];
    listboxStr{i} = str;
end
set( handles.listRegs, 'String', listboxStr );


function buttonOK_Callback(hObject, eventdata, handles)
% apply changes

if(handles.doubleClick(1)) % listImages
    handles.doubleClick(1) = false; % restore double click mode
    % apply changes
    iInd = get(handles.listImages, 'Value');
    handles.loadedImg{iInd}.dVoxelsize(1) = str2double(get(handles.voxel_x,'String'));
    handles.loadedImg{iInd}.dVoxelsize(2) = str2double(get(handles.voxel_y,'String'));
    handles.loadedImg{iInd}.dVoxelsize(3) = str2double(get(handles.voxel_z,'String'));
    contents = cellstr(get(handles.pm_corsagtrans,'String'));
    SagCorTrans = contents{get(handles.pm_corsagtrans,'value')};
    handles.loadedImg{iInd}.orientation = SagCorTrans;
    % update show
    if(length(handles.loadedImg{iInd}.sShowname) < handles.iMinShow)
        handles.loadedImg{iInd}.sShow = sprintf('%s - %d x %d x %d x %d (%.2f x %.2f x %.2f) - %s', handles.loadedImg{iInd}.sShowname(1:min([handles.iMinShow,length(handles.loadedImg{iInd}.sShowname)])), handles.loadedImg{iInd}.size(1), handles.loadedImg{iInd}.size(2), handles.loadedImg{iInd}.size(3), handles.loadedImg{iInd}.i4D, handles.loadedImg{iInd}.dVoxelsize(1), handles.loadedImg{iInd}.dVoxelsize(2), handles.loadedImg{iInd}.dVoxelsize(3), handles.loadedImg{iInd}.orientation);
    else
        handles.loadedImg{iInd}.sShow = sprintf('%s... - %d x %d x %d x %d (%.2f x %.2f x %.2f) - %s', handles.loadedImg{iInd}.sShowname(1:min([handles.iMinShow,length(handles.loadedImg{iInd}.sShowname)])), handles.loadedImg{iInd}.size(1), handles.loadedImg{iInd}.size(2), handles.loadedImg{iInd}.size(3), handles.loadedImg{iInd}.i4D, handles.loadedImg{iInd}.dVoxelsize(1), handles.loadedImg{iInd}.dVoxelsize(2), handles.loadedImg{iInd}.dVoxelsize(3), handles.loadedImg{iInd}.orientation);
    end
    
    % make everything black again
    pre = '<HTML><FONT color="';
    post = '</FONT></HTML>'; 
    listboxStr = cell(length(handles.loadedImg),1);
    for i = 1:length(handles.loadedImg)      
        str = [pre '000000">' handles.loadedImg{i}.sShow post]; 
        listboxStr{i} = str;
    end
    set( handles.listImages, 'String', listboxStr );
end
if(handles.doubleClick(2)) % listReg
    handles.doubleClick(2) = false; % restore double click mode
    % apply changes
    iInd = get(handles.listRegs, 'Value');
    lChange = true;
    if(handles.allReg{iInd,2}.lDone)
        sAnswer = questdlg('Registration already performed. Previous results will be deleted!', 'Registration already performed', 'OK', 'Abort', 'OK');
        if(strcmp(sAnswer,'Abort'))
            lChange = false;
        end
    end

    if(lChange)
        handles.allReg{iInd,2}.sRegMethod = handles.sRegMethod;
        handles.allReg{iInd,2}.sParafile = handles.sParFile;
        if(length(handles.allReg{iInd,1}{1}.size) == 2) % just 2D image
            handles.allReg{iInd,2}.iDim = 1;
        else
            handles.allReg{iInd,2}.iDim = get(handles.popupDim,'value'); % 1=2D, 2=3D
        end
        % update show
        sTmp = {'2D','3D'};
        [~,sTmpParafile] = fileparts(handles.allReg{iInd,2}.sParafile);
        sShownames = cellfun(@(x) x.sShowname, handles.allReg{iInd,1},'UniformOutput',false);
        if(length(unique(sShownames)) < length(sShownames)) % same datasets are included
            if(length(handles.allReg{iInd,1}) > 4)
                handles.allReg{iInd,2}.sShow = sprintf('reg %02d: %s %s (%s) %s (%02d) <-- %s (%02d), %s (%02d), %s (%02d), ...', iInd, handles.allReg{iInd,2}.sRegMethod, sTmp{handles.allReg{iInd,2}.iDim}, sTmpParafile, handles.allReg{iInd,1}{1}.sShowname(1:min([length(handles.allReg{iInd,1}{1}.sShowname) handles.iMinShow])), handles.allReg{iInd,1}{1}.i4D, handles.allReg{iInd,1}{2}.sShowname(1:min([length(handles.allReg{iInd,1}{2}.sShowname) handles.iMinShow])), handles.allReg{iInd,1}{2}.i4D, handles.allReg{iInd,1}{3}.sShowname(1:min([length(handles.allReg{iInd,1}{3}.sShowname) handles.iMinShow])), handles.allReg{iInd,1}{3}.i4D, handles.allReg{iInd,1}{4}.sShowname(1:min([length(handles.allReg{iInd,1}{4}.sShowname) handles.iMinShow])), handles.allReg{iInd,1}{4}.i4D);
            elseif(length(handles.allReg{iInd,1}) == 4)
                handles.allReg{iInd,2}.sShow = sprintf('reg %02d: %s %s (%s) %s (%02d) <-- %s (%02d), %s (%02d), %s (%02d)', iInd, handles.allReg{iInd,2}.sRegMethod, sTmp{handles.allReg{iInd,2}.iDim}, sTmpParafile, handles.allReg{iInd,1}{1}.sShowname(1:min([length(handles.allReg{iInd,1}{1}.sShowname) handles.iMinShow])), handles.allReg{iInd,1}{1}.i4D, handles.allReg{iInd,1}{2}.sShowname(1:min([length(handles.allReg{iInd,1}{2}.sShowname) handles.iMinShow])), handles.allReg{iInd,1}{2}.i4D, handles.allReg{iInd,1}{3}.sShowname(1:min([length(handles.allReg{iInd,1}{3}.sShowname) handles.iMinShow])), handles.allReg{iInd,1}{3}.i4D, handles.allReg{iInd,1}{4}.sShowname(1:min([length(handles.allReg{iInd,1}{4}.sShowname) handles.iMinShow])), handles.allReg{iInd,1}{4}.i4D);
            elseif(length(handles.allReg{iInd,1}) == 3)
                handles.allReg{iInd,2}.sShow = sprintf('reg %02d: %s %s (%s) %s (%02d) <-- %s (%02d), %s (%02d)',iInd, handles.allReg{iInd,2}.sRegMethod, sTmp{handles.allReg{iInd,2}.iDim}, sTmpParafile, handles.allReg{iInd,1}{1}.sShowname(1:min([length(handles.allReg{iInd,1}{1}.sShowname) handles.iMinShow])), handles.allReg{iInd,1}{1}.i4D, handles.allReg{iInd,1}{2}.sShowname(1:min([length(handles.allReg{iInd,1}{2}.sShowname) handles.iMinShow])), handles.allReg{iInd,1}{2}.i4D, handles.allReg{iInd,1}{3}.sShowname(1:min([length(handles.allReg{iInd,1}{3}.sShowname) handles.iMinShow])), handles.allReg{iInd,1}{3}.i4D);
            else
                handles.allReg{iInd,2}.sShow = sprintf('reg %02d: %s %s (%s) %s (%02d) <-- %s (%02d)',iInd, handles.allReg{iInd,2}.sRegMethod, sTmp{handles.allReg{iInd,2}.iDim}, sTmpParafile, handles.allReg{iInd,1}{1}.sShowname(1:min([length(handles.allReg{iInd,1}{1}.sShowname) handles.iMinShow])), handles.allReg{iInd,1}{1}.i4D, handles.allReg{iInd,1}{2}.sShowname(1:min([length(handles.allReg{iInd,1}{2}.sShowname) handles.iMinShow])), handles.allReg{iInd,1}{2}.i4D);
            end
        else       
            if(length(handles.allReg{iInd,1}) > 4)
                handles.allReg{iInd,2}.sShow = sprintf('reg %02d: %s %s (%s) %s <-- %s, %s, %s, ...',iInd, handles.allReg{iInd,2}.sRegMethod, sTmp{handles.allReg{iInd,2}.iDim}, sTmpParafile, handles.allReg{iInd,1}{1}.sShowname(1:min([length(handles.allReg{iInd,1}{1}.sShowname) handles.iMinShow])), handles.allReg{iInd,1}{2}.sShowname(1:min([length(handles.allReg{iInd,1}{2}.sShowname) handles.iMinShow])), handles.allReg{iInd,1}{3}.sShowname(1:min([length(handles.allReg{iInd,1}{3}.sShowname) handles.iMinShow])), handles.allReg{iInd,1}{4}.sShowname(1:min([length(handles.allReg{iInd,1}{4}.sShowname) handles.iMinShow])));
            elseif(length(handles.allReg{iInd,1}) == 4)
                handles.allReg{iInd,2}.sShow = sprintf('reg %02d: %s %s (%s) %s <-- %s, %s, %s',iInd, handles.allReg{iInd,2}.sRegMethod, sTmp{handles.allReg{iInd,2}.iDim}, sTmpParafile, handles.allReg{iInd,1}{1}.sShowname(1:min([length(handles.allReg{iInd,1}{1}.sShowname) handles.iMinShow])), handles.allReg{iInd,1}{2}.sShowname(1:min([length(handles.allReg{iInd,1}{2}.sShowname) handles.iMinShow])), handles.allReg{iInd,1}{3}.sShowname(1:min([length(handles.allReg{iInd,1}{3}.sShowname) handles.iMinShow])), handles.allReg{iInd,1}{4}.sShowname(1:min([length(handles.allReg{iInd,1}{4}.sShowname) handles.iMinShow])));
            elseif(length(handles.allReg{iInd,1}) == 3)
                handles.allReg{iInd,2}.sShow = sprintf('reg %02d: %s %s (%s) %s <-- %s, %s',iInd, handles.allReg{iInd,2}.sRegMethod, sTmp{handles.allReg{iInd,2}.iDim}, sTmpParafile, handles.allReg{iInd,1}{1}.sShowname(1:min([length(handles.allReg{iInd,1}{1}.sShowname) handles.iMinShow])), handles.allReg{iInd,1}{2}.sShowname(1:min([length(handles.allReg{iInd,1}{2}.sShowname) handles.iMinShow])), handles.allReg{iInd,1}{3}.sShowname(1:min([length(handles.allReg{iInd,1}{3}.sShowname) handles.iMinShow])));
            else
                handles.allReg{iInd,2}.sShow = sprintf('reg %02d: %s %s (%s) %s <-- %s',iInd, handles.allReg{iInd,2}.sRegMethod, sTmp{handles.allReg{iInd,2}.iDim}, sTmpParafile, handles.allReg{iInd,1}{1}.sShowname(1:min([length(handles.allReg{iInd,1}{1}.sShowname) handles.iMinShow])), handles.allReg{iInd,1}{2}.sShowname(1:min([length(handles.allReg{iInd,1}{2}.sShowname) handles.iMinShow])));
            end
        end
    end
    % make everything black again
    bColors = {'FFFFFF','0d9f1e'};
    pre1 = '<HTML><BODY bgcolor="';
    pre2 = '"><FONT color="';
    post = '</FONT></BODY></HTML>'; 
    listboxStr = cell(size(handles.allReg,1),1);
    for i = 1:size(handles.allReg,1)
        str = [pre1 bColors{handles.allReg{i,2}.lDone+1} pre2 '000000">' handles.allReg{i,2}.sShow post];
        listboxStr{i} = str;
    end
    set( handles.listRegs, 'String', listboxStr );
end

set(handles.buttonOK,'Visible','off');
set(handles.buttonCancel,'Visible','off');
guidata(hObject,handles);


function buttonCancel_Callback(hObject, eventdata, handles)
% abort updating image/registration information

if(handles.doubleClick(1)) % listImages
    handles.doubleClick(1) = false; % restore double click mode
        
    % make everything black again
    pre = '<HTML><FONT color="';
    post = '</FONT></HTML>'; 
    listboxStr = cell(length(handles.loadedImg),1);
    for i = 1:length(handles.loadedImg)      
        str = [pre '000000">' handles.loadedImg{i}.sShow post]; 
        listboxStr{i} = str;
    end
    set( handles.listImages, 'String', listboxStr );
end
if(handles.doubleClick(2)) % listReg
    handles.doubleClick(2) = false; % restore double click mode
       
    % make everything black again
    bColors = {'FFFFFF','0d9f1e'};
    pre1 = '<HTML><BODY bgcolor="';
    pre2 = '"><FONT color="';
    post = '</FONT></BODY></HTML>'; 
    listboxStr = cell(size(handles.allReg,1),1);
    for i = 1:size(handles.allReg,1)
        str = [pre1 bColors{handles.allReg{i,2}.lDone+1} pre2 '000000">' handles.allReg{i,2}.sShow post];
        listboxStr{i} = str;
    end
    set( handles.listRegs, 'String', listboxStr );
end

set(handles.buttonOK,'Visible','off');
set(handles.buttonCancel,'Visible','off');
guidata(hObject,handles);


function handles = plotImageAndDefField(handles, iReg)
% plot images and deformation field of performed registration

% check if solutions available for current regi
if(~handles.allReg{iReg,2}.lDone)
    return;
end
% clear all first
handles = fClearAllAxes(handles);    
for iJ=1:4 % show images on axis    
    if(iJ > length(handles.allReg{iReg,1}))
        continue;
    end
	iOpen = find(handles.lOpen,1,'first');
	handles.lOpen(iOpen) = false;
    handles.iIndGlobal(iJ) = -iReg;
	dImg = scaleImg(handles.allReg{iReg,1}{iJ}.dImg,[0 1]);
	handles.slice(:,iOpen) = [round(size(dImg,3)/2); size(dImg,3)];
	eval(sprintf('set(handles.hI%d, ''CData'', dImg(:,:,round(size(dImg,3)/2)));',iOpen));
    if(isfield(handles.allReg{iReg,3},'SDeform') && ~isempty(handles.allReg{iReg,3}.SDeform(iJ).dFx))
        dFx = handles.allReg{iReg,3}.SDeform(iJ).dFx;
        dFy = handles.allReg{iReg,3}.SDeform(iJ).dFy;
        iX = 1:handles.quiverFac:size(dFx, 2);
        iY = 1:handles.quiverFac:size(dFy, 1);
        eval(sprintf('set(handles.hQ%d, ''XData'', iX, ''YData'', iY, ''UData'', dFx(1:handles.quiverFac:end, 1:handles.quiverFac:end, handles.slice(1,iJ)), ''VData'', dFy(1:handles.quiverFac:end, 1:handles.quiverFac:end, handles.slice(1,iJ)), ''Color'', ''y'', ''Visible'', ''on'');',iJ));
    end    
	eval(sprintf('set(handles.sliceNo%d, ''String'', ''%02d/%02d'');', iOpen, handles.slice(1,iOpen), handles.slice(2,iOpen)));
	eval(sprintf('set(handles.showName%d, ''String'', ''%s - %02d'');', iOpen, handles.allReg{iReg,1}{iJ}.sShowname, handles.allReg{iReg,1}{iJ}.i4D));
end


function handles = fClearAllAxes(handles)
% clear all axes
set(handles.hI1, 'CData', zeros(256,256));
set(handles.hI2, 'CData', zeros(256,256));
set(handles.hI3, 'CData', zeros(256,256));
set(handles.hI4, 'CData', zeros(256,256));
dFx = zeros(256,256,1);
iX = 1:handles.quiverFac:size(dFx, 2);
iY = 1:handles.quiverFac:size(dFx, 1);
set(handles.hQ2, 'XData', iX, 'YData', iY, 'UData', dFx(1:handles.quiverFac:end, 1:handles.quiverFac:end, 1), 'VData', dFx(1:handles.quiverFac:end, 1:handles.quiverFac:end, 1), 'Visible', 'off');
set(handles.hQ3, 'XData', iX, 'YData', iY, 'UData', dFx(1:handles.quiverFac:end, 1:handles.quiverFac:end, 1), 'VData', dFx(1:handles.quiverFac:end, 1:handles.quiverFac:end, 1), 'Visible', 'off');
set(handles.hQ4, 'XData', iX, 'YData', iY, 'UData', dFx(1:handles.quiverFac:end, 1:handles.quiverFac:end, 1), 'VData', dFx(1:handles.quiverFac:end, 1:handles.quiverFac:end, 1), 'Visible', 'off');

handles.lOpen = true(1,4);
handles.slice = ones(2,4);
handles.iIndGlobal = zeros(1,4);

set(handles.sliceNo1, 'String', '');
set(handles.sliceNo2, 'String', '');
set(handles.sliceNo3, 'String', '');
set(handles.sliceNo4, 'String', '');
set(handles.showName1, 'String', '');
set(handles.showName2, 'String', '');
set(handles.showName3, 'String', '');
set(handles.showName4, 'String', '');


function pb_editParam_Callback(hObject, eventdata, handles)
% edit parameter file
if(~strcmp(handles.sParFile,'NA'))
    edit(handles.sParFile);
end
