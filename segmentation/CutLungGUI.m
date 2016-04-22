function varargout = CutLungGUI(varargin)
% CUTLUNGGUI MATLAB code for CutLungGUI.fig
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% -------------------------------------------------------------------------
% (c) 2015: Thomas Kuestner, Verena Neumann
% -------------------------------------------------------------------------

% Last Modified by GUIDE v2.5 23-Jul-2015 11:00:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CutLungGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @CutLungGUI_OutputFcn, ...
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


% --- Executes just before CutLungGUI is made visible.
function CutLungGUI_OpeningFcn(hObject, eventdata, handles, varargin)

if(nargin > 3)
    for index = 1:2:(nargin-3),
        if nargin-3==index, break, end
        switch lower(varargin{index})
            case 'inarg'
                handles.hData = varargin{index+1};
        end
    end
end

% Choose default command line output for CutLungGUI
handles.output = hObject;

axes(handles.ax_lung)
dImg = handles.hData.dImg(:,:,:,handles.hData.nGate);
handles.slice = floor(size(dImg,3)/2);
lung3d = handles.hData.lung3d{1};
rgbImg = (dImg(:,:,:,[1 1 1])); % make the first a grey-scale image with three channels so it will not be affected by the colormap later on
rgbImg_r = (lung3d(:,:,:)).*(rgbImg(:,:,:,1)+ones(size(lung3d,1),size(lung3d,2),size(lung3d,3))/2);
handles.rgbImg_pm = rgbImg;
handles.rgbImg_pm(:,:,:,2) = rgbImg(:,:,:,2) + rgbImg_r;
handles.rbgImg_pm(handles.rgbImg_pm > 1) = 1;
handles.hImg = imshow(reshape(handles.rgbImg_pm(:,:,handles.slice,:),size(dImg,1),size(dImg,2),3));
hold on;
xX = 1:size(dImg,1); yY = 100*ones(size(dImg,2),1);
handles.hLine = plot(xX,yY,'r');

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes CutLungGUI wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CutLungGUI_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = -round(get(handles.sl_Line,'Value'));
% The figure can be deleted now
delete(handles.figure1);


% --- Executes on slider movement.
function sl_Line_Callback(hObject, eventdata, handles)
% hObject    handle to sl_Line (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
get(hObject,'Value');
dImg = handles.hData.dImg(:,:,:,handles.hData.nGate);
% size(dImg,1)
% size(dImg,2)
set(handles.hLine,'YData',-round(get(hObject,'Value')*ones(size(dImg,2),1)));
guidata(hObject, handles)


function sl_Line_CreateFcn(hObject, eventdata, handles)
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pb_continue.
function pb_continue_Callback(hObject, eventdata, handles)

% The GUI is still in UIWAIT, us UIRESUME
uiresume(gcbf);
guidata(hObject, handles);
return;


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.sl_Line,'Value',-1);
if isequal(get(hObject, 'waitstatus'), 'waiting')
% The GUI is still in UIWAIT, us UIRESUME
uiresume(hObject);
else
% The GUI is no longer waiting, just close it
delete(hObject);
end




% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonUpFcn(hObject, eventdata, handles)


% --- Executes on scroll wheel click while the figure is in focus.
function figure1_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)
% get scroll events
if eventdata.VerticalScrollCount < 0
    handles.slice = max([1 handles.slice - 1]);
else
    handles.slice = min([size(handles.hData.dImg, 3) handles.slice + 1]);
end

dImg = handles.hData.dImg(:,:,:,handles.hData.nGate);
set(handles.hImg, 'CData', reshape(handles.rgbImg_pm(:,:,handles.slice,:),size(dImg,1),size(dImg,2),3));
guidata(hObject, handles)
