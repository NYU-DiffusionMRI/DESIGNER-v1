function varargout = IMGUI(varargin)
% IMGUI(IMG,SCALE,TIME)
% IMG=[4D Time Series] (It won't crash for 3D time series either) 
% SCALE = [min(IMG(:) max(IMG(:))] Colorbar scale 
% TIME=[size(IMG,4),1] Time series with the same length as the 4th column
% of IMG
%
% Arrow Key Commands
% The arrow keys can be used to traverse the each image
% -/+ will descend/ascend through the time series
%
% IMGUI MATLAB code for IMGUI.fig
%      IMGUI, by itself, creates a new IMGUI or raises the existing
%      singleton*.
%
%      H = IMGUI returns the handle to a new IMGUI or the handle to
%      the existing singleton*.
%
%      IMGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMGUI.M with the given input arguments.
%
%      IMGUI('Property','Value',...) creates a new IMGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before IMGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to IMGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help IMGUI

% Last Modified by GUIDE v2.5 30-Nov-2016 18:30:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @IMGUI_OpeningFcn, ...
    'gui_OutputFcn',  @IMGUI_OutputFcn, ...
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


% --- Executes just before IMGUI is made visible.
function IMGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to IMGUI (see VARARGIN)

% Choose default command line output for IMGUI
handles.output = hObject;

if isempty(varargin)
    I=repmat(phantom,[1,1,10,10]);
    for i=1:10
        I(:,:,:,i)=I(:,:,:,i).*exp(-i./10);
    end
    R=[0 max(I(:))];
    DT=1:size(I,4);
else
    I = varargin{1};
    
    if size(varargin,2)==1
        R=[0 max(I(:))];
        DT=1:size(I,4);
    else
        R = varargin{2};
    end
    
    if size(varargin,2)==2
        DT=1:size(I,4);
    end
    if size(varargin,2)==3
        DT=varargin{3};
    end
end
%if (ndims(volume) ~= 3 && ndims(volume) ~= 4)
%    error('Input volume must have 3 or 4 dimensions.');
%end;
handles.I = I;
handles.R = R;
handles.DT= DT;
% axes(handles.XY)
% imagesc(I(:,:,1,1),R)
% colormap gray
% set(gca,'xtick',[],'ytick',[])
% title('XY')
%
% axes(handles.XZ)
% imagesc(squeeze(I(:,round(size(I,2)/2),:,1)),R)
% colormap gray
% set(gca,'xtick',[],'ytick',[])
% title('XZ')
%
% axes(handles.YZ)
% imagesc(squeeze(I(round(size(I,1)/2),:,:,1)),R)
% colormap gray
% set(gca,'xtick',[],'ytick',[])
% title('XZ')


% init 4d pointer
vol_sz = size(I);
if (ndims(I) == 3)
    vol_sz(4) = 1;
end;
pointer4dt = floor(vol_sz/2)+1;

%pointer4dt = round(size(I,1)/2)+1;
handles.vol_sz= vol_sz;
handles.pointer4dt = pointer4dt;

plot6sub(hObject, handles);
% stores ID of last axis window
% (0 means that no axis was clicked yet)
handles.last_axis_id = 0;

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes IMGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = IMGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function XY_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Subplot1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% This object contains the XY slice

pt=get(gca,'currentpoint');
xpos=round(pt(1,2)); ypos=round(pt(1,1));
zpos = handles.pointer4dt(3);
tpos = handles.pointer4dt(4);
handles.pointer4dt = [xpos ypos zpos tpos];
handles.pointer4dt = clipointer4d(handles.pointer4dt,handles.vol_sz);
plot6sub(hObject, handles);
% store this axis as last clicked region
handles.last_axis_id = 1;
% Update handles structure
guidata(hObject, handles);

function XZ_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Subplot1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% This object contains the XZ slice

pt=get(gca,'currentpoint');
xpos=round(pt(1,2)); zpos=round(pt(1,1));
ypos = handles.pointer4dt(2);
tpos = handles.pointer4dt(4);
handles.pointer4dt = [xpos ypos zpos tpos];
handles.pointer4dt = clipointer4d(handles.pointer4dt,handles.vol_sz);
plot6sub(hObject, handles);
% store this axis as last clicked region
handles.last_axis_id = 2;
% Update handles structure
guidata(hObject, handles);

function YZ_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Subplot1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% This object contains the YZ slice

pt=get(gca,'currentpoint');
ypos=round(pt(1,2)); zpos=round(pt(1,1));
xpos = handles.pointer4dt(1);
tpos = handles.pointer4dt(4);
handles.pointer4dt = [xpos ypos zpos tpos];
handles.pointer4dt = clipointer4d(handles.pointer4dt,handles.vol_sz);
plot6sub(hObject, handles);
% store this axis as last clicked region
handles.last_axis_id = 3;
% Update handles structure
guidata(hObject, handles);

function XT_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Subplot1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% This object contains the XT slice

pt=get(gca,'currentpoint');
xpos=round(pt(1,2)); tpos=round(pt(1,1));
ypos = handles.pointer4dt(2);
zpos = handles.pointer4dt(3);
handles.pointer4dt = [xpos ypos zpos tpos];
handles.pointer4dt = clipointer4d(handles.pointer4dt,handles.vol_sz);
plot6sub(hObject, handles);
% store this axis as last clicked region
handles.last_axis_id = 4;
% Update handles structure
guidata(hObject, handles);

function YT_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Subplot1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% This object contains the YT slice

pt=get(gca,'currentpoint');
ypos=round(pt(1,2)); tpos=round(pt(1,1));
xpos = handles.pointer4dt(1);
zpos = handles.pointer4dt(3);
handles.pointer4dt = [xpos ypos zpos tpos];
handles.pointer4dt = clipointer4d(handles.pointer4dt,handles.vol_sz);
plot6sub(hObject, handles);
% store this axis as last clicked region
handles.last_axis_id = 5;
% Update handles structure
guidata(hObject, handles);

function ZT_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Subplot1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% This object contains the ZT slice

pt=get(gca,'currentpoint');
zpos=round(pt(1,2)); tpos=round(pt(1,1));
xpos = handles.pointer4dt(1);
ypos = handles.pointer4dt(2);
handles.pointer4dt = [xpos ypos zpos tpos];
handles.pointer4dt = clipointer4d(handles.pointer4dt,handles.vol_sz);
plot6sub(hObject, handles);
% store this axis as last clicked region
handles.last_axis_id = 6;
% Update handles structure
guidata(hObject, handles);

% --- Plots all 6 slices XY, YZ, XZ, XT, YT, ZT into 6 subplots
function plot6sub(hObject, handles)
% pointer4d     4d coordinates in volume matrix (integers)

value4dt = handles.I(handles.pointer4dt(1), handles.pointer4dt(2), handles.pointer4dt(3), handles.pointer4dt(4));
clims=handles.R;
DT=handles.DT;
text_str = ['[X:' int2str(handles.pointer4dt(2)) ...
    ', Y:' int2str(handles.pointer4dt(1)) ...
    ', Z:' int2str(handles.pointer4dt(3)) ...
    ', Time:' int2str(handles.pointer4dt(4)) '/' int2str(handles.vol_sz(4)) ...
    '], value:' num2str(value4dt)];

set(handles.Point, 'String', text_str);
guidata(hObject, handles);

sliceXY = squeeze(handles.I(:,:,handles.pointer4dt(3),handles.pointer4dt(4)));
sliceYZ = squeeze(handles.I(handles.pointer4dt(1),:,:,handles.pointer4dt(4)));
sliceXZ = squeeze(handles.I(:,handles.pointer4dt(2),:,handles.pointer4dt(4)));
%sliceZY = squeeze(permute(sliceYZ, [2 1 3]));

sliceXT = squeeze(handles.I(:,handles.pointer4dt(2),handles.pointer4dt(3),:));
sliceYT = squeeze(handles.I(handles.pointer4dt(1),:,handles.pointer4dt(3),:));
sliceZT = squeeze(handles.I(handles.pointer4dt(1),handles.pointer4dt(2),:,:));

NT=squeeze(handles.I((handles.pointer4dt(1)-1):(handles.pointer4dt(1)+1),(handles.pointer4dt(2)-1):(handles.pointer4dt(2)+1),handles.pointer4dt(3),:));
NT=squeeze(mean(mean(NT,2),1));
axes(handles.XY);
imagesc(sliceXY, clims);
title('Slice XY');
line([handles.pointer4dt(2) handles.pointer4dt(2)], [0 size(handles.I,1)]);
line([0 size(handles.I,2)], [handles.pointer4dt(1) handles.pointer4dt(1)]);
%set(allchild(gca),'ButtonDownFcn',@Subplot1_ButtonDownFcn);
set(allchild(gca),'ButtonDownFcn','IMGUI(''XY_ButtonDownFcn'',gca,[],guidata(gcbo))');
set(gca,'xtick',[],'ytick',[])
colormap gray

axes(handles.XZ);
imagesc(sliceXZ, clims);
title('Slice YZ');
line([handles.pointer4dt(3) handles.pointer4dt(3)], [0 size(handles.I,1)]);
line([0 size(handles.I,3)], [handles.pointer4dt(1) handles.pointer4dt(1)]);
%set(allchild(gca),'ButtonDownFcn',@Subplot2_ButtonDownFcn);
set(allchild(gca),'ButtonDownFcn','IMGUI(''XZ_ButtonDownFcn'',gca,[],guidata(gcbo))');
set(gca,'xtick',[],'ytick',[])
colormap gray

axes(handles.YZ);
imagesc(sliceYZ, clims);
title('Slice XZ');
line([handles.pointer4dt(3) handles.pointer4dt(3)], [0 size(handles.I,2)]);
line([0 size(handles.I,3)], [handles.pointer4dt(2) handles.pointer4dt(2)]);
set(gca,'xtick',[],'ytick',[])
%set(allchild(gca),'ButtonDownFcn',@Subplot3_ButtonDownFcn);
set(allchild(gca),'ButtonDownFcn','IMGUI(''YZ_ButtonDownFcn'',gca,[],guidata(gcbo))');
colormap gray

axes(handles.YT);
imagesc(sliceXT, clims);
title('Slice YT');
line([handles.pointer4dt(4) handles.pointer4dt(4)], [0 size(handles.I,1)]);
line([0 size(handles.I,4)], [handles.pointer4dt(1) handles.pointer4dt(1)]);
%set(allchild(gca),'ButtonDownFcn',@Subplot1_ButtonDownFcn);
set(allchild(gca),'ButtonDownFcn','IMGUI(''XT_ButtonDownFcn'',gca,[],guidata(gcbo))');
set(gca,'xtick',[],'ytick',[])
colormap gray

axes(handles.XT);
imagesc(sliceYT, clims);
title('Slice XT');
line([handles.pointer4dt(4) handles.pointer4dt(4)], [0 size(handles.I,2)]);
line([0 size(handles.I,4)], [handles.pointer4dt(2) handles.pointer4dt(2)]);
%set(allchild(gca),'ButtonDownFcn',@Subplot2_ButtonDownFcn);
set(allchild(gca),'ButtonDownFcn','IMGUI(''YT_ButtonDownFcn'',gca,[],guidata(gcbo))');
set(gca,'xtick',[],'ytick',[])
colormap gray

axes(handles.ZT);
imagesc(sliceZT, clims);
title('Slice ZT');
line([handles.pointer4dt(4) handles.pointer4dt(4)], [0 size(handles.I,3)]);
line([0 size(handles.I,4)], [handles.pointer4dt(3) handles.pointer4dt(3)]);
set(gca,'xtick',[],'ytick',[])
%set(allchild(gca),'ButtonDownFcn',@Subplot3_ButtonDownFcn);
set(allchild(gca),'ButtonDownFcn','IMGUI(''ZT_ButtonDownFcn'',gca,[],guidata(gcbo))');
colormap gray

%DT(handles.pointer4dt(4));

axes(handles.Curves)
cla(handles.Curves,'reset')
plot(DT,NT,'.-','linewidth',2); hold on;
plot([DT(handles.pointer4dt(4)) DT(handles.pointer4dt(4))],[min(NT(:))*0.99, max(NT(:))*1.01],'r--')
if get(handles.xSemilog,'Value')
set(gca,'xscale','log');
end
if get(handles.ySemilog,'Value')
set(gca,'yscale','log');
end
%set(gca,'xtick',[])


function pointer4d_out = clipointer4d(pointer4d_in,vol_size)
pointer4d_out = pointer4d_in;
for p_id=1:4
    if (pointer4d_in(p_id) > vol_size(p_id))
        pointer4d_out(p_id) = vol_size(p_id);
    end;
    if (pointer4d_in(p_id) < 1)
        pointer4d_out(p_id) = 1;
    end;
end;

function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to SliceBrowserFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%disp('KeyPressFcn');
curr_char = int8(get(gcf,'CurrentCharacter'));


if isempty(curr_char)
    return;
end;

xpos = handles.pointer4dt(1);
ypos = handles.pointer4dt(2);
zpos = handles.pointer4dt(3);
tpos = handles.pointer4dt(4);
% Keys:
% - up:   30
% - down:   31
% - left:   28
% - right:   29
% - '1': 49
% - '2': 50
% - '3': 51
% - 'e': 101
% - plus:  43
% - minus:  45
switch curr_char
    case 30 %UP
        switch handles.last_axis_id
            case 1
                xpos = xpos -1;
            case 2
                xpos = xpos -1;
            case 3
                ypos = ypos -1;
            case 4
                xpos = xpos -1;
            case 5
                ypos = ypos -1;
            case 6
                zpos = zpos -1;
            case 0
        end;
    case 31 %DOWN
        switch handles.last_axis_id
            case 1
                xpos = xpos +1;
            case 2
                xpos = xpos +1;
            case 3
                ypos = ypos +1;
            case 4
                xpos = xpos +1;
            case 5
                ypos = ypos +1;
            case 6
                zpos = zpos +1;
            case 0
        end;
    case 28 %LEFT
        switch handles.last_axis_id
            case 1
                ypos = ypos -1;
            case 2
                zpos = zpos -1;
            case 3
                zpos = zpos -1;
            case 4
                tpos = tpos -1;
            case 5
                tpos = tpos -1;
            case 6
                tpos = tpos -1;
            case 0
        end;
    case 29 %RIGHT
        switch handles.last_axis_id
            case 1
                ypos = ypos +1;
            case 2
                zpos = zpos +1;
            case 3
                zpos = zpos +1;
            case 4
                tpos = tpos +1;
            case 5
                tpos = tpos +1;
            case 6
                tpos = tpos +1;
            case 0
        end;
    case 43
        % plus key
        tpos = tpos+1;
    case 45
        % minus key
        tpos = tpos-1;
    case 49
        % key 1
        handles.last_axis_id = 1;
    case 50
        % key 2
        handles.last_axis_id = 2;
    case 51
        % key 3
        handles.last_axis_id = 3;
    otherwise
        return
end;
handles.pointer4dt = [xpos ypos zpos tpos];
handles.pointer4dt = clipointer4d(handles.pointer4dt,handles.vol_sz);
plot6sub(hObject, handles);
% Update handles structure
guidata(hObject, handles);


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in xSemilog.
function xSemilog_Callback(hObject, eventdata, handles)
% hObject    handle to xSemilog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of xSemilog


% --- Executes on button press in ySemilog.
function ySemilog_Callback(hObject, eventdata, handles)
% hObject    handle to ySemilog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ySemilog
