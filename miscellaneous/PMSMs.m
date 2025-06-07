function varargout = PMSMs(varargin)
% PMSMS MATLAB code for PMSMs.fig
%      PMSMS, by itself, creates a new PMSMS or raises the existing
%      singleton*.
%
%      H = PMSMS returns the handle to a new PMSMS or the handle to
%      the existing singleton*.
%
%      PMSMS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PMSMS.M with the given input arguments.
%
%      PMSMS('Property','Value',...) creates a new PMSMS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PMSMs_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PMSMs_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PMSMs

% Last Modified by GUIDE v2.5 03-Jun-2017 04:43:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PMSMs_OpeningFcn, ...
                   'gui_OutputFcn',  @PMSMs_OutputFcn, ...
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


% --- Executes just before PMSMs is made visible.
function PMSMs_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PMSMs (see VARARGIN)

% Choose default command line output for PMSMs
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PMSMs wait for user response (see UIRESUME)
% uiwait(handles.main);

handles.dimt.ColumnEditable = [false,true,false,false]; 
handles.dimt.ColumnName = {'Parameter','Value','Unit','Description'}; 

index = 0;

index = index + 1;
handles.dimt.Data{index,1} = 'Ns';
handles.dimt.Data{index,2} = 18;
handles.dimt.Data{index,3} = '';
handles.dimt.Data{index,4} = 'Number of Stator Slots';

index = index + 1;
handles.dimt.Data{index,1} = 'Nm';
handles.dimt.Data{index,2} = 8;
handles.dimt.Data{index,3} = '';
handles.dimt.Data{index,4} = 'Number of Magnets';

index = index + 1;
handles.dimt.Data{index,1} = 'm';
handles.dimt.Data{index,2} = 3;
handles.dimt.Data{index,3} = '';
handles.dimt.Data{index,4} = 'Number of Phases';

index = index + 1;
handles.dimt.Data{index,1} = 'Rso';
handles.dimt.Data{index,2} = 50;
handles.dimt.Data{index,3} = '[mm]';
handles.dimt.Data{index,4} = 'Outer Radius of Stator Lamination';

index = index + 1;
handles.dimt.Data{index,1} = 'Rro';
handles.dimt.Data{index,2} = 27.5;
handles.dimt.Data{index,3} = '[mm]';
handles.dimt.Data{index,4} = 'Outer Rotor Radius';

index = index + 1;
handles.dimt.Data{index,1} = 'wsy';
handles.dimt.Data{index,2} = 6.9;
handles.dimt.Data{index,3} = '[mm]';
handles.dimt.Data{index,4} = 'Width of Stator Yoke';

index = index + 1;
handles.dimt.Data{index,1} = 'wry';
handles.dimt.Data{index,2} = 6.9;
handles.dimt.Data{index,3} = '[mm]';
handles.dimt.Data{index,4} = 'Width of Rotor Yoke';

index = index + 1;
handles.dimt.Data{index,1} = 'wtb';
handles.dimt.Data{index,2} = 5.1;
handles.dimt.Data{index,3} = '[mm]';
handles.dimt.Data{index,4} = 'Width of Stator Teeth';

index = index + 1;
handles.dimt.Data{index,1} = 'g';
handles.dimt.Data{index,2} = 1;
handles.dimt.Data{index,3} = '[mm]';
handles.dimt.Data{index,4} = 'Air Gap Length';

index = index + 1;
handles.dimt.Data{index,1} = 'Lst';
handles.dimt.Data{index,2} = 100;
handles.dimt.Data{index,3} = '[mm]';
handles.dimt.Data{index,4} = 'Stack Length';

index = index + 1;
handles.dimt.Data{index,1} = 'hm';
handles.dimt.Data{index,2} = 4;
handles.dimt.Data{index,3} = '[mm]';
handles.dimt.Data{index,4} = 'Height of Magnets';

index = index + 1;
handles.dimt.Data{index,1} = 'betam';
handles.dimt.Data{index,2} = 0.89;
handles.dimt.Data{index,3} = '';
handles.dimt.Data{index,4} = 'Ration of Magnet Arc to Rotor Pole Pitch';

index = index + 1;
handles.dimt.Data{index,1} = 'alphas';
handles.dimt.Data{index,2} = 0.3;
handles.dimt.Data{index,3} = '';

index = index + 1;
handles.dimt.Data{index,1} = 'hss1';
handles.dimt.Data{index,2} = 1.5;
handles.dimt.Data{index,3} = '[mm]';

index = index + 1;
handles.dimt.Data{index,1} = 'hss2';
handles.dimt.Data{index,2} = 0.5;
handles.dimt.Data{index,3} = '[mm]';

% set winding table

handles.wint.ColumnName = {'Parameter','Value','Unit','Description'};

index = 0;

index = index + 1;
handles.wint.Data{index,1} = 'a';
handles.wint.Data{index,2} = 1;
handles.wint.Data{index,3} = '';
handles.wint.Data{index,4} = 'Number of parallel path';

index = index + 1;
handles.wint.Data{index,1} = 'Kf';
handles.wint.Data{index,2} = 0.45;
handles.wint.Data{index,3} = '';
handles.wint.Data{index,4} = 'Fill factor';

index = index + 1;
handles.wint.Data{index,1} = 'Ncoil';
handles.wint.Data{index,2} = 27;
handles.wint.Data{index,3} = '';
handles.wint.Data{index,4} = 'Number of coil turns';

index = index + 1;
handles.wint.Data{index,1} = 'sigma';
handles.wint.Data{index,2} = 57e6;
handles.wint.Data{index,3} = '[S/m]';
handles.wint.Data{index,4} = 'Electrical conductivity of winding material at 25[C]';

index = index + 1;
handles.wint.Data{index,1} = 'dwb';
handles.wint.Data{index,2} = 0.5;
handles.wint.Data{index,3} = '[mm]';
handles.wint.Data{index,4} = 'Bared wire diameter';

index = index + 1;
handles.wint.Data{index,1} = 'nd';
handles.wint.Data{index,2} = 20;
handles.wint.Data{index,3} = '';
handles.wint.Data{index,4} = '';

index = index + 1;
handles.wint.Data{index,1} = 'winT';
handles.wint.Data{index,2} = 60;
handles.wint.Data{index,3} = '[C]';
handles.wint.Data{index,4} = 'Average winding termerature';

index = index + 1;
handles.wint.Data{index,1} = 'alphaT';
handles.wint.Data{index,2} = 3.81e-3;
handles.wint.Data{index,3} = '[1/K]';
handles.wint.Data{index,4} = 'Temperature coefficient of resistivity';

% set loss table

handles.losst.ColumnEditable = [false,true,false,false];
handles.losst.ColumnName = {'Parameter','Value','Unit','Description'}; 

index = 0;

index = index + 1;
handles.losst.Data{index,1} = 'Kh';
handles.losst.Data{index,2} = 6.673865e-3;
handles.losst.Data{index,3} = '[W/Kg]';
handles.losst.Data{index,4} = 'Hysteresis Loss Coefficient';

index = index + 1;
handles.losst.Data{index,1} = 'alpha';
handles.losst.Data{index,2} = 1.2916;
handles.losst.Data{index,3} = '';
handles.losst.Data{index,4} = 'Power of Frequency in Hysteresis Term of Iron Loss Model';

index = index + 1;
handles.losst.Data{index,1} = 'Ke';
handles.losst.Data{index,2} = 66.156534e-6;
handles.losst.Data{index,3} = '[W/Kg]';
handles.losst.Data{index,4} = 'Eddy Current Loss Coefficient';

index = index + 1;
handles.losst.Data{index,1} = 'Density';
handles.losst.Data{index,2} = 7650;
handles.losst.Data{index,3} = '[Kg/m^3]';
handles.losst.Data{index,4} = 'Density of Electrical Steel';

% set datasheet table

handles.cont.ColumnEditable = [false,true,false,false];
handles.cont.ColumnName = {'Parameter','Value','Unit','Description'}; 

index = 0;

index = index + 1;
handles.cont.Data{index,1} = 'Pout';
handles.cont.Data{index,2} = 0.5;
handles.cont.Data{index,3} = '[KW]';
handles.cont.Data{index,4} = 'Rated Output Power';

index = index + 1;
handles.cont.Data{index,1} = 'rpm';
handles.cont.Data{index,2} = 1500;
handles.cont.Data{index,3} = '[rpm]';
handles.cont.Data{index,4} = 'Rated Speed';

index = index + 1;
handles.cont.Data{index,1} = 'Vrms';
handles.cont.Data{index,2} = 220;
handles.cont.Data{index,3} = '[V]';
handles.cont.Data{index,4} = 'Phase Voltage';

index = index + 1;
handles.cont.Data{index,1} = 'eff';
handles.cont.Data{index,2} = 0.8;
handles.cont.Data{index,3} = '';
handles.cont.Data{index,4} = 'Efficiency';

index = index + 1;
handles.cont.Data{index,1} = 'pf';
handles.cont.Data{index,2} = 0.95;
handles.cont.Data{index,3} = '';
handles.cont.Data{index,4} = 'Input Power Factor';

% set mst table
index = 0;

index = index + 1;
handles.mst.Data{index,1} = 'thetag';
handles.mst.Data{index,2} = 1.5;
handles.mst.Data{index,3} = '[Degree]';

index = index + 1;
handles.mst.Data{index,1} = 'thetar';
handles.mst.Data{index,2} = 5;
handles.mst.Data{index,3} = '[Degree]';

index = index + 1;
handles.mst.Data{index,1} = 'thetas';
handles.mst.Data{index,2} = 3;
handles.mst.Data{index,3} = '[Degree]';

handles.mst.Data = handles.mst.Data(1:index,:);

% set sim table
index = 0;

index = index + 1;
handles.simt.Data{index,1} = 'Nsim';
handles.simt.Data{index,2} = 8;
handles.simt.Data{index,3} = '';

index = index + 1;
handles.simt.Data{index,1} = 'Hc';
handles.simt.Data{index,2} = 1.1629e6;
handles.simt.Data{index,3} = '[A/m]';

index = index + 1;
handles.simt.Data{index,1} = 'RelErr';
handles.simt.Data{index,2} = 0.01;
handles.simt.Data{index,3} = '';

index = index + 1;
handles.simt.Data{index,1} = 'MaxIter';
handles.simt.Data{index,2} = 20;
handles.simt.Data{index,3} = '';

handles.simt.Data = handles.simt.Data(1:index,:);

% set cloos table
index = 0;

index = index + 1;
handles.closst.Data{index,1} = 'Stator Hysteresis';
handles.closst.Data{index,2} = 0;
handles.closst.Data{index,3} = '[W]';

index = index + 1;
handles.closst.Data{index,1} = 'Rotor Hysteresis';
handles.closst.Data{index,2} = 0;
handles.closst.Data{index,3} = '[W]';

index = index + 1;
handles.closst.Data{index,1} = 'Stator Eddy Current';
handles.closst.Data{index,2} = 0;
handles.closst.Data{index,3} = '[W]';

index = index + 1;
handles.closst.Data{index,1} = 'Rotor Eddy Current';
handles.closst.Data{index,2} = 0;
handles.closst.Data{index,3} = '[W]';

index = index + 1;
handles.closst.Data{index,1} = 'Copper';
handles.closst.Data{index,2} = 0;
handles.closst.Data{index,3} = '[W]';

index = index + 1;
handles.closst.Data{index,1} = 'Mechanical';
handles.closst.Data{index,2} = 0;
handles.closst.Data{index,3} = '[W]';

index = index + 1;
handles.closst.Data{index,1} = 'Magnet';
handles.closst.Data{index,2} = 0;
handles.closst.Data{index,3} = '[W]';

handles.main.UserData.states.isWinSet = false;
handles.main.UserData.states.isMGSet = false;
handles.main.UserData.states.isMGSplit = true;
handles.main.UserData.states.isRunnedNoLoad = false;
handles.main.UserData.states.isRunnedStatic = false;
handles.main.UserData.states.isRunnedInductance = false;
handles.main.UserData.states.isStatorLam = true;
handles.main.UserData.states.isRotorLam = true;
handles.main.UserData.states.isMagnetSeg = false;

ImportDataInt(handles)
handles.outt.Data = cell(1,3);

ih = 0.02;
iw = 0.03;
Hmsg = 0.15;
Hwinp = 0.35;
Hdimp = 1-Hwinp-Hmsg;
Hmsp = 0.25;
Hconp = 0.35;
Hlossp = 1-Hmsp-Hconp-Hmsg;
Houtp = 0.3;
Hclossp = 0.3;
Hsimp = 1-Hclossp-Houtp-Hmsg;
Wp = (1-4*iw)/3;

handles.msg.Units = 'normalized';
x1 = iw;
y1 = ih;
x2 = 1-2*iw;
y2 = Hmsg-3*ih/2;
handles.msg.Position = [x1,y1,x2,y2];

handles.winp.Units = 'normalized';
x1 = iw;
y1 = Hmsg+ih/2;
x2 = Wp;
y2 = Hwinp-ih;
handles.winp.Position = [x1,y1,x2,y2];

handles.dimp.Units = 'normalized';
x1 = iw;
y1 = Hmsg+Hwinp+ih/2;
x2 = Wp;
y2 = Hdimp-3*ih/2;
handles.dimp.Position = [x1,y1,x2,y2];

handles.msp.Units = 'normalized';
x1 = Wp+2*iw;
y1 = Hmsg+ih/2;
x2 = Wp;
y2 = Hmsp-ih;
handles.msp.Position = [x1,y1,x2,y2];

handles.conp.Units = 'normalized';
x1 = Wp+2*iw;
y1 = Hmsg+Hmsp+ih/2;
x2 = Wp;
y2 = Hconp-ih;
handles.conp.Position = [x1,y1,x2,y2];

handles.lossp.Units = 'normalized';
x1 = Wp+2*iw;
y1 = Hmsg+Hconp+Hmsp+ih/2;
x2 = Wp;
y2 = Hlossp-3*ih/2;
handles.lossp.Position = [x1,y1,x2,y2];


handles.outp.Units = 'normalized';
x1 = 2*Wp+3*iw;
y1 = Hmsg+ih/2;
x2 = Wp;
y2 = Houtp-ih;
handles.outp.Position = [x1,y1,x2,y2];

handles.clossp.Units = 'normalized';
x1 = 2*Wp+3*iw;
y1 = Hmsg+Houtp+ih/2;
x2 = Wp;
y2 = Hclossp-ih;
handles.clossp.Position = [x1,y1,x2,y2];

handles.simp.Units = 'normalized';
x1 = 2*Wp+3*iw;
y1 = Hmsg+Houtp+Hclossp+ih/2;
x2 = Wp;
y2 = Hsimp-3*ih/2;
handles.simp.Position = [x1,y1,x2,y2];

iWt = 0.05;
iHt = 0.05;
handles.dimt.Units = 'normalized';
handles.dimt.Position = [iWt,iHt,1-2*iWt,1-2*iHt];

handles.wint.Units = 'normalized';
handles.wint.Position = [iWt,iHt,1-2*iWt,1-2*iHt];

handles.losst.Units = 'normalized';
handles.losst.Position = [iWt,iHt,1-2*iWt,1-2*iHt];

handles.cont.Units = 'normalized';
handles.cont.Position = [iWt,iHt,1-2*iWt,1-2*iHt];

handles.mst.Units = 'normalized';
handles.mst.Position = [iWt,iHt,1-2*iWt,1-2*iHt];

handles.simt.Units = 'normalized';
handles.simt.Position = [iWt,iHt,1-2*iWt,1-2*iHt];

handles.closst.Units = 'normalized';
handles.closst.Position = [iWt,iHt,1-2*iWt,1-2*iHt];

handles.outt.Units = 'normalized';
handles.outt.Position = [iWt,iHt,1-2*iWt,1-2*iHt];


% --- Outputs from this function are returned to the command line.
function varargout = PMSMs_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;
varargout{1} = handles.main.UserData;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_4_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_5_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_13_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RunNoLoad(handles)

% --------------------------------------------------------------------
function m_cgm_Callback(hObject, eventdata, handles)
% hObject    handle to m_cgm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GenerateMesh(hObject, eventdata, handles)

% --------------------------------------------------------------------
function m_shgm_Callback(hObject, eventdata, handles)
% hObject    handle to m_shgm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ShowGlobalMesh(handles)

function ShowGlobalMesh(handles)
if ~handles.main.UserData.states.isMGSet
    WriteMsg(handles,0,'Geometry and Mesh does not created ...')
    return
end
handles.main.UserData.m.showmesh;

% --------------------------------------------------------------------
function m_shwf_Callback(hObject, eventdata, handles)
% hObject    handle to m_shwf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Showmwf(handles)

% --------------------------------------------------------------------
function m_shfb_Callback(hObject, eventdata, handles)
% hObject    handle to m_shfb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function m_gQ_Callback(hObject, eventdata, handles)
% hObject    handle to m_gQ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GetMeshQuality(handles)

% --------------------------------------------------------------------
function Untitled_6_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SetWinding(hObject, eventdata, handles)


% --------------------------------------------------------------------
function Untitled_7_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_14_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiopen([cd,'\Report.pdf'],1)

% --------------------------------------------------------------------
function Untitled_15_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgbox({'An application for loss analysis of surface permanent',...
    'magnet synchronous motors.',...
        ' _________________________________________ ',...
        ' @Developed by Ali Jamalifard.'},...
    'About','modal')
% --------------------------------------------------------------------
function Untitled_16_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_17_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_18_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_19_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_23_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
InitialDesign;

% --------------------------------------------------------------------
function Untitled_20_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ImportData(handles)

% --------------------------------------------------------------------
function Untitled_21_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ExportData(handles)



% --------------------------------------------------------------------
function Untitled_22_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ExportDataInt(handles)
delete(handles.main)
% --------------------------------------------------------------------
function Untitled_25_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function slot1_Callback(hObject, eventdata, handles)
% hObject    handle to slot1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function slot2_Callback(hObject, eventdata, handles)
% hObject    handle to slot2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_24_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in msg.
function msg_Callback(hObject, eventdata, handles)
% hObject    handle to msg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns msg contents as cell array
%        contents{get(hObject,'Value')} returns selected item from msg


% --- Executes during object creation, after setting all properties.
function msg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to msg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function m_shmzs_Callback(hObject, eventdata, handles)
% hObject    handle to m_shmzs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Showmzs(handles)

function ReadData(handles)


for i = 1:size(handles.dimt.Data,1)
    handles.main.UserData.dim.(handles.dimt.Data{i,1}) = ...
        handles.dimt.Data{i,2};
end


for i = 1:size(handles.losst.Data,1)
    handles.main.UserData.loss.(handles.losst.Data{i,1}) = ...
        handles.losst.Data{i,2};  
end


for i = 1:size(handles.wint.Data,1)
    handles.main.UserData.win.(handles.wint.Data{i,1}) = ...
        handles.wint.Data{i,2};  
end


for i = 1:size(handles.cont.Data,1)
    handles.main.UserData.con.(handles.cont.Data{i,1}) = ...
        handles.cont.Data{i,2};  
end


for i = 1:size(handles.simt.Data,1)
    handles.main.UserData.sim.(handles.simt.Data{i,1}) = ...
        handles.simt.Data{i,2};  
end


for i = 1:size(handles.mst.Data,1)
    handles.main.UserData.ms.(handles.mst.Data{i,1}) = ...
        handles.mst.Data{i,2};  
end


% --------------------------------------------------------------------
function Untitled_29_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function uipushtool7_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Untitled_20_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function uipushtool6_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Untitled_21_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function uipushtool8_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GenerateMesh(hObject, eventdata, handles)


% --------------------------------------------------------------------
function Untitled_31_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_32_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_34_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_35_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_33_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_37_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RunNoLoad(handles)

% --------------------------------------------------------------------
function Untitled_38_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RunInductance(handles)

% --------------------------------------------------------------------
function simc_Callback(hObject, eventdata, handles)
% hObject    handle to simc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_40_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SetWinding(hObject, eventdata, handles)

% --------------------------------------------------------------------
function Untitled_41_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_43_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function winc_Callback(hObject, eventdata, handles)
% hObject    handle to winc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_44_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function GenerateMesh(hObject, eventdata, handles)

WriteMsg(handles,3,'Creation of geometry and mesh generation ...')

materialdir = [cd,'\MaterialsData'];
WaitFig = waitbar(1/8,'Calling Geometry Editor','WindowStyle','modal',...
    'Name','MG0');
g = GDBC2D;

tmp = struct;
for i = 1:size(handles.mst.Data,1)
    tmp.(handles.mst.Data{i,1}) = handles.mst.Data{i,2};
end
thetar = tmp.thetar;
thetas = tmp.thetas;
thetag = tmp.thetag;
clear tmp;

handles.MsgBox(handles);
pm = handles.main.UserData.pm;

%% creation of stator lamination
p1 = [pm.Rro+pm.g,0];
p2 = [pm.Rso,0];
p3 = pmove(p2,'theta',pm.Ts/2);
p4 = pmove(p3,'r',-pm.wsy);
p5 = pmove(p4,'theta',-pm.Tso2/2);
p6 = pmove(p5,'x',-pm.hss3);
p8 = pmove(p1,'theta',(pm.Ts-pm.Tso1)/2);
p7 = pmove(p8,'p',pm.hss1,[cos(pm.Ts/2),sin(pm.Ts/2)]);
tmp1 = 2*(pm.Rro+pm.g)*sin(thetag*pi/180/2);
tmp2 = 2*(pm.Rso)*sin(thetas*pi/180/2);
tmp3 = 2*(pm.Rso-pm.wsy)*sin(thetas*pi/180/2);
[g,L1] = g.newdlinewdkps(p1,p2,'l1',tmp1,'l2',tmp2);
[g,A1] = g.newdarccppwdkps([0,0],p2,p3,'maxDegree',thetas);
[g,L2] = g.newdlinewdkps(p3,p4,'l1',tmp2,'l2',tmp3);
[g,A2] = g.newdarccppwdkps([0,0],p4,p5,'direction',-1,'maxDegree',thetas);
[g,L3] = g.newdlinewdkps(p6,p5,'l1',tmp1,'l2',tmp3);
[g,L4] = g.newdlinewdkps(p6,p7,'maxLength',tmp1);
[g,L5] = g.newdlinewdkps(p7,p8,'maxLength',tmp1);
[g,A3] = g.newdarccppwdkps([0,0],p8,p1,'direction',-1,'maxDegree',thetag);
g = g.newcb('s1',L1,1,A1,1,L2,1,A2,1,L3,-1,L4,1,L5,1,A3,1);
%% cration of coil
p9 = pmove(p1,'theta',pm.Ts/2);
p9 = pmove(p9,'r',pm.hss1);
g = g.newdlinewdkps(p7,p9,'maxLength',tmp1);
g = g.newdlinewdkps(p4,p9,'l1',tmp3,'l2',tmp1);
g = g.newcb('c11','L3',1,'A2',-1,'L7',1,'L6',-1,'L4',-1);
%% creation of slot air
g = g.cmirrorbd('L5',[cos(pm.Ts/2),sin(pm.Ts/2)]);
g = g.cmirrorbd('L6',[cos(pm.Ts/2),sin(pm.Ts/2)]);
g = g.newdarccpp('kp8','kp11',[0,0],'maxDegree',thetag);
g = g.newcb('sap1','L5',-1,'L6',1,'L9',-1,'L8',1,'A4',-1);
%% magnet
p1 = [pm.rsh+pm.wry,0];
p2 = p1+[pm.hm,0];
p3 = pmove(p2,'theta',pm.betam*pi/pm.Nm);
p4 = pmove(p3,'r',-pm.hm);
p5 = pmove(p4,'theta',(1-pm.betam)*pi/pm.Nm);
p7 = [pm.rsh,0];
p6 = pmove(p7,'theta',pi/pm.Nm);
g = g.newdlinewdkps(p1,p2,'maxLength',tmp1);
g = g.newdarccppwdkps([0,0],p2,p3,'maxDegree',thetag);
g = g.newdlinewdkps(p4,p3,'maxLength',tmp1);
g = g.newdarccppwdkps([0,0],p4,p1,'direction',-1,'maxDegree',thetag);
g = g.newcb('m1','L10',1,'A5',1,'L11',-1,'A6',1);
%% rotor lamination
tmp2 = 2*(pm.rsh)*sin(thetar*pi/180/2);
g = g.newdarccppwdkps([0,0],p4,p5,'maxDegree',thetag);
g = g.newdlinewdkps(p5,p6,'l1',tmp1,'l2',tmp2);
g = g.newdarccppwdkps([0,0],p6,p7,'direction',-1,'maxDegree',thetar);
g = g.newdlinewdkps(p1,p7,'l1',tmp1,'l2',tmp2);
g = g.newcb('r1','L13',-1,'A6',-1,'A7',1,'L12',1,'A8',1);
%% rotor air pocket
g = g.newdlinewdkps(p5,pmove(p2,'theta',pi/pm.Nm),'maxLength',tmp1);
g = g.newdarccpp('kp14','kp19',[0,0],'maxDegree',thetag);
g = g.newcb('rap1','L11',1,'A9',1,'L14',-1,'A7',-1);
%% creation of domains
waitbar(2/8,WaitFig,'MG for Stator');
g = g.newdDM('s1','s1');
waitbar(3/8,WaitFig,'MG for Coil');
g = g.newdDM('c11','c11');
waitbar(4/8,WaitFig,'MG for Stator Air Pocket');
g = g.newdDM('sap1','sap1');
waitbar(5/8,WaitFig,'MG for Magnet');
g = g.newdDM('m1','m1');
waitbar(6/8,WaitFig,'MG for Rotor');
g = g.newdDM('r1','r1');
waitbar(7/8,WaitFig,'MG for Rotor Air Pocket');
g = g.newdDM('rap11','rap1');
waitbar(1,WaitFig,'Mapping Mesh Zones and Creation of Global Mesh');
%% creation of mesh
m = MDBCT(g);clear g;
% adding needed materials
m = m.addMaterial(materialdir,'air');
m = m.addMaterial(materialdir,'m19_24ga');
m = m.addMaterial(materialdir,'copper');
% setting materials
m = m.setMaterial('s1','m19_24ga');
m = m.setMaterial('c11','copper');
m = m.setMaterial('r1','m19_24ga');
% setting mesh zone colors
m = m.setmzColor('s1',[0 191 255]/255);
m = m.setmzColor('r1',[0 191 255]/255);
m = m.setmzColor('c11',[255 140 0]/255);
m = m.setmzColor('sap1','w');
m = m.setmzColor('rap11','w');
m = m.setmzColor('m1','m');
% creation of other mesh zones
m = m.cmirrormz('s2','s1',[cos(pm.Ts/2),sin(pm.Ts/2)]);
m = m.cmirrormz('c21','c11',[cos(pm.Ts/2),sin(pm.Ts/2)]);
m = m.setmzColor('c21','r');
m = m.cmirrormz('rap21','rap11',[1,0]);

if strcmp(handles.mSplit.Checked,'on')
    sNsplit = pm.Ns/pm.Nsplit;
    rNsplit = pm.Nm/pm.Nsplit;
else
    sNsplit = pm.Ns;
    rNsplit = pm.Nm;
end

for i = 1:2:2*(sNsplit-1)
    m = m.crotatemz(['s',num2str(i+2)],['s',num2str(i)],pm.Ts);
    m = m.crotatemz(['s',num2str(i+3)],['s',num2str(i+1)],pm.Ts);
end
m = m.cmirrormz('r2','r1',[1,0]);
m = m.cmirrormz('m2','m1',[1,0]);
m = m.joinmzs('magnet1','m1','m2');
for i = 1:(sNsplit-1)
    m = m.crotatemz(['c1',num2str(i+1)],['c1',num2str(i)],pm.Ts);
    m = m.crotatemz(['c2',num2str(i+1)],['c2',num2str(i)],pm.Ts);
    m = m.crotatemz(['sap',num2str(i+1)],['sap',num2str(i)],pm.Ts);
end
for i = 1:(rNsplit-1)
    m = m.crotatemz(['magnet',num2str(i+1)],['magnet',num2str(i)],pm.Tr);
    m = m.crotatemz(['rap1',num2str(i+1)],['rap1',num2str(i)],pm.Tr);
    m = m.crotatemz(['rap2',num2str(i+1)],['rap2',num2str(i)],pm.Tr);
end
for i = 1:2:2*(rNsplit-1)
    m = m.crotatemz(['r',num2str(i+2)],['r',num2str(i)],pm.Tr);
    m = m.crotatemz(['r',num2str(i+3)],['r',num2str(i+1)],pm.Tr);
end
% joining stator mesh zones
temp = getlist('s',1:2*sNsplit);
m = m.joinmzs('stator',temp{:});
% joining rotor mesh zones
temp = getlist('r',1:2*rNsplit);
m = m.joinmzs('rotor',temp{:});
% smoothing stator mesh
mz = TMZPC(m.mzs.stator.t,m.mzs.stator.p);
mz = mz.moveNodes;
m.mzs.stator.p = mz.nodes;
% smoothing rotor mesh
mz = TMZPC(m.mzs.rotor.t,m.mzs.rotor.p);
mz = mz.moveNodes;
m.mzs.rotor.p = mz.nodes;
m = m.ggmesh;
m = m.setdata;

WriteMsg(handles,1,'Creation of geometry and mesh generation compeleted!')
waitbar(1,WaitFig,'Completed');
close(WaitFig)
handles.main.UserData.m = m;
handles.main.UserData.p.isMG = true;

handles.main.UserData.states.isMGSet = true;
handles.main.UserData.states.isRunnedNoLoad = false;

% --------------------------------------------------------------------
function Untitled_46_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GenerateMesh(hObject, eventdata, handles)

% --------------------------------------------------------------------
function meshc_Callback(hObject, eventdata, handles)
% hObject    handle to meshc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function SetWinding(hObject, eventdata, handles)

ReadData(handles);

[FileName,PathName] = uigetfile('*.win','Select .win File','win');
if ~FileName
    return
end
f = fopen([PathName,FileName],'r');
win = fscanf(f,'%d');
fclose(f);

Ns = win(1);
m = win(2);
S = win(3);
L = win(4);

if Ns ~= handles.main.UserData.dim.Ns
    WriteMsg(handles,0,['Windings data does not set properly',...
        ' ... Number of coils must be equal to the number of slots'])
    return
end

if L ~= 2
    WriteMsg(handles,0,'Number of layers should be 2 ...')
    return
end
 
win = reshape(win(5:end),2,[]);
win = win';
win = reshape(win,[],2*m);

for i = 1:m
    phase = win(:,[i,i+m]);
    if any(bitor(phase(:,1)>Ns,phase(:,1)<1))
        WriteMsg(handles,0,['Windings data does not set properly',...
            ' ... Some coil index is grater than Ns'])
        return
    end
    handles.main.UserData.win.(['Phase',num2str(i)]) = phase;
end

handles.main.UserData.win.S = S;
handles.main.UserData.win.m = m;
handles.main.UserData.win.L = L;

coils = (1:Ns)';
coils = [coils,circshift(coils,-S)];
handles.main.UserData.win.coils = coils;

WriteMsg(handles,1,'Winding configuration completed!')
handles.main.UserData.states.isWinSet = true;


% --------------------------------------------------------------------
function Untitled_47_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function EX_NL_Callback(hObject, eventdata, handles)
% hObject    handle to EX_NL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.EX_NL.Checked = 'on';
handles.EX_FL.Checked = 'off';
handles.EX_HA.Checked = 'off';

% --------------------------------------------------------------------
function EX_FL_Callback(hObject, eventdata, handles)
% hObject    handle to EX_FL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.EX_NL.Checked = 'off';
handles.EX_FL.Checked = 'on';
handles.EX_HA.Checked = 'off';

% --------------------------------------------------------------------
function EX_HA_Callback(hObject, eventdata, handles)
% hObject    handle to EX_HA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.EX_NL.Checked = 'off';
handles.EX_FL.Checked = 'off';
handles.EX_HA.Checked = 'on';


% --------------------------------------------------------------------
function Untitled_48_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function Untitled_49_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_50_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function handles.MsgBox(handles)
pm = SD_PMSM0;
for i = 1:size(handles.dimt.Data,1)
    pm.(handles.dimt.Data{i,1}) = handles.dimt.Data{i,2};
end
for i = 1:size(handles.wint.Data,1)
    pm.(handles.wint.Data{i,1}) = handles.wint.Data{i,2};
end
for i = 1:size(handles.cont.Data,1)
    pm.(handles.cont.Data{i,1}) = handles.cont.Data{i,2};
end
handles.main.UserData.pm = pm;


% --------------------------------------------------------------------
function Untitled_52_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_53_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ShowMeshFrames(handles)

% --------------------------------------------------------------------
function Untitled_54_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ShowBmagFrames(handles)

function ShowMeshFrames(handles)

% figure
% for j = 1:length(handles.main.UserData.mF)
%     imshow(handles.main.UserData.mF(j).cdata)
%     pause(0.01);
% end
implay(handles.main.UserData.mF(1:end-1));

function ShowBmagFrames(handles)

figure('Name','Movie Player','WindowStyle','modal','NumberTitle','off'...
,'Color','k')
for i = 1:100
    for j = 1:length(handles.main.UserData.bF)-1
        imshow(handles.main.UserData.bF(j).cdata)
        pause(0.01);
    end
end
% hold all
% axis off equal
% movie(handles.main.UserData.bF,10000,20)
% implay(handles.main.UserData.bF(1:end-1));

% --------------------------------------------------------------------
function Untitled_55_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_56_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_56 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_61_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_61 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GetMeshQuality(handles)

% --------------------------------------------------------------------
function Untitled_57_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_57 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Showmzs(handles)

% --------------------------------------------------------------------
function Untitled_58_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_58 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ShowGlobalMesh(handles)

% --------------------------------------------------------------------
function Untitled_59_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_59 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Showmwf(handles)

% --------------------------------------------------------------------
function Untitled_60_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_60 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function slot3_Callback(hObject, eventdata, handles)
% hObject    handle to slot3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function uipushtool10_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RunNoLoad(handles)

function WriteMsg(handles,type,txt)

handles.msg.String = '';

switch type
    case 0
        handles.msg.ForegroundColor = 'r';
        handles.msg.String{1} = 'Error:';
    case 1
        handles.msg.ForegroundColor = 'b';
        handles.msg.String{1} = 'Done!';
    case 2
        handles.msg.ForegroundColor = [249,124,0]/255;
        handles.msg.String{1} = 'Warning:';
    case 3
        handles.msg.ForegroundColor = 'r';
        handles.msg.String{1} = 'Busy:';
    case 4
        handles.msg.String{1} = 'Results:';
        handles.msg.ForegroundColor = 'k';
end
 
if isa(txt,'cell')
    handles.msg.String = {handles.msg.String(:),txt(:)};
else
    handles.msg.String{2} = txt;
end

function Showmzs(handles)

if ~handles.main.UserData.states.isMGSet
    WriteMsg(handles,0,'Geometry and Mesh does not created ...')
    return
end

handles.main.UserData.m.showmzs;
set(gcf,'WindowStyle','modal')

function Showmwf(handles)

if ~handles.main.UserData.states.isMGSet
    WriteMsg(handles,0,'Geometry and Mesh does not created ...')
    return
end

figure
handles.main.UserData.m.plotwf;
set(gcf,'WindowStyle','modal','NumberTitle','off','Name','Wire Frame Mesh')

% --- Executes when entered data in editable cell(s) in dimt.
function dimt_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to dimt (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
if  eventdata.NewData ~= eventdata.PreviousData
    WriteMsg(handles,2,'Dimensions Changed, you muste create new geometry and mesh ...')
    handles.main.UserData.states.isMGSet = false;
end

handles.outt.Data = cell(1,3);


function ShowHyst(handles)
% hObject    handle to ShowHyst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~handles.main.UserData.states.isRunnedNoLoad
    WriteMsg(handles,0,'No Load simulation is not done')
    return
end

pm = handles.main.UserData.pm;
%% Calculation of Hysteresis Loss Using GSE Method
% angularSpeed
angularSpeed = pm.wm;
%% lamination data
% hysteresis loss coefficient
Kh =  handles.main.UserData.loss.Kh;
% power of frequency
alpha = handles.main.UserData.loss.alpha;
% power of peak of magnetic flux density
beta = 2;
% calculation of KGSE
temp = linspace(0,2*pi,1000);
KGSE = Kh/((2*pi)^(alpha-1))/(sum((abs(cos(temp))).^alpha .* ...
    (abs(sin(temp))).^(beta-alpha))*temp(2));
% zone
ZoneIndex = handles.main.UserData.m.mzs.stator.zi;
bx = handles.main.UserData.results.NoLoad.Bxs;
by = handles.main.UserData.results.NoLoad.Bys;

dbxdt = [bx(:,2)-bx(:,1),...
    (bx(:,3:end)-bx(:,1:end-2))/2,...
    bx(:,end)-bx(:,end-1)];

dbydt = [by(:,2)-by(:,1),...
    (by(:,3:end)-by(:,1:end-2))/2,...
    by(:,end)-by(:,end-1)];

phystx = (abs(dbxdt).^alpha) .* (abs(bx).^(beta-alpha));
physty = (abs(dbydt).^alpha) .* (abs(by).^(beta-alpha));

phystx(:,1) = phystx(:,1)/2;
phystx(:,end) = phystx(:,end)/2;

physty(:,1) = physty(:,1)/2;
physty(:,end) = physty(:,end)/2;

phystx = sum(phystx ,2);
physty = sum(physty ,2);
physt = phystx+physty;

xrotorAngle = handles.main.UserData.results.NoLoad.xrotorAngle;
simulationAngle = handles.main.UserData.results.NoLoad.simulationAngle;

physt = physt*7650*KGSE* angularSpeed^alpha * length(xrotorAngle)^(alpha-1)/(simulationAngle)^alpha;

figure('Name','Hysteresis Loss Density','NumberTitle','off',...
    'WindowStyle','modal')
hold all
t = handles.main.UserData.m.t(handles.main.UserData.m.tzi(:,ZoneIndex),1:3);
tt = 1:3*size(t,1);
t = t';
t = t(:);
tt = reshape(tt,3,[]);
c = colorbar;
c.Color = 'w';
c.Ticks = linspace(min(physt),max(physt),8);
c.Limits = [min(physt),max(physt)];
c.FontSize = 20;
physt = repmat(physt,1,3);
physt = physt';
trisurf(tt',handles.main.UserData.m.p(t,1),...
    handles.main.UserData.m.p(t,2),physt(:),'edgecolor','none');
axis off equal;
view([0,0,1]);
colormap jet;
set(gcf,'Color','k')
title('Hysteresis Loss Density [W/m^3] (GSE Method)',...
    'color','w','fontsize',15);
set(gcf,'Renderer','opengl');
zoom on;

function ShowRotorHyst(handles)
% hObject    handle to ShowHyst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~handles.main.UserData.states.isRunnedNoLoad
    WriteMsg(handles,0,'No Load simulation is not done')
    return
end

pm = handles.main.UserData.pm;
%% Calculation of Hysteresis Loss Using GSE Method
% angularSpeed
angularSpeed = pm.wm;
%% lamination data
% hysteresis loss coefficient
Kh =  handles.main.UserData.loss.Kh;
% power of frequency
alpha = handles.main.UserData.loss.alpha;
% power of peak of magnetic flux density
beta = 2;
% calculation of KGSE
temp = linspace(0,2*pi,1000);
KGSE = Kh/((2*pi)^(alpha-1))/(sum((abs(cos(temp))).^alpha .* ...
    (abs(sin(temp))).^(beta-alpha))*temp(2));
% zone
ZoneIndex = handles.main.UserData.m.mzs.rotor.zi;
bx = handles.main.UserData.results.NoLoad.Bxr;
by = handles.main.UserData.results.NoLoad.Byr;

dbxdt = [bx(:,2)-bx(:,1),...
    (bx(:,3:end)-bx(:,1:end-2))/2,...
    bx(:,end)-bx(:,end-1)];

dbydt = [by(:,2)-by(:,1),...
    (by(:,3:end)-by(:,1:end-2))/2,...
    by(:,end)-by(:,end-1)];

phystx = (abs(dbxdt).^alpha) .* (abs(bx).^(beta-alpha));
physty = (abs(dbydt).^alpha) .* (abs(by).^(beta-alpha));

phystx(:,1) = phystx(:,1)/2;
phystx(:,end) = phystx(:,end)/2;

physty(:,1) = physty(:,1)/2;
physty(:,end) = physty(:,end)/2;

phystx = sum(phystx ,2);
physty = sum(physty ,2);
physt = phystx+physty;

xrotorAngle = handles.main.UserData.results.NoLoad.xrotorAngle;
simulationAngle = handles.main.UserData.results.NoLoad.simulationAngle;

physt = physt*7650*KGSE* angularSpeed^alpha * length(xrotorAngle)^(alpha-1)/(simulationAngle)^alpha;

figure('Name','Hysteresis Loss Density','NumberTitle','off',...
    'WindowStyle','modal')
hold all
t = handles.main.UserData.m.t(handles.main.UserData.m.tzi(:,ZoneIndex),1:3);
tt = 1:3*size(t,1);
t = t';
t = t(:);
tt = reshape(tt,3,[]);
c = colorbar;
c.Color = 'w';
c.Ticks = linspace(min(physt),max(physt),8);
c.Limits = [min(physt),max(physt)];
c.FontSize = 20;
physt = repmat(physt,1,3);
physt = physt';
trisurf(tt',handles.main.UserData.m.p(t,1),...
    handles.main.UserData.m.p(t,2),physt(:),'edgecolor','none');
axis off equal;
view([0,0,1]);
colormap jet;
set(gcf,'Color','k')
title('Hysteresis Loss Density [W/m^3] (GSE Method)',...
    'color','w','fontsize',15);
set(gcf,'Renderer','opengl');
zoom on;
% --------------------------------------------------------------------
function Untitled_65_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_65 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_66_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_66 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_63_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_63 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function Untitled_64_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_64 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function ShowBEMFs(handles)

if ~handles.main.UserData.states.isRunnedNoLoad
    WriteMsg(handles,0,'No Load simulation is not done')
    return
end

pm = handles.main.UserData.pm;
figure('Name','Phase to Phase and Line to Line BEMFs','NumberTitle','off',...
    'WindowStyle','modal')
subplot(211)
hold all
for i = 1:handles.main.UserData.win.m
    e.(['e',num2str(i)]) = ppval(handles.main.UserData.results.NoLoad.(['e',num2str(i)]),...
        handles.main.UserData.results.NoLoad.xi);
    plot(handles.main.UserData.results.NoLoad.xi*180/pi * pm.Nm/2,e.(['e',num2str(i)]),'LineWidth',1.2)
end

tmp = cell(1,handles.main.UserData.win.m);
for i = 1:handles.main.UserData.win.m
    tmp{i} = ['Phase ',num2str(i)];
end

legend(tmp)
xlabel('Electrical Angle')
ylabel('Back EMF [V]')
set(gca,'XLim',[0,360],'XTick',[0,90,180,270,360],'XTickLabel',...
    {'0','\pi/2','\pi','3\pi/2','2\pi'},'FontSize',12,'GridLineStyle',':',...
    'GridColor','k','GridAlpha',1)
grid on
zoom on
title('Phase BEMFs');

subplot(212)
hold all
tmp1 = [1:handles.main.UserData.win.m,1];
for i = 1:handles.main.UserData.win.m
    plot(handles.main.UserData.results.NoLoad.xi*180/pi * pm.Nm/2,...
        (e.(['e',num2str(tmp1(i))])-e.(['e',num2str(tmp1(i+1))])),'LineWidth',1.2);
end
for i = 1:handles.main.UserData.win.m
    tmp{i} = ['E_{',num2str(tmp1(i)),num2str(tmp1(i+1)),'}'];
end
legend(tmp)
xlabel('Electrical Angle')
ylabel('Back EMF [V]')
set(gca,'XLim',[0,360],'XTick',[0,90,180,270,360],'XTickLabel',...
    {'0','\pi/2','\pi','3\pi/2','2\pi'},'FontSize',12,'GridLineStyle',':',...
    'GridColor','k','GridAlpha',1)
grid on
zoom on
title('Line to Line BEMFs');


% --------------------------------------------------------------------
function Untitled_67_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_67 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ShowBEMFs(handles)

function ShowEddy(handles)

if ~handles.main.UserData.states.isRunnedNoLoad
    WriteMsg(handles,0,'No Load simulation is not done')
    return
end

pm = handles.main.UserData.pm;
angularSpeed = pm.wm;

ZoneIndex = handles.main.UserData.m.mzs.stator.zi;
bx = handles.main.UserData.results.NoLoad.Bxs;
by = handles.main.UserData.results.NoLoad.Bys;

dbxdt = [bx(:,2)-bx(:,1),...
    (bx(:,3:end)-bx(:,1:end-2))/2,...
    bx(:,end)-bx(:,end-1)];

dbydt = [by(:,2)-by(:,1),...
    (by(:,3:end)-by(:,1:end-2))/2,...
    by(:,end)-by(:,end-1)];

peddy = (dbxdt.^2+dbydt.^2);

peddy(:,1) = peddy(:,1) / 2;
peddy(:,end) = peddy(:,end) / 2;

peddy = sum(peddy,2);

simulationAngle = handles.main.UserData.results.NoLoad.simulationAngle;
xrotorAngle = handles.main.UserData.results.NoLoad.xrotorAngle;

Ke = handles.main.UserData.loss.Ke;
peddy =  (Ke/2/pi^2)*peddy*7650*angularSpeed^2/simulationAngle/xrotorAngle(2);

figure('Name','Eddy Current Loss Density','NumberTitle','off',...
    'WindowStyle','modal')
hold all

t = handles.main.UserData.m.t(handles.main.UserData.m.tzi(:,ZoneIndex),1:3);
tt = 1:3*size(t,1);
t = t';
t = t(:);
tt = reshape(tt,3,[]);

c = colorbar;
c.Color = 'w';
c.Ticks = linspace(min(peddy),max(peddy),8);
c.Limits = [min(peddy),max(peddy)];
c.FontSize = 20;

peddy = repmat(peddy,1,3);
peddy = peddy';

trisurf(tt',handles.main.UserData.m.p(t,1),...
    handles.main.UserData.m.p(t,2),peddy(:),'edgecolor','none');

axis off equal;
view([0,0,1]);
colormap jet;
set(gcf,'Color','k')
title('Eddy Current Loss Density [W/m^3]',...
    'color','w','fontsize',15);
set(gcf,'Renderer','opengl');
zoom on;

function ShowRotorEddy(handles)

if ~handles.main.UserData.states.isRunnedNoLoad
    WriteMsg(handles,0,'No Load simulation is not done')
    return
end

pm = handles.main.UserData.pm;
angularSpeed = pm.wm;

ZoneIndex = handles.main.UserData.m.mzs.rotor.zi;
bx = handles.main.UserData.results.NoLoad.Bxr;
by = handles.main.UserData.results.NoLoad.Byr;

dbxdt = [bx(:,2)-bx(:,1),...
    (bx(:,3:end)-bx(:,1:end-2))/2,...
    bx(:,end)-bx(:,end-1)];

dbydt = [by(:,2)-by(:,1),...
    (by(:,3:end)-by(:,1:end-2))/2,...
    by(:,end)-by(:,end-1)];

peddy = (dbxdt.^2+dbydt.^2);

peddy(:,1) = peddy(:,1) / 2;
peddy(:,end) = peddy(:,end) / 2;

peddy = sum(peddy,2);

simulationAngle = handles.main.UserData.results.NoLoad.simulationAngle;
xrotorAngle = handles.main.UserData.results.NoLoad.xrotorAngle;

Ke = handles.main.UserData.loss.Ke;
peddy =  (Ke/2/pi^2)*peddy*7650*angularSpeed^2/simulationAngle/xrotorAngle(2);

figure('Name','Eddy Current Loss Density','NumberTitle','off',...
    'WindowStyle','modal')
hold all

t = handles.main.UserData.m.t(handles.main.UserData.m.tzi(:,ZoneIndex),1:3);
tt = 1:3*size(t,1);
t = t';
t = t(:);
tt = reshape(tt,3,[]);

c = colorbar;
c.Color = 'w';
c.Ticks = linspace(min(peddy),max(peddy),8);
c.Limits = [min(peddy),max(peddy)];
c.FontSize = 20;

peddy = repmat(peddy,1,3);
peddy = peddy';

trisurf(tt',handles.main.UserData.m.p(t,1),...
    handles.main.UserData.m.p(t,2),peddy(:),'edgecolor','none');

axis off equal;
view([0,0,1]);
colormap jet;
set(gcf,'Color','k')
title('Eddy Current Loss Density [W/m^3]',...
    'color','w','fontsize',15);
set(gcf,'Renderer','opengl');
zoom on;

function EvalStatorEddyLoss(handles)

pm = handles.main.UserData.pm;
angularSpeed = pm.wm;

ZoneIndex = handles.main.UserData.m.mzs.stator.zi;
bx = handles.main.UserData.results.NoLoad.Bxs;
by = handles.main.UserData.results.NoLoad.Bys;

dbxdt = [bx(:,2)-bx(:,1),...
    (bx(:,3:end)-bx(:,1:end-2))/2,...
    bx(:,end)-bx(:,end-1)];

dbydt = [by(:,2)-by(:,1),...
    (by(:,3:end)-by(:,1:end-2))/2,...
    by(:,end)-by(:,end-1)];

peddy = (dbxdt.^2+dbydt.^2);

peddy(:,1) = peddy(:,1) / 2;
peddy(:,end) = peddy(:,end) / 2;

peddy = sum(peddy,2);

simulationAngle = handles.main.UserData.results.NoLoad.simulationAngle;
xrotorAngle = handles.main.UserData.results.NoLoad.xrotorAngle;

Ke = handles.main.UserData.loss.Ke;
Density =  handles.main.UserData.loss.Density;
peddy =  (Ke/2/pi^2)*peddy*Density*angularSpeed^2/simulationAngle/xrotorAngle(2);

TotalEddy = handles.main.UserData.m.gta(handles.main.UserData.m.tzi(:,ZoneIndex)) * peddy;

if strcmp(handles.mSplit.Checked,'on')
    Nsplit = pm.Nsplit;
else
    Nsplit = 1;
end

TotalEddy = TotalEddy * 1e-9 * pm.Lst * Nsplit;

handles.main.UserData.LossData.StatorEddy = TotalEddy;

function EvalRotorEddyLoss(handles)

pm = handles.main.UserData.pm;
angularSpeed = pm.wm;

ZoneIndex = handles.main.UserData.m.mzs.rotor.zi;
bx = handles.main.UserData.results.NoLoad.Bxr;
by = handles.main.UserData.results.NoLoad.Byr;

dbxdt = [bx(:,2)-bx(:,1),...
    (bx(:,3:end)-bx(:,1:end-2))/2,...
    bx(:,end)-bx(:,end-1)];

dbydt = [by(:,2)-by(:,1),...
    (by(:,3:end)-by(:,1:end-2))/2,...
    by(:,end)-by(:,end-1)];

peddy = (dbxdt.^2+dbydt.^2);

peddy(:,1) = peddy(:,1) / 2;
peddy(:,end) = peddy(:,end) / 2;

peddy = sum(peddy,2);

simulationAngle = handles.main.UserData.results.NoLoad.simulationAngle;
xrotorAngle = handles.main.UserData.results.NoLoad.xrotorAngle;

Ke = handles.main.UserData.loss.Ke;
Density =  handles.main.UserData.loss.Density;
peddy =  (Ke/2/pi^2)*peddy*Density*angularSpeed^2/simulationAngle/xrotorAngle(2);

TotalEddy = handles.main.UserData.m.gta(handles.main.UserData.m.tzi(:,ZoneIndex)) * peddy;

if strcmp(handles.mSplit.Checked,'on')
    Nsplit = pm.Nsplit;
else
    Nsplit = 1;
end

TotalEddy = TotalEddy * 1e-9 * pm.Lst * Nsplit;

handles.main.UserData.LossData.RotorEddy = TotalEddy;

function EvalStatorHystLoss(handles)

pm = handles.main.UserData.pm;
angularSpeed = pm.wm;
% lamination data
Kh =  handles.main.UserData.loss.Kh;
alpha = handles.main.UserData.loss.alpha;
Density =  handles.main.UserData.loss.Density;
beta = 2;
% calculation of KGSE
temp = linspace(0,2*pi,1000);
KGSE = Kh/((2*pi)^(alpha-1))/(sum((abs(cos(temp))).^alpha .* ...
    (abs(sin(temp))).^(beta-alpha))*temp(2));
% zone
ZoneIndex = handles.main.UserData.m.mzs.stator.zi;
bx = handles.main.UserData.results.NoLoad.Bxs;
by = handles.main.UserData.results.NoLoad.Bys;

dbxdt = [bx(:,2)-bx(:,1),...
    (bx(:,3:end)-bx(:,1:end-2))/2,...
    bx(:,end)-bx(:,end-1)];

dbydt = [by(:,2)-by(:,1),...
    (by(:,3:end)-by(:,1:end-2))/2,...
    by(:,end)-by(:,end-1)];

phystx = (abs(dbxdt).^alpha) .* (abs(bx).^(beta-alpha));
physty = (abs(dbydt).^alpha) .* (abs(by).^(beta-alpha));

phystx(:,1) = phystx(:,1)/2;
phystx(:,end) = phystx(:,end)/2;

physty(:,1) = physty(:,1)/2;
physty(:,end) = physty(:,end)/2;

phystx = sum(phystx ,2);
physty = sum(physty ,2);
physt = phystx+physty;

xrotorAngle = handles.main.UserData.results.NoLoad.xrotorAngle;
simulationAngle = handles.main.UserData.results.NoLoad.simulationAngle;

physt = physt*Density*KGSE* angularSpeed^alpha * length(xrotorAngle)^(alpha-1)/(simulationAngle)^alpha;

TotalHyst = handles.main.UserData.m.gta(handles.main.UserData.m.tzi(:,ZoneIndex))* physt;

if strcmp(handles.mSplit.Checked,'on')
    Nsplit = pm.Nsplit;
else
    Nsplit = 1;
end

TotalHyst = TotalHyst * 1e-9 * Nsplit * pm.Lst;

handles.main.UserData.LossData.StatorHyst = TotalHyst;

function EvalRotorHystLoss(handles)
% get motor data
pm = handles.main.UserData.pm;
angularSpeed = pm.wm;
% lamination data
Kh =  handles.main.UserData.loss.Kh;
alpha = handles.main.UserData.loss.alpha;
Density =  handles.main.UserData.loss.Density;
beta = 2;
% calculation of KGSE
temp = linspace(0,2*pi,1000);
KGSE = Kh/((2*pi)^(alpha-1))/(sum((abs(cos(temp))).^alpha .* ...
    (abs(sin(temp))).^(beta-alpha))*temp(2));
% getting zone data
ZoneIndex = handles.main.UserData.m.mzs.rotor.zi;
bx = handles.main.UserData.results.NoLoad.Bxr;
by = handles.main.UserData.results.NoLoad.Byr;
% evaluation of derivative of flux density
% bx
dbxdt = [bx(:,2)-bx(:,1),...
    (bx(:,3:end)-bx(:,1:end-2))/2,...
    bx(:,end)-bx(:,end-1)];
% by
dbydt = [by(:,2)-by(:,1),...
    (by(:,3:end)-by(:,1:end-2))/2,...
    by(:,end)-by(:,end-1)];

phystx = (abs(dbxdt).^alpha) .* (abs(bx).^(beta-alpha));
physty = (abs(dbydt).^alpha) .* (abs(by).^(beta-alpha));

phystx(:,1) = phystx(:,1)/2;
phystx(:,end) = phystx(:,end)/2;

physty(:,1) = physty(:,1)/2;
physty(:,end) = physty(:,end)/2;

phystx = sum(phystx ,2);
physty = sum(physty ,2);
physt = phystx+physty;

xrotorAngle = handles.main.UserData.results.NoLoad.xrotorAngle;
simulationAngle = handles.main.UserData.results.NoLoad.simulationAngle;

physt = physt*Density*KGSE* angularSpeed^alpha * length(xrotorAngle)^(alpha-1)/(simulationAngle)^alpha;

TotalHyst = handles.main.UserData.m.gta(handles.main.UserData.m.tzi(:,ZoneIndex))* physt;

if strcmp(handles.mSplit.Checked,'on')
    Nsplit = pm.Nsplit;
else
    Nsplit = 1;
end

TotalHyst = TotalHyst * 1e-9 * Nsplit * pm.Lst;

handles.main.UserData.LossData.RotorHyst = TotalHyst;

function EvalCopperLoss(handles)

Pin = handles.main.UserData.con.Pout/...
    handles.main.UserData.con.eff;

Irms = Pin*1e3/handles.main.UserData.dim.m/...
    handles.main.UserData.con.Vrms/...
    handles.main.UserData.con.pf;

handles.main.UserData.LossData.Copper = handles.main.UserData.dim.m*Irms^2*...
   handles.main.UserData.p.Rdc; 

function EvalLosses(handles)

EvalStatorHystLoss(handles)
EvalRotorHystLoss(handles)
EvalStatorEddyLoss(handles)
EvalRotorEddyLoss(handles)
EvalCopperLoss(handles)

handles.closst.Data{1,2} = handles.main.UserData.LossData.StatorHyst;
handles.closst.Data{2,2} = handles.main.UserData.LossData.RotorHyst;
handles.closst.Data{3,2} = handles.main.UserData.LossData.StatorEddy;
handles.closst.Data{4,2} = handles.main.UserData.LossData.RotorEddy;
handles.closst.Data{5,2} = handles.main.UserData.LossData.Copper;


% --------------------------------------------------------------------
function uipushtool12_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
InitialDesign;

% --- Executes when user attempts to close main.
function main_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
ExportDataInt(handles)
delete(hObject);

function WriteData(handles)


for i = 1:size(handles.dimt.Data,1)
    handles.dimt.Data{i,2} = ...
        handles.main.UserData.pm.(handles.dimt.Data{i,1});
end


% --------------------------------------------------------------------
function Untitled_69_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_69 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ImportData(handles)

% --------------------------------------------------------------------
function Untitled_70_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_70 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ExportData(handles)

% --------------------------------------------------------------------
function dimc_Callback(hObject, eventdata, handles)
% hObject    handle to dimc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function ImportData(handles)
[FileName,PathName] = uigetfile('*.pmsm');
if ~FileName
    return
end

f = fopen([PathName,FileName],'r');
tables = {'dimt','wint','losst','simt','mst','cont'};

while true
    str = fgets(f);
    if str==-1
        break
    end
    str = strsplit(str);
    for i = 1:length(tables)
        table = tables{i};
        [check,index] = ismember(str{1},handles.(table).Data(:,1));
        if check
            handles.(table).Data{index,1} = str{1};
            handles.(table).Data{index,2} = str2double(str{2});
            break
        end
%         if i == length(tables)
%             handles.msg.ForegroundColor = 'r';	
%             handles.msg.String = '';
%             handles.msg.String{1} = 'Error:';
%             handles.msg.String{2} = 'File format is wrong ...';
%             return
%         end
    end
end

WriteMsg(handles,1,'Done!')
handles.main.UserData.states.isWinSet = false;
handles.main.UserData.states.isMGSet = false;

function ExportData(handles)

[FileName,PathName] = uiputfile('*.pmsm','File Selector',[cd,'\project.pmsm']);

if ~FileName
    return
end

f = fopen([PathName,FileName],'w');

tables = {'dimt','wint','losst','simt','mst','cont'};
for i = 1:length(tables)
    table = tables{i};
    for j = 1:size(handles.(table).Data,1)
        fprintf(f,'%s\t%f\n',handles.(table).Data{j,1},handles.(table).Data{j,2});
    end
end

fclose(f);
WriteMsg(handles,1,'Done!')


function RunInductance(handles)
WriteMsg(handles,3,'Runnig Inductance Calculation ...')
% getting windings data
if ~handles.main.UserData.states.isWinSet
    WriteMsg(handles,0,'Winding is not set ...')
    return
else
    win = handles.main.UserData.win;
end
% get mesh data
if ~handles.main.UserData.states.isMGSet
    WriteMsg(handles,0,'Geometry and Mesh is not created ...')
    return
else
    % definition of solver
    s = IHNLNRMSTL3old(handles.main.UserData.m);
end
% getting motor data
pm = handles.main.UserData.pm;
% reading tables data
ReadData(handles)
% starting simulation 
simTime = tic;
% setting solver units
% length [mm]
s.scs.l = 1e-3;
% current density to [J/mm^2]
s.scs.f = 1e6;
% set sSplit and rSplit
if strcmp(handles.mSplit.Checked,'on')
    sNsplit = pm.Ns/pm.Nsplit;
    rNsplit = pm.Nm/pm.Nsplit;
else
    sNsplit = pm.Ns;
    rNsplit = pm.Nm;
end
% vector for saving linkage fluxes
LFlux = zeros(win.m,win.m);
% phase currents
currs = eye(win.m,win.m)/pm.a;
% relative error and maximum iterations
RelErr = handles.main.UserData.sim.RelErr;
MaxIter = handles.main.UserData.sim.MaxIter;
% loop for sequence of simulations
for i = 1:win.m
    % air gap remesh
    s.m = s.m.ggmesh;
    k1 = s.m.getIndexOnCircle([0,0],pm.Rro);
    k2 = s.m.getIndexOnCircle([0,0],pm.Rro+pm.g);
    if (pm.Nsplit == 1) || (~handles.main.UserData.states.isMGSplit)
        s.m = s.m.makeAG(k1,k2);
    else
        s.m = s.m.makeAAG(k1,k2);
    end
    % setting boundary conditions
    s = s.clearallbcs;
    k0 = [s.m.getIndexOnCircle([0,0],pm.rsh)
        s.m.getIndexOnCircle([0,0],pm.Rso)];
    s = s.setdbc(k0,0);
    if rNsplit<pm.Nm
        k = s.m.getfb;
        k = setdiff(unique(k(:)),k0);
        [km,ks] = s.m.splitPeriodic(k,rNsplit*2*pi/pm.Nm);
        if rem(rNsplit,2) == 0
            s = s.setepbc(km,ks);
        else
            s = s.setopbc(km,ks);
        end
    end
    % setting excitations
    for j = 1:win.m
        Phase = win.(['Phase',num2str(j)]);
        % loop over each coils
        for k = 1:size(Phase,1)
            % phase k
            coil = win.coils(Phase(k,1),:);
            % arm1
            if coil(1)<=sNsplit
                if Phase(k,2)>0
                    s = s.setExcitation(['c1',num2str(coil(1))],currs(i,j)*pm.Ncoil,'C');
                else
                    s = s.setExcitation(['c1',num2str(coil(1))],-currs(i,j)*pm.Ncoil,'C');
                end
            end
            % arm2
            if coil(2)<=sNsplit
                if Phase(k,2)>0
                    s = s.setExcitation(['c2',num2str(coil(2))],-currs(i,j)*pm.Ncoil,'C');
                else
                    s = s.setExcitation(['c2',num2str(coil(2))],currs(i,j)*pm.Ncoil,'C');
                end
            end
        end
    end
    % runnig solver
    s.m = s.m.evalKeFeC('TL3');
    s = s.assignEdata;
    s = s.solve(RelErr,MaxIter);
    % evaluation of linkages fluxes
    % loop over phase
    for j = 1:win.m
        Phase = win.(['Phase',num2str(j)]);
        % loop over each coils
        tmpflux = 0;
        for k = 1:size(Phase,1)
            % phase k
            coil = win.coils(Phase(k,1),:);
            % arm1
            if coil(1)<=sNsplit
                if Phase(k,2)>0
                    tmpflux = tmpflux + s.evalLF(['c1',num2str(coil(1))]);
                else
                    tmpflux = tmpflux - s.evalLF(['c1',num2str(coil(1))]);
                end
            end
            % arm2
            if coil(2)<=sNsplit
                if Phase(k,2)>0
                    tmpflux = tmpflux - s.evalLF(['c2',num2str(coil(2))]);
                else
                    tmpflux = tmpflux + s.evalLF(['c2',num2str(coil(2))]);
                end
            end
        end
        LFlux(i,j) = tmpflux;
    end
    handles.msg.String{2} = ['Total Sim = ',num2str(win.m),...
        ', Sim ',num2str(i),' Compeleted ...'];
end
% evaluation of linkage fluxes
if handles.main.UserData.states.isMGSplit
    Nsplit = pm.Nsplit;
else
    Nsplit = 1;
end
LFlux = LFlux*pm.Ncoil*pm.Lst*Nsplit/1000/pm.a;
% writng runnig time
WriteMsg(handles,1,['Inductance calculation finished, Total time = ',...
    num2str(toc(simTime)),' [Sec]'])
% Inuctance calculation state
handles.main.UserData.states.isRunnedInductance = true;
% saving results
handles.main.UserData.results.Inductance.LFlux = round(LFlux*1e3,4);
% displaying results
WriteInductances(handles)

function WriteInductances(handles)

if ~handles.main.UserData.states.isRunnedInductance
    WriteMsg(handles,0,'Inductance calculation is not done ...');
    return
end

LFlux = handles.main.UserData.results.Inductance.LFlux;

index = 0;
handles.outt.Data = cell(size(LFlux,1)*size(LFlux,2),3);

for i = 1:size(LFlux,1)
    for j = 1:size(LFlux,2)
        index = index + 1;
        handles.outt.Data{index,1} = ['L',num2str(i),num2str(j)];
        handles.outt.Data{index,2} = LFlux(i,j);
        handles.outt.Data{index,3} = '[mH]';
    end
end

% --------------------------------------------------------------------
function Untitled_71_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_71 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RunInductance(handles)


% --------------------------------------------------------------------
function Untitled_72_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_72 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_73_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_73 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_74_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_74 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_75_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_75 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_76_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_76 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RunInductance(handles)

% --------------------------------------------------------------------
function Untitled_77_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_77 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_78_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_78 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_79_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_79 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function RunNoLoad(handles)

WriteMsg(handles,3,'Runnig Magnetostatic Simulations ...')

% getting windings data
if ~CheckWinExistance(handles)
    return
end
win = handles.main.UserData.win;

% get mesh data
if ~CheckMeshExistance(handles)
    return
end
% definition of solver
s = IHNLNRMSTL3old(handles.main.UserData.m);

% getting motor data
pm = handles.main.UserData.pm;
% getting needed data for loss calculation
szi = s.m.mzs.stator.zi;
rzi = s.m.mzs.rotor.zi;
Nts = sum(s.m.tzi(:,szi));
Ntr = sum(s.m.tzi(:,rzi));
% reading tables data
ReadData(handles)
% number of needed magnetostatics simulations
Nsim = handles.main.UserData.sim.Nsim;
% starting no load simulation
simTime = tic;
% setting solvere units
% length [mm]
s.scs.l = 1e-3;
% current density to [J/mm^2]
s.scs.f = 1e6;
% set sSplit and rSplit
if strcmp(handles.mSplit.Checked,'on')
    sNsplit = pm.Ns/pm.Nsplit;
    rNsplit = pm.Nm/pm.Nsplit;
else
    sNsplit = pm.Ns;
    rNsplit = pm.Nm;
end
% setting magnets magnetizations
Hc = handles.main.UserData.sim.Hc;
for i = 1:rNsplit
    if rem(i,2) == 0
        s = s.setMagnetization(['magnet',num2str(i)],...
            getRotate(msMagnetization(Hc,'r'),i*pm.Tr));
    else
        s = s.setMagnetization(['magnet',num2str(i)],...
            getRotate(msMagnetization(-Hc,'r'),i*pm.Tr));
    end
end
% setting initial rotor angle
thetao = 0;
% one pole rotation
simulationAngle = 2*pi/pm.Nm;
% rotor positions
xrotorAngle = linspace(0,simulationAngle,Nsim);
% linkage fluxes
for i = 1:win.m
    handles.main.UserData.results.NoLoad.(['lf',num2str(i)]) = zeros(Nsim,1);
end
% matrices for saving field's data
Bxs = zeros(Nts,Nsim);
Bys = zeros(Nts,Nsim);
Bxr = zeros(Ntr,Nsim);
Byr = zeros(Ntr,Nsim);
% relative error and maximum iteration
RelErr = handles.main.UserData.sim.RelErr;
MaxIter = handles.main.UserData.sim.MaxIter;
% loop for sequence of simulations
WaitFig = waitbar(0,'Running','WindowStyle','modal',...
    'Name','IHNLNRMSTL3');
for i = 1:Nsim
    % setting rotor at new positions
    thetan = xrotorAngle(i);
    s.m = s.m.rotatemz('rotor',thetan-thetao);
    for j = 1:rNsplit
        s.m = s.m.rotatemz(['rap1',num2str(j)],thetan-thetao);
        s.m = s.m.rotatemz(['rap2',num2str(j)],thetan-thetao);
        s.m = s.m.rotatemz(['magnet',num2str(j)],thetan-thetao);
    end
    thetao = thetan;
    % air gap remesh
    s.m = s.m.ggmesh;
    k1 = s.m.getIndexOnCircle([0,0],pm.Rro);
    k2 = s.m.getIndexOnCircle([0,0],pm.Rro+pm.g);
    if (pm.Nsplit == 1) || (~handles.main.UserData.states.isMGSplit)
        s.m = s.m.makeAG(k1,k2);
    else
        s.m = s.m.makeAAG(k1,k2);
    end
    % setting boundary conditions
    s = s.clearallbcs;
    k0 = [s.m.getIndexOnCircle([0,0],pm.rsh)
        s.m.getIndexOnCircle([0,0],pm.Rso)];
    s = s.setdbc(k0,0);
    if rNsplit<pm.Nm
        k = s.m.getfb;
        k = setdiff(unique(k(:)),k0);
        [km,ks] = s.m.splitPeriodic(k,rNsplit*2*pi/pm.Nm);
        if rem(rNsplit,2) == 0
            s = s.setepbc(km,ks);
        else
            s = s.setopbc(km,ks);
        end
    end
    % runnig solver
    s.m = s.m.evalKeFeC('TL3');
    s = s.assignEdata;
    s = s.solve(RelErr,MaxIter);
    % evaluation of linkages flux
    % loop over phase
    for j = 1:win.m
        Phase = win.(['Phase',num2str(j)]);
        % loop over each coils
        LFlux = 0;
        for k = 1:size(Phase,1)
            % phase k
            coil = win.coils(Phase(k,1),:);
            % arm1
            if coil(1)<=sNsplit
                if Phase(k,2)>0
                    LFlux = LFlux - s.evalLF(['c1',num2str(coil(1))]);
                else
                    LFlux = LFlux + s.evalLF(['c1',num2str(coil(1))]);
                end
            end
            % arm2
            if coil(2)<=sNsplit
                if Phase(k,2)>0
                    LFlux = LFlux + s.evalLF(['c2',num2str(coil(2))]);
                else
                    LFlux = LFlux - s.evalLF(['c2',num2str(coil(2))]);
                end
            end
        end
        handles.main.UserData.results.NoLoad.(['lf',num2str(j)])(i) = LFlux;
    end
    % getting field's data
    Bxs(:,i) = s.B(s.m.tzi(:,szi),1);
    Bys(:,i) = s.B(s.m.tzi(:,szi),2);
    Bxr(:,i) = s.B(s.m.tzi(:,rzi),1);
    Byr(:,i) = s.B(s.m.tzi(:,rzi),2);
    handles.msg.String{2} = ['Total Sim = ',num2str(Nsim),...
        ', Sim ',num2str(i),' Compeleted ...'];
    waitbar(i/Nsim,WaitFig,['Total Sim = ',num2str(Nsim),...
        ', Sim ',num2str(i),' Compeleted ...']);
end
% evalution of linkage fluxes
if handles.main.UserData.states.isMGSplit
    Nsplit = pm.Nsplit;
else
    Nsplit = 1;
end
for i = 1:win.m
    handles.main.UserData.results.NoLoad.(['lf',num2str(i)]) = ...
        handles.main.UserData.results.NoLoad.(['lf',num2str(i)]) * ...
        pm.Ncoil*pm.Lst*Nsplit/1000/pm.a;
end
% construction of half period
xi = xrotorAngle';
for i = 1:win.m
    handles.main.UserData.results.NoLoad.(['lf',num2str(i)]) = ...
        [handles.main.UserData.results.NoLoad.(['lf',num2str(i)]); ...
        -handles.main.UserData.results.NoLoad.(['lf',num2str(i)])(2:end)];
end
xi = [xi;xi(2:end)+xi(end)];
% evaluation od BEMFs
for i = 1:win.m
    handles.main.UserData.results.NoLoad.(['e',num2str(i)]) = ...
        spline(xi,handles.main.UserData.results.NoLoad.(['lf',num2str(i)]));
    handles.main.UserData.results.NoLoad.(['e',num2str(i)]).coefs = ...
        pm.wm*handles.main.UserData.results.NoLoad.(['e',num2str(i)]).coefs * ...
        diag(3:-1:1,1);
end
% saving no load data
handles.main.UserData.results.NoLoad.xi = xi;
handles.main.UserData.results.NoLoad.xrotorAngle = xrotorAngle;
handles.main.UserData.results.NoLoad.simulationAngle = simulationAngle;
handles.main.UserData.results.NoLoad.Bxs = Bxs;
handles.main.UserData.results.NoLoad.Bys = Bys;
handles.main.UserData.results.NoLoad.Bxr = Bxr;
handles.main.UserData.results.NoLoad.Byr = Byr;
handles.main.UserData.m.gta = s.m.gta;
% evaluation phase DC resistance
EvalRdc(handles)

% writing runnig time
waitbar(i/Nsim,WaitFig,'Simulation Completed.');
WriteMsg(handles,1,['No Load simulation finished, Total time = ',...
    num2str(toc(simTime)),' [Sec]'])
waitbar(i/Nsim,WaitFig,'Evaluation of Losses');
EvalLosses(handles)
close(WaitFig)
handles.main.UserData.states.isRunnedNoLoad = true;


% --------------------------------------------------------------------
function mSplit_Callback(hObject, eventdata, handles)
% hObject    handle to mSplit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.main.UserData.states.isMGSplit = true;
handles.mSplit.Checked = 'on';
handles.mTotal.Checked = 'off';

% --------------------------------------------------------------------
function mTotal_Callback(hObject, eventdata, handles)
% hObject    handle to mTotal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.main.UserData.states.isMGSplit = false;
handles.mSplit.Checked = 'off';
handles.mTotal.Checked = 'on';


% --------------------------------------------------------------------
function Untitled_80_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_80 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_81_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_81 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_82_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_82 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_83_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_83 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ShowLinkageFluxes(handles)

% --------------------------------------------------------------------
function Untitled_84_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_84 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_85_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_85 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_86_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_86 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_87_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_87 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function ShowLinkageFluxes(handles)

if ~handles.main.UserData.states.isRunnedNoLoad
    WriteMsg(handles,0,'No Load simulation is not done')
    return
end

pm = handles.main.UserData.pm;
figure('Name','Phase Linkage Flux','NumberTitle','off','WindowStyle','modal')
hold all
for i = 1:handles.main.UserData.win.m
    plot(handles.main.UserData.results.NoLoad.xi*180/pi * pm.Nm/2,...
        handles.main.UserData.results.NoLoad.(['lf',num2str(i)]),'LineWidth',1.2)
end

tmp = cell(1,handles.main.UserData.win.m);
for i = 1:handles.main.UserData.win.m
    tmp{i} = ['Phase ',num2str(i)];
end

legend(tmp)

xlabel('Electrical Angle')
ylabel('Linkage Flux [wb]')
set(gca,'XLim',[0,360],'XTick',[0,90,180,270,360],'XTickLabel',...
    {'0','\pi/2','\pi','3\pi/2','2\pi'},'FontSize',12,'GridLineStyle',':',...
    'GridColor','k','GridAlpha',1)
title('Phase Linkage Flux')
grid on
zoom on


% --------------------------------------------------------------------
function Untitled_89_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_89 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
EvalTotalLoss(handles)

% --------------------------------------------------------------------
function Untitled_90_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_90 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ironLoss = handles.closst.Data{1,2}+handles.closst.Data{2,2}+...
    handles.closst.Data{3,2}+handles.closst.Data{4,2};
copperLoss = handles.closst.Data{5,2};

WriteMsg(handles,1,['Absolute difference between Copper and Iron Loss = ',...
    num2str(abs(ironLoss-copperLoss)),' [W]'])

% --------------------------------------------------------------------
function clossc_Callback(hObject, eventdata, handles)
% hObject    handle to clossc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function EvalTotalLoss(handles)

Loss = 0;

for i = 1:7
    Loss = Loss + handles.closst.Data{i,2};
end

WriteMsg(handles,1,['Total Loss = ',num2str(Loss),' [W]'])

function EvalTotalIronLoss(handles)

Loss = 0;

for i = 1:4
    Loss = Loss + handles.closst.Data{i,2};
end

WriteMsg(handles,1,['Total Iron Loss = ',num2str(Loss),' [W]'])


function GetMeshQuality(handles)

if ~CheckMeshExistance(handles)
    return
end

[y1,y2,y3] = handles.main.UserData.m.getQuality;
WriteMsg(handles,1,['Mean = ',num2str(y1),...
    ', Min = ',num2str(y2),', Max = ',num2str(y3)])


% --------------------------------------------------------------------
function Untitled_91_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_91 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ShowBEMFs(handles)

% --------------------------------------------------------------------
function Untitled_92_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_92 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ShowLinkageFluxes(handles)

% --------------------------------------------------------------------
function Untitled_93_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_93 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_94_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_94 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_95_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_95 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RunStatic(handles)


% --------------------------------------------------------------------
function Untitled_96_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_96 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_97_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_97 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.main.UserData.states.isRunnedStatic
    handles.main.UserData.results.static.plotBmag;
else
    WriteMsg(handles,0,'Static Simulation is not done ...')
end

function RunStatic(handles)
% getting windings data
if ~handles.main.UserData.states.isWinSet
    WriteMsg(handles,0,'Winding is not set ...')
    return
else
    win = handles.main.UserData.win;
end
% get mesh data
if ~handles.main.UserData.states.isMGSet
    WriteMsg(handles,0,'Geometry and Mesh is not created ...')
    return
else
    % definition of solver
    s = IHNLNRMSTL3old(handles.main.UserData.m);
end
% reading tables data
ReadData(handles)
% dialog for getting phase current and magnetization
Hc = handles.main.UserData.sim.Hc;
prompt = cell(1,win.m+1);
defaultans = cell(1,win.m+1);
for i = 1:win.m
    prompt{i} = ['i',num2str(i),' [A]'];
    defaultans{i} = '0';
end
prompt{end} = 'Hc [A/m]';
defaultans{end} = num2str(Hc);
name = 'Run Static Simulation';
answer = inputdlg(prompt,name,[1 50],defaultans);
% check for cancel
if isempty(answer)
    return
end
% str 2 numeric
currs = zeros(1,win.m);
for i = 1:win.m
    currs(i) = str2double(answer{i});
end
Hc = str2double(answer{end});
% getting motor data
pm = handles.main.UserData.pm;
% starting simulation
WriteMsg(handles,3,'Runnig Static Simulation ...')
simTime = tic;
% setting solver units
% length [mm]
s.scs.l = 1e-3;
% current density to [J/mm^2]
s.scs.f = 1e6;
% set sSplit and rSplit
if strcmp(handles.mSplit.Checked,'on')
    sNsplit = pm.Ns/pm.Nsplit;
    rNsplit = pm.Nm/pm.Nsplit;
    s.m = s.m.rotatemz('rotor',pi/pm.Nm);
    for j = 1:rNsplit
        s.m = s.m.rotatemz(['rap1',num2str(j)],pi/pm.Nm);
        s.m = s.m.rotatemz(['rap2',num2str(j)],pi/pm.Nm);
        s.m = s.m.rotatemz(['magnet',num2str(j)],pi/pm.Nm);
    end
    s.m = s.m.ggmesh;
    s.m = s.m.setdata;
else
    sNsplit = pm.Ns;
    rNsplit = pm.Nm;
end
% setting magnets magnetizations
for i = 1:rNsplit
    if rem(i,2) == 0
        s = s.setMagnetization(['magnet',num2str(i)],...
            getRotate(msMagnetization(Hc,'r'),i*pm.Tr));
    else
        s = s.setMagnetization(['magnet',num2str(i)],...
            getRotate(msMagnetization(-Hc,'r'),i*pm.Tr));
    end
end
% relative error and maximum iterations
RelErr = handles.main.UserData.sim.RelErr;
MaxIter = handles.main.UserData.sim.MaxIter;
% air gap remesh
k1 = s.m.getIndexOnCircle([0,0],pm.Rro);
k2 = s.m.getIndexOnCircle([0,0],pm.Rro+pm.g);
if (pm.Nsplit == 1) || (~handles.main.UserData.states.isMGSplit)
    s.m = s.m.makeAG(k1,k2);
else
    s.m = s.m.makeAAG(k1,k2);
end
% setting boundary conditions
s = s.clearallbcs;
k0 = [s.m.getIndexOnCircle([0,0],pm.rsh)
    s.m.getIndexOnCircle([0,0],pm.Rso)];
s = s.setdbc(k0,0);
if rNsplit<pm.Nm
    k = s.m.getfb;
    k = setdiff(unique(k(:)),k0);
    [km,ks] = s.m.splitPeriodic(k,rNsplit*2*pi/pm.Nm);
    if rem(rNsplit,2) == 0
        s = s.setepbc(km,ks);
    else
        s = s.setopbc(km,ks);
    end
end
% setting excitations
for j = 1:win.m
    Phase = win.(['Phase',num2str(j)]);
    % loop over each coils
    for k = 1:size(Phase,1)
        % phase k
        coil = win.coils(Phase(k,1),:);
        % arm1
        if coil(1)<=sNsplit
            if Phase(k,2)>0
                s = s.setExcitation(['c1',num2str(coil(1))],currs(j)*pm.Ncoil,'C');
            else
                s = s.setExcitation(['c1',num2str(coil(1))],-currs(j)*pm.Ncoil,'C');
            end
        end
        % arm2
        if coil(2)<=sNsplit
            if Phase(k,2)>0
                s = s.setExcitation(['c2',num2str(coil(2))],-currs(j)*pm.Ncoil,'C');
            else
                s = s.setExcitation(['c2',num2str(coil(2))],currs(j)*pm.Ncoil,'C');
            end
        end
    end
end
% runnig solver
s.m = s.m.evalKeFeC('TL3');
s = s.assignEdata;
s = s.solve(RelErr,MaxIter);
% show simulation run time
WriteMsg(handles,1,['Static simulation finished, Total time = ',...
    num2str(toc(simTime)),' [Sec]'])
handles.main.UserData.states.isRunnedStatic = true;
handles.main.UserData.results.static = s;


% --------------------------------------------------------------------
function Untitled_98_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_98 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.main.UserData.states.isRunnedStatic
    handles.main.UserData.results.static.plotBvec;
else
    WriteMsg(handles,0,'Static Simulation is not done ...')
end


% --------------------------------------------------------------------
function Untitled_99_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_99 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.main.UserData.states.isRunnedStatic
    handles.main.UserData.results.static.plotAmag;
else
    WriteMsg(handles,0,'Static Simulation is not done ...')
end


% --------------------------------------------------------------------
function Untitled_100_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_100 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_101_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_101 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_102_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_102 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GetBEMFFFT(handles)

function GetBEMFFFT(handles)

if ~handles.main.UserData.states.isRunnedNoLoad
    WriteMsg(handles,0,'No Load simulation is not done')
end

ea = ppval(handles.main.UserData.results.NoLoad.ea,handles.main.UserData.results.NoLoad.xi);
eb = ppval(handles.main.UserData.results.NoLoad.eb,handles.main.UserData.results.NoLoad.xi);
eab = ea-eb;

ea = abs(fft(ea));
eab = abs(fft(eab));


figure('Name','Harmonic Content of Phase to Phase and Line to Line BEMFs','NumberTitle','off',...
    'WindowStyle','modal')

subplot(211)
N = floor((length(ea)-1)/2);
ea = ea(2:N);
stem(2:N-1,ea(2:end)*100/ea(1),'LineWidth',1.2)
ylabel('|E_n|/|E_1|,%')
set(gca,'XLim',[1,N],'FontSize',12)
zoom on
title('Phase BEMF Harmonic Amplitude Relative to Fundamental');


% --------------------------------------------------------------------
function Untitled_103_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_103 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function slot4_Callback(hObject, eventdata, handles)
% hObject    handle to slot4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function slot5_Callback(hObject, eventdata, handles)
% hObject    handle to slot5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function ExportDataInt(handles)


f = fopen([cd,'\Int.pmsm'],'w');

tables = {'dimt','wint','losst','simt','mst','cont'};
for i = 1:length(tables)
    table = tables{i};
    for j = 1:size(handles.(table).Data,1)
        fprintf(f,'%s\t%f\n',handles.(table).Data{j,1},handles.(table).Data{j,2});
    end
end

fclose(f);

function ImportDataInt(handles)

f = fopen([cd,'\Int.pmsm'],'r');

if f<0
    return
end

tables = {'dimt','wint','losst','simt','mst','cont'};
while true
    str = fgets(f);
    if str==-1
        break
    end
    str = strsplit(str);
    for i = 1:length(tables)
        table = tables{i};
        [check,index] = ismember(str{1},handles.(table).Data(:,1));
        if check
            handles.(table).Data{index,1} = str{1};
            handles.(table).Data{index,2} = str2double(str{2});
            break
        end
    end
end
fclose(f);


% --------------------------------------------------------------------
function Untitled_104_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_104 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% get mesh data
EvalRdc(handles)

function EvalRdc(handles)

% getting windings data
if ~handles.main.UserData.states.isWinSet
    WriteMsg(handles,0,'Winding is not set ...')
    return
end

% getting mesh data
if ~handles.main.UserData.states.isMGSet
    WriteMsg(handles,0,'Geometry and Mesh is not created ...')
    return
end

handles.MsgBox(handles)
pm = handles.main.UserData.pm;

handles.main.UserData.p.Rdc = pm.getRph(...
    handles.main.UserData.m.g.ds.c11.value.area*2,...
    handles.main.UserData.win.S);

handles.outt.Data = cell(2,3);
handles.outt.Data{1,1} = 'Rdc,25';
handles.outt.Data{1,2} = handles.main.UserData.p.Rdc;
handles.outt.Data{1,3} = '[Ohm]';

handles.main.UserData.p.Rdc = handles.main.UserData.p.Rdc*(1+pm.alphaT*(pm.winT-25));

handles.outt.Data{2,1} = ['Rdc,',num2str(pm.winT)];
handles.outt.Data{2,2} = handles.main.UserData.p.Rdc;
handles.outt.Data{2,3} = '[Ohm]';

handles.main.UserData.p.Rdc = handles.main.UserData.p.Rdc*pm.getKR;

handles.outt.Data{3,1} = ['Rac,',num2str(pm.winT)];
handles.outt.Data{3,2} = handles.main.UserData.p.Rdc;
handles.outt.Data{3,3} = '[Ohm]';

% --------------------------------------------------------------------
function Untitled_105_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_105 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
EvalRdc(handles)


% --------------------------------------------------------------------
function uipushtool13_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SetWinding(hObject, eventdata, handles)


% --------------------------------------------------------------------
function Untitled_106_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_106 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.main.UserData.states.isRunnedStatic
    handles.main.UserData.results.static.plotBmagSmooth;
else
    WriteMsg(handles,0,'Static Simulation is not done ...')
end


% --------------------------------------------------------------------
function Untitled_107_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_107 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.MsgBox(handles);
pm = handles.main.UserData.pm;
if handles.main.UserData.states.isRunnedStatic
    handles.main.UserData.results.static.plotBrBtOnCircle(pm.Rro+pm.g/2);
else
    WriteMsg(handles,0,'Static Simulation is not done ...')
end


% --------------------------------------------------------------------
function Untitled_108_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_108 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_109_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_109 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_110_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_110 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% displaying results
WriteInductances(handles)


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --------------------------------------------------------------------
function Untitled_111_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_111 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_112_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_112 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_113_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_113 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_114_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_114 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_115_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_115 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function rLam_Callback(hObject, eventdata, handles)
% hObject    handle to rLam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.main.UserData.states.isRotorLam = true;
handles.rLam.Checked = 'on';
handles.rSolid.Checked = 'off';

% --------------------------------------------------------------------
function rSolid_Callback(hObject, eventdata, handles)
% hObject    handle to rSolid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.main.UserData.states.isRotorLam = false;
handles.rLam.Checked = 'off';
handles.rSolid.Checked = 'on';

% --------------------------------------------------------------------
function Untitled_118_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_118 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function sLam_Callback(hObject, eventdata, handles)
% hObject    handle to sLam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.main.UserData.states.isStatorLam = true;
handles.sLam.Checked = 'on';
handles.sSolid.Checked = 'off';

% --------------------------------------------------------------------
function sSolid_Callback(hObject, eventdata, handles)
% hObject    handle to sSolid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.main.UserData.states.isStatorLam = false;
handles.sLam.Checked = 'off';
handles.sSolid.Checked = 'on';

% --------------------------------------------------------------------
function Untitled_121_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_121 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function munseg_Callback(hObject, eventdata, handles)
% hObject    handle to munseg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mseg_Callback(hObject, eventdata, handles)
% hObject    handle to mseg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_122_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_122 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_123_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_123 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_124_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_124 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ShowHyst(handles)

% --------------------------------------------------------------------
function Untitled_125_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_125 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ShowRotorHyst(handles)

% --------------------------------------------------------------------
function Untitled_126_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_126 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ShowEddy(handles)

% --------------------------------------------------------------------
function Untitled_127_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_127 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ShowRotorEddy(handles)

% --------------------------------------------------------------------
function Untitled_130_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_130 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ShowEddy(handles)

% --------------------------------------------------------------------
function Untitled_131_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_131 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ShowRotorEddy(handles)

% --------------------------------------------------------------------
function Untitled_128_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_128 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ShowHyst(handles)

% --------------------------------------------------------------------
function Untitled_129_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_129 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ShowRotorHyst(handles)


% --------------------------------------------------------------------
function Untitled_132_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_132 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% get mesh data
EvalBWD(handles)


% --------------------------------------------------------------------
function Untitled_133_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_133 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Evalnd(handles)


% --------------------------------------------------------------------
function Untitled_134_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_134 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
EvalBWD(handles)

% --------------------------------------------------------------------
function Untitled_135_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_135 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Evalnd(handles)


% --------------------------------------------------------------------
function Untitled_136_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_136 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
EvalTotalIronLoss(handles)


function isExist = CheckMeshExistance(handles)

if handles.main.UserData.states.isMGSet
    isExist = true;
    return
else
    WriteMsg(handles,0,'Geometry and Mesh are not created ...')
    isExist = false;
end

function isExist = CheckWinExistance(handles)

if handles.main.UserData.states.isWinSet
    isExist = true;
    return
else
    WriteMsg(handles,0,'Winding is not set ...')
    isExist = false;
end

function EvalBWD(handles)

if ~CheckMeshExistance(handles)
    return
end

pm = handles.main.UserData.pm;
As = handles.main.UserData.m.g.ds.c11.value.area;
Awb = pm.Kf*As/pm.Ncoil/2;
dwb = sqrt(4*Awb/pi);

WriteMsg(handles,1,['Bared Wire Diameter [mm] =',num2str(dwb)]);

for i = 1:size(handles.wint.Data,1)
    if strcmpi(handles.wint.Data{i,1},'dwb')
        handles.wint.Data{i,2} = dwb;
        break
    end
end

function Evalnd(handles)

handles.MsgBox(handles);
pm = handles.main.UserData.pm;

nd = pm.hss3/pm.dwb;
nd = floor(nd)-1;

WriteMsg(handles,1,['nd = ',num2str(nd)])

for i = 1:size(handles.wint.Data,1)
    if strcmpi(handles.wint.Data{i,1},'nd')
        handles.wint.Data{i,2} = nd;
        break
    end
end
