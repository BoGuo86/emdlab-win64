function varargout = PMSM_MOTOR(varargin)
% PMSM_MOTOR MATLAB code for PMSM_MOTOR.fig
%      PMSM_MOTOR, by itself, creates a new PMSM_MOTOR or raises the existing
%      singleton*.
%
%      H = PMSM_MOTOR returns the handle to a new PMSM_MOTOR or the handle to
%      the existing singleton*.
%
%      PMSM_MOTOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PMSM_MOTOR.CGM with the given input arguments.
%
%      PMSM_MOTOR('Property','Value',...) creates a new PMSM_MOTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PMSM_MOTOR_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PMSM_MOTOR_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE'RunSimulation Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PMSM_MOTOR

% Last Modified by GUIDE v2.5 27-Apr-2017 00:36:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PMSM_MOTOR_OpeningFcn, ...
                   'gui_OutputFcn',  @PMSM_MOTOR_OutputFcn, ...
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


% --- Executes just before PMSM_MOTOR is made visible.
function PMSM_MOTOR_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PMSM_MOTOR (see VARARGIN)

% Choose default command line output for PMSM_MOTOR
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

handles.EX.SelectedObject = handles.NL;
% UIWAIT makes PMSM_MOTOR wait for user response (see UIRESUME)
% uiwait(handles.main);


% --- Outputs from this function are returned to the command line.
function varargout = PMSM_MOTOR_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function Ns_Callback(hObject, eventdata, handles)
% hObject    handle to Ns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ns as text
%        str2double(get(hObject,'String')) returns contents of Ns as a double

handles.cgm.Enable = 'off';
OffButtuns(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function Ns_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Nm_Callback(hObject, eventdata, handles)
% hObject    handle to Nm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Nm as text
%        str2double(get(hObject,'String')) returns contents of Nm as a double

handles.cgm.Enable = 'off';
OffButtuns(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function Nm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Nm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Rso_Callback(hObject, eventdata, handles)
% hObject    handle to Rso (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Rso as text
%        str2double(get(hObject,'String')) returns contents of Rso as a double
OffButtuns(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function Rso_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rso (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Rro_Callback(hObject, eventdata, handles)
% hObject    handle to Rro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Rro as text
%        str2double(get(hObject,'String')) returns contents of Rro as a double
OffButtuns(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function Rro_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wsy_Callback(hObject, eventdata, handles)
% hObject    handle to wsy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wsy as text
%        str2double(get(hObject,'String')) returns contents of wsy as a double
OffButtuns(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function wsy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wsy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wry_Callback(hObject, eventdata, handles)
% hObject    handle to wry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wry as text
%        str2double(get(hObject,'String')) returns contents of wry as a double
OffButtuns(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function wry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wtb_Callback(hObject, eventdata, handles)
% hObject    handle to wtb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wtb as text
%        str2double(get(hObject,'String')) returns contents of wtb as a double
OffButtuns(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function wtb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wtb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hm_Callback(hObject, eventdata, handles)
% hObject    handle to hm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hm as text
%        str2double(get(hObject,'String')) returns contents of hm as a double
OffButtuns(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function hm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function g_Callback(hObject, eventdata, handles)
% hObject    handle to g (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of g as text
%        str2double(get(hObject,'String')) returns contents of g as a double
OffButtuns(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function g_CreateFcn(hObject, eventdata, handles)
% hObject    handle to g (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cgm.
function cgm_Callback(hObject, eventdata, handles)

% hObject    handle to cgm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Msg.ForegroundColor = 'r';
handles.Msg.String = 'Busy: Creation of Geometry And Mesh Generation ...';
materialdir = [cd,'\MaterialsData'];

% read data
ReadingData(hObject, eventdata, handles)
pm = handles.main.UserData.pm;

% maximum gap degree
thetag = handles.main.UserData.p.thetag;
angle1 = handles.main.UserData.p.angle1;
angle2 = handles.main.UserData.p.angle2;

g = GDBC2D;
%% creation of stator lamination
p1 = [pm.Rro+pm.g,0];
p2 = [pm.Rso,0];
p3 = pmove(p2,'theta',pm.Ts/2);
p4 = pmove(p3,'r',-pm.wsy);
p5 = pmove(p4,'theta',-pm.Tso2/2);
p6 = pmove(p5,'x',-pm.hss3);
p8 = pmove(p1,'theta',(pm.Ts-pm.Tso1)/2);
p7 = pmove(p8,'r',pm.hss1);
tmp1 = 2*(pm.Rro+pm.g)*sin(thetag*pi/180/2);
tmp3 = angle2;
tmp2 = 2*(pm.Rso)*sin(tmp3*pi/180/2);
tmp5 = angle2;
tmp4 = 2*(pm.Rso-pm.wsy)*sin(tmp5*pi/180/2);
[g,L1] = g.newdlinewdkps(p1,p2,'l1',tmp1,'l2',tmp2);
[g,A1] = g.newdarccppwdkps([0,0],p2,p3,'maxDegree',tmp3);
[g,L2] = g.newdlinewdkps(p3,p4,'l1',tmp2,'l2',tmp4);
[g,A2] = g.newdarccppwdkps([0,0],p4,p5,'direction',-1,'maxDegree',tmp5);
[g,L3] = g.newdlinewdkps(p6,p5,'l1',tmp1,'l2',tmp4);
[g,L4] = g.newdlinewdkps(p6,p7,'maxLength',tmp1);
[g,L5] = g.newdlinewdkps(p7,p8,'maxLength',tmp1);
[g,A3] = g.newdarccppwdkps([0,0],p8,p1,'direction',-1,'maxDegree',thetag);
g = g.newcb('s1',L1,1,A1,1,L2,1,A2,1,L3,-1,L4,1,L5,1,A3,1);
%% cration of coil
p9 = pmove(p1,'theta',pm.Ts/2);
p9 = pmove(p9,'r',pm.hss1);
g = g.newdlinewdkps(p7,p9,'maxLength',tmp1);
g = g.newdlinewdkps(p4,p9,'l1',tmp4,'l2',tmp1);
g = g.newcb('c11','L3',1,'A2',-1,'L7',1,'L6',-1,'L4',-1);
%% creation of slot air
g = g.cmirrorbd('L5',[cos(pm.Ts/2),sin(pm.Ts/2)]);
g = g.cmirrorbd('L6',[cos(pm.Ts/2),sin(pm.Ts/2)]);
g = g.newdarccpp('kp8','kp11',[0,0],'maxDegree',thetag);
g = g.newcb('sap1','L5',-1,'L6',1,'L9',-1,'L8',1,'A4',-1);
%% magnet
p1 = [pm.rsh+pm.wry,0];
p2 = p1+[pm.hm,0];
p3 = pmove(p2,'theta',pm.beta_m*pi/pm.Nm);
p4 = pmove(p3,'r',-pm.hm);
p5 = pmove(p4,'theta',(1-pm.beta_m)*pi/pm.Nm);
p7 = [pm.rsh,0];
p6 = pmove(p7,'theta',pi/pm.Nm);
g = g.newdlinewdkps(p1,p2,'maxLength',tmp1);
g = g.newdarccppwdkps([0,0],p2,p3,'maxDegree',thetag);
g = g.newdlinewdkps(p4,p3,'maxLength',tmp1);
g = g.newdarccppwdkps([0,0],p4,p1,'direction',-1,'maxDegree',thetag);
g = g.newcb('m1','L10',1,'A5',1,'L11',-1,'A6',1);
%% rotor lamination
tmp2 = 2*(pm.rsh)*sin(angle1*pi/180/2);
g = g.newdarccppwdkps([0,0],p4,p5,'maxDegree',thetag);
g = g.newdlinewdkps(p5,p6,'l1',tmp1,'l2',tmp2);
g = g.newdarccppwdkps([0,0],p6,p7,'direction',-1,'maxDegree',angle1);
g = g.newdlinewdkps(p1,p7,'l1',tmp1,'l2',tmp2);
g = g.newcb('r1','L13',-1,'A6',-1,'A7',1,'L12',1,'A8',1);
%% rotor air pocket
g = g.newdlinewdkps(p5,pmove(p2,'theta',pi/pm.Nm),'maxLength',tmp1);
g = g.newdarccpp('kp14','kp19',[0,0],'maxDegree',thetag);
g = g.newcb('rap1','L11',1,'A9',1,'L14',-1,'A7',-1);
%% creation of domains
g = g.newdDM('s1','s1');
close(gcf);
g = g.newdDM('r1','r1');
close(gcf);
g = g.newdDM('rap11','rap1');
close(gcf);
g = g.newdDM('m1','m1');
close(gcf);
g = g.newdDM('sap1','sap1');
close(gcf);
g = g.newdDM('c11','c11');
close(gcf);
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
for i = 1:2:2*(pm.Nsplit*pm.Ns/pm.Nm-1)
    m = m.crotatemz(['s',num2str(i+2)],['s',num2str(i)],pm.Ts);
    m = m.crotatemz(['s',num2str(i+3)],['s',num2str(i+1)],pm.Ts);
end
m = m.cmirrormz('r2','r1',[1,0]);
m = m.cmirrormz('m2','m1',[1,0]);
m = m.joinmzs('magnet1','m1','m2');
for i = 1:(pm.Nsplit*pm.Ns/pm.Nm-1)
    m = m.crotatemz(['c1',num2str(i+1)],['c1',num2str(i)],pm.Ts);
    m = m.crotatemz(['c2',num2str(i+1)],['c2',num2str(i)],pm.Ts);
    m = m.crotatemz(['sap',num2str(i+1)],['sap',num2str(i)],pm.Ts);
end
for i = 1:(pm.Nsplit-1)
    m = m.crotatemz(['magnet',num2str(i+1)],['magnet',num2str(i)],pm.Tr);
    m = m.crotatemz(['rap1',num2str(i+1)],['rap1',num2str(i)],pm.Tr);
    m = m.crotatemz(['rap2',num2str(i+1)],['rap2',num2str(i)],pm.Tr);
end
for i = 1:2:2*(pm.Nsplit-1)
    m = m.crotatemz(['r',num2str(i+2)],['r',num2str(i)],pm.Tr);
    m = m.crotatemz(['r',num2str(i+3)],['r',num2str(i+1)],pm.Tr);
end
% joining stator mesh zones
temp = getlist('s',1:2*pm.Nsplit*pm.Ns/pm.Nm);
m = m.joinmzs('stator',temp{:});
% joining rotor mesh zones
temp = getlist('r',1:2*pm.Nsplit);
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

handles.main.UserData.m = m;
handles.main.UserData.pm = pm;
handles.Msg.ForegroundColor = 'b';
handles.Msg.String = 'Creation of Geometry And Mesh Generation Compeleted !';

handles.Showmzs.Enable = 'on';
handles.Showgm.Enable = 'on';
handles.RunSimulation.Enable = 'on';
handles.MQ.Enable = 'on';

% --- Executes on button press in Showmzs.
function Showmzs_Callback(hObject, eventdata, handles)
% hObject    handle to Showmzs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.main.UserData.m.showmzs;



function gmsize_Callback(hObject, eventdata, handles)
% hObject    handle to gmsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gmsize as text
%        str2double(get(hObject,'String')) returns contents of gmsize as a double


% --- Executes during object creation, after setting all properties.
function gmsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gmsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function thetag_Callback(hObject, eventdata, handles)
% hObject    handle to thetag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thetag as text
%        str2double(get(hObject,'String')) returns contents of thetag as a double


% --- Executes during object creation, after setting all properties.
function thetag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thetag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function angle1_Callback(hObject, eventdata, handles)
% hObject    handle to angle1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of angle1 as text
%        str2double(get(hObject,'String')) returns contents of angle1 as a double


% --- Executes during object creation, after setting all properties.
function angle1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angle1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function angle2_Callback(hObject, eventdata, handles)
% hObject    handle to angle2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of angle2 as text
%        str2double(get(hObject,'String')) returns contents of angle2 as a double


% --- Executes during object creation, after setting all properties.
function angle2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angle2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Nsim_Callback(hObject, eventdata, handles)
% hObject    handle to Nsim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Nsim as text
%        str2double(get(hObject,'String')) returns contents of Nsim as a double


% --- Executes during object creation, after setting all properties.
function Nsim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Nsim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Showgm.
function Showgm_Callback(hObject, eventdata, handles)
% hObject    handle to Showgm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ff = handles.main.UserData.m.showmesh;
ff.Visible = 'on';

% --- Executes on button press in RunSimulation.
function RunSimulation_Callback(hObject, eventdata, handles)
% hObject    handle to RunSimulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Msg.ForegroundColor = 'r';
handles.Msg.String = 'Busy: Runnig Magnetostatic Simulations ...';

pm = handles.main.UserData.pm;
m = handles.main.UserData.m;

szi = m.mzs.stator.zi;
rzi = m.mzs.rotor.zi;
Nts = sum(m.tzi(:,szi));
Ntr = sum(m.tzi(:,rzi));

% number of needed magnetostatics simulations
Nsim = str2double(handles.Nsim.String);

s = IHNLNRMSTL3(m);clear m;
% setting units
% length [mm]
s.scs.l = 1e-3;
% current density to [J/mm^2]
s.scs.f = 1e6;
% getting windings data
coils = handles.main.UserData.win.coils;
phaseA = handles.main.UserData.win.phaseA;
phaseB = handles.main.UserData.win.phaseB;
phaseC = handles.main.UserData.win.phaseC;

if isempty(coils)
    handles.Msg.String = 'Error: Winiding configurations does not set.';
    return
end

% setting magnets magnetizations
Hc = str2double(handles.Hc.String);
for i = 1:pm.Nsplit
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
%% solver setting
% one pole rotation
simulationAngle = 2*pi/pm.Nm;

% rotor positions
xrotorAngle = linspace(0,simulationAngle,Nsim);
Ntheta = length(xrotorAngle);
F(Ntheta) = struct('cdata',[],'colormap',[]);
M(Ntheta) = struct('cdata',[],'colormap',[]);

aLFlux = zeros(Ntheta,1);
bLFlux = zeros(Ntheta,1);
cLFlux = zeros(Ntheta,1);

Bxs = zeros(Nts,Ntheta);
Bys = zeros(Nts,Ntheta);
Bxr = zeros(Ntr,Ntheta);
Byr = zeros(Ntr,Ntheta);

switch handles.EX.SelectedObject.Tag
    case 'NL'
        ia = zeros(Ntheta,1);
        ib = zeros(Ntheta,1);
        ic = zeros(Ntheta,1);
    case 'FL'
        ia = ones(Ntheta,1);
        ib = zeros(Ntheta,1);
        ic = zeros(Ntheta,1);
    case 'HA' 
        ia = zeros(Ntheta,1);
        ib = zeros(Ntheta,1);
        ic = zeros(Ntheta,1);
end

% loop for sequence of simulations
for i = 1:Ntheta
    % setting rotor at new positions
    thetan = xrotorAngle(i);
    s.m = s.m.rotatemz('rotor',thetan-thetao);
    for j = 1:pm.Nsplit
        s.m = s.m.rotatemz(['rap1',num2str(j)],thetan-thetao);
        s.m = s.m.rotatemz(['rap2',num2str(j)],thetan-thetao);
        s.m = s.m.rotatemz(['magnet',num2str(j)],thetan-thetao);
    end
    thetao = thetan;
    % air gap remesh
    s.m = s.m.ggmesh;
    k1 = s.m.getIndexOnCircle([0,0],pm.Rro);
    k2 = s.m.getIndexOnCircle([0,0],pm.Rro+pm.g);
    if pm.Nsplit<pm.Nm
        s.m = s.m.makeAAG(k1,k2);
    else
        s.m = s.m.makeAG(k1,k2,4);
    end
    % setting boundary conditions
    s = s.clearallbcs;
    k0 = [s.m.getIndexOnCircle([0,0],pm.rsh)
        s.m.getIndexOnCircle([0,0],pm.Rso)];
    s = s.setdbc(k0,0);
    if pm.Nsplit<pm.Nm
        k = s.m.getfb;
        k = setdiff(unique(k(:)),k0);
        [km,ks] = s.m.splitPeriodic(k,pm.Nsplit*2*pi/pm.Nm);
        if rem(pm.Nsplit,2) == 0
            s = s.setepbc(km,ks);
        else
            s = s.setopbc(km,ks);
        end
    end
    % saving rotor mesh plot
    ff = s.m.showmesh;
    M(i) = getframe(ff);
    close(ff)
    
    % setting stator winding currents excitations
    iai = interp1(xrotorAngle,ia,thetan);
    ibi = interp1(xrotorAngle,ib,thetan);
    ici = interp1(xrotorAngle,ic,thetan);
    
    for j = 1:size(phaseA,1)
        % phase A
        coil = coils(phaseA(j,1),:);
        % arm1
        if coil(1)<=pm.Ns*pm.Nsplit/pm.Nm
            if phaseA(j,2)>0
                s = s.setExcitation(['c1',num2str(coil(1))],iai*pm.Ncoil,'C');
            else
                s = s.setExcitation(['c1',num2str(coil(1))],-iai*pm.Ncoil,'C');
            end
        end
        % arm2
        if coil(2)<=pm.Ns*pm.Nsplit/pm.Nm
            if phaseA(j,2)>0
                s = s.setExcitation(['c2',num2str(coil(2))],-iai*pm.Ncoil,'C');
            else
                s = s.setExcitation(['c2',num2str(coil(2))],iai*pm.Ncoil,'C');
            end
        end
        % phase B
        coil = coils(phaseB(j,1),:);
        % arm1
        if coil(1)<=pm.Ns*pm.Nsplit/pm.Nm
            if phaseB(j,2)>0
                s = s.setExcitation(['c1',num2str(coil(1))],ibi*pm.Ncoil,'C');
            else
                s = s.setExcitation(['c1',num2str(coil(1))],-ibi*pm.Ncoil,'C');
            end
        end
        % arm2
        if coil(2)<=pm.Ns*pm.Nsplit/pm.Nm
            if phaseB(j,2)>0
                s = s.setExcitation(['c2',num2str(coil(2))],-ibi*pm.Ncoil,'C');
            else
                s = s.setExcitation(['c2',num2str(coil(2))],ibi*pm.Ncoil,'C');
            end
        end
        % phase C
        coil = coils(phaseC(j,1),:);
        % arm1
        if coil(1)<=pm.Ns*pm.Nsplit/pm.Nm
            if phaseC(j,2)>0
                s = s.setExcitation(['c1',num2str(coil(1))],ici*pm.Ncoil,'C');
            else
                s = s.setExcitation(['c1',num2str(coil(1))],-ici*pm.Ncoil,'C');
            end
        end
        % arm2
        if coil(2)<=pm.Ns*pm.Nsplit/pm.Nm
            if phaseC(j,2)>0
                s = s.setExcitation(['c2',num2str(coil(2))],-ici*pm.Ncoil,'C');
            else
                s = s.setExcitation(['c2',num2str(coil(2))],ici*pm.Ncoil,'C');
            end
        end
    end
    
    % runnig solver
    s.m = s.m.evalKeFeC('TL3');
    s = s.assignEdata;
    s = s.solve(1e-6,20);
    % evaluation of linkages flux
    for j = 1:size(phaseA,1)
        % phase A
        coil = coils(phaseA(j,1),:);
        % arm1
        if coil(1)<=pm.Ns*pm.Nsplit/pm.Nm
            if phaseA(j,2)>0
                aLFlux(i) = aLFlux(i) - s.evalLF(['c1',num2str(coil(1))]);
            else
                aLFlux(i) = aLFlux(i) + s.evalLF(['c1',num2str(coil(1))]);
            end
        end
        % arm2
        if coil(2)<=pm.Ns*pm.Nsplit/pm.Nm
            if phaseA(j,2)>0
                aLFlux(i) = aLFlux(i) + s.evalLF(['c2',num2str(coil(2))]);
            else
                aLFlux(i) = aLFlux(i) - s.evalLF(['c2',num2str(coil(2))]);
            end
        end
        % phase B
        coil = coils(phaseB(j,1),:);
        % arm1
        if coil(1)<=pm.Ns*pm.Nsplit/pm.Nm
            if phaseB(j,2)>0
                bLFlux(i) = bLFlux(i) - s.evalLF(['c1',num2str(coil(1))]);
            else
                bLFlux(i) = bLFlux(i) + s.evalLF(['c1',num2str(coil(1))]);
            end
        end
        % arm2
        if coil(2)<=pm.Ns*pm.Nsplit/pm.Nm
            if phaseB(j,2)>0
                bLFlux(i) = bLFlux(i) + s.evalLF(['c2',num2str(coil(2))]);
            else
                bLFlux(i) = bLFlux(i) - s.evalLF(['c2',num2str(coil(2))]);
            end
        end
        % phase C
        coil = coils(phaseC(j,1),:);
        % arm1
        if coil(1)<=pm.Ns*pm.Nsplit/pm.Nm
            if phaseC(j,2)>0
                cLFlux(i) = cLFlux(i) - s.evalLF(['c1',num2str(coil(1))]);
            else
                cLFlux(i) = cLFlux(i) + s.evalLF(['c1',num2str(coil(1))]);
            end
        end
        % arm2
        if coil(2)<=pm.Ns*pm.Nsplit/pm.Nm
            if phaseC(j,2)>0
                cLFlux(i) = cLFlux(i) + s.evalLF(['c2',num2str(coil(2))]);
            else
                cLFlux(i) = cLFlux(i) - s.evalLF(['c2',num2str(coil(2))]);
            end
        end
    end

    % get Bmag frame
    a = getlist('magnet',1:pm.Nsplit);
    
    ff = s.plotBmagw('stator','rotor',a{:});
    F(i) = getframe(ff);
    close(ff)

    Bxs(:,i) = s.B(s.m.tzi(:,szi),1);
    Bys(:,i) = s.B(s.m.tzi(:,szi),2);
    Bxr(:,i) = s.B(s.m.tzi(:,rzi),1);
    Byr(:,i) = s.B(s.m.tzi(:,rzi),2);
    handles.Msg.String = ['Busy: Total Sim = ',num2str(Ntheta),...
        ', Sim ',num2str(i),' Compeleted ...'];
end


aLFlux = aLFlux*pm.Ncoil*pm.Lst*pm.Nm/1000/pm.Nsplit;
bLFlux = bLFlux*pm.Ncoil*pm.Lst*pm.Nm/1000/pm.Nsplit;
cLFlux = cLFlux*pm.Ncoil*pm.Lst*pm.Nm/1000/pm.Nsplit;
xi = xrotorAngle';

aLFlux = [aLFlux;-aLFlux(2:end)];
bLFlux = [bLFlux;-bLFlux(2:end)];
cLFlux = [cLFlux;-cLFlux(2:end)];
xi = [xi;xi(2:end)+xi(end)];

ea = spline(xi,aLFlux);
ea.coefs = pm.wm*ea.coefs*diag(3:-1:1,1);
eb = spline(xi,bLFlux);
eb.coefs = pm.wm*eb.coefs*diag(3:-1:1,1);
ec = spline(xi,cLFlux);
ec.coefs = pm.wm*ec.coefs*diag(3:-1:1,1);

handles.main.UserData.aLFlux = aLFlux;
handles.main.UserData.bLFlux = bLFlux;
handles.main.UserData.cLFlux = cLFlux;
handles.main.UserData.ea = ea;
handles.main.UserData.eb = eb;
handles.main.UserData.ec = ec;
handles.main.UserData.xi = xi;


handles.main.UserData.p.xrotorAngle = xrotorAngle;
handles.main.UserData.p.simulationAngle = simulationAngle;
handles.main.UserData.p.Nsim = Nsim;
handles.main.UserData.mF = M;
handles.main.UserData.bF = F;
handles.main.UserData.Bxs = Bxs;
handles.main.UserData.Bys = Bys;
handles.main.UserData.Bxr = Bxr;
handles.main.UserData.Byr = Byr;
handles.main.UserData.m.gta = s.m.gta;

handles.main.UserData.p.Rdc = pm.getRph(s.m.mzs.c11.area*2,...
    handles.main.UserData.win.span);

handles.Msg.ForegroundColor = 'b';
handles.Msg.String = 'Finish!';


UpdateLoss(hObject, eventdata, handles);

handles.ShowSimMesh.Enable = 'on';
handles.ShowSimBmag.Enable = 'on';
handles.UpdateLosses.Enable = 'on';
handles.ShowHyst.Enable = 'on';
handles.ShowEddy.Enable = 'on';
handles.Showbemfs.Enable = 'on';
handles.Showlf.Enable = 'on';
handles.Showbemfs.Enable = 'on';

function Kh_Callback(hObject, eventdata, handles)
% hObject    handle to Kh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Kh as text
%        str2double(get(hObject,'String')) returns contents of Kh as a double


% --- Executes during object creation, after setting all properties.
function Kh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Kh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alpha_Callback(hObject, eventdata, handles)
% hObject    handle to alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha as text
%        str2double(get(hObject,'String')) returns contents of alpha as a double


% --- Executes during object creation, after setting all properties.
function alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ke_Callback(hObject, eventdata, handles)
% hObject    handle to Ke (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ke as text
%        str2double(get(hObject,'String')) returns contents of Ke as a double


% --- Executes during object creation, after setting all properties.
function Ke_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ke (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ShowSimMesh.
function ShowSimMesh_Callback(hObject, eventdata, handles)
% hObject    handle to ShowSimMesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure
for j = 1:length(handles.main.UserData.mF)
    imshow(handles.main.UserData.mF(j).cdata)
    pause(0.01);
end



% --- Executes on button press in ShowSimBmag.
function ShowSimBmag_Callback(hObject, eventdata, handles)
% hObject    handle to ShowSimBmag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure
for j = 1:length(handles.main.UserData.bF)
    imshow(handles.main.UserData.bF(j).cdata)
    pause(0.01);
end




function edit31_Callback(hObject, eventdata, handles)
% hObject    handle to edit31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit31 as text
%        str2double(get(hObject,'String')) returns contents of edit31 as a double


% --- Executes during object creation, after setting all properties.
function edit31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit30_Callback(hObject, eventdata, handles)
% hObject    handle to edit30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit30 as text
%        str2double(get(hObject,'String')) returns contents of edit30 as a double


% --- Executes during object creation, after setting all properties.
function edit30_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit29_Callback(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit29 as text
%        str2double(get(hObject,'String')) returns contents of edit29 as a double


% --- Executes during object creation, after setting all properties.
function edit29_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit28_Callback(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit28 as text
%        str2double(get(hObject,'String')) returns contents of edit28 as a double


% --- Executes during object creation, after setting all properties.
function edit28_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ShowHyst.
function ShowHyst_Callback(hObject, eventdata, handles)
% hObject    handle to ShowHyst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pm = handles.main.UserData.pm;
%% Calculation of Hysteresis Loss Using GSE Method
% angularSpeed
angularSpeed = str2double(handles.rpm.String)*pi/30;
%% lamination data
% hysteresis loss coefficient
Kh =  str2double(handles.Kh.String);
% power of frequency
alpha = str2double(handles.alpha.String);
% power of peak of magnetic flux density
beta = 2;
% calculation of KGSE
temp = linspace(0,2*pi,1000);
KGSE = Kh/((2*pi)^(alpha-1))/(sum((abs(cos(temp))).^alpha .* ...
    (abs(sin(temp))).^(beta-alpha))*temp(2));
% zone
ZoneIndex = handles.main.UserData.m.mzs.stator.zi;
bx = handles.main.UserData.Bxs;
by = handles.main.UserData.Bys;

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

xrotorAngle = handles.main.UserData.p.xrotorAngle;
simulationAngle = handles.main.UserData.p.simulationAngle;

physt = physt*7650*KGSE* angularSpeed^alpha * length(xrotorAngle)^(alpha-1)/(simulationAngle)^alpha;

TotalHyst = handles.main.UserData.m.gta(handles.main.UserData.m.tzi(:,ZoneIndex))* physt;

TotalHyst = TotalHyst * 1e-9 * (pm.Nm/pm.Nsplit) * pm.Lst;

handles.Hyst.String = num2str(TotalHyst);

figure
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


% --- Executes on button press in ShowEddy.
function ShowEddy_Callback(hObject, eventdata, handles)
% hObject    handle to ShowEddy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pm = handles.main.UserData.pm;
%% Calculation of Hysteresis Loss Using GSE Method

angularSpeed = str2double(handles.rpm.String)*pi/30;

%% Calculation of Eddy Current Loss in Electrical Lamination
% zone index

ZoneIndex = handles.main.UserData.m.mzs.stator.zi;
bx = handles.main.UserData.Bxs;
by = handles.main.UserData.Bys;

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

simulationAngle = handles.main.UserData.p.simulationAngle;
xrotorAngle = handles.main.UserData.p.xrotorAngle;

Ke = str2double(handles.Ke.String);
peddy =  (Ke/2/pi^2)*peddy*7650*angularSpeed^2/simulationAngle/xrotorAngle(2);

TotalEddy = handles.main.UserData.m.gta(handles.main.UserData.m.tzi(:,ZoneIndex)) * peddy;

TotalEddy = TotalEddy * 1e-9 * pm.Lst * (pm.Nm/pm.Nsplit);

handles.Eddy.String = num2str(TotalEddy);

figure
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

% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileName,PathName] = uiputfile('*.pmsm','File Selector',[cd,'\project.pmsm']);

if ~FileName
    return
end

f = fopen([PathName,FileName],'w');
fprintf(f,['Ns\t',handles.Ns.String,'\n']);
fprintf(f,['Nm\t',handles.Nm.String,'\n']);
fprintf(f,['Rso\t',handles.Rso.String,'\n']);
fprintf(f,['Rro\t',handles.Rro.String,'\n']);
fprintf(f,['wsy\t',handles.wsy.String,'\n']);
fprintf(f,['wry\t',handles.wry.String,'\n']);
fprintf(f,['wtb\t',handles.wtb.String,'\n']);
fprintf(f,['g\t',handles.g.String,'\n']);
fprintf(f,['hm\t',handles.hm.String,'\n']);
fprintf(f,['Lst\t',handles.Lst.String,'\n']);
fprintf(f,['beta_m\t',handles.beta_m.String,'\n']);

fprintf(f,['a\t',handles.a.String,'\n']);
fprintf(f,['Ncoil\t',handles.Ncoil.String,'\n']);
fprintf(f,['Kf\t',handles.Kf.String,'\n']);
fprintf(f,['con\t',handles.con.String,'\n']);

fprintf(f,['Nsim\t',handles.Nsim.String,'\n']);
fprintf(f,['rpm\t',handles.rpm.String,'\n']);
fprintf(f,['Pin\t',handles.Pin.String,'\n']);
fprintf(f,['Vrms\t',handles.Vrms.String,'\n']);
fprintf(f,['Hc\t',handles.Hc.String,'\n']);
fprintf(f,['pf\t',handles.pf.String,'\n']);

fprintf(f,['Kh\t',handles.Kh.String,'\n']);
fprintf(f,['alpha\t',handles.alpha.String,'\n']);
fprintf(f,['Ke\t',handles.Ke.String,'\n']);
fprintf(f,['Density\t',handles.Density.String,'\n']);

fprintf(f,['thetag\t',handles.thetag.String,'\n']);
fprintf(f,['angle1\t',handles.angle1.String,'\n']);
fprintf(f,['angle2\t',handles.angle2.String,'\n']);

fclose(f);


function edit34_Callback(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit34 as text
%        str2double(get(hObject,'String')) returns contents of edit34 as a double


% --- Executes during object creation, after setting all properties.
function edit34_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit33_Callback(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit33 as text
%        str2double(get(hObject,'String')) returns contents of edit33 as a double


% --- Executes during object creation, after setting all properties.
function edit33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit32_Callback(hObject, eventdata, handles)
% hObject    handle to edit32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit32 as text
%        str2double(get(hObject,'String')) returns contents of edit32 as a double


% --- Executes during object creation, after setting all properties.
function edit32_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Lst_Callback(hObject, eventdata, handles)
% hObject    handle to Lst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Lst as text
%        str2double(get(hObject,'String')) returns contents of Lst as a double

% --- Executes during object creation, after setting all properties.
function Lst_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Lst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rpm_Callback(hObject, eventdata, handles)
% hObject    handle to rpm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rpm as text
%        str2double(get(hObject,'String')) returns contents of rpm as a double


% --- Executes during object creation, after setting all properties.
function rpm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rpm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Pin_Callback(hObject, eventdata, handles)
% hObject    handle to Pin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pin as text
%        str2double(get(hObject,'String')) returns contents of Pin as a double


% --- Executes during object creation, after setting all properties.
function Pin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Vrms_Callback(hObject, eventdata, handles)
% hObject    handle to Vrms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Vrms as text
%        str2double(get(hObject,'String')) returns contents of Vrms as a double


% --- Executes during object creation, after setting all properties.
function Vrms_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Vrms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Density_Callback(hObject, eventdata, handles)
% hObject    handle to Density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Density as text
%        str2double(get(hObject,'String')) returns contents of Density as a double


% --- Executes during object creation, after setting all properties.
function Density_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in win.
function win_Callback(hObject, eventdata, handles)
% hObject    handle to win (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileName,PathName] = uigetfile('*.win','Select .win File','win');
if ~FileName
    return
end
handles.main.UserData.p.windir = [PathName,FileName];

f = fopen(handles.main.UserData.p.windir,'r');
win = fscanf(f,'%d');
fclose(f);

if rem((length(win)-1),6)
    handles.Msg.ForegroundColor = 'r';
    handles.Msg.String = ['Error: Windings data does not set properly',...
        ' ... For balanced winding number of coils must be an multiplier of 3'];
    return
end

Ns = str2double(handles.Ns.String);
if (length(win)-1)/2 ~= Ns
    handles.Msg.ForegroundColor = 'r';
    handles.Msg.String = ['Error: Windings data does not set properly',...
        ' ... Number of coils must be equal to the number of slots'];
    return
end

span = win(1);

coils = (1:str2double(handles.Ns.String))';
coils = [coils,circshift(coils,-span)];

win = reshape(win(2:end),2,[]);
win = win';
win = reshape(win,[],6);

phaseA = win(:,[1,4]);
phaseB = win(:,[2,5]);
phaseC = win(:,[3,6]);

if any(bitor(phaseA(:,1)>Ns,phaseA(:,1)<1))
    handles.Msg.ForegroundColor = 'r';
    handles.Msg.String = ['Error: Windings data does not set properly',...
        ' ... Some coil index is grater than Ns'];
    return
end

if any(bitor(phaseB(:,1)>Ns,phaseB(:,1)<1))
    handles.Msg.ForegroundColor = 'r';
    handles.Msg.String = ['Error: Windings data does not set properly',...
        ' ... Some coil index is grater than Ns'];
    return
end

if any(bitor(phaseC(:,1)>Ns,phaseC(:,1)<1))
    handles.Msg.ForegroundColor = 'r';
    handles.Msg.String = ['Error: Windings data does not set properly',...
        ' ... Some coil index is grater than Ns'];
    return
end

handles.main.UserData.win.span = span;
handles.main.UserData.win.coils = coils;
handles.main.UserData.win.phaseA = phaseA ;
handles.main.UserData.win.phaseB = phaseB ;
handles.main.UserData.win.phaseC = phaseC ;

handles.cgm.Enable = 'on';
handles.Msg.ForegroundColor = 'b';
handles.Msg.String = 'Windings are configurated completed!';

% --- Executes on button press in UpdateLosses.
function UpdateLosses_Callback(hObject, eventdata, handles)
% hObject    handle to UpdateLosses (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UpdateLoss(hObject, eventdata, handles)
% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile('*.pmsm');

if ~FileName
    return
end

f = fopen([PathName,FileName],'r');
while true
    str = fgets(f);
    if str==-1
        break
    end
    str = strsplit(str);
    handles.(str{1}).String = str{2};
end
fclose(f);

function y = evalEddyLoss(hObject, eventdata, handles)


% Calculation of Hysteresis Loss Using GSE Method
pm = handles.main.UserData.pm;
angularSpeed = str2double(handles.rpm.String)*pi/30;
xrotorAngle = handles.main.UserData.p.xrotorAngle;
simulationAngle = handles.main.UserData.p.simulationAngle;

% stator zone 
ZoneIndex = handles.main.UserData.m.mzs.stator.zi;
bx = handles.main.UserData.Bxs;
by = handles.main.UserData.Bys;

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


Ke = str2double(handles.Ke.String);
peddy =  (Ke/2/pi^2)*peddy*7650*angularSpeed^2/simulationAngle/xrotorAngle(2);

TotalEddy = handles.main.UserData.m.gta(handles.main.UserData.m.tzi(:,ZoneIndex)) * peddy;

y = TotalEddy * 1e-9 * pm.Lst * (pm.Nm/pm.Nsplit);

function y = evalHystLoss(hObject, eventdata, handles)

% Calculation of Hysteresis Loss Using GSE Method
pm = handles.main.UserData.pm;
angularSpeed = str2double(handles.rpm.String)*pi/30;
xrotorAngle = handles.main.UserData.p.xrotorAngle;
simulationAngle = handles.main.UserData.p.simulationAngle;

% lamination data
% hysteresis loss coefficient
Kh =  str2double(handles.Kh.String);
% power of frequency
alpha = str2double(handles.alpha.String);
% power of peak of magnetic flux density
beta = 2;
% calculation of KGSE
temp = linspace(0,2*pi,1000);
KGSE = Kh/((2*pi)^(alpha-1))/(sum((abs(cos(temp))).^alpha .* ...
    (abs(sin(temp))).^(beta-alpha))*temp(2));

% stator zone
ZoneIndex = handles.main.UserData.m.mzs.stator.zi;
bx = handles.main.UserData.Bxs;
by = handles.main.UserData.Bys;

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

physt = physt*7650*KGSE* angularSpeed^alpha * ...
    length(xrotorAngle)^(alpha-1)/(simulationAngle)^alpha;
TotalHyst = handles.main.UserData.m.gta(handles.main.UserData.m.tzi(:,ZoneIndex))* physt;
y = TotalHyst * 1e-9 * (pm.Nm/pm.Nsplit) * pm.Lst;

function UpdateLoss(hObject, eventdata, handles)
handles.Eddy.String = num2str(evalEddyLoss(hObject, eventdata, handles));
handles.Hyst.String = num2str(evalHystLoss(hObject, eventdata, handles));

handles.Ohmic.String = num2str(evalOhmicLoss(hObject, eventdata, handles));
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
function Untitled_6_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_7_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox({'Loss Analysis of PMSM','Developed by Ali Jamalifard'});


% --------------------------------------------------------------------
function Untitled_8_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_9_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_10_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pushbutton9_Callback(hObject, eventdata, handles)



function a_Callback(hObject, eventdata, handles)
% hObject    handle to a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a as text
%        str2double(get(hObject,'String')) returns contents of a as a double

% --- Executes during object creation, after setting all properties.
function a_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ncoil_Callback(hObject, eventdata, handles)
% hObject    handle to Ncoil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ncoil as text
%        str2double(get(hObject,'String')) returns contents of Ncoil as a double

% --- Executes during object creation, after setting all properties.
function Ncoil_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ncoil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Kf_Callback(hObject, eventdata, handles)
% hObject    handle to Kf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Kf as text
%        str2double(get(hObject,'String')) returns contents of Kf as a double

% --- Executes during object creation, after setting all properties.
function Kf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Kf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function con_Callback(hObject, eventdata, handles)
% hObject    handle to con (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of con as text
%        str2double(get(hObject,'String')) returns contents of con as a double

% --- Executes during object creation, after setting all properties.
function con_CreateFcn(hObject, eventdata, handles)
% hObject    handle to con (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_11_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pushbutton12_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function Untitled_12_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_13_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function Untitled_14_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function pf_Callback(hObject, eventdata, handles)
% hObject    handle to pf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pf as text
%        str2double(get(hObject,'String')) returns contents of pf as a double


% --- Executes during object creation, after setting all properties.
function pf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function y = evalOhmicLoss(hObject, eventdata, handles)

Pin = str2double(handles.Pin.String);
Vrms = str2double(handles.Vrms.String);
pf = str2double(handles.pf.String);

Irms = Pin*1e3/sqrt(3)/Vrms/pf;

y = 3*handles.main.UserData.p.Rdc*Irms^2;



function span_Callback(hObject, eventdata, handles)
% hObject    handle to span (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of span as text
%        str2double(get(hObject,'String')) returns contents of span as a double


% --- Executes during object creation, after setting all properties.
function span_CreateFcn(hObject, eventdata, handles)
% hObject    handle to span (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Hc_Callback(hObject, eventdata, handles)
% hObject    handle to Hc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Hc as text
%        str2double(get(hObject,'String')) returns contents of Hc as a double


% --- Executes during object creation, after setting all properties.
function Hc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Hc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when main is resized.
function main_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit57_Callback(hObject, eventdata, handles)
% hObject    handle to edit57 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit57 as text
%        str2double(get(hObject,'String')) returns contents of edit57 as a double


% --- Executes during object creation, after setting all properties.
function edit57_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit57 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Hyst_Callback(hObject, eventdata, handles)
% hObject    handle to Hyst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Hyst as text
%        str2double(get(hObject,'String')) returns contents of Hyst as a double


% --- Executes during object creation, after setting all properties.
function Hyst_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Hyst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Eddy_Callback(hObject, eventdata, handles)
% hObject    handle to Eddy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Eddy as text
%        str2double(get(hObject,'String')) returns contents of Eddy as a double


% --- Executes during object creation, after setting all properties.
function Eddy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Eddy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ohmic_Callback(hObject, eventdata, handles)
% hObject    handle to Ohmic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ohmic as text
%        str2double(get(hObject,'String')) returns contents of Ohmic as a double


% --- Executes during object creation, after setting all properties.
function Ohmic_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ohmic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit61_Callback(hObject, eventdata, handles)
% hObject    handle to edit61 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit61 as text
%        str2double(get(hObject,'String')) returns contents of edit61 as a double


% --- Executes during object creation, after setting all properties.
function edit61_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit61 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit62_Callback(hObject, eventdata, handles)
% hObject    handle to edit62 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit62 as text
%        str2double(get(hObject,'String')) returns contents of edit62 as a double


% --- Executes during object creation, after setting all properties.
function edit62_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit62 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function beta_m_Callback(hObject, eventdata, handles)
% hObject    handle to beta_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of beta_m as text
%        str2double(get(hObject,'String')) returns contents of beta_m as a double
OffButtuns(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function beta_m_CreateFcn(hObject, eventdata, handles)
% hObject    handle to beta_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ReadingData(hObject, eventdata, handles)

% motor data
handles.main.UserData.pm = PMSM;
handles.main.UserData.pm.Ns = str2double(handles.Ns.String);
handles.main.UserData.pm.Nm = str2double(handles.Nm.String);
handles.main.UserData.pm.Rso = str2double(handles.Rso.String);
handles.main.UserData.pm.Rro = str2double(handles.Rro.String);
handles.main.UserData.pm.wry = str2double(handles.wry.String);
handles.main.UserData.pm.wsy = str2double(handles.wsy.String);
handles.main.UserData.pm.wtb = str2double(handles.wtb.String);
handles.main.UserData.pm.hm = str2double(handles.hm.String);
handles.main.UserData.pm.Lst = str2double(handles.Lst.String);
handles.main.UserData.pm.a = str2double(handles.a.String);
handles.main.UserData.pm.con = str2double(handles.con.String);
handles.main.UserData.pm.Kf = str2double(handles.Kf.String);
handles.main.UserData.pm.Ncoil = str2double(handles.Ncoil.String);
handles.main.UserData.pm.beta_m = str2double(handles.beta_m.String);

% mesh sizes
handles.main.UserData.p = struct;
handles.main.UserData.p.thetag = str2double(handles.thetag.String);
handles.main.UserData.p.angle1 = str2double(handles.angle1.String);
handles.main.UserData.p.angle2 = str2double(handles.angle2.String);

% Loss Data
handles.main.UserData.p.Kh = str2double(handles.Kh.String);
handles.main.UserData.p.alpha = str2double(handles.alpha.String);
handles.main.UserData.p.Ke = str2double(handles.Ke.String);
handles.main.UserData.p.Density = str2double(handles.Density.String);


% --- Executes during object creation, after setting all properties.
function main_CreateFcn(hObject, eventdata, handles)
% hObject    handle to main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

function OffButtuns(hObject, eventdata, handles)

handles.ShowSimMesh.Enable = 'off';
handles.ShowSimBmag.Enable = 'off';
handles.UpdateLosses.Enable = 'off';
handles.ShowHyst.Enable = 'off';
handles.ShowEddy.Enable = 'off';
handles.Showmzs.Enable = 'off';
handles.Showgm.Enable = 'off';
handles.RunSimulation.Enable = 'off';
handles.Showlf.Enable = 'of';
handles.Showbemfs.Enable = 'of';
handles.MQ.Enable = 'off';



function edit64_Callback(hObject, eventdata, handles)
% hObject    handle to edit64 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit64 as text
%        str2double(get(hObject,'String')) returns contents of edit64 as a double


% --- Executes during object creation, after setting all properties.
function edit64_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit64 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Showbemfs.
function Showbemfs_Callback(hObject, eventdata, handles)
% hObject    handle to Showbemfs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pm = handles.main.UserData.pm;
ea = ppval(handles.main.UserData.ea,handles.main.UserData.xi);
eb = ppval(handles.main.UserData.eb,handles.main.UserData.xi);
ec = ppval(handles.main.UserData.ec,handles.main.UserData.xi);

figure
subplot(211)
hold all
plot(handles.main.UserData.xi*180/pi * pm.Nm/2,ea);
plot(handles.main.UserData.xi*180/pi * pm.Nm/2,eb);
plot(handles.main.UserData.xi*180/pi * pm.Nm/2,ec);
legend('phase A','phase B','phase C')
xlabel('Electrical Angle')
ylabel('Back EMF [V]')
set(gca,'XLim',[0,360])
title('Phase BEMFs');

subplot(212)
hold all
plot(handles.main.UserData.xi*180/pi * pm.Nm/2,ea-eb);
plot(handles.main.UserData.xi*180/pi * pm.Nm/2,eb-ec);
plot(handles.main.UserData.xi*180/pi * pm.Nm/2,ec-ea);
legend('E_{AB}','E_{BC}','E_{CA}')
xlabel('Electrical Angle')
ylabel('Back EMF [V]')
set(gca,'XLim',[0,360])
title('Line to Line BEMFs');




% --- Executes on button press in Showlf.
function Showlf_Callback(hObject, eventdata, handles)
% hObject    handle to Showlf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pm = handles.main.UserData.pm;
figure
hold all
plot(handles.main.UserData.xi*180/pi * pm.Nm/2,handles.main.UserData.aLFlux)
plot(handles.main.UserData.xi*180/pi * pm.Nm/2,handles.main.UserData.bLFlux)
plot(handles.main.UserData.xi*180/pi * pm.Nm/2,handles.main.UserData.cLFlux)
legend('phase A','phase B','phase C')
xlabel('Electrical Angle')
ylabel('Linkage Flux [wb]')
set(gca,'XLim',[0,360])


% --- Executes when user attempts to close main.
function main_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


% --------------------------------------------------------------------
function Untitled_15_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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


% --- Executes on button press in MQ.
function MQ_Callback(hObject, eventdata, handles)
% hObject    handle to MQ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[y1,y2,y3] = handles.main.UserData.m.getQuality;

handles.Msg.String = [' ave = ',num2str(y1),...
    '   *   min = ',num2str(y2),'   *   max = ',num2str(y3)];

% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in listbox4.
function listbox4_Callback(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox4


% --- Executes during object creation, after setting all properties.
function listbox4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox3.
function listbox3_Callback(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox3


% --- Executes during object creation, after setting all properties.
function listbox3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in NL.
function NL_Callback(hObject, eventdata, handles)
% hObject    handle to NL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of NL

