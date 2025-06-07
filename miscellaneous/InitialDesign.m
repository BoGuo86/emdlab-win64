function varargout = InitialDesign(varargin)
% INITIALDESIGN MATLAB code for InitialDesign.fig
%      INITIALDESIGN, by itself, creates a new INITIALDESIGN or raises the existing
%      singleton*.
%
%      H = INITIALDESIGN returns the handle to a new INITIALDESIGN or the handle to
%      the existing singleton*.
%
%      INITIALDESIGN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INITIALDESIGN.M with the given input arguments.
%
%      INITIALDESIGN('Property','Value',...) creates a new INITIALDESIGN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before InitialDesign_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to InitialDesign_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help InitialDesign

% Last Modified by GUIDE v2.5 08-May-2017 16:39:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @InitialDesign_OpeningFcn, ...
                   'gui_OutputFcn',  @InitialDesign_OutputFcn, ...
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


% --- Executes just before InitialDesign is made visible.
function InitialDesign_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to InitialDesign (see VARARGIN)

% Choose default command line output for InitialDesign
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes InitialDesign wait for user response (see UIRESUME)
% uiwait(handles.main);

handles.msg.Units = 'normalized';
handles.msg.Position = [0.05,0.05,0.9,0.15];

handles.inputs.Units = 'normalized';
handles.inputs.Position = [0.05,0.25,0.425,0.7];
handles.outputs.Units = 'normalized';
handles.outputs.Position = [0.525,0.25,0.425,0.7];

handles.it.Units = 'normalized';
handles.it.Position = [0.05,0.05,0.9,0.9];
handles.ot.Units = 'normalized';
handles.ot.Position = [0.05,0.05,0.9,0.9];



index = 0;

index = index + 1;
handles.it.Data{index,1} = 'Pout';
handles.it.Data{index,2} = 0.5;
handles.it.Data{index,3} = '[KW]';

index = index + 1;
handles.it.Data{index,1} = 'rpm';
handles.it.Data{index,2} = 1500;
handles.it.Data{index,3} = '[rpm]';

index = index + 1;
handles.it.Data{index,1} = 'Vrms';
handles.it.Data{index,2} = 127;
handles.it.Data{index,3} = '[V]';

index = index + 1;
handles.it.Data{index,1} = 'Rso';
handles.it.Data{index,2} = 50;
handles.it.Data{index,3} = '';

index = index + 1;
handles.it.Data{index,1} = 'D';
handles.it.Data{index,2} = 0.6;
handles.it.Data{index,3} = '';

index = index + 1;
handles.it.Data{index,1} = 'Ns';
handles.it.Data{index,2} = 18;
handles.it.Data{index,3} = '';

index = index + 1;
handles.it.Data{index,1} = 'Nm';
handles.it.Data{index,2} = 8;
handles.it.Data{index,3} = '';

index = index + 1;
handles.it.Data{index,1} = 'm';
handles.it.Data{index,2} = 3;
handles.it.Data{index,3} = '';

index = index + 1;
handles.it.Data{index,1} = 'g';
handles.it.Data{index,2} = 1;
handles.it.Data{index,3} = '[mm]';

index = index + 1;
handles.it.Data{index,1} = 'Lst';
handles.it.Data{index,2} = 100;
handles.it.Data{index,3} = '[mm]';

index = index + 1;
handles.it.Data{index,1} = 'hm';
handles.it.Data{index,2} = 4;
handles.it.Data{index,3} = '[mm]';

index = index + 1;
handles.it.Data{index,1} = 'betam';
handles.it.Data{index,2} = 0.9;
handles.it.Data{index,3} = '';

index = index + 1;
handles.it.Data{index,1} = 'Bg';
handles.it.Data{index,2} = 0.87;
handles.it.Data{index,3} = '[Tesla]';

index = index + 1;
handles.it.Data{index,1} = 'Bt';
handles.it.Data{index,2} = 1.4;
handles.it.Data{index,3} = '[Tesla]';

index = index + 1;
handles.it.Data{index,1} = 'Bsy';
handles.it.Data{index,2} = 1.4;
handles.it.Data{index,3} = '[Tesla]';

index = index + 1;
handles.it.Data{index,1} = 'Bry';
handles.it.Data{index,2} = 1.4;
handles.it.Data{index,3} = '[Tesla]';

index = index + 1;
handles.it.Data{index,1} = 'alphas';
handles.it.Data{index,2} = 0.1;
handles.it.Data{index,3} = '';

index = index + 1;
handles.it.Data{index,1} = 'hss1';
handles.it.Data{index,2} = 1;
handles.it.Data{index,3} = '[mm]';

index = index + 1;
handles.it.Data{index,1} = 'hss2';
handles.it.Data{index,2} = 0.5;
handles.it.Data{index,3} = '[mm]';

index = index + 1;
handles.it.Data{index,1} = 'a';
handles.it.Data{index,2} = 1;
handles.it.Data{index,3} = '';

index = index + 1;
handles.it.Data{index,1} = 'sigma';
handles.it.Data{index,2} = 57e6;
handles.it.Data{index,3} = '[Ms]';

index = index + 1;
handles.it.Data{index,1} = 'Kf';
handles.it.Data{index,2} = 0.45;
handles.it.Data{index,3} = '';


handles.main.UserData.pm = ID_PMSM0;
ReadData(handles)

index = 0;

index = index + 1;
handles.ot.Data{index,1} = 'wsy';
handles.ot.Data{index,2} = 0;
handles.ot.Data{index,3} = '[mm]';

index = index + 1;
handles.ot.Data{index,1} = 'wry';
handles.ot.Data{index,2} = 0;
handles.ot.Data{index,3} = '[mm]';

index = index + 1;
handles.ot.Data{index,1} = 'wtb';
handles.ot.Data{index,2} = 0;
handles.ot.Data{index,3} = '[mm]';

index = index + 1;
handles.ot.Data{index,1} = 'Rro';
handles.ot.Data{index,2} = 0;
handles.ot.Data{index,3} = '[mm]';

index = index + 1;
handles.ot.Data{index,1} = 'rsh';
handles.ot.Data{index,2} = 0;
handles.ot.Data{index,3} = '[mm]';

index = index + 1;
handles.ot.Data{index,1} = 'Ncoil';
handles.ot.Data{index,2} = 0;
handles.ot.Data{index,3} = '';

WriteData(handles)

% --- Outputs from this function are returned to the command line.
function varargout = InitialDesign_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.main.UserData.pm;


% --- Executes when entered data in editable cell(s) in it.
function it_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to it (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

ReadData(handles)
WriteData(handles)

function ReadData(handles)


for i = 1:size(handles.it.Data,1)
    switch handles.it.Data{i,1}
        case 'Pout'
            handles.main.UserData.pm.(handles.it.Data{i,1}) = ...
                handles.it.Data{i,2}*1e3;
        otherwise
            handles.main.UserData.pm.(handles.it.Data{i,1}) = ...
                handles.it.Data{i,2};
    end
end

function WriteData(handles)

for i = 1:size(handles.ot.Data,1)
        handles.ot.Data{i,2} = ...
        handles.main.UserData.pm.(handles.ot.Data{i,1});
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

[FileName,PathName] = uiputfile('*.pmsm','File Selector',[cd,'\project.pmsm']);

if ~FileName
    return
end

f = fopen([PathName,FileName],'w');

tables = {'it','ot'};
for i = 1:length(tables)
    table = tables{i};
    for j = 1:size(handles.(table).Data,1)
        fprintf(f,'%s\t%f\n',handles.(table).Data{j,1},handles.(table).Data{j,2});
    end
end

fclose(f);

% --- Executes when user attempts to close main.
function y = main_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


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

options.Interpreter = 'tex';
prompt={'\chi = L/D','42500 < \sigma <22500'};
name = 'D2L: LD_sigma';
defaultans = {'1','33500'};
answer = inputdlg(prompt,name,[1 40],defaultans,options);

if isempty(answer)
    return
end

T = handles.main.UserData.pm.T;
LD = str2double(answer{1});
TRV = 2*str2double(answer{2});
V = T*1e9/TRV;
D = nthroot(4*V/pi/LD,3);

str{1} = ['Rro[mm] = ',num2str(D/2)];
str{2} = ['L[mm] = ',num2str(LD*D)];
WriteMsg(handles,1,str)



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


% --------------------------------------------------------------------
function Untitled_10_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_11_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_15_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

g = handles.main.UserData.pm.g;
prompt={'Bg [Tesla]','g [mm]','Hc [A/m]'};
name = 'Find hm';
defaultans = {'0.87',num2str(g),'9.8524e+05'};
answer = inputdlg(prompt,name,[1 40],defaultans);
hm = 1.2*str2double(answer{1})*str2double(answer{2})/...
    (str2double(answer{3})*4*pi*1e-7);
msgbox(['Estimated Value: hm[mm]=',num2str(hm)]);

% --------------------------------------------------------------------
function Untitled_16_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

options.Interpreter = 'tex';
Bg = handles.main.UserData.pm.Bg;
C = handles.main.UserData.pm.betam;
g = handles.main.UserData.pm.g;
prompt={'B_g [Tesla]','g [mm]','B_r [A/m]','K_r','K_l','\mu_R','C_{\phi}'};
name = 'Find hm';
defaultans = {num2str(Bg),num2str(g),'1.3','1.2','0.95','1.05',num2str(C)};
answer = inputdlg(prompt,name,[1 40],defaultans,options);

if isempty(answer)
    return
end

Bg = str2double(answer{1});
g = str2double(answer{2});
Br = str2double(answer{3});
Kr = str2double(answer{4});
Kl = str2double(answer{5});
mur = str2double(answer{6});
C = str2double(answer{7});

c = Kl*C*Br/Bg;
b = Kr*mur*C*g;
hm = b/(c-1);
WriteMsg(handles,1,['Estimated Value: hm[mm] = ',num2str(hm)]);


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

y = handles.main.UserData.pm.getD;

WriteMsg(handles,1,['D = ',num2str(y)])


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
    handles.msg.String = [handles.msg.String(:);txt(:)];
else
    handles.msg.String{2} = txt;
end
