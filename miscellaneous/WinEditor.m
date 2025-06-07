function varargout = WinEditor(varargin)
% WINEDITOR MATLAB code for WinEditor.fig
%      WINEDITOR, by itself, creates a new WINEDITOR or raises the existing
%      singleton*.
%
%      H = WINEDITOR returns the handle to a new WINEDITOR or the handle to
%      the existing singleton*.
%
%      WINEDITOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WINEDITOR.M with the given input arguments.
%
%      WINEDITOR('Property','Value',...) creates a new WINEDITOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before WinEditor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to WinEditor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help WinEditor

% Last Modified by GUIDE v2.5 01-May-2017 23:35:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @WinEditor_OpeningFcn, ...
                   'gui_OutputFcn',  @WinEditor_OutputFcn, ...
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


% --- Executes just before WinEditor is made visible.
function WinEditor_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to WinEditor (see VARARGIN)

% Choose default command line output for WinEditor
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes WinEditor wait for user response (see UIRESUME)
% uiwait(handles.figure1);

iw = 0.025;
ih = 0.02;
Wpc = 0.35;
Wv = 1-Wpc;
Hp = 0.3;
Hc = 1-Hp;

handles.cp.Units = 'normalized';
handles.cp.Position = [iw,ih,Wpc-3*iw/2,Hc-3*ih/2];

handles.pp.Units = 'normalized';
handles.pp.Position = [iw,Hc+ih/2,Wpc-3*iw/2,Hp-3*ih/2];

handles.vp.Units = 'normalized';
handles.vp.Position = [Wpc+iw/2,ih,Wv-3*iw/2,1-2*ih];

Wt = 0.025;
Ht = 0.05;

handles.pt.Units = 'normalized';
handles.pt.Position = [Wt,Ht,1-2*Wt,1-2*Ht];

handles.ct.Units = 'normalized';
handles.ct.Position = [Wt,Ht,1-2*Wt,1-2*Ht];

index = 0;

index = index + 1;
handles.pt.Data{index,1} = 'm';
handles.pt.Data{index,2} = 3;
handles.pt.Data{index,3} = 'Number of phases';

index = index + 1;
handles.pt.Data{index,1} = 'Ns';
handles.pt.Data{index,2} = 18;
handles.pt.Data{index,3} = 'Number of slots';

index = index + 1;
handles.pt.Data{index,1} = 'S';
handles.pt.Data{index,2} = 2;
handles.pt.Data{index,3} = 'Coil span';

index = index + 1;
handles.pt.Data{index,1} = 'L';
handles.pt.Data{index,2} = 2;
handles.pt.Data{index,3} = 'Number of layers';

Setct(handles)

% --- Outputs from this function are returned to the command line.
function varargout = WinEditor_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


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


function Setct(handles)

m = handles.pt.Data{1,2};
Ns = handles.pt.Data{2,2};
S = handles.pt.Data{3,2};
L = handles.pt.Data{4,2};

handles.ct.Data = cell(Ns,2*m);
Name = cell(1,2*m);
Format = cell(1,2*m);
for i = 1:2:2*m
    Name{i} = ['Ph',num2str((i+1)/2)];
    Name{i+1} = ['d',num2str((i+1)/2)];
    Format{i} = 'numeric';
    Format{i+1} = 'numeric';
end
handles.ct.ColumnName = Name; 
handles.ct.ColumnFormat = Format;
handles.ct.ColumnEditable = true(1,2*m);
handles.ct.ColumnWidth = {35};
% --- Executes when entered data in editable cell(s) in pt.
function pt_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to pt (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
Setct(handles)

function ExportWin(handles)

[FileName,PathName] = uiputfile('*.win','File Selector',[cd,'\win.win']);

if ~FileName
    return
end

f = fopen([PathName,FileName],'w');

S = handles.pt.Data{3,2};
fprintf(f,'%d',S);
t = handles.ct.Data;
for i = 1:2:size(t,2)
    for j = 1:size(t,1)
        if ~isempty(t{j,i})
            fprintf('\n');
            fprintf(f,'%d',j);
            fprintf('\n');
            fprintf(f,'%d',t{j,i+1});
        end
    end
end
    
fclose(f);


% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ExportWin(handles)


% --------------------------------------------------------------------
function Untitled_4_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
