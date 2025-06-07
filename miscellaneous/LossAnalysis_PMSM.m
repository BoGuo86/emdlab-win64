function varargout = LossAnalysis_PMSM(varargin)
% LOSSANALYSIS_PMSM MATLAB code for LossAnalysis_PMSM.fig
%      LOSSANALYSIS_PMSM, by itself, creates a new LOSSANALYSIS_PMSM or raises the existing
%      singleton*.
%
%      H = LOSSANALYSIS_PMSM returns the handle to a new LOSSANALYSIS_PMSM or the handle to
%      the existing singleton*.
%
%      LOSSANALYSIS_PMSM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOSSANALYSIS_PMSM.M with the given input arguments.
%
%      LOSSANALYSIS_PMSM('Property','Value',...) creates a new LOSSANALYSIS_PMSM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LossAnalysis_PMSM_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LossAnalysis_PMSM_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LossAnalysis_PMSM

% Last Modified by GUIDE v2.5 03-Aug-2017 11:29:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LossAnalysis_PMSM_OpeningFcn, ...
                   'gui_OutputFcn',  @LossAnalysis_PMSM_OutputFcn, ...
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


% --- Executes just before LossAnalysis_PMSM is made visible.
function LossAnalysis_PMSM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LossAnalysis_PMSM (see VARARGIN)

% UIWAIT makes LossAnalysis_PMSM wait for user response (see UIRESUME)
% uiwait(handles.main);

hObject.UserData.states.isSetMesh = false;
hObject.UserData.states.isSetWin = false;
hObject.UserData.sdclass = pmsm0_design;

%% tab group definition
tabgp = uitabgroup(hObject);
tabgp.TabLocation = 'left';
tabgp.Units = 'normalized';
tabgp.Position = [0,0.15,1,0.85];
tabgp.Tag = 'tabgp';

%% control tab
conTab = uitab(tabgp,'Title','Data Sheet');
conTab.Tag = 'conTab';
coniTable = uitable(conTab);
coniTable.Tag = 'coniTable';

% independent
coniTable.Units = 'normalized';
coniTable.Position = [0,0,0.5,1];
coniTable.FontSize = 13;
coniTable.ColumnEditable = [false,true,false,false]; 
coniTable.ColumnName = {'Parameter','Value','Unit','Description'}; 

index = 0;

index = index + 1;
coniTable.Data{index,1} = 'Pout';
coniTable.Data{index,2} = hObject.UserData.sdclass.(coniTable.Data{index,1});
coniTable.Data{index,3} = '[KW]';
coniTable.Data{index,4} = 'Rated Output Power';

index = index + 1;
coniTable.Data{index,1} = 'rpm';
coniTable.Data{index,2} = hObject.UserData.sdclass.(coniTable.Data{index,1});
coniTable.Data{index,3} = '[rpm]';
coniTable.Data{index,4} = 'Rated Speed';

index = index + 1;
coniTable.Data{index,1} = 'Vrms';
coniTable.Data{index,2} = hObject.UserData.sdclass.(coniTable.Data{index,1});
coniTable.Data{index,3} = '[V]';
coniTable.Data{index,4} = 'Phase Voltage';

index = index + 1;
coniTable.Data{index,1} = 'm';
coniTable.Data{index,2} = hObject.UserData.sdclass.(coniTable.Data{index,1});
coniTable.Data{index,3} = '';
coniTable.Data{index,4} = 'Number of Phases';

index = index + 1;
coniTable.Data{index,1} = 'eff';
coniTable.Data{index,2} = hObject.UserData.sdclass.(coniTable.Data{index,1});
coniTable.Data{index,3} = '';
coniTable.Data{index,4} = 'Efficiency';

index = index + 1;
coniTable.Data{index,1} = 'pf';
coniTable.Data{index,2} = hObject.UserData.sdclass.(coniTable.Data{index,1});
coniTable.Data{index,3} = '';
coniTable.Data{index,4} = 'Input Power Factor';

condTable = uitable(conTab);
condTable.Tag = 'condTable';

% dependent
condTable.Units = 'normalized';
condTable.Position = [0.5,0,0.5,1];
condTable.FontSize = 13;
condTable.ColumnEditable = [false,false,false,false]; 
condTable.ColumnName = {'Parameter','Value','Unit','Description'}; 

index = 0;

index = index + 1;
condTable.Data{index,1} = 'Irms';
condTable.Data{index,2} = hObject.UserData.sdclass.(condTable.Data{index,1});
condTable.Data{index,3} = '[A]';
condTable.Data{index,4} = 'RMS of phasecurent';

index = index + 1;
condTable.Data{index,1} = 'T';
condTable.Data{index,2} = hObject.UserData.sdclass.(condTable.Data{index,1});
condTable.Data{index,3} = '[N.m]';
condTable.Data{index,4} = 'Desired Output Torque';

%% dimension tab
dimTab = uitab(tabgp,'Title','Dimensions');
dimTab.Tag = 'dimTab';

% independent
dimiTable = uitable(dimTab);
dimiTable.Tag = 'dimiTable';
dimiTable.Units = 'normalized';
dimiTable.Position = [0,0,0.5,1];
dimiTable.FontSize = 13;
dimiTable.ColumnEditable = [false,true,false,false]; 
dimiTable.ColumnName = {'Parameter','Value','Unit','Description'};

% setting data
index = 0;

index = index + 1;
dimiTable.Data{index,1} = 'Ns';
dimiTable.Data{index,2} = hObject.UserData.sdclass.(dimiTable.Data{index,1});
dimiTable.Data{index,3} = '';
dimiTable.Data{index,4} = 'Number of Stator Slots';

index = index + 1;
dimiTable.Data{index,1} = 'Nm';
dimiTable.Data{index,2} = hObject.UserData.sdclass.(dimiTable.Data{index,1});
dimiTable.Data{index,3} = '';
dimiTable.Data{index,4} = 'Number of Magnets';

index = index + 1;
dimiTable.Data{index,1} = 'Rso';
dimiTable.Data{index,2} = hObject.UserData.sdclass.(dimiTable.Data{index,1});
dimiTable.Data{index,3} = '[mm]';
dimiTable.Data{index,4} = 'Outer Radius of Stator Lamination';

index = index + 1;
dimiTable.Data{index,1} = 'D';
dimiTable.Data{index,2} = hObject.UserData.sdclass.(dimiTable.Data{index,1});
dimiTable.Data{index,3} = '';
dimiTable.Data{index,4} = 'Ration of Rro respect to Rso';

index = index + 1;
dimiTable.Data{index,1} = 'wsy';
dimiTable.Data{index,2} = hObject.UserData.sdclass.(dimiTable.Data{index,1});
dimiTable.Data{index,3} = '[mm]';
dimiTable.Data{index,4} = 'Width of Stator Yoke';

index = index + 1;
dimiTable.Data{index,1} = 'wry';
dimiTable.Data{index,2} = hObject.UserData.sdclass.(dimiTable.Data{index,1});
dimiTable.Data{index,3} = '[mm]';
dimiTable.Data{index,4} = 'Width of Rotor Yoke';

index = index + 1;
dimiTable.Data{index,1} = 'wtb';
dimiTable.Data{index,2} = hObject.UserData.sdclass.(dimiTable.Data{index,1});
dimiTable.Data{index,3} = '[mm]';
dimiTable.Data{index,4} = 'Width of Stator Teeth';

index = index + 1;
dimiTable.Data{index,1} = 'g';
dimiTable.Data{index,2} = hObject.UserData.sdclass.(dimiTable.Data{index,1});
dimiTable.Data{index,3} = '[mm]';
dimiTable.Data{index,4} = 'Air Gap Length';

index = index + 1;
dimiTable.Data{index,1} = 'Lst';
dimiTable.Data{index,2} = hObject.UserData.sdclass.(dimiTable.Data{index,1});
dimiTable.Data{index,3} = '[mm]';
dimiTable.Data{index,4} = 'Stack Length';

index = index + 1;
dimiTable.Data{index,1} = 'Kst';
dimiTable.Data{index,2} = hObject.UserData.sdclass.(dimiTable.Data{index,1});
dimiTable.Data{index,3} = '';
dimiTable.Data{index,4} = 'Stacking Factor of Lamination';

index = index + 1;
dimiTable.Data{index,1} = 'hm';
dimiTable.Data{index,2} = hObject.UserData.sdclass.(dimiTable.Data{index,1});
dimiTable.Data{index,3} = '[mm]';
dimiTable.Data{index,4} = 'Height of Magnets';

index = index + 1;
dimiTable.Data{index,1} = 'betam';
dimiTable.Data{index,2} = hObject.UserData.sdclass.(dimiTable.Data{index,1});
dimiTable.Data{index,3} = '';
dimiTable.Data{index,4} = 'Ration of Magnet Arc to Rotor Pole Pitch';

index = index + 1;
dimiTable.Data{index,1} = 'bss1';
dimiTable.Data{index,2} = hObject.UserData.sdclass.(dimiTable.Data{index,1});
dimiTable.Data{index,3} = '[mm]';

index = index + 1;
dimiTable.Data{index,1} = 'hss1';
dimiTable.Data{index,2} = hObject.UserData.sdclass.(dimiTable.Data{index,1});
dimiTable.Data{index,3} = '[mm]';

index = index + 1;
dimiTable.Data{index,1} = 'Bg';
dimiTable.Data{index,2} = hObject.UserData.sdclass.(dimiTable.Data{index,1});
dimiTable.Data{index,3} = '[Tesla]';

index = index + 1;
dimiTable.Data{index,1} = 'Bt';
dimiTable.Data{index,2} = hObject.UserData.sdclass.(dimiTable.Data{index,1});
dimiTable.Data{index,3} = '[Tesla]';

index = index + 1;
dimiTable.Data{index,1} = 'Bsy';
dimiTable.Data{index,2} = hObject.UserData.sdclass.(dimiTable.Data{index,1});
dimiTable.Data{index,3} = '[Tesla]';

index = index + 1;
dimiTable.Data{index,1} = 'Bry';
dimiTable.Data{index,2} = hObject.UserData.sdclass.(dimiTable.Data{index,1});
dimiTable.Data{index,3} = '[Tesla]';

% dependent
dimdTable = uitable(dimTab);
dimdTable.Tag = 'dimdTable';
dimdTable.Units = 'normalized';
dimdTable.Position = [0.5,0,0.5,1];
dimdTable.FontSize = 13;
dimdTable.ColumnEditable = [false,false,false,false]; 
dimdTable.ColumnName = {'Parameter','Value','Unit','Description'}; 

% setting data
index = 0;

index = index + 1;
dimdTable.Data{index,1} = 'Rro';
dimdTable.Data{index,2} = hObject.UserData.sdclass.(dimdTable.Data{index,1});
dimdTable.Data{index,3} = '[mm]';
dimdTable.Data{index,4} = 'Outer Rotor Radius';

index = index + 1;
dimdTable.Data{index,1} = 'rsh';
dimdTable.Data{index,2} = hObject.UserData.sdclass.(dimdTable.Data{index,1});
dimdTable.Data{index,3} = '[mm]';
dimdTable.Data{index,4} = 'Shaft Radius';

index = index + 1;
dimdTable.Data{index,1} = 'Ts';
dimdTable.Data{index,2} = hObject.UserData.sdclass.(dimdTable.Data{index,1});
dimdTable.Data{index,3} = '[rad]';
dimdTable.Data{index,4} = 'Stator Pole Pitch';

index = index + 1;
dimdTable.Data{index,1} = 'Tr';
dimdTable.Data{index,2} = hObject.UserData.sdclass.(dimdTable.Data{index,1});
dimdTable.Data{index,3} = '[rad]';
dimdTable.Data{index,4} = 'Rotor Pole Pitch';

index = index + 1;
dimdTable.Data{index,1} = 'Tso1';
dimdTable.Data{index,2} = hObject.UserData.sdclass.(dimdTable.Data{index,1});
dimdTable.Data{index,3} = '[rad]';
dimdTable.Data{index,4} = 'Slot Openning';

index = index + 1;
dimdTable.Data{index,1} = 'rss1';
dimdTable.Data{index,2} = hObject.UserData.sdclass.(dimdTable.Data{index,1});
dimdTable.Data{index,3} = '[mm]';
dimdTable.Data{index,4} = 'Shaft Radius';

index = index + 1;
dimdTable.Data{index,1} = 'rss2';
dimdTable.Data{index,2} = hObject.UserData.sdclass.(dimdTable.Data{index,1});
dimdTable.Data{index,3} = '[mm]';
dimdTable.Data{index,4} = 'Shaft Radius';

index = index + 1;
dimdTable.Data{index,1} = 'phi_tot';
dimdTable.Data{index,2} = hObject.UserData.sdclass.(dimdTable.Data{index,1});
dimdTable.Data{index,3} = '[wb]';
dimdTable.Data{index,4} = 'Total air gap flux';

index = index + 1;
dimdTable.Data{index,1} = 'phi_p';
dimdTable.Data{index,2} = hObject.UserData.sdclass.(dimdTable.Data{index,1});
dimdTable.Data{index,3} = '[wb]';
dimdTable.Data{index,4} = 'Pole flux';

index = index + 1;
dimdTable.Data{index,1} = 'Aslot';
dimdTable.Data{index,2} = hObject.UserData.sdclass.(dimdTable.Data{index,1});
dimdTable.Data{index,3} = '[mm^2]';
dimdTable.Data{index,4} = 'Slot Area';

%% windings tab
winTab = uitab(tabgp,'Title','Winding');
winTab.Tag = 'winTab';

% independent
winiTable = uitable(winTab);
winiTable.Tag = 'winiTable';
winiTable.Units = 'normalized';
winiTable.Position = [0,0,0.5,1];
winiTable.FontSize = 13;
winiTable.ColumnEditable = [false,true,false,false]; 
winiTable.ColumnName = {'Parameter','Value','Unit','Description'}; 

% setting table
index = 0;

index = index + 1;
winiTable.Data{index,1} = 'a';
winiTable.Data{index,2} = hObject.UserData.sdclass.(winiTable.Data{index,1});
winiTable.Data{index,3} = '';
winiTable.Data{index,4} = 'Number of parallel path';

index = index + 1;
winiTable.Data{index,1} = 'Ncoil';
winiTable.Data{index,2} = hObject.UserData.sdclass.(winiTable.Data{index,1});
winiTable.Data{index,3} = '';
winiTable.Data{index,4} = 'Number of coil turns';

index = index + 1;
winiTable.Data{index,1} = 'Span';
winiTable.Data{index,2} = hObject.UserData.sdclass.(winiTable.Data{index,1});
winiTable.Data{index,3} = '';
winiTable.Data{index,4} = 'Coil Span';

index = index + 1;
winiTable.Data{index,1} = 'Kf';
winiTable.Data{index,2} = hObject.UserData.sdclass.(winiTable.Data{index,1});
winiTable.Data{index,3} = '';
winiTable.Data{index,4} = 'Fill factor';

index = index + 1;
winiTable.Data{index,1} = 'sigma';
winiTable.Data{index,2} = hObject.UserData.sdclass.(winiTable.Data{index,1});
winiTable.Data{index,3} = '[S/m]';
winiTable.Data{index,4} = 'Electrical conductivity of winding material at 25[C]';

index = index + 1;
winiTable.Data{index,1} = 'Lend';
winiTable.Data{index,2} = hObject.UserData.sdclass.(winiTable.Data{index,1});
winiTable.Data{index,3} = '[mm]';
winiTable.Data{index,4} = '';


index = index + 1;
winiTable.Data{index,1} = 'winT';
winiTable.Data{index,2} = hObject.UserData.sdclass.(winiTable.Data{index,1});
winiTable.Data{index,3} = '[C]';
winiTable.Data{index,4} = 'Average winding termerature';

index = index + 1;
winiTable.Data{index,1} = 'alphaT';
winiTable.Data{index,2} = hObject.UserData.sdclass.(winiTable.Data{index,1});
winiTable.Data{index,3} = '[1/K]';
winiTable.Data{index,4} = 'Temperature coefficient of resistivity';

index = index + 1;
winiTable.Data{index,1} = 'kR';
winiTable.Data{index,2} = hObject.UserData.sdclass.(winiTable.Data{index,1});
winiTable.Data{index,3} = '';
winiTable.Data{index,4} = 'AC Resistance Factor';

% dependent
windTable = uitable(winTab);
windTable.Tag = 'windTable';
windTable.Units = 'normalized';
windTable.Position = [0.5,0,0.5,1];
windTable.FontSize = 13;
windTable.ColumnEditable = [false,false,false,false]; 
windTable.ColumnName = {'Parameter','Value','Unit','Description'}; 

index = 0;

index = index + 1;
windTable.Data{index,1} = 'Rdc_coil';
windTable.Data{index,2} = hObject.UserData.sdclass.(windTable.Data{index,1});
windTable.Data{index,3} = '[Ohm]';
windTable.Data{index,4} = 'Coil Resistance';

index = index + 1;
windTable.Data{index,1} = 'Rdc_coil_T';
windTable.Data{index,2} = hObject.UserData.sdclass.(windTable.Data{index,1});
windTable.Data{index,3} = '[Ohm]';
windTable.Data{index,4} = 'Coil Resistance At Nominal Temperature';

index = index + 1;
windTable.Data{index,1} = 'Rac_coil_T';
windTable.Data{index,2} = hObject.UserData.sdclass.(windTable.Data{index,1});
windTable.Data{index,3} = '[Ohm]';
windTable.Data{index,4} = 'Coil Resistance at Nominal Temperature and Frequency';

index = index + 1;
windTable.Data{index,1} = 'Rph';
windTable.Data{index,2} = hObject.UserData.sdclass.(windTable.Data{index,1});
windTable.Data{index,3} = '[Ohm]';
windTable.Data{index,4} = 'Phase Resistance';

index = index + 1;
windTable.Data{index,1} = 'CopperLoss';
windTable.Data{index,2} = hObject.UserData.sdclass.(windTable.Data{index,1});
windTable.Data{index,3} = '[W]';
windTable.Data{index,4} = 'Copper Loss';

%% iron loss tab
ilossTab = uitab(tabgp,'Title','Iron Loss');
ilossTab.Tag = 'ilossTab';

ilossTable = uitable(ilossTab);
ilossTable.Tag = 'ilossTable';
ilossTable.Units = 'normalized';
ilossTable.Position = [0,0,1,1];
ilossTable.FontSize = 13;
ilossTable.ColumnEditable = [false,true,false,false]; 
ilossTable.ColumnName = {'Parameter','Value','Unit','Description'}; 

index = 0;

index = index + 1;
ilossTable.Data{index,1} = 'Kh';
ilossTable.Data{index,2} = 6.673865e-3;
ilossTable.Data{index,3} = '[W/Kg]';
ilossTable.Data{index,4} = 'Hysteresis Loss Coefficient';

index = index + 1;
ilossTable.Data{index,1} = 'alpha';
ilossTable.Data{index,2} = 1.2916;
ilossTable.Data{index,3} = '';
ilossTable.Data{index,4} = 'Power of Frequency in Hysteresis Term of Iron Loss Model';

index = index + 1;
ilossTable.Data{index,1} = 'Ke';
ilossTable.Data{index,2} = 66.156534e-6;
ilossTable.Data{index,3} = '[W/Kg]';
ilossTable.Data{index,4} = 'Eddy Current Loss Coefficient';

index = index + 1;
ilossTable.Data{index,1} = 'Density';
ilossTable.Data{index,2} = 7650;
ilossTable.Data{index,3} = '[Kg/m^3]';
ilossTable.Data{index,4} = 'Density of Electrical Steel';

%% mechanical loss tab
mlossTab = uitab(tabgp,'Title','Mechanical Loss');
mlossTab.Tag = 'mlossTab';

mlossTable = uitable(mlossTab);
mlossTable.Tag = 'mlossTable';
mlossTable.Units = 'normalized';
mlossTable.Position = [0,0,1,1];
mlossTable.FontSize = 13;
mlossTable.ColumnEditable = [false,true,false,false]; 
mlossTable.ColumnName = {'Parameter','Value','Unit','Description'}; 

%% simulation parameters tab
simTab = uitab(tabgp,'Title','Simulation Parameters');
simTab.Tag = 'simTab';

simTabgp = uitabgroup(simTab);

solTab = uitab(simTabgp,'Title','Solver');
solTab.Tag = 'solTab';
meshTab = uitab(simTabgp,'Title','Mesh');
meshTab.Tag = 'meshTab';

solTable = uitable(solTab);
solTable.Tag = 'solTable';

% setting table
solTable.Units = 'normalized';
solTable.Position = [0,0,1,1];
solTable.FontSize = 13;
solTable.ColumnEditable = [false,true,false,false]; 
solTable.ColumnName = {'Parameter','Value','Unit','Description'}; 

index = 0;

index = index + 1;
solTable.Data{index,1} = 'Nsim';
solTable.Data{index,2} = 20;
solTable.Data{index,3} = '';
solTable.Data{index,4} = 'Number of magneto-static simulations';

meshTable = uitable(meshTab);
meshTable.Tag = 'meshTable';

% setting table
meshTable.Units = 'normalized';
meshTable.Position = [0,0,1,1];
meshTable.FontSize = 13;
meshTable.ColumnEditable = [false,true,false,false]; 
meshTable.ColumnName = {'Parameter','Value','Unit','Description'}; 

index = 0;

index = index + 1;
meshTable.Data{index,1} = 'thetar';
meshTable.Data{index,2} = 4;
meshTable.Data{index,3} = '[Degree]';
meshTable.Data{index,4} = 'Mesh angle at inner surface of rotor';

index = index + 1;
meshTable.Data{index,1} = 'thetag';
meshTable.Data{index,2} = 1;
meshTable.Data{index,3} = '[Degree]';
meshTable.Data{index,4} = 'Mesh angle at air gap';

index = index + 1;
meshTable.Data{index,1} = 'thetas';
meshTable.Data{index,2} = 3;
meshTable.Data{index,3} = '[Degree]';
meshTable.Data{index,4} = 'Mesh angle at outer surface of stator';

%% output data tab
outTab = uitab(tabgp,'Title','Outputs');
outTab.Tag = 'outTab';
outTable = uitable(outTab);
outTable.Units = 'normalized';
outTable.Position = [0,0,1,1];

%% veiw tab
viewTab = uitab(tabgp,'Title','View');
viewTab.Tag = 'viewTab';
viewAxis = axes(viewTab);
viewAxis.Tag = 'viewAxis';

%% Msg Box
MsgBox = uicontrol(hObject,'Style','listbox');
MsgBox.Tag = 'MsgBox';
MsgBox.Units = 'normalized';
MsgBox.Position = [0.0025,0,0.995,0.15];
MsgBox.FontSize = 13;
MsgBox.String{1} = 'EMDLAB:';
MsgBox.String{2} = 'An application for loss analysis of PMSMs';

dimiTable.CellEditCallback = @dimiTable_CellEditCallback;
coniTable.CellEditCallback = @coniTable_CellEditCallback; 
winiTable.CellEditCallback = @winiTable_CellEditCallback;
meshTable.CellEditCallback = @meshTable_CellEditCallback;

% Choose default command line output for LossAnalysis_PMSM
% Update handles structure
handles = guihandles(hObject);
handles.output = hObject;
guidata(hObject, handles);

%% setting cmenus

dimiTable.UIContextMenu = handles.dim_cmenu;
dimdTable.UIContextMenu = handles.dim_cmenu;

windTable.UIContextMenu = handles.win_cmenu;
winiTable.UIContextMenu = handles.win_cmenu;

%% setting data
ImportDataOnOpen(handles)
Read_sdclass(handles)
Write_sdclass(handles)


% --- Outputs from this function are returned to the command line.
function varargout = LossAnalysis_PMSM_OutputFcn(hObject, eventdata, handles) 
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
function Untitled_9_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_10_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_10 (see GCBO)
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


% --------------------------------------------------------------------
function Untitled_8_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_8 (see GCBO)
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
function dim_cmenu_Callback(hObject, eventdata, handles)
% hObject    handle to dim_cmenu (see GCBO)
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
function tlocation_Top_Callback(hObject, eventdata, handles)
% hObject    handle to tlocation_Top (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.tlocation_Top.Checked = 'on';
handles.tlocation_Left.Checked = 'off';
handles.tlocation_Down.Checked = 'off';
handles.tabgp.TabLocation = 'top';

% --------------------------------------------------------------------
function tlocation_Left_Callback(hObject, eventdata, handles)
% hObject    handle to tlocation_Left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.tlocation_Top.Checked = 'off';
handles.tlocation_Left.Checked = 'on';
handles.tlocation_Down.Checked = 'off';
handles.tabgp.TabLocation = 'left';

% --------------------------------------------------------------------
function tlocation_Down_Callback(hObject, eventdata, handles)
% hObject    handle to tlocation_Down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.tlocation_Top.Checked = 'off';
handles.tlocation_Left.Checked = 'off';
handles.tlocation_Down.Checked = 'on';
handles.tabgp.TabLocation = 'bottom';

% --------------------------------------------------------------------
function Untitled_22_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CreateGeom(handles)

% --------------------------------------------------------------------
function Untitled_23_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isSetGeom(handles)
    cla
    handles.main.UserData.Geometry.plotfs;
    handles.tabgp.SelectedTab = handles.viewTab;
end

function WriteMsg(handles,type,txt)

handles.MsgBoxBox.String = '';

switch type
    case 0
        handles.MsgBox.ForegroundColor = 'r';
        handles.MsgBox.String{1} = 'Error:';
    case 1
        handles.MsgBox.ForegroundColor = 'b';
        handles.MsgBox.String{1} = 'Done!';
    case 2
        handles.MsgBox.ForegroundColor = [249,124,0]/255;
        handles.MsgBox.String{1} = 'Warning:';
    case 3
        handles.MsgBox.ForegroundColor = 'r';
        handles.MsgBox.String{1} = 'Busy:';
    case 4
        handles.MsgBox.String{1} = 'Results:';
        handles.MsgBox.ForegroundColor = 'k';
    case 5
        handles.MsgBox.String{1} = 'Message:';
        handles.MsgBox.ForegroundColor = 'b';
end
 
if isa(txt,'cell')
    handles.MsgBox.String = {handles.MsgBox.String(:),txt(:)};
else
    handles.MsgBox.String{2} = txt;
end

function ExportData(handles)

[FileName,PathName] = uiputfile('*.pmsm','Enter Your File Name',[cd,'\project.pmsm']);

if ~FileName
    return
end

f = fopen([PathName,FileName],'w');

tables = {'dimiTable','dimdTable','winiTable','windTable','ilossTable','coniTable','meshTable'};
for i = 1:length(tables)
    table = tables{i};
    for j = 1:size(handles.(table).Data,1)
        fprintf(f,'%s\t%f\n',handles.(table).Data{j,1},handles.(table).Data{j,2});
    end
end

fclose(f);
WriteMsg(handles,1,'Finish')

function ExportDataOnClose(handles)

f = fopen([cd,'\Internal.pmsm'],'w');

tables = {'dimiTable','dimdTable','winiTable','windTable','ilossTable','coniTable','meshTable'};
for i = 1:length(tables)
    table = tables{i};
    for j = 1:size(handles.(table).Data,1)
        fprintf(f,'%s\t%f\n',handles.(table).Data{j,1},handles.(table).Data{j,2});
    end
end

fclose(f);
WriteMsg(handles,1,'Finish')

function ImportData(handles)

[FileName,PathName] = uigetfile('*.pmsm');
if ~FileName
    return
end

f = fopen([PathName,FileName],'r');
tables = {'dimiTable','winiTable','ilossTable','coniTable','meshTable'};

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
fclose(f);

WriteMsg(handles,1,'Finish')
hObject.UserData.states.isSetMesh = false;
hObject.UserData.states.isSetWin = false;

function ImportDataOnOpen(handles)
f = fopen([cd,'\LastDesignData.txt'],'r');
tables = {'dimiTable','winiTable','ilossTable','coniTable','meshTable'};

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
fclose(f);

WriteMsg(handles,1,'Finish')
hObject.UserData.states.isSetMesh = false;
hObject.UserData.states.isSetWin = false;

% --------------------------------------------------------------------
function Untitled_24_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ImportData(handles)

% --------------------------------------------------------------------
function Untitled_25_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ExportData(handles)


function y = isSetMesh(handles)

if handles.main.UserData.states.isSetMesh
     y = true;
else
    y = false;
    WriteMsg(handles,0,'Mesh does not generated ...');
end
    

% --------------------------------------------------------------------
function Untitled_26_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isSetGeom(handles)
    cla
    handles.main.UserData.Geometry.plotwf;
    handles.tabgp.SelectedTab = handles.viewTab;
end

% --------------------------------------------------------------------
function Untitled_27_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isSetGeom(handles)
    cla
    handles.main.UserData.Geometry.plotim;
    handles.tabgp.SelectedTab = handles.viewTab;
end

function GenerateMesh(handles)

if handles.main.UserData.sdclass.checkIterdependencies(handles.MsgBox)
    return;
end

WaitFig = waitbar(0/8,'Creation of Geometry And Mesh',...
    'WindowStyle','modal',...
    'Name','MG0');

materialdir = [cd,'\MaterialsData'];
p = handles.main.UserData.sdclass;

waitbar(1/8,WaitFig,'Writing Par File');
p.writeParFile;

waitbar(2/8,WaitFig,'Creation of Geometry');
geom_pmsm0;

m = TMDBC;
waitbar(3/8,WaitFig,'Calling MG0');
m.read_g2d_bin('geom.g2d');

% adding needed materials
m.addMaterial(materialdir,'air');
m.addMaterial(materialdir,'m19_24ga');
m.addMaterial(materialdir,'copper');

% setting materials
m.setMaterial('s1','m19_24ga');
m.setMaterial('c11','copper');
m.setMaterial('r1','m19_24ga');
m.setmzColor('rap1','w');
m.setmzColor('sap1','w');
m.setmzColor('s1',[133, 193, 233]/255);
m.setmzColor('r1',[133, 193, 233]/255);

waitbar(4/8,WaitFig,'Creation of Stator Mesh Zone');
m.cmirrormz('s2','s1',[1,0]);
for i = 1:2:2*(p.Ns-1)
    m.crotatemz(['s',num2str(i+2)],['s',num2str(i)],p.Ts);
    m.crotatemz(['s',num2str(i+3)],['s',num2str(i+1)],p.Ts);
end
temp = getlist('s',1:2*p.Ns);
m.joinmzs('stator',temp{:});

waitbar(5/8,WaitFig,'Creation of Coil Mesh Zones');
m.cmirrormz('c21','c11',[1,0]);
m.setmzColor('c11',[241,195,15]/255);
m.setmzColor('c21','c');
for i = 1:p.Ns-1
    m.crotatemz(['c1',num2str(i+1)],['c1',num2str(i)],p.Ts);
    m.crotatemz(['c2',num2str(i+1)],['c2',num2str(i)],p.Ts);
    m.crotatemz(['sap',num2str(i+1)],['sap',num2str(i)],p.Ts);
end

waitbar(6/8,WaitFig,'Creation of Rotor Mesh Zone');
m.cmirrormz('r2','r1',[1,0]);
for i = 1:2:2*(p.Nm-1)
    m.crotatemz(['r',num2str(i+2)],['r',num2str(i)],p.Tr);
    m.crotatemz(['r',num2str(i+3)],['r',num2str(i+1)],p.Tr);
end
temp = getlist('r',1:2*p.Nm);
m.joinmzs('rotor',temp{:});

waitbar(7/8,WaitFig,'Creation of Magnets Mesh Zones');
m.cmirrormz('m2','m1',[1,0]);
m.joinmzs('magnet1','m1','m2');
m.setmzColor('magnet1',[195, 155, 211]/255)
for i = 1:p.Nm-1
    m.crotatemz(['magnet',num2str(i+1)],['magnet',num2str(i)],p.Tr);
    m.crotatemz(['rap',num2str(i+1)],['rap',num2str(i)],p.Tr);
end

waitbar(1,WaitFig,'Setting Data For Global Mesh');
m.ggmesh;
close(WaitFig);

handles.main.UserData.Mesh = m;
WriteMsg(handles,1,'Mesh Generation Completed!')
handles.main.UserData.states.isSetMesh = true;


% --------------------------------------------------------------------
function Untitled_28_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GenerateMesh(handles);

% --------------------------------------------------------------------
function Untitled_29_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isSetMesh(handles)
    cla
    handles.main.UserData.Mesh.showmzsf;
    handles.tabgp.SelectedTab = handles.viewTab;
    set(gcf,'color',[70, 130, 180]/255);
end


function dimiTable_CellEditCallback(hObject, eventdata, handles)
handles = guidata(hObject);
if handles.main.UserData.sdclass.setVariable(...
        hObject.Data{eventdata.Indices(1),1},...
        eventdata.NewData,handles.MsgBox)
    hObject.Data{eventdata.Indices(1),2} = eventdata.PreviousData;
end
Write_sdclass(handles)

function coniTable_CellEditCallback(hObject, eventdata, handles)
handles = guidata(hObject);
if handles.main.UserData.sdclass.setVariable(...
        hObject.Data{eventdata.Indices(1),1},...
        eventdata.NewData,handles.MsgBox)
    hObject.Data{eventdata.Indices(1),2} = eventdata.PreviousData;
end
Write_sdclass(handles)

function winiTable_CellEditCallback(hObject, eventdata, handles)
handles = guidata(hObject);
if handles.main.UserData.sdclass.setVariable(...
        hObject.Data{eventdata.Indices(1),1},...
        eventdata.NewData,handles.MsgBox)
    hObject.Data{eventdata.Indices(1),2} = eventdata.PreviousData;
end
Write_sdclass(handles)

function meshTable_CellEditCallback(hObject, eventdata, handles)
handles = guidata(hObject);
if handles.main.UserData.sdclass.setVariable(...
        hObject.Data{eventdata.Indices(1),1},...
        eventdata.NewData,handles.MsgBox)
    hObject.Data{eventdata.Indices(1),2} = eventdata.PreviousData;
end

function Read_sdclass(handles)
tables = {'dimiTable','coniTable','winiTable'};
for t = tables
    for i = 1:size(handles.(t{1}).Data,1)
        handles.main.UserData.sdclass.setVariable(handles.(t{1}).Data{i,1},...
            handles.(t{1}).Data{i,2});
    end
end

function Write_sdclass(handles)
tables = {'dimdTable','dimiTable','condTable','coniTable','winiTable','windTable'};
for t = tables
    for i = 1:size(handles.(t{1}).Data,1)
        handles.(t{1}).Data{i,2} = ...
            handles.main.UserData.sdclass.(handles.(t{1}).Data{i,1});
    end
end


% --------------------------------------------------------------------
function Untitled_30_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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
function Untitled_36_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_36 (see GCBO)
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
function Set_Loss_Data_For_Callback(hObject, eventdata, handles)
% hObject    handle to Set_Loss_Data_For (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_37_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Read_sdclass(handles)
handles.main.UserData.sdclass.update_wtb_For_Bt;
Write_sdclass(handles)

% --------------------------------------------------------------------
function Untitled_38_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Read_sdclass(handles)
handles.main.UserData.sdclass.update_wsy_For_Bsy;
Write_sdclass(handles)
% --------------------------------------------------------------------
function Untitled_39_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Read_sdclass(handles)
handles.main.UserData.sdclass.update_wry_For_Bry;
Write_sdclass(handles)


% --------------------------------------------------------------------
function Untitled_40_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Read_sdclass(handles)
handles.main.UserData.sdclass.update_wtb_For_Bt;
handles.main.UserData.sdclass.update_wsy_For_Bsy;
handles.main.UserData.sdclass.update_wry_For_Bry;
handles.main.UserData.sdclass.update_geom;
Write_sdclass(handles)


% --- Executes when user attempts to close main.
function main_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ExportDataOnClose(handles)
% Hint: delete(hObject) closes the figure
delete(hObject);


% --------------------------------------------------------------------
function Untitled_42_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_43_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Read_sdclass(handles)
handles.main.UserData.sdclass.update_Ncoil_For_BEMF;
Write_sdclass(handles)

% --------------------------------------------------------------------
function win_cmenu_Callback(hObject, eventdata, handles)
% hObject    handle to win_cmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function Untitled_44_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function Untitled_45_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isSetMesh(handles)
    cla
    handles.main.UserData.Mesh.showmzs;
    handles.tabgp.SelectedTab = handles.viewTab;
    set(gcf,'color',[70, 130, 180]/255);
end
