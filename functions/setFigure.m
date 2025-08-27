function CurrentAxisHandle = setFigure(FigureTitle, tmpbg)

if ~nargin
  FigureTitle = '[@EMDLab]';
end

f = figure('Visible', 'off');
set(f, 'position', [0,0,800,600]);
movegui(f, 'center');
set(f, 'Menu', 'none');
set(f, 'ToolBar', 'none');
set(f, 'Tag', 'main');
set(f, 'Name', FigureTitle);
set(f, 'NumberTitle', 'off');
set(f, 'Clipping', 'off');
set(f, 'color', [0.2,0.2,0.2]);
set(f,'Renderer','opengl');
set(f, 'Units', 'pixel');

f.UserData.isMouseHoled = false;
f.UserData.isMouseRight = false;
f.UserData.isMouseLeft = false;

h = 0.05;
w = 1/15;


axes(f);
set(gca, 'Tag', 'ca');
cor = axes(f, 'Tag', 'cor');
cor.Units = 'pixels';
cor.Position = [20,20,100,100];
axis(gca, 'off', 'vis3d')
hold all
% quiver3(0,0,0,1,0,0, 'LineWidth', 2);
% quiver3(0,0,0,0,1,0, 'LineWidth', 2);
% quiver3(0,0,0,0,0,1, 'LineWidth', 2);
emdlab_flib_arrow3d([0,1],[0,0],[0,0],0.85,0.045,0.1,'r');
emdlab_flib_arrow3d([0,0],[0,1],[0,0],0.85,0.045,0.1,'b');
emdlab_flib_arrow3d([0,0],[0,0],[0,1],0.85,0.045,0.1,'g');
patch('faces',[1,2;3,4;5,6],'vertices',[-1,0,0;1,0,0;0,-1,0;0,1,0;0,0,-1;0,0,1],'EdgeColor','none')
text(cor,1.1,0,0,'X', 'color', 'r', 'fontsize', 10);
text(cor,0,1.1,0,'Y', 'color', 'b', 'fontsize', 10);
text(cor,0,0,1.1,'Z', 'color', [0,0.6,0], 'fontsize', 10);

if nargin < 2
axes(f,'Tag','bg');
end
h = guihandles(f);
% setting figure background
guidata(f, h);
if nargin < 2
  SetBackGround(h)
end


% figure callbacks
f.WindowButtonDownFcn = @fButtonDownFcn;
f.WindowButtonMotionFcn = @fButtonMotionFcn;
f.WindowButtonUpFcn = @fButtonUpFcn;
f.SizeChangedFcn = @resizeFcn;
f.WindowKeyPressFcn = @fKeyFcn;
f.WindowScrollWheelFcn = @scrollFcn;



% f.HandleVisibility = 'off';
cor.HandleVisibility = 'off';

% setting current axis
camzoom(h.ca, 0.8);
h.ca.Clipping = 'off';
axis(h.ca,'off','vis3d');

h.ca.UserData.Npz = 0;
h.ca.UserData.Nnz = 0;
h.ca.UserData.Nlm = 0;
h.ca.UserData.Nrm = 0;
h.ca.UserData.Num = 0;
h.ca.UserData.Ndm = 0;

% output
h.ca.Units = 'normalized';
h.ca.Position = [0,0,1,1];
CurrentAxisHandle = h.ca;

function fButtonDownFcn(source, e)
handles = guidata(source);
handles.main.UserData.isMouseHoled = true;
handles.main.UserData.cp = handles.main.CurrentPoint;
if strcmpi(get(source,'selectiontype'),'alt')

 
handles.main.UserData.isMouseRight = true;
  handles.main.UserData.isMouseLeft = false;

handles.main.Pointer = 'fleur';

elseif strcmpi(get(source,'selectiontype'),'normal')
  handles.main.UserData.isMouseRight = false;
  handles.main.UserData.isMouseLeft = true;
else
  handles.main.UserData.isMouseRight = false;
  handles.main.UserData.isMouseLeft = false;
end

function fButtonMotionFcn(source,eventdata)
handles = guidata(source);


if handles.main.UserData.isMouseHoled
if handles.main.UserData.isMouseRight
  tmp = handles.main.CurrentPoint - handles.main.UserData.cp;
  frame = handles.main.Position;
  camdolly(handles.ca,-tmp(1)/(frame(3)-frame(1)),-tmp(2)/(frame(4)-frame(2)), 0);
  drawnow;
  handles.main.UserData.cp = handles.main.CurrentPoint ; 
elseif handles.main.UserData.isMouseLeft
  
   
else
%   tmp = handles.main.CurrentPoint - handles.main.UserData.cp;
%   frame = handles.main.Position;
%   camorbit(handles.ca,-tmp(1)*20/(frame(3)-frame(1)),-tmp(2)*20/(frame(4)-frame(2)));
%   drawnow;
%   handles.main.UserData.cp = handles.main.CurrentPoint ; 
end
end


function fButtonUpFcn(source,eventdata)
handles = guidata(source);
handles.main.UserData.isMouseHoled = false;
handles.main.Pointer = 'arrow';

% tmp = handles.main.CurrentPoint - handles.main.UserData.cp;
% %   tmp1 = tmp./norm(tmp);
% 
%   camorbit(handles.ca, -tmp(1),-tmp(2),'camera');
%   camorbit(handles.cor, -tmp(1),-tmp(2),'camera');
%   drawnow;
% %   handles.main.UserData.cp = tmp;

function SetBackGround(h)

bg = h.bg;
bg.Units = 'normalized';
bg.Position = [0,0,1,1];
uistack(bg,'bottom');

n = 250;
d = n/3;

v1 = 1;
v2 = 0.8;
v3 = 1;
c1 = [1,1,1;d^2,d,1;n^2,n,1]\[v1;v2;v3];

v1 = 1;
v2 = 0.81;
v3 = 1;
c2 = [1,1,1;d^2,d,1;n^2,n,1]\[v1;v2;v3];

v1 = 1;
v2 = 0.82;
v3 = 1;
c3 = [1,1,1;d^2,d,1;n^2,n,1]\[v1;v2;v3];

m = zeros(n,n,3);
for i = 1:n
  m(i,:,1) = c1(1)*i^2+c1(2)*i+c1(3);
  m(i,:,2) = c2(1)*i^2+c2(2)*i+c2(3);
  m(i,:,3) = c3(1)*i^2+c3(2)*i+c3(3); 
end
imagesc(bg,m);
% pcolor(bg, m);
set(bg,'handlevisibility','off')
set(bg,'visible','off')
% freezeColors(bg)
colormap jet
tmp = linspace(0,1,10)';
colormap(repmat(tmp,1,3))

function fKeyFcn(source,eventdata)
handles = guidata(source);
m = eventdata.Modifier;
nMod = numel(m);

if strcmpi(eventdata.Key,'rightarrow')&&~nMod
  camorbit(handles.ca, -10,0,'camera');
  camorbit(handles.cor, -10,0,'camera');
  return
end
if strcmpi(eventdata.Key,'leftarrow')&&~nMod
  camorbit(handles.ca, 10,0,'camera');
  camorbit(handles.cor, 10,0,'camera');
  return
end
if strcmpi(eventdata.Key,'uparrow')&&~nMod
  camorbit(handles.ca, 0,-10,'camera');
  camorbit(handles.cor, 0,-10,'camera');
  return
end
if strcmpi(eventdata.Key,'downarrow')&&~nMod
  camorbit(handles.ca, 0,10,'camera');
  camorbit(handles.cor, 0,10,'camera');
  return
end
if strcmpi(eventdata.Key,'uparrow')&&nMod==1&&strcmpi(m{1},'shift')
  camzoom(handles.ca, 1.1);
  handles.ca.UserData.Npz = handles.ca.UserData.Npz + 1;
  return
end
if strcmpi(eventdata.Key,'downarrow')&&nMod==1&&strcmpi(m{1},'shift')
  camzoom(handles.ca, 0.9);
  handles.ca.UserData.Nnz = handles.ca.UserData.Nnz + 1;
  return
end
if strcmpi(eventdata.Key,'rightarrow')&&nMod==1&&strcmpi(m{1},'shift')
  camroll(handles.ca, -5);
  camroll(handles.cor, -5);
  return
end
if strcmpi(eventdata.Key,'leftarrow')&&nMod==1&&strcmpi(m{1},'shift')
  camroll(handles.ca, 5);
  camroll(handles.cor, 5);
  return
end
if strcmpi(eventdata.Key,'rightarrow')&&nMod==1&&strcmpi(m{1},'control')
  handles.ca.UserData.Nrm = handles.ca.UserData.Nrm + 1;
  camdolly(handles.ca,-0.1,0,0)
  return
end
if strcmpi(eventdata.Key,'leftarrow')&&nMod==1&&strcmpi(m{1},'control')
  handles.ca.UserData.Nlm = handles.ca.UserData.Nlm + 1;
  camdolly(handles.ca,0.1,0,0)
  return
end
if strcmpi(eventdata.Key,'uparrow')&&nMod==1&&strcmpi(m{1},'control')
  handles.ca.UserData.Num = handles.ca.UserData.Num + 1;
  camdolly(handles.ca,0,-0.1,0);
  return
end
if strcmpi(eventdata.Key,'downarrow')&&nMod==1&&strcmpi(m{1},'control')
  handles.ca.UserData.Ndm = handles.ca.UserData.Ndm + 1;
  camdolly(handles.ca,0,0.1,0)
  return
end
if nMod==0&&strcmpi(eventdata.Key,'z')
  view(handles.ca,[0,0,1]);
  view(handles.cor,[0,0,1]);
  return
end
if nMod==0&&strcmpi(eventdata.Key,'x')
  view(handles.ca,[1,0,0]);
  view(handles.cor,[1,0,0]);
  return
end
if nMod==0&&strcmpi(eventdata.Key,'y')
  view(handles.ca,[0,1,0]);
  view(handles.cor,[0,1,0]);
  camroll(handles.ca,-90)
  camroll(handles.cor,-90)
  return
end
if nMod==0&&strcmpi(eventdata.Key,'i')
  camorbit(handles.ca, 30,30,'camera');
  camorbit(handles.cor, 30,30,'camera');
  return
end
if nMod==0&&strcmpi(eventdata.Key,'space')
  camroll(handles.ca,90)
  camroll(handles.cor,90)
% fitView(handles)
  return
end
if strcmpi(eventdata.Key,'s')&&nMod==1&&strcmpi(m{1},'control')
   handles.bg.Children(1).Visible = 'off';
  saveas(handles.ca,'view.jpg');
   handles.bg.Children(1).Visible = 'on';
  return
end

function resizeFcn(source,eventdata)
handles = guidata(source);
% % cb = handles.main
% handles.ColorBar.Units = 'pixels';
% handles.ColorBar.Position(1) = handles.ColorBar.UserData.x1;
% handles.ColorBar.Position(3) = handles.ColorBar.UserData.x2;
% handles.ColorBar.Units = 'normalized';

function fitView(handles)

% camdolly(handles.ca,0.1*handles.ca.UserData.Nrm,0,0)
% handles.ca.UserData.Nrm = 0;
% 
% camdolly(handles.ca,-0.1*handles.ca.UserData.Nlm,0,0)
% handles.ca.UserData.Nlm = 0;
% 
% camdolly(handles.ca,0,0.1*handles.ca.UserData.Num,0)
% handles.ca.UserData.Num = 0;
% 
% camdolly(handles.ca,0,-0.1*handles.ca.UserData.Ndm,0)
% handles.ca.UserData.Ndm = 0;

zoomFactor = (1/0.9)^(handles.ca.UserData.Nnz) * (1/1.1)^(handles.ca.UserData.Npz);
camzoom(handles.ca,zoomFactor)
handles.ca.UserData.Npz = 0;
handles.ca.UserData.Nnz = 0;

function scrollFcn(source,eventdata)
handles = guidata(source);
if eventdata.VerticalScrollCount<0
  camzoom(handles.ca, 0.8);
else
  camzoom(handles.ca, 1.2);
end

