function f = GraphicWindowUPD()

% main figure
f = figure('Visible', 'off');
set(f, 'position', [0,0,800,600]);
movegui(f, 'center');
set(f, 'Menu', 'none');
set(f, 'ToolBar', 'none');
set(f, 'Tag', 'main');
set(f, 'Name', 'EMDLab');
set(f, 'NumberTitle', 'off');
set(f, 'Clipping', 'off');
set(f, 'color', [0.2,0.2,0.2]);
set(f, 'Renderer', 'opengl');
set(f, 'Units', 'pixel');
f.UserData.cp = [0,0];
f.UserData.isMouseHoled = false;
f.UserData.isMouseRight = false;
f.UserData.isMouseLeft = false;

% view axes
va = axes(f, 'Tag', 'va');
va.Units = 'normalized';
va.Position = [0,0,1,1];
camzoom(va, 0.8);
va.Clipping = 'off';
axis(va,'off','vis3d');

% setting background
bg = axes(f, 'Tag', 'bg');
bg.Units = 'normalized';
bg.Position = [0,0,1,1];
bg.HandleVisibility = 'off';
bg.Clipping = 'off';
uistack(bg,'bottom');
SetBackGround()

% coordinate axes
ca = axes(f, 'Tag', 'ca');
ca.Units = 'pixels';
ca.Position = [20,20,100,100];
axis(ca, 'off', 'vis3d')
ca.NextPlot = 'add';
arrow3d([0,1],[0,0],[0,0],0.85,0.045,0.1,'r');
arrow3d([0,0],[0,1],[0,0],0.85,0.045,0.1,'b');
arrow3d([0,0],[0,0],[0,1],0.85,0.045,0.1,'g');
patch(ca,'faces',[1,2;3,4;5,6],'vertices',[-1,0,0;1,0,0;0,-1,0;0,1,0;0,0,-1;0,0,1],'EdgeColor','none')
text(ca,1.1,0,0,'X', 'color', 'r', 'fontsize', 10);
text(ca,0,1.1,0,'Y', 'color', 'b', 'fontsize', 10);
text(ca,0,0,1.1,'Z', 'color', [0,0.6,0], 'fontsize', 10);

% figure callbacks
f.WindowButtonDownFcn = @fButtonDownFcn;
f.WindowButtonMotionFcn = @fButtonMotionFcn;
f.WindowButtonUpFcn = @fButtonUpFcn;
f.SizeChangedFcn = @fResizeFcn;
f.WindowKeyPressFcn = @fKeyFcn;
f.WindowScrollWheelFcn = @fScrollFcn;

% saving gui components
set(f, 'HandleVisibility', 'off');
guidata(f, guihandles(f));
f.Visible = 'on';

%%%%%%%%%% figure callback
  function fButtonDownFcn(~,~)
    f.UserData.isMouseHoled = true;
    f.UserData.cp = f.CurrentPoint;
    if strcmpi(f.SelectionType,'alt')
      f.UserData.isMouseRight = true;
      f.UserData.isMouseLeft = false;
      f.Pointer = 'fleur';
    elseif strcmpi(f.SelectionType,'normal')
      f.UserData.isMouseRight = false;
      f.UserData.isMouseLeft = true;
    else
      f.UserData.isMouseRight = false;
      f.UserData.isMouseLeft = false;
    end
  end

%%%%%%%%%% figure callback
  function fButtonMotionFcn(~,~)
    if f.UserData.isMouseHoled
      if f.UserData.isMouseRight
        tmp = f.CurrentPoint - f.UserData.cp;
        frame = f.Position;
        camdolly(va, -tmp(1)/(frame(3)-frame(1)),-tmp(2)/(frame(4)-frame(2)), 0);
        drawnow;
        f.UserData.cp = f.CurrentPoint;
      elseif f.UserData.isMouseLeft
        %       else
      end
    end
  end

%%%%%%%%%% figure callback
  function fButtonUpFcn(~,~)
    f.UserData.isMouseHoled = false;
    if f.UserData.isMouseRight
      f.Pointer = 'arrow';
    elseif f.UserData.isMouseLeft
    else
    end
  end

%%%%%%%%%% figure callback
  function fResizeFcn(~,~)
  end

%%%%%%%%%% figure callback
  function fScrollFcn(~, eventData)
    if eventData.VerticalScrollCount<0
      camzoom(va, 0.8);
    else
      camzoom(va, 1.2);
    end
  end

%%%%%%%%%% figure callback
  function fKeyFcn(~,eventData)
    % number of modifires
    m = eventData.Modifier;
    nMod = numel(m);
    % right rotation
    if strcmpi(eventData.Key,'rightarrow') && ~nMod
      camorbit(ca, -10,0,'camera');
      camorbit(va, -10,0,'camera');
      return
    end
    % left rotation
    if strcmpi(eventData.Key,'leftarrow') && ~nMod
      camorbit(ca, 10,0,'camera');
      camorbit(va, 10,0,'camera');
      return
    end
    % up rotation
    if strcmpi(eventData.Key,'uparrow') && ~nMod
      camorbit(ca, 0,-10,'camera');
      camorbit(va, 0,-10,'camera');
      return
    end
    % down rotation
    if strcmpi(eventData.Key,'downarrow') && ~nMod
      camorbit(ca, 0,10,'camera');
      camorbit(va, 0,10,'camera');
      return
    end
    % zoom in
    if strcmpi(eventData.Key,'uparrow') && nMod==1 && strcmpi(m{1},'shift')
      camzoom(va, 1.1);
      return
    end
    % zoom out
    if strcmpi(eventData.Key,'downarrow') && nMod==1 && strcmpi(m{1},'shift')
      camzoom(va, 0.9);
      return
    end
    % anti-clockwise rotate
    if strcmpi(eventData.Key,'rightarrow') && nMod==1 && strcmpi(m{1},'shift')
      camroll(va, -5);
      camroll(ca, -5);
      return
    end
    % clockwise rotate
    if strcmpi(eventData.Key,'leftarrow')&&nMod==1&&strcmpi(m{1},'shift')
      camroll(va, 5);
      camroll(ca, 5);
      return
    end
    % right shift
    if strcmpi(eventData.Key,'rightarrow')&&nMod==1&&strcmpi(m{1},'control')
      camdolly(va,-0.1,0,0)
      return
    end
    % left shift
    if strcmpi(eventData.Key,'leftarrow')&&nMod==1&&strcmpi(m{1},'control')
      camdolly(va,0.1,0,0)
      return
    end
    % up shift
    if strcmpi(eventData.Key,'uparrow')&&nMod==1&&strcmpi(m{1},'control')
      camdolly(va,0,-0.1,0);
      return
    end
    % down shift
    if strcmpi(eventData.Key,'downarrow')&&nMod==1&&strcmpi(m{1},'control')
      camdolly(va,0,0.1,0)
      return
    end
    % z axis view
    if nMod==0 && strcmpi(eventData.Key,'z')
      view(va,[0,0,1]);
      view(ca,[0,0,1]);
      return
    end
    % x axis view
    if nMod==0 && strcmpi(eventData.Key,'x')
      view(va,[1,0,0]);
      view(ca,[1,0,0]);
      return
    end
    % y axis view
    if nMod==0 && strcmpi(eventData.Key,'y')
      view(va,[0,1,0]);
      view(ca,[0,1,0]);
      return
    end
    % isotropic rotation
    if nMod==0 && strcmpi(eventData.Key,'i')
      camorbit(va, 30,30,'camera');
      camorbit(ca, 30,30,'camera');
      return
    end
    % 90 degree rotation
    if nMod==0 && strcmpi(eventData.Key,'space')
      camroll(va,90)
      camroll(ca,90)
      return
    end
    % saving figure
    if strcmpi(eventData.Key,'s')&&nMod==1&&strcmpi(m{1},'control')
      bg.Visible = 'off';
      saveas(va,'view.jpg');
      bg.Visible = 'on';
      return
    end
  end

%%%%%%%%%% setting background image
  function SetBackGround()
    % number of vertical strip
    n = 250;
    % dividing figure into three parts
    d = n/3;
    % interpolation with a quadratic function
    % first coefficient
    rgb = [1;0.8;1];
    c1 = [1,1,1;d^2,d,1;n^2,n,1]\rgb;
    % second coefficient
    rgb = [1;0.81;1];
    c2 = [1,1,1;d^2,d,1;n^2,n,1]\rgb;
    % third coefficient
    rgb = [1;0.82;1];
    c3 = [1,1,1;d^2,d,1;n^2,n,1]\rgb;
    % image pixels
    m = zeros(n,n,3);
    for i = 1:n
      m(i,:,1) = c1(1)*i^2+c1(2)*i+c1(3);
      m(i,:,2) = c2(1)*i^2+c2(2)*i+c2(3);
      m(i,:,3) = c3(1)*i^2+c3(2)*i+c3(3);
    end
    imagesc(bg, m);
    axis(bg, 'off');
    bg.Tag = 'bg';
  end

end
