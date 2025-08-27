function [f, ax] = emdlab_r3d_mesh()

% generate a new figure
f = figure('Visible', 'off');
set(f, 'position', [0,0,800,600]);
movegui(f, 'center');
set(f, 'Tag', 'main');
set(f, 'Name', 'EMDLAB: 3D RENDERING');
set(f, 'clipping', 'off');
set(f, 'color', [0.9,0.9,0.9]);
set(f, 'Units', 'pixel');
set(f, 'Renderer', 'opengl');
movegui(f, 'center');

% saving gui components
guidata(f, guihandles(f));

% view axes
ax = axes(f, 'Tag', 'va');
ax.Units = 'normalized';
ax.Position = [0,0,1,1];
camzoom(ax, 0.8);
ax.Clipping = 'off';
axis(ax,'off','vis3d');

% coordinate axes
ca = axes(f, 'Tag', 'ca');
ca.Units = 'pixels';
ca.Position = [20,20,100,100];
axis(ca, 'off', 'vis3d')
ca.NextPlot = 'add';
set(ca,'handlevisibility','off');
% arrow3d([0,1],[0,0],[0,0],0.85,0.045,0.1,'r',ca);
% arrow3d([0,0],[0,1],[0,0],0.85,0.045,0.1,'b',ca);
% arrow3d([0,0],[0,0],[0,1],0.85,0.045,0.1,'g',ca);
% patch(ca,'faces',[1,2;3,4;5,6],'vertices',[-1,0,0;1,0,0;0,-1,0;0,1,0;0,0,-1;0,0,1],'EdgeColor','none')
% text(ca,1.4,0,0,'X', 'color', 'r', 'fontsize', 10);
% text(ca,0,1.4,0,'Y', 'color', 'b', 'fontsize', 10);
% text(ca,0,0,1.4,'Z', 'color', [0,0.6,0], 'fontsize', 10);

f.UserData.cp = [0,0];
f.UserData.isMouseHoled = false;
f.UserData.isMouseRight = false;
f.UserData.isMouseLeft = false;
f.UserData.zf = 1;

% figure callbacks
f.WindowButtonDownFcn = @fButtonDownFcn;
f.WindowButtonMotionFcn = @fButtonMotionFcn;
f.WindowButtonUpFcn = @fButtonUpFcn;
f.SizeChangedFcn = @fResizeFcn;
f.WindowKeyPressFcn = @fKeyFcn;
f.WindowScrollWheelFcn = @fScrollFcn;


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
                if tmp(1) > 0 && tmp(2) > 0
                    camdolly(ax, -10,-10, 0);
                elseif tmp(1) > 0 && tmp(2) < 0
                    camdolly(ax, -10,-10, 0);
                elseif tmp(1) < 0 && tmp(2) > 0
                    camdolly(ax, -10,-10, 0);
                else
                    camdolly(ax, -10,-10, 0);
                end
                drawnow;
                f.UserData.cp = f.CurrentPoint;
            elseif f.UserData.isMouseLeft
                tmp = f.CurrentPoint - f.UserData.cp;
                frame = f.Position;
                camorbit(ax, -50*tmp(1)/(frame(3)-frame(1)),-50*tmp(2)/(frame(4)-frame(2)), 0);
                camorbit(ca, -50*tmp(1)/(frame(3)-frame(1)),-50*tmp(2)/(frame(4)-frame(2)), 0);
                drawnow;
                f.UserData.cp = f.CurrentPoint;
            end
        end
        %     import java.awt.Robot;
        %
        %     import java.awt.event.*;
        %     mouse = Robot;
        %                 mouse.mousePress(InputEvent.BUTTON1_MASK);
        %                 mouse.mouseRelease(InputEvent.BUTTON1_MASK)

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
        handles = guihandles(f);
        if ismember('cba', fieldnames(handles))
            handles.cba.Units = 'normalized';
            handles.cba.Position(2) = 0.2;
            handles.cba.Position(4) = 0.7;
            handles.cba.Units = 'pixels';
            handles.cba.Position(1) = 30;
            handles.cba.Position(3) = 20;
        end
    end

%%%%%%%%%% figure callback
    function fScrollFcn(~, eventData)
        if eventData.VerticalScrollCount<0
            camzoom(ax, 0.8);
            f.UserData.zf = f.UserData.zf*0.8;
        else
            camzoom(ax, 1.2);
            f.UserData.zf = f.UserData.zf*1.2;
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
            camorbit(ax, -10,0,'camera');
            return
        end
        % left rotation
        if strcmpi(eventData.Key,'leftarrow') && ~nMod
            camorbit(ca, 10,0,'camera');
            camorbit(ax, 10,0,'camera');
            return
        end
        % up rotation
        if strcmpi(eventData.Key,'uparrow') && ~nMod
            camorbit(ca, 0,-10,'camera');
            camorbit(ax, 0,-10,'camera');
            return
        end
        % down rotation
        if strcmpi(eventData.Key,'downarrow') && ~nMod
            camorbit(ca, 0,10,'camera');
            camorbit(ax, 0,10,'camera');
            return
        end
        % zoom in
        if strcmpi(eventData.Key,'uparrow') && nMod==1 && strcmpi(m{1},'shift')
            camzoom(ax, 1.1);
            f.UserData.zf = f.UserData.zf*1.1;
            return
        end
        % zoom out
        if strcmpi(eventData.Key,'downarrow') && nMod==1 && strcmpi(m{1},'shift')
            camzoom(ax, 0.9);
            f.UserData.zf = f.UserData.zf*0.9;
            return
        end
        % anti-clockwise rotate
        if strcmpi(eventData.Key,'rightarrow') && nMod==1 && strcmpi(m{1},'shift')
            camroll(ax, -5);
            camroll(ca, -5);
            return
        end
        % clockwise rotate
        if strcmpi(eventData.Key,'leftarrow')&&nMod==1&&strcmpi(m{1},'shift')
            camroll(ax, 5);
            camroll(ca, 5);
            return
        end
        % right shift
        if strcmpi(eventData.Key,'rightarrow')&&nMod==1&&strcmpi(m{1},'control')
            camdolly(ax,-0.1,0,0)
            return
        end
        % left shift
        if strcmpi(eventData.Key,'leftarrow')&&nMod==1&&strcmpi(m{1},'control')
            camdolly(ax,0.1,0,0)
            return
        end
        % up shift
        if strcmpi(eventData.Key,'uparrow')&&nMod==1&&strcmpi(m{1},'control')
            camdolly(ax,0,-0.1,0);
            return
        end
        % down shift
        if strcmpi(eventData.Key,'downarrow')&&nMod==1&&strcmpi(m{1},'control')
            camdolly(ax,0,0.1,0)
            return
        end
        % z axis view
        if nMod==0 && strcmpi(eventData.Key,'z')
            view(ax,[0,0,1]);
            view(ca,[0,0,1]);
            return
        end
        % x axis view
        if nMod==0 && strcmpi(eventData.Key,'x')
            view(ax,[1,0,0]);
            view(ca,[1,0,0]);
            return
        end
        % y axis view
        if nMod==0 && strcmpi(eventData.Key,'y')
            view(ax,[0,1,0]);
            view(ca,[0,1,0]);
            return
        end
        % isotropic rotation
        if nMod==0 && strcmpi(eventData.Key,'i')
            camorbit(ax, 30,30,'camera');
            camorbit(ca, 30,30,'camera');
            return
        end
        % 90 degree rotation
        if nMod==0 && strcmpi(eventData.Key,'space')
            camroll(ax,90)
            camroll(ca,90)
            return
        end
        % saving figure
        if strcmpi(eventData.Key,'s')&&nMod==1&&strcmpi(m{1},'control')
            bg.Visible = 'off';
            f.Renderer = 'painters';
            saveas(ax,'view.eps');
            bg.Visible = 'on';
            return
        end
        % zoom fit
        if strcmpi(eventData.Key,'f')
            camzoom(ax, 1/f.UserData.zf);
            f.UserData.zf = 1;
            return
        end
    end




end