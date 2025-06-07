function setggg

f = gcf;
f.MenuBar = 'none';

a = gca;
a.Clipping = 'off';

f.UserData.z = zoom;
f.UserData.r = rotate3d;


set(gcf,'WindowButtonDownFcn',@RotateOn)
set(gcf,'WindowButtonUpFcn',@RotateOff) 


set(gcf,'WindowScrollWheelFcn',@ZoomIn)


% mm = uimenu(gcf,'Label','Zoom');
% uimenu(mm,'Label','Zoom On','CallBack','zoom on','Accelerator','z');
% uimenu(mm,'Label','Zoom Off','CallBack','zoom off','Accelerator','c');
% uimenu(mm,'Label','Set View','CallBack','zoom out','Accelerator','x');
% set(gcf,'UserData',inf)

end



function RotateOn(src,evn)
f = gcf;
f.UserData.r.Enable = 'on';
end

function RotateOff(src,evn)
f = gcf;
f.UserData.r.Enable = 'off';
end

function ZoomIn(src,evn)

rotate3d off;
z = zoom;
z.Enable = 'on';

if evn.VerticalScrollCount>0
    z.Direction = 'In';
else
    z.Direction = 'Out';
end

pause(0.5)
z.Enable = 'off';
rotate3d on;
end