
function TestFigure()

f = figure('MenuBar','none');

p1 = uipanel(f);
p2 = uipanel(f);
p3 = uipanel(f);

f.Units = 'pixels';
p1.Units = 'pixels';
p2.Units = 'pixels';
p3.Units = 'pixels';
p = f.Position;
    p1.Position = [0,0,100,p(4)];
    p2.Position = [100,0,5,p(4)];
    p3.Position = [105,0,p(3)-105,p(4)];
    
f.SizeChangedFcn = @fResizeFcn;
f.WindowButtonMotionFcn = @fMotionFcn;
f.UserData.isMouseHoled = false;
f.UserData.isMouseOverSliderHoled = false;
f.WindowButtonDownFcn = @fButtonDownFcn;
f.WindowButtonUpFcn = @fButtonUpFcn;

p1.BorderType = 'none';
p2.BorderType = 'none';
p3.BorderType = 'none';

% p1.BackgroundColor = 'r';
p2.BackgroundColor = [0.8,0.8,0.8];
% p3.BackgroundColor = 'b';


%     tgp = uitabgroup(p1);
%     tgp.Units = 'normalized';
%     tgp.Position = [0,0,1,1];
%     
%     uitab(tgp,'Title', 'Stesdgfs'); 
%     uitab(tgp,'Title', 'Stesss'); 

a = axes(p1);
a.Position = [0,0,1,1];
x = 0:0.01:2*pi;
plot(a,x,sin(x));
    
a = axes(p3);
a.Units = 'normalized';
a.Position = [0,0,1,1];
x = 0:0.01:2*pi;
plot(a,x,sin(x));

  function fResizeFcn(~,~)
    p = f.Position;
    if p(3)-p1.Position(3)-p2.Position(3)<1
      return
    end
    p1.Position(4) = p(4);
    p2.Position(4) = p(4);
    p3.Position(3) = p(3)-p1.Position(3)-p2.Position(3);
    p3.Position(4) = p(4);
  end

  function fMotionFcn(~,~)
    if f.UserData.isMouseOverSliderHoled
      p = f.Position;
      if f.CurrentPoint(1)>1 && f.CurrentPoint(1)<p(3)-7
      p1.Position = [0,0,f.CurrentPoint(1),p(4)];
      p2.Position = [f.CurrentPoint(1),0,5,p(4)];
      p3.Position = [f.CurrentPoint(1)+5,0,p(3)-f.CurrentPoint(1)-5,p(4)];
      end
      return
    end
    p = p2.Position;
    if f.CurrentPoint(1)>p(1) && f.CurrentPoint(1)<p(1)+p(3)
      f.Pointer = 'right'; 
    else
      f.Pointer = 'arrow';
    end
  end

  function fButtonDownFcn(~,~)
    f.UserData.isMouseHoled = true;
    p = p2.Position;
    if f.CurrentPoint(1)>p(1) && f.CurrentPoint(1)<p(1)+p(3)
      f.UserData.isMouseOverSliderHoled = true;
    end
  end

 function fButtonUpFcn(~,~)
    f.UserData.isMouseHoled = false;
    f.UserData.isMouseOverSliderHoled = false;
  end
end
