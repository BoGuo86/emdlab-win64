
function test_java


f = figure('menu','figure','numbertitle', 'off', 'name', 'NASR Company', 'visible', 'off', 'handlevisibility', 'off');
% uimenu(f, 'text', 'File');
% uimenu(f, 'text', 'Help');
drawnow;

bp = javax.swing.JPanel;
[~, wbg] = javacomponent(bp, [], f);
set(wbg, 'units', 'normalized', 'position', [0,0,1,1]);
uistack(wbg, 'bottom')
% w.BackgroundColor = 'r';
drawnow;

top_p = javax.swing.JPanel;
mid_p = javax.swing.JPanel;
bottom_p = javax.swing.JPanel;

top_p.setPreferredSize(java.awt.Dimension(10,60))
bottom_p.setPreferredSize(java.awt.Dimension(10,20))

bp.setLayout(java.awt.BorderLayout);
bp.add(top_p, java.awt.BorderLayout.NORTH);
bp.add(mid_p, java.awt.BorderLayout.CENTER);
bp.add(bottom_p, java.awt.BorderLayout.SOUTH);



p1 = javax.swing.JPanel;
p2 = javax.swing.JPanel;

vp = uipanel(f,'BorderType', 'none', 'clipping', 'off');
vp.Units = 'pixels'; vp.Position = [0,0,600,400];
% vp.Units = 'norm'; vp.Position = [0,0,1,1];
% tgp = uitabgroup(vp);
% ta = uitab(tgp);
drawnow

% p3 = vp.JavaFrame.getGUIDEView.getParent;
% p3.setSize(java.awt.Dimension(600,400));
p3 = javax.swing.JPanel;
% sp3 = javax.swing.JScrollPane(p3);
pp3 = javax.swing.JPanel;
% p3.setMinimumSize(java.awt.Dimension(200,200));
% p3.setMaximumSize(java.awt.Dimension(600,400));
% sp3.setMaximumSize(java.awt.Dimension(600,400));

ss3 = javax.swing.JSplitPane(javax.swing.JSplitPane.HORIZONTAL_SPLIT, p3, pp3);

% ss3.getLeftComponent.setMinimumSize(java.awt.Dimension(300,300));
% ss3.getLeftComponent.setMaximumSize(java.awt.Dimension(600,400));

p4 = javax.swing.JPanel;

p5 = javax.swing.JPanel;
p6 = javax.swing.JPanel;

s1 = javax.swing.JSplitPane(javax.swing.JSplitPane.HORIZONTAL_SPLIT, p1, p2);
s2 = javax.swing.JSplitPane(javax.swing.JSplitPane.VERTICAL_SPLIT, ss3, p4);
s3 = javax.swing.JSplitPane(javax.swing.JSplitPane.HORIZONTAL_SPLIT, p5, p6);

s1.setContinuousLayout(true)
s2.setContinuousLayout(true)
s3.setContinuousLayout(true)
ss3.setContinuousLayout(true)

p2.setLayout(java.awt.BorderLayout);
p2.add(s2, java.awt.BorderLayout.CENTER);
p11 = javax.swing.JPanel;
p12 = javax.swing.JPanel;

t = javax.swing.JTree;
sp = javax.swing.JScrollPane(t);
p11.setLayout(java.awt.BorderLayout);
p11.add(sp, java.awt.BorderLayout.CENTER);

s11 = javax.swing.JSplitPane(javax.swing.JSplitPane.VERTICAL_SPLIT, p11, p12);

s11.setContinuousLayout(true)
p1.setLayout(java.awt.BorderLayout);
p1.add(s11, java.awt.BorderLayout.CENTER);

mid_p.setLayout(java.awt.BorderLayout);
mid_p.add(s1, java.awt.BorderLayout.CENTER);

p4.setLayout(java.awt.BorderLayout);
p4.add(s3, java.awt.BorderLayout.CENTER);

p = javax.swing.JProgressBar;

p6.setLayout(java.awt.BorderLayout);
p6.add(p, java.awt.BorderLayout.NORTH);

b = javax.swing.JButton('RUN Simulation');

top_p.setLayout(java.awt.GridLayout(1,10));
top_p.add(b);
set(b, 'MouseClickedCallback', @my_f);
set(b, 'MouseEnteredCallback', @(src,evn) b.setForeground(java.awt.Color.RED));

for i = 1:9
    top_p.add(javax.swing.JButton('RUN Simulation'))
end

tx = javax.swing.JTextArea;
stx = javax.swing.JScrollPane(tx);
p5.setLayout(java.awt.BorderLayout);
p5.add(stx, java.awt.BorderLayout.CENTER);
tx.setText('salam ali ..');
tx.setEditable(false);

ptx = javax.swing.JTextArea;
sptx = javax.swing.JScrollPane(ptx);

p6.add(sptx, java.awt.BorderLayout.CENTER);


p.setStringPainted(true)
s1.setResizeWeight(0.2);
s2.setResizeWeight(0.8);
s3.setResizeWeight(0.5);
s11.setResizeWeight(0.5);
ss3.setResizeWeight(0.7);



ff = GraphicWindow;
hh = guidata(ff);
hh.bg.Parent = vp;
hh.bg.HandleVisibility = 'off';
aaa = hh.va;
aaa.Parent = vp;
delete(ff);
% aaa = axes(vp, 'units', 'pixels', 'position', [0,0,500,500]);
% aaa = axes(vp, 'units', 'norm');
sphere(aaa)
% surf(aaa,peaks)
aaa.Clipping = 'off';
% addlistener(aaa, 'SizeChanged', @(src,evn) axis(aaa,'equal'))
% camzoom(aaa, 0.8);
aaa.Units = 'pixels';
aaa.Position(3) = 500;
aaa.Position(4) = 500;
axis(aaa, 'off', 'equal')
aaa.Toolbar.Visible = 'on';
set(p3, 'ComponentResizedCallback',  @ali_ff)

% set(aaa,'units', 'norm', 'position', [0,0,1,1]);



% bp.setBackground(java.awt.Color(0,0,0,0))
% bottom_p.setBackground(java.awt.Color(0,0,0,0))
f.Visible = 'on';

% uistack([wbg, vp, aaa])
    function my_f(~,e)
        
        wf = figure('menu', 'none', 'windowstyle', 'modal', 'position',[0,0,200,150]);
        movegui(wf, 'center');
        for i = 0:100
            p.setValue(i)
            ptx.setText(['percent ', num2str(i), ' %']);
            pause(0.005)
        end
        ptx.setText('Completed!');
        delete(wf);
        
        
    end

    function ali_ff(~,~)
%         tmp = getpixelposition(vp);
%         aaa.Position(1) = tmp(3)/2-250;
%         aaa.Position(2) = tmp(4)/2-250;
%         
% p3.getSize.getWidth
%         tmp = getpixelposition(vp)
% xdim = p3.getSize.getWidth;
% ydim = p3.getSize.getHeight;
% x = 0.95*min(xdim,ydim);
%         tmpx = x/xdim;
%         tmpy = x/ydim;
%         aaa.Position = [0.5-tmpx/2, 0.5-tmpy/2, tmpx, tmpy];
%         getpixelposition(aaa);
%         axis(aaa, 'equal')
% %         drawnow
%         vp.Position
%         p3.getBounds
vp.Position(1) = t.getSize.getWidth + 13;
vp.Position(2) = tx.getSize.getHeight + 33;
        vp.Position(3) = p3.getSize.getWidth;
        vp.Position(4) = p3.getSize.getHeight;
aaa.Position(1) = (vp.Position(3) - aaa.Position(3))/2;
aaa.Position(2) = (vp.Position(4) - aaa.Position(4))/2;
    end

end