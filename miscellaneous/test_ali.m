

function test_ali

% f = figure('menu', 'figure', 'color', 'g');
% a = axes(f);
% sphere(a);
% axis(a, 'off', 'equal');
% f = GraphicWindow;
%% Initialization

materialdir = [cd,'\MaterialsData'];
%% mesh
m = TMDBC;
glib_srm_st1_3dt(48.18,78.4,89.8, 8, 20.2*pi/180, 5, 5, 6);
m.read_g2d_bin('geom.g2d', 'MG0', 10);
ew = m.mzs.ew.getCopy;
m.removemz('ew');
% extursion vectors
z1 = linspace(0,75,16);
z2 = [z1, 80];
z3 = [z1, 80, 85, 90, 95, 100];
tm = TTMDBC;
m.ggmesh;
m.getMeshZoneExtrude(tm, 'stator', z1);
m.getMeshZoneExtrude(tm, 'coil', z2);
m.getMeshZoneExtrude(tm, 'housing', z3);
m.addmz('ew', ew);
m.removemz('stator');
m.jmzs('ewin', 'ew', 'coil');
m.ggmesh
m.getMeshZoneExtrude(tm, 'ewin', [80, 85, 90]);
% assiging named selections
tm.ggmesh;
tm.addFacetNamedSelection('ins',[tm.getfbfiop([0,0,0],[0,0,1]);tm.getfbfiop([0,0,0],[0,1,0]);tm.getfbfiop([0,0,0],[-sin(pi/8),cos(pi/8),0])]);
tm.addFacetNamedSelection('oc',tm.getbfioic([0,0,0],[0,0,1],89.8+5));
tm.addFacetNamedSelection('ic',tm.getbfioc([0,0,0],[0,0,1],48.18 ,0, 75));
% assigning materials
tm.addMaterial(materialdir,'aluminium');
tm.setMaterial('housing','aluminium');
tm.addMaterial(materialdir,'m290_50a');
tm.setMaterial('stator','m290_50a');
tm.addMaterial(materialdir,'ew');
tm.setMaterial('ewin','ew');
tm.addMaterial(materialdir,'eqwin');
tm.setMaterial('coil','eqwin');
%% solver
s = NIHLTHSTTL4(tm);
s.setUnit('length', 'mm');
s.setUnit('volumeLossDensity', 'w/mm^3');
s.setUnit('power', 'w');
%% process
% boundary conditions
s.bcs.clearAllBCs;
s.setnfConvectionBC('oc', 5, 25);
s.setnfConvectionBC('ic', 60, 25);
s.setnfConvectionBC('none', 5, 25);
% internal losses
s.setInternalLoss('stator',5,'l');
s.setInternalLoss('coil',2,'l');
% calling solver
s.solve;
f=s.plotTemperature('w');


%% definition of panels
[jTopPanel, wTopPanel] = javacomponent(javax.swing.JPanel, [], f);
[jLeftPanel, wLeftPanel] = javacomponent(javax.swing.JPanel, [], f);
[jRightPanel, wRightPanel] = javacomponent(javax.swing.JPanel, [], f);
[jBottomPanel, wBottomPanel] = javacomponent(javax.swing.JPanel, [], f);
[jHintPanel, wHintPanel] = javacomponent(javax.swing.JPanel, [], f);
[jLeftDivider, wLeftDivider] = javacomponent(javax.swing.JPanel, [], f);
[jRightDivider, wRightDivider] = javacomponent(javax.swing.JPanel, [], f);
[jBottomDivider, wBottomDivider] = javacomponent(javax.swing.JPanel, [], f);

%% setting initial sizes of panels
wLeftPanel.Position(3) = 200;
wTopPanel.Position(4) = 80;
wHintPanel.Position(4) = 20;
wRightPanel.Position(3) = 200;
wBottomPanel.Position(4) = 100;
wBottomDivider.Position(3) = 5;
wRightDivider.Position(3) = 5;
wRightDivider.Position(4) = 5;

setLayout;
f.SizeChangedFcn = @setLayout;

set(jLeftDivider, 'MouseEnteredCallback', @leftDividerMouseEneter_cbk);
set(jLeftDivider, 'MouseExitedCallback', @leftDividerMouseExit_cbk);
set(jLeftDivider, 'MouseReleasedCallback', @leftDividerMouseReleased_cbk);
set(jLeftDivider, 'MouseDraggedCallback', @leftDividerMouseDragged_cbk);

set(jRightDivider, 'MouseEnteredCallback', @leftDividerMouseEneter_cbk);
set(jRightDivider, 'MouseExitedCallback', @leftDividerMouseExit_cbk);
set(jRightDivider, 'MouseReleasedCallback', @rightDividerMouseReleased_cbk);
set(jRightDivider, 'MouseDraggedCallback', @rightDividerMouseDragged_cbk);

set(jBottomDivider, 'MouseEnteredCallback', @bottomDividerMouseEneter_cbk);
set(jBottomDivider, 'MouseExitedCallback', @bottomDividerMouseExit_cbk);
set(jBottomDivider, 'MouseReleasedCallback', @bottomDividerMouseReleased_cbk);
set(jBottomDivider, 'MouseDraggedCallback', @bottomDividerMouseDragged_cbk);

%% setting of left panel
t = javax.swing.JTree;
st = javax.swing.JScrollPane(t);
p = javax.swing.JPanel;
s = javax.swing.JSplitPane(javax.swing.JSplitPane.VERTICAL_SPLIT, st, p);
s.setResizeWeight(0.1);
s.setContinuousLayout(true);
jLeftPanel.setLayout(java.awt.BorderLayout);
jLeftPanel.add(s, java.awt.BorderLayout.CENTER);


%% setting of right panel
msgtxa = javax.swing.JTextArea;
sp_msgtxa = javax.swing.JScrollPane(msgtxa);

progp = javax.swing.JPanel;
progtxa = javax.swing.JTextArea;
sp_progtxa = javax.swing.JScrollPane(progtxa);
progbar = javax.swing.JProgressBar;
progbar.setStringPainted(true);
progp.setLayout(java.awt.BorderLayout);
progp.add(progbar, java.awt.BorderLayout.NORTH);
progp.add(sp_progtxa, java.awt.BorderLayout.CENTER);

s = javax.swing.JSplitPane(javax.swing.JSplitPane.HORIZONTAL_SPLIT, sp_msgtxa, progp);
s.setResizeWeight(0.6);
s.setContinuousLayout(true);
jBottomPanel.setLayout(java.awt.BorderLayout);
jBottomPanel.add(s, java.awt.BorderLayout.CENTER);

%% setting right panel
t = javax.swing.JTree;
jRightPanel.setLayout(java.awt.BorderLayout);
jRightPanel.add(t, java.awt.BorderLayout.CENTER);

    function setLayout(~,~)
        
        f_p = getpixelposition(f);
        
        wTopPanel.Position = [0, f_p(4)-wTopPanel.Position(4), f_p(3), wTopPanel.Position(4)];
        
        wHintPanel.Position = [0, 0, f_p(3), wHintPanel.Position(4)];
        
        wLeftPanel.Position = [0, wHintPanel.Position(4), wLeftPanel.Position(3), f_p(4)-wTopPanel.Position(4)-wHintPanel.Position(4)];
        
        wLeftDivider.Position = [wLeftPanel.Position(3), wHintPanel.Position(4), 5, f_p(4)-wTopPanel.Position(4)-wHintPanel.Position(4)];
        
        wBottomPanel.Position = [wLeftPanel.Position(3)+5, wHintPanel.Position(4), f_p(3)-wLeftPanel.Position(3)-5, wBottomPanel.Position(4)];
        
        wBottomDivider.Position = [wLeftPanel.Position(3)+5, wHintPanel.Position(4)+wBottomPanel.Position(4), f_p(3)-wLeftPanel.Position(3)-5, 5];
        
        wRightPanel.Position = [f_p(3)-wRightPanel.Position(3), wBottomPanel.Position(4)+wHintPanel.Position(4)+5, wRightPanel.Position(3), f_p(4)-wTopPanel.Position(4)-wBottomPanel.Position(4)-wHintPanel.Position(4)-5];
        
        wRightDivider.Position = [f_p(3)-wRightPanel.Position(3)-5, wBottomPanel.Position(4)+wHintPanel.Position(4)+5, 5, f_p(4)-wTopPanel.Position(4)-wBottomPanel.Position(4)-wHintPanel.Position(4)-5];
        
        h = guidata(f);
        h.ca.Position(1) = wLeftPanel.Position(3) + 5;
        h.ca.Position(2) = wBottomPanel.Position(4) + wHintPanel.Position(4) + 5;
        
    end

    function leftDividerMouseEneter_cbk(~,~)
        f.Pointer = 'right';
    end

    function leftDividerMouseExit_cbk(~,~)
        f.Pointer = 'arrow';
    end

    function leftDividerMouseDragged_cbk(~,eventData)
        f.Pointer = 'right';
        wLeftPanel.Position(3) =eventData.getXOnScreen - f.Position(1);
        setLayout;
    end

    function leftDividerMouseReleased_cbk(~,eventData)
        f.Pointer = 'arrow';
        wLeftPanel.Position(3) =eventData.getXOnScreen - f.Position(1);
        setLayout;
    end

    function rightDividerMouseDragged_cbk(~,eventData)
        f.Pointer = 'right';
        wRightPanel.Position(3) = f.Position(3) - (eventData.getXOnScreen - f.Position(1)) - 5;
        setLayout;
    end

    function rightDividerMouseReleased_cbk(~,eventData)
        f.Pointer = 'arrow';
        wRightPanel.Position(3) = f.Position(3) - (eventData.getXOnScreen - f.Position(1)) - 5;
        setLayout;
    end

    function bottomDividerMouseEneter_cbk(~,~)
        f.Pointer = 'top';
    end

    function bottomDividerMouseExit_cbk(~,~)
        f.Pointer = 'arrow';
    end

    function bottomDividerMouseDragged_cbk(~,eventData)
        f.Pointer = 'top';
        tmp = get(0, 'screensize');
        tmp = tmp(4) - eventData.getYOnScreen - wHintPanel.Position(4) - f.Position(2);
        wBottomPanel.Position(4) = tmp;
        setLayout;
    end

    function bottomDividerMouseReleased_cbk(~,eventData)
        f.Pointer = 'arrow';
        tmp = get(0, 'screensize');
        tmp = tmp(4) - eventData.getYOnScreen - wHintPanel.Position(4) - f.Position(2);
        wBottomPanel.Position(4) = tmp;
        setLayout;
    end
end