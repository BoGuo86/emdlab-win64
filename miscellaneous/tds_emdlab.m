

function tds_main

f = GraphicWindow;

%% definition of panels
[jTopPanel, wTopPanel] = javacomponent(javax.swing.JPanel, [], f);
[jLeftPanel, wLeftPanel] = javacomponent(javax.swing.JPanel, [], f);
[jRightPanel, wRightPanel] = javacomponent(javax.swing.JPanel, [], f);
[jBottomPanel, wBottomPanel] = javacomponent(javax.swing.JPanel, [], f);
[jHintPanel, wHintPanel] = javacomponent(javax.swing.JPanel, [], f);
[jLeftDivider, wLeftDivider] = javacomponent(javax.swing.JPanel, [], f);
[jRightDivider, wRightDivider] = javacomponent(javax.swing.JPanel, [], f);
[jBottomDivider, wBottomDivider] = javacomponent(javax.swing.JPanel, [], f);
viewPanel = uipanel(f, 'units', 'pixels');

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
t = tds_getTree();
st = javax.swing.JScrollPane(t);
p = javax.swing.JPanel;
s = javax.swing.JSplitPane(javax.swing.JSplitPane.VERTICAL_SPLIT, st, p);
s.setResizeWeight(0.1);
s.setContinuousLayout(true);
jLeftPanel.setLayout(java.awt.BorderLayout);
jLeftPanel.add(s, java.awt.BorderLayout.CENTER);


%% setting of bottom panel

% command text area
cmdtxa = javax.swing.JTextArea;
sp_cmdtxa = javax.swing.JScrollPane(cmdtxa);

% message text area
msgtxa = javax.swing.JTextArea;
sp_msgtxa = javax.swing.JScrollPane(msgtxa);

% progress panel and progress text area
progp = javax.swing.JPanel;
progtxa = javax.swing.JTextArea;
sp_progtxa = javax.swing.JScrollPane(progtxa);
progbar = javax.swing.JProgressBar;
progbar.setStringPainted(true);
progp.setLayout(java.awt.BorderLayout);
progp.add(progbar, java.awt.BorderLayout.NORTH);
progp.add(sp_progtxa, java.awt.BorderLayout.CENTER);

s = javax.swing.JSplitPane(javax.swing.JSplitPane.HORIZONTAL_SPLIT, sp_cmdtxa, sp_msgtxa);
s.setResizeWeight(0.5);
s.setContinuousLayout(true);
s = javax.swing.JSplitPane(javax.swing.JSplitPane.HORIZONTAL_SPLIT, s, progp);
s.setResizeWeight(0.8);
s.setContinuousLayout(true);

jBottomPanel.setLayout(java.awt.BorderLayout);
jBottomPanel.add(s, java.awt.BorderLayout.CENTER);

%% setting right panel

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
        
        leftMargin = wLeftPanel.Position(3) + wLeftDivider.Position(3);
        rightMargin = wRightPanel.Position(3) + wRightDivider.Position(3);
        topMargin = wTopPanel.Position(4);
        bottomMargin = wHintPanel.Position(4) + wBottomPanel.Position(4) + wBottomDivider.Position(4);
        viewPanel.Position = [leftMargin, bottomMargin, f_p(3)-(leftMargin+rightMargin), f_p(4)-(bottomMargin+topMargin) ];
        
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