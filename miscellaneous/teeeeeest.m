

f = figure('menubar','none');
a = axes(f, 'tag', 'a');
p=patch(a,'faces',[1,2,3],'vertices',[0,0;1,0;0,1],'tag','sg');
patch(a,'faces',[1,2,3],'vertices',[1,0;1,1;0,1],'tag','sdg');


% a.ButtonDownFcn = @selectPatchCallback;
% f.ButtonDownFcn = @(h,e) get(h,'selectiontype');
% f.WindowButtonDownFcn = @(h,e) get(h,'selectiontype');

f.WindowButtonMotionFcn = @highlightPatchCallback;

h = guihandles(f);
guidata(f,h);

axis off equal
