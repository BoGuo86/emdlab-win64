


function highlightPatchCallback(f, e)

% f.ButtonDownFcn( f, []);
% f.WindowButtonDownFcn( f, []);
h = guihandles(f);
% 
% feval(h.a.ButtonDownFcn, h.a, []);

names = fieldnames(h);

for i = 1:numel(names)
  
  optr = h.(names{i});
  
  if isa(optr, 'matlab.graphics.primitive.Patch')
    optr.FaceColor = 'c';
  end
end

if isa(f.CurrentObject, 'matlab.graphics.primitive.Patch')
f.CurrentObject.FaceColor = 'r';
end
