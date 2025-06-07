function f = emdlab_r2d_mesh()

% generate a new figure
f = figure('Visible', 'off');
set(f, 'position', [0,0,800,600]);
movegui(f, 'center');
set(f, 'Tag', 'main');
set(f, 'Name', 'EMDLAB: 2D RENDERING');
set(f, 'clipping', 'off');
set(f, 'color', [0.9,0.9,0.9]);
set(f, 'Units', 'pixel');

% saving gui components
guidata(f, guihandles(f));
f.Visible = 'on';
zoom on;

end
