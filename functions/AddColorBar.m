function AddColorBar(f, minValue, maxValue, cbTitle, nContour)

if nargin < 5
    nContour = 10;
end

% colorbar axes
cba = axes(f);
cba.Units = 'normalized';
cba.Position(2) = 0.2;
cba.Position(4) = 0.7;
cba.Units = 'pixels';
cba.Position(1) = 30;
cba.Position(3) = 20;

% colorbar image
m = linspace(minValue, maxValue, 10);
im = imagesc(cba, flipud(m'));

tmp = cell(length(m),1);
for i = 1:numel(tmp)
  tmp{i} = sprintf('%.2e',m(i));
end

cba.Tag = 'cba';
cba.XTick = [];
cba.YTickLabel = flipud(tmp);
cba.YAxisLocation = 'right';
set(cba, 'fontsize', 14);
set(cba, 'FontName', 'Consolas');
if nargin == 4 || nargin == 5
  cba.Title.HorizontalAlignment = 'left';
  cba.Title.String = cbTitle;
end
va = f.findobj('tag','va');
va.Colormap = repmat(linspace(0, 1, 5)',1,3);

colormap(va, jet(nContour));
colormap(cba, jet(nContour));
im.Tag = 'colorbar';

guidata(f, guihandles(f));

end
