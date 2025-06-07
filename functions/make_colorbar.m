function make_colorbar(x,y,varargin)
x=0;
y=1;
f = gcf;
a = axes(f, 'Tag', 'ColorBar');
n = 20;

minV = x;
maxV = y;
m = linspace(minV,maxV,n)';
m = flipud(m);
p= imagesc(m);
colormap jet

a.Units = 'normalized';
a.Position = [0.05,0.3,0.05,0.6];
% axis(a, 'off')
a.UserData.x1 = 50;
a.UserData.x2 = 15;

a.Units = 'pixels';
a.Position(1) = a.UserData.x1;
a.Position(3) = a.UserData.x2;
a.Units = 'normalized';
% 
% a.YColorMode = 'manual';
% a.YColor = 'w';
% 
if numel(varargin)
title(varargin{1},'color','k');
end

set(a,'ylim',[1,n]);
tmp = linspace(1, n,8);
yticks(a, tmp);
tmp = linspace(minV, maxV, 8);
tmp = fliplr(tmp);

% tmp = round(,4);

tmpN = cell(1,numel(tmp));
for i = 1:numel(tmp)
  tmpN{i} = sprintf('%.2e',tmp(i));
end
yticklabels(a, tmpN);

tmp = [];
xticks(a, tmp);
xticklabels(a, tmp);

set(a, 'fontsize', 14)
set(a, 'FontName', 'Consolas')
% set(a, 'color', 'w')
a.YAxisLocation = 'right';


end
