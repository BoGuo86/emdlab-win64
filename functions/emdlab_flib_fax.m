% f: is figure
% a: is axes

function [f,ax] = emdlab_flib_fax(varargin)
switch numel(varargin)
    case 0
        f = figure('color', [0.9,0.9,0.9], 'position', [0,0,800,600], 'Visible','off');
        movegui(f, 'center');
        ax = axes(f, 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'normal');
        set(f, 'Visible','on');

    case 1

        if isa(varargin{1}, 'matlab.graphics.axis.Axes')
            ax = varargin{1};
            f = ax.Parent;
        else
            f = figure('color', [0.9,0.9,0.9], 'position', [0,0,800,600], 'Visible','off');
            movegui(f, 'center');
            ax = axes(f, 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'normal');
            set(f, 'Visible','on');
        end

    otherwise
        error('Too many input arguments');
end
end