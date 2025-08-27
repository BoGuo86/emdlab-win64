% f: is figure
% a: is axes
% N: is default number of contours

function [f,ax,N] = emdlab_flib_faxN(varargin)

switch numel(varargin)
    case 0
        f = figure('color', [0.9,0.9,0.9], 'position', [0,0,800,600], 'Visible','off');
        movegui(f, 'center');
        ax = axes(f, 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'normal');
        set(f, 'Visible','on');
        N = 12;

    case 1

        if isa(varargin{1}, 'matlab.graphics.axis.Axes')
            ax = varargin{1};
            f = ax.Parent;
            N = 12;

        elseif isscalar(varargin{1}) && isinteger(varargin{1}) && varargin{1}>0
            N = varargin{1};
            f = figure('color', [0.9,0.9,0.9], 'position', [0,0,800,600], 'Visible','off');
            movegui(f, 'center');
            ax = axes(f, 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'normal');
            set(f, 'Visible','on');

        else
            error('Wrong input type');
        end

    case 2

        if isa(varargin{1}, 'matlab.graphics.axis.Axes')
            ax = varargin{1};
            f = ax.Parent;
        else
            error('Wrong input type');
        end

        if isscalar(varargin{2}) && isinteger(varargin{2}) && varargin{2}>0
            N = varargin{2};
        else
            error('Wrong input type');
        end

    otherwise
        error('Too many input arguments');
end
end