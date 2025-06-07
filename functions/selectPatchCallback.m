function selectPatchCallback(h, e)

    if e.Button == 1

        if isequal('FaceAlpha', 1)
            return;
        end

        % view axis
        a = h.Parent;
        % patch
        p = a.findobj('FaceAlpha', 1);

        if ~ isempty(p)
            set(p, 'FaceColor', p.UserData.color, 'FaceAlpha', 0.8);
        end

        set(h, 'FaceColor', 'r', 'FaceAlpha', 1);
    end

end
