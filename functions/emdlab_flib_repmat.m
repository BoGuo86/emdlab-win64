function out = emdlab_flib_repmat(obj, m, n)
    % Deep copy version of repmat for handle objects
    % Works like MATLAB's repmat for value classes, 
    % but makes independent copies for Copyable handle classes

    % Detect if it's a handle object
    if isa(obj, 'handle')
        % Detect if it supports copy()
        if ismethod(obj, 'copy')
            % Deep copy for each position
            copies = arrayfun(@(x) obj.copy(), 1:(m*n), 'UniformOutput', false);
            out = reshape([copies{:}], m, n);
        else
            error('Object is a handle but does not implement copy() method.');
        end
    else
        % For value classes or basic types, just use MATLAB's repmat
        out = repmat(obj, m, n);
    end
end
