function names = emdlab_flib_varargin2StringList(varargin)
if all(cellfun(@ischar, varargin))
    names = string(varargin);
else
    names = string([varargin{:}]);
end
end