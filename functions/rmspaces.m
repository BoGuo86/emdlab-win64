function y = rmspaces(x)
% removing extra spaces
if ~ischar(x)
    error('Input must be string.');
end
y = strjoin(strsplit(x,' '),'');
end