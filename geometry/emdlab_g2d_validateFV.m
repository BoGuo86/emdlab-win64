function emdlab_g2d_validateFV(f,v)

f = f(:);
if any(f<0) || max(f)>size(v,1)
    error('illegal g2d file: incorrect edge.');
end
f = sort(f);
f = reshape(f,2,[]);
f = sum(diff(f));

% if f~= 0
%   error('illegal g2d file.');
% end
end
