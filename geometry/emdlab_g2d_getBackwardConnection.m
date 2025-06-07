function y = emdlab_g2d_getBackwardConnection(v)

y = (size(v,1):-1:1)';
y = [y,[y(2:end);y(1)]];

end