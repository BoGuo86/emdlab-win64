function y = emdlab_g2d_getForwardConnection(v)

y = (1:size(v,1))';
y = [y,[y(2:end);y(1)]];

end