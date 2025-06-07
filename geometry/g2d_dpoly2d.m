function y = g2d_dpoly2d(f, v, p)

y = ext_dpoly2dnew(p,f,v);
[~,index] = min(abs(y),[],2);

end