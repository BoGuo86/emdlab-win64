function p = sortPointsOnLineU(p, u)
    tmp = p*u';
    [~, index] = sort(tmp);
    p = p(index,:);
end