function p = sortPointsOnLineX(p)
    isPointsOnLine(p);
    [~, index] = sort(p(:,1));
    p = p(index,:);
end