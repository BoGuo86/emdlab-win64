function p = sortPointsOnLineY(p)
    isPointsOnLine(p);
    [~, index] = sort(p(:,2));
    p = p(index,:);
end